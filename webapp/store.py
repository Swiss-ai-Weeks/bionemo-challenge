"""Small SQLite store. Every operation owns a short-lived connection."""
import json
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat()


def uid():
    return uuid.uuid4().hex


class Store:
    def __init__(self, folder):
        self.folder = Path(folder)
        self.folder.mkdir(parents=True, exist_ok=True)
        self.path = self.folder / 'workspace.sqlite3'
        with self.connect() as db:
            db.executescript('''
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS targets(id TEXT PRIMARY KEY,name TEXT,sequence TEXT,accession TEXT);
            CREATE TABLE IF NOT EXISTS variants(id TEXT PRIMARY KEY,target_id TEXT,name TEXT,sequence TEXT,changes TEXT);
            CREATE TABLE IF NOT EXISTS ligands(id TEXT PRIMARY KEY,name TEXT,smiles TEXT UNIQUE,svg TEXT,mw REAL);
            CREATE TABLE IF NOT EXISTS measurements(target_id TEXT,ligand_id TEXT,endpoint TEXT,value TEXT,unit TEXT,source TEXT,
                UNIQUE(target_id,ligand_id,endpoint,value,source));
            CREATE TABLE IF NOT EXISTS runs(id TEXT PRIMARY KEY,name TEXT,created TEXT,cancelled INTEGER DEFAULT 0);
            CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY,run_id TEXT,snapshot TEXT,status TEXT,attempt INTEGER DEFAULT 0,
                created TEXT,started TEXT,finished TEXT,score TEXT,error TEXT);
            CREATE TABLE IF NOT EXISTS attempts(job_id TEXT,number INTEGER,status TEXT,started TEXT,finished TEXT,error TEXT,
                PRIMARY KEY(job_id,number));
            CREATE TABLE IF NOT EXISTS benchmarks(run_id TEXT PRIMARY KEY,config TEXT,state TEXT,sampling TEXT,error TEXT,created TEXT);
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);
            CREATE INDEX IF NOT EXISTS jobs_queue ON jobs(status,created);
            CREATE INDEX IF NOT EXISTS jobs_run ON jobs(run_id);
            ''')

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=20)
        db.row_factory = sqlite3.Row
        try:
            yield db
            db.commit()
        except BaseException:
            db.rollback()
            raise
        finally:
            db.close()

    def all(self, sql, params=()):
        with self.connect() as db:
            return [dict(row) for row in db.execute(sql, params).fetchall()]

    def one(self, sql, params=()):
        rows = self.all(sql, params)
        return rows[0] if rows else None

    def meta(self, key, default=None):
        row = self.one('SELECT value FROM meta WHERE key=?', (key,))
        return json.loads(row['value']) if row else default

    def set_meta(self, key, value):
        with self.connect() as db:
            db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)', (key, json.dumps(value)))

    def recover(self):
        with self.connect() as db:
            active = db.execute("SELECT id,attempt FROM jobs WHERE status='running'").fetchall()
            for job in active:
                db.execute("UPDATE jobs SET status='interrupted',error=?,finished=? WHERE id=?",
                           ('Application stopped during inference; backend completion is uncertain.', now(), job['id']))
                db.execute("UPDATE attempts SET status='interrupted',finished=? WHERE job_id=? AND number=?",
                           (now(), job['id'], job['attempt']))
            if active:
                db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)', ('paused', json.dumps(True)))

    def claim(self):
        with self.connect() as db:
            db.execute('BEGIN IMMEDIATE')
            paused = db.execute("SELECT value FROM meta WHERE key='paused'").fetchone()
            if paused and json.loads(paused['value']):
                return None
            if db.execute("SELECT 1 FROM jobs WHERE status='running' LIMIT 1").fetchone():
                return None
            job = db.execute("SELECT jobs.* FROM jobs JOIN runs ON runs.id=jobs.run_id LEFT JOIN benchmarks ON benchmarks.run_id=runs.id WHERE jobs.status='queued' AND runs.cancelled=0 AND (benchmarks.run_id IS NULL OR benchmarks.state='ready') ORDER BY CASE WHEN benchmarks.run_id IS NULL THEN 0 ELSE 1 END, jobs.created, jobs.rowid LIMIT 1").fetchone()
            if not job:
                return None
            result = dict(job)
            result['attempt'] += 1
            result['started'] = now()
            db.execute("UPDATE jobs SET status='running',attempt=?,started=?,finished=NULL,error=NULL WHERE id=?",
                       (result['attempt'],result['started'],result['id']))
            db.execute('INSERT INTO attempts(job_id,number,status,started) VALUES (?,?,?,?)',
                       (result['id'],result['attempt'],'running',result['started']))
            return result

    def finish(self, job, status, score=None, error=None, pause=False):
        with self.connect() as db:
            db.execute('UPDATE jobs SET status=?,score=?,error=?,finished=? WHERE id=?',
                       (status,json.dumps(score) if score else None,error,now(),job['id']))
            db.execute('UPDATE attempts SET status=?,finished=?,error=? WHERE job_id=? AND number=?',
                       (status,now(),error,job['id'],job['attempt']))
            if pause:
                db.execute('INSERT OR REPLACE INTO meta VALUES (?,?)', ('paused', 'true'))

    def run(self, run_id):
        run = self.one('SELECT * FROM runs WHERE id=?', (run_id,))
        if not run:
            return None
        jobs = self.all('SELECT * FROM jobs WHERE run_id=? ORDER BY created,rowid', (run_id,))
        for job in jobs:
            job['snapshot'] = json.loads(job['snapshot'])
            job['score'] = json.loads(job['score']) if job['score'] else None
            job['rank'] = None
            job['has_structure'] = (self.folder/'artifacts'/job['id']/str(job['attempt'])/'complex.cif').is_file()
        variants = {j['snapshot']['variant_id'] for j in jobs}
        for v in variants:
            group = sorted([j for j in jobs if j['snapshot']['variant_id']==v and j['status']=='succeeded'],
                           key=lambda j: j['score']['raw'])
            for j in group:
                positions=[i+1 for i,k in enumerate(group) if k['score']['raw']==j['score']['raw']]
                j['rank']=sum(positions)/len(positions)
        counts = {state: sum(j['status']==state for j in jobs) for state in
                  ['queued','running','succeeded','failed','cancelled','interrupted']}
        run['counts'] = counts
        run['total'] = len(jobs)
        run['status'] = ('running' if counts['running'] else 'queued' if counts['queued'] else
                         'cancelled' if run['cancelled'] else 'completed with errors' if counts['failed'] or counts['interrupted'] else 'complete')
        run['jobs'] = jobs
        return run
