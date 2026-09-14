"""Independent inference worker; web app restarts do not interrupt GPU jobs."""
import os
import fcntl
from pathlib import Path
import signal
import threading
from webapp.service import Service
from webapp.store import Store


def main():
    store=Store(os.getenv('APP_DATA_DIR','app-data'))
    lock=(store.folder/'worker.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
    service=Service(store,os.getenv('BOLTZ_URL','http://boltz2:8000'),bootstrap=False)
    stopped=threading.Event()
    def shutdown(*_):stopped.set()
    signal.signal(signal.SIGTERM,shutdown)
    signal.signal(signal.SIGINT,shutdown)
    service.start(worker=True,prepare=False)
    try:
        stopped.wait()
    finally:
        service.close()


if __name__=='__main__':main()
