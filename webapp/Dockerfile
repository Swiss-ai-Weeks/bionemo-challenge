FROM python:3.12-slim
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 APP_DATA_DIR=/workspace-data
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends libxrender1 libxext6 libexpat1 && rm -rf /var/lib/apt/lists/*
COPY webapp/requirements.txt /app/webapp/requirements.txt
RUN pip install --no-cache-dir -r webapp/requirements.txt
COPY webapp /app/webapp
COPY scripts/download_bindingdb.py /app/scripts/download_bindingdb.py
COPY assets/vendor /app/assets/vendor
RUN useradd -u 1001 -m app && mkdir -p /workspace-data && chown app:app /workspace-data
USER app
EXPOSE 8080
CMD ["uvicorn", "webapp.main:app", "--host", "0.0.0.0", "--port", "8080", "--workers", "1"]
