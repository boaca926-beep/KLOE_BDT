# Dockerfile
FROM python:3.12-slim # slim = smaller size, fewer pre-installed packages

# Set working directory
WORKDIR /app # Creates and switches to /app directory inside container

# Installs system packages needed for Python MySQL libraries
# gcc: C compiler (compiles MySQL Python package)
# default-libmysqlclient-dev: MySQL client library files
# pkg-config: Helps find MySQL libraries
# curl: Used for health checks
# && and rm -rf cleans up to keep image smaller
RUN apt-get update && apt-get install -y \
    gcc \
    default-libmysqlclient-dev \
    pkg-config \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first (for better caching)
# --no-cache-dir = don't save downloaded packages
# Typical packages: uvicorn, fastapi, mysql-connector-python, etc.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code, Python files into the container
COPY api_mysql.py .
COPY mysql_db.py .
COPY init_bdt.py .
COPY models/ ./models/

# Expose port. Documents that the container listens on port 8000
# Actual publishing happens in docker-compose.yml: ports: - "8000:8000"
EXPOSE 8000

# Health check
# Tells Docker how to check if the container is healthy
# Check every 30 seconds; Wait max 3s for response; Give 5s before first check; Fail after 3 failed checks
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

# Run the application
# Default command when container starts
# In docker-compose.yml, the init-db service overrides this with command: python init_bdt.py
CMD ["uvicorn", "api_mysql:app", "--host", "0.0.0.0", "--port", "8000"]
