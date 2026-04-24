#!/bin/bash
# start_api.sh - Simple version

echo "=========================================="
echo "KLOE BDT API - Starting..."
echo "=========================================="

# Step 1: Initialize database
echo -e "\n📦 Initializing database..."
echo "DROP" | uv run init_bdt.py 2>/dev/null
echo "YES" | uv run init_bdt.py
echo "✅ Database ready"

# Step 2: Stop existing API if running
echo -e "\n🛑 Stopping existing API..."
pkill -f "uvicorn api_mysql:app" 2>/dev/null
rm -f api.pid
sleep 2

# Step 3: Start API
echo -e "\n🚀 Starting API server..."
nohup uv run uvicorn api_mysql:app --host 0.0.0.0 --port 8000 > api.log 2>&1 & # Runing API in background
echo $! > api.pid
echo "✅ API started with PID: $(cat api.pid)"

# Step 4: Wait for API to be ready
echo -e "\n⏳ Waiting for API to be ready..."
sleep 5

# Step 5: Verify API is running
if curl -s http://localhost:8000/health > /dev/null 2>&1; then
    echo -e "\n✅ API IS RUNNING!"
    echo -e "\n📊 API Status:"
    curl -s http://localhost:8000/health | python3 -m json.tool
else
    echo -e "\n❌ API failed to start!"
    echo "Check api.log for errors:"
    tail -20 api.log
    exit 1
fi

echo -e "\n=========================================="
echo "📍 API URL: http://localhost:8000"
echo "📚 Swagger: http://localhost:8000/docs"
echo "📝 Logs: tail -f api.log"
echo "=========================================="
