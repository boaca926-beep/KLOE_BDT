#!/bin/bash
# stop_api.sh - Stop API server

if [ -f api.pid ]; then
    PID=$(cat api.pid)
    echo "Stopping API (PID: $PID)..."
    kill $PID
    rm api.pid
    echo "API stopped"
else
    echo "No PID file found, trying pkill..."
    pkill -f uvicorn
fi