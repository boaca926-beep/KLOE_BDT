#!/bin/bash
# make_pred.sh - Updated with 10 features

# Stop any running API
#pkill -f uvicorn

# Start API
#uvicorn api_mysql:app --reload --host 0.0.0.0 --port 8000 &
./start_api.sh

# Wait a moment
sleep 3

curl -X POST "http://localhost:8000/predict-and-save" \
  -H "Content-Type: application/json" \
  -d '{
    "run_number": 12345,
    "event_number": 67890,
    "photon_pairs": [
      {
        "m_gg": 135.2,
        "opening_angle": 0.85,
        "cos_theta": 0.5,
        "E_asym": 0.12,
        "e_min_x_angle": 42.1,
        "E1": 50.0,
        "E2": 48.0,
        "E3": 37.2,
        "asym_x_angle": 0.10,
        "E_diff": 15.3
      }
    ]
  }' | python3 -m json.tool

echo -e "\n✅ Prediction complete! Server is still running for more requests."
echo "To stop the server, run: ./stop_api.sh"