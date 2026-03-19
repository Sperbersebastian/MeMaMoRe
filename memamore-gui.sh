#!/usr/bin/env bash
# memamore-gui.sh
# Launcher script for the MeMaMoRe Web Dashboard
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
PORT="${1:-5000}"
HOST="127.0.0.1"
URL="http://$HOST:$PORT"

echo "============================================================"
echo "                   🧬 MeMaMoRe GUI 🧬                     "
echo "============================================================"
echo "Starting local web server on $URL"
echo "Please wait a moment..."

# Start the Flask server in the background
"$ROOT/bin/main.sh" gui --port "$PORT" --host "$HOST" &
SERVER_PID=$!

# Wait a few seconds for the server to spin up
sleep 3

# Try to open the browser automatically
echo "Opening browser to $URL ..."
if command -v xdg-open >/dev/null 2>&1; then
  xdg-open "$URL" >/dev/null 2>&1 &
elif command -v open >/dev/null 2>&1; then
  open "$URL" >/dev/null 2>&1 &
elif command -v python3 >/dev/null 2>&1; then
  python3 -m webbrowser "$URL" >/dev/null 2>&1 &
else
  echo "Could not detect web browser launcher."
  echo "Please open your browser manually and go to: $URL"
fi

echo ""
echo "Press Ctrl+C to stop the GUI server and exit."
echo "============================================================"

# Wait for the server process (keeps the script running so Ctrl+C works)
wait $SERVER_PID
