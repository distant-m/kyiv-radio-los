#!/bin/bash
# Double-click on macOS to install requirements (first time) and start the app.
cd "$(dirname "$0")"
if ! command -v python3 >/dev/null; then
  echo "Python 3 is not installed. Get it from https://www.python.org/downloads/ and run this again."
  read -n 1 -s -r -p "Press any key to close"; exit 1
fi
python3 -c "import rasterio, scipy, shapely, pyproj" 2>/dev/null || python3 -m pip install --user -r requirements.txt
python3 app.py
