@echo off
REM Double-click on Windows to install requirements (first time) and start the app.
cd /d "%~dp0"
where py >nul 2>nul && (set PY=py -3) || (set PY=python)
%PY% -c "import rasterio, scipy, shapely, pyproj" 2>nul || %PY% -m pip install --user -r requirements.txt
%PY% app.py
pause
