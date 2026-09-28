@echo off
setlocal
where python >nul 2>nul
if errorlevel 1 (
  echo ERROR: python not found on PATH
  exit /b 1
)
python scripts\smoke_test.py
if errorlevel 1 exit /b 1
pause
