@echo off
setlocal
cd /d "%~dp0"
where py >nul 2>nul
if %errorlevel%==0 (set PY=py -3) else (set PY=python)
%PY% scripts\thesis_validation.py
if errorlevel 1 (echo VALIDATION_FAILED& exit /b 1)
echo VALIDATION_PASS
