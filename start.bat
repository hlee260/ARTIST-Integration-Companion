@echo off
echo =======================================================================================
echo                    ARTISTIC: ARTIST Integration Companion
echo    Software by Heonjoon Lee, Samuel W. Schaffter, and Rebecca Schulman
echo                 Schulman Lab at the Johns Hopkins University
echo    Cellular Engineering Group at the National Institute of Standards and Technology
echo =======================================================================================
echo.

python --version >nul 2>&1
if errorlevel 1 (
    echo Python is not installed or not on PATH.m
    echo Please install Python 3.10+ from https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during install.
    pause
    exit /b 1
)

python -c "import RNA" >nul 2>&1
if errorlevel 1 (
    echo ViennaRNA not found. Please install it with Python bindings:
    echo   https://www.tbi.univie.ac.at/RNA/#download
    pause
    exit /b 1
)

echo Installing/updating dependencies...
pip install -r requirements.txt --quiet

echo Starting ARTISTIC...
python application.py
