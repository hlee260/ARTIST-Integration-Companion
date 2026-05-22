#!/bin/bash
echo "================================================================"
echo  "ARTISTIC: ARTIST Integration Companion"
echo  "Software by Heonjoon Lee, Byunghwa Kang, Samuel W. Schaffter"
echo  "Schulman Lab at the Johns Hopkins University"
echo  "Cellular Engineering Group at National Institute of Standards and Technology"
echo "================================================================  "
echo

if ! command -v python3 &>/dev/null; then
    echo "Python 3 is not installed."
    echo "Install it from https://www.python.org/downloads/ or via your package manager."
    exit 1
fi

if ! python3 -c "import RNA" &>/dev/null; then
    echo "ViennaRNA not found. Please install it with Python bindings:"
    echo "  https://www.tbi.univie.ac.at/RNA/#download"
    exit 1
fi

echo "Installing/updating dependencies..."
pip3 install -r requirements.txt --quiet

echo "Starting ARTISTIC..."
python3 application.py
