ARTISTIC: ARTIST Integration Companion
Software by Heonjoon Lee, Samuel W. Schaffter, and Rebecca Schulman

Schulman Lab at the Johns Hopkins University

Cellular Engineering Group at the National Institute of Standards and Technology

Overview
ARTISTIC is a Python-based graphical application for designing and simulating RNA-based molecular sensors using the ARTIST (Aptamer-Regulated dART-Initiated Synthesis Translation) platform. The software provides tools for:

dART Design: Generate optimized DNA templates for aptamer-regulated txn based on input aptamer 
Kd,apparent Fitting: Normalize RFU vs time kinetic plots and extract dissociation constants from kinetic plot 
Analog:Simulate dose-response curve and prescribe experimenta for analog sensor
Digital: Simulate dose-response curve and prescribe experiment for digital sensor

Operating System
Windows: Native support with start.bat launcher
macOS/Linux: Run via python application.py

Python Version
Python 3.10 or higher (required)
Ensure Python is added to your system PATH

Dependencies
The following Python packages are required (listed in requirements.txt):
PySide6          # GUI framework
numpy            # Numerical computing
scipy            # Scientific computing & optimization
ViennaRNA        # RNA secondary structure prediction
matplotlib       # Plotting
draw_rna         # RNA structure visualization (optional)
pandas           # Data manipulation
openpyxl         # Excel file handling

Installation
Step 1: Install Python
Download Python 3.10+ from python.org/downloads
Important: Check "Add Python to PATH" during installation

Step 3: Install Python Dependencies
Navigate to the ARTISTIC directory and run:
pip install -r requirements.txt

Running the Program

Windows:
Double-click start.bat in the ARTISTIC directory. The batch script will:
Check for Python installation
Verify ViennaRNA is available
Install/update dependencies automatically
Launch the application

macOS/Linux:
Open a terminal in the ARTISTIC directory and run:
bash
python application.py

Features & Usage
ARTISTIC provides four main tabs for different workflows:

1. dART Design Tab
Purpose: Generate optimized dART templates from aptamer sequences.

Inputs:
Aptamer sequence: DNA or RNA sequence (only A, T, C, G, U accepted)
Reported Kd: Dissociation constant in nanomolar (nM)
Salt/buffer condition: Reported Salt/buffer composition (e.g., 150 mM NaCl, 5 mM MgCl2, PBS)
Output domain: Select O1, O2, or O3 domains

Program tests multiple insulation sequences for optimal RNA folding, evaluates secondary structure stability using ViennaRNA, and ranks designs by minimum free energy (MFE). Program then outputs to the user best dART design with sequences, RNA secondary structure visualization, and strands to order (template strand, non-template strands (nt; e.g., Prom-nt, Output-nt).

User is required to order the sequences as prescribed by the software (e.g., via Integrated DNA Technologies) as standard desalting, PAGE, or HPLC purified oligo DNA or RNA strands (if RNA aptamer is used).
Upon arrival, strands are recommended for dilution at 100 µM using MilliQ water.
To anneal strands into dARTs, each strand is added into a 200 µL PCR tube at concentrations of 1 µL per strand in 1X NEB RNAPol reaction buffer. 
Mixture is heated to 90 oC, incubated for 5 min, then cooled to 20 oC at a rate of 1 oC min−1.
dARTs whose aptamers are predicted to form a G-quadruplex* are supplemented with 100 mM of KCl during annealing.

*G-quadruplex predictions are available under: https://bioinformatics.ramapo.edu/QGRS/index.php (Kikin and Bagga et al., Nuc. Acid Res. 34, 676-682 (2006).)

Detailed methods are available in the main text of the paper.

2. Kd Fitting Tab
Purpose: Normalizes raw time-lapse fluorescence data to reacted reporter kinetics and extract Kd,apparent from fit.

Required Excel Format (from application.py):

Column 1: Time points (minutes)
Columns 2+: RFU (Relative Fluorescence Units) for each ligand concentration
Rows: Must match Kd multipliers (0, 1×, 2×, 2.5×, 5×, 10×, etc.)
First row: 0× Kd, no ligand added

Parameters:
dART conc. (nM): dART template concentration
Reporter conc. (nM): Total reporter concentration
Slope time point (min): Time at which to compute slope (default: 20 min)
Kd multipliers: Comma-separated list matching Excel columns (e.g., 0, 1, 2, 2.5, 5, 10, 25, 50, 100)
Process:

Normalizes RFU data
Fits basal transcription rate from 0 nM ligand condition
Computes experimental slopes at specified time point
Optimizes Kd to minimize residuals between experimental and simulated slopes
Outputs:

Fitted Kd: Best-fit dissociation constant ± error
Fitted k_txn: Transcription rate
Slope Comparison Table: Experimental vs. simulated slopes for each concentration
Plots:
Normalized RFU time courses (overlaid)
Experimental vs. simulated slope comparison

Detailed methods are available in the main text of the paper.

3. Analog Sensor Tab
Purpose: Simulate dose-response curves for analog sensors.

Parameters:
Kd (nM): Dissociation constant (> 0 nM)
k_txn: Transcription rate (range: 0.0001 – 1.0)
dART (nM): dART template concentration (range: 10 – 200 nM)

Modes:
Single Curve: Simulate one dose-response curve at specified dART concentration
Sweep: Generate family of curves across multiple dART concentrations
Outputs:

Metrics Table:
L10: Ligand concentration at 10% max response
L90: Ligand concentration at 90% max response
Inflection point (EC50)
Sigmoid slope (Sensitivity)
Maximum signal (maxY)
Ligand Concentration Series: Recommended experimental concentrations
Dose-Response Plot: Reporter signal vs. ligand concentration

4. Digital Sensor Tab
Purpose: Model threshold-based digital (bistable) RNA sensors.

Parameters:
Kd (nM): Aptamer-ligand dissociation constant
k_txn Ref template: Reference template transcription rate
k_txn Inverter dART: Inverter dART transcription rate
Inverter dART (nM): Fixed at 50 nM. User may toggle if desired.

Modes:
Single Curve: Simulate at one Reference template concentration (10 – 100 nM slider)
Sweep: Generate curves across multiple Ref dART values
Outputs:

Threshold concentration: Ligand level triggering switch
Maximum signal (maxY)
Ligand titration series: 7-point series around threshold (0, 3x below threshold, 3x above threshold)
Threshold plot: Reporter signal vs. ligand (log scale)

ARTISTIC/
│
├── application.py       # Main GUI application
├── ARTISTIC.py          # Core simulation & design functions
├── requirements.txt     # Python dependencies
├── start.bat            # Windows launcher script
└── README.md            # This file

Troubleshooting
Common Issues
"Python is not installed or not on PATH"

Reinstall Python and check "Add Python to PATH"
Verify with: python --version
"ViennaRNA not found"

Install ViennaRNA with Python bindings
Test: python -c "import RNA"
"draw_rna not available"

RNA structures won't display visually
Install: pip install draw_rna
The program will still function without it
Excel file not loading

Ensure .xlsx format (not .xls)
First column must contain numeric time values
Check that Kd multipliers match number of data columns

Citation
If you use ARTISTIC in your research, please cite:

[Citation information to be provided by authors]

Contact & Support
For questions, bug reports, or feature requests:

Schulman Lab @ Johns Hopkins University : Heonjoon Lee (hlee260@jhmi.edu) and Rebecca Schulman (rschulm3@jhu.edu)
Cellular Engineering Group @ NIST: Samuel W. Schaffter ()
GitHub: https://github.com/hlee260/ARTISTIC
License
[License information to be provided by authors]
