# ARTISTIC: ARTIST Integration Companion
**Software by Heonjoon Lee, Samuel W. Schaffter, and Rebecca Schulman**

Schulman Lab at the Johns Hopkins University

Cellular Engineering Group at the National Institute of Standards and Technology

## Overview
ARTISTIC is a Python-based software designed to methodize molecular signal transduction using the ARTIST (Aptamer-Regulated Transcription for In-vitro Sensing and Transduction) platform. The software provides tools for:

- dART Design: Generate DNA templates for aptamer-regulated transcription (dARTs) based on an input aptamer, reported Kd, and reported salt/buffer conditions
- Kd,apparent Fit: Normalize RFU vs time kinetic plots and extract dissociation constants from kinetic plot 
- Analog: Prescribe experimentation for analog response based on desired detection range and simulate dose-response curves
- Digital: Prescribe experimentation for digital response based on desired threshold and simulate dose-response curves

### Operating System:
- Windows: Native support with start.bat launcher
- macOS/Linux: Run via python application.py

### Python Version:
- Python 3.10 or higher (required)
- Ensure Python is added to your system PATH

### Dependencies:
The following Python packages are required (listed in requirements.txt):
- PySide6          GUI framework
- numpy            Numerical computing
- scipy            Scientific computing & optimization
- ViennaRNA        RNA secondary structure prediction
- matplotlib       Plotting
- draw_rna         RNA structure visualization (optional)
- pandas           Experimental data upload and handling
- openpyxl         Excel file handling

### Installation:
Step 1: Install Python
- Download Python 3.10+ from python.org/downloads
- Important: Check "Add Python to PATH" during installation

Step 2: Install Python Dependencies
- Navigate to the ARTISTIC directory and run:
- pip install -r requirements.txt

### Running the Program
Windows:
- Double-click start.bat in the ARTISTIC directory. The batch script will:
- Check for Python installation
- Verify ViennaRNA is available
- Install/update dependencies automatically
- Launch the application

macOS/Linux:
- Open a terminal in the ARTISTIC directory and run:
  
  cd ~/ARTIST-Integration-Companion
  
  chmod +x start.sh
  
  ./start.sh

### Features & Usage
ARTISTIC provides four main tabs for different workflows:

## 1. dART Design
### Purpose:
Generate DNA templates for aptamer-regulated transcription (dARTs) based on an input aptamer, reported Kd, and reported salt/buffer conditions.

### Inputs:
- Aptamer sequence: DNA or RNA sequence (only A, T, C, G, U accepted)
- Reported Kd: Reported dissociation constant in nanomolar (nM)
- Salt/buffer condition: Reported Salt/buffer composition (e.g., 150 mM NaCl, 5 mM MgCl2, PBS)
- Output domain: Select O1, O2, or O3 domains
- Tests multiple insulation sequences for optimal RNA folding, evaluates secondary structure stability using ViennaRNA, and ranks designs by minimum free energy (MFE).

### Outputs:
- Best dART design with sequences
- RNA secondary structure visualization
- Strands to order (template strand, non-template strands (nt; e.g., Prom-nt, Output-nt), reporter strands corresponding to output (O1, O2, and O3).
- o1_f, o2_f, and o3_f ordered standard desalted in the study. o1_q, o2_q, and o3_q must be ordered HPLC purified. 

- User is prescribed sequences to order by the software (e.g., via Integrated DNA Technologies) as standard desalting, PAGE, or HPLC purified oligo DNA or RNA strands (if RNA aptamer is used).
  
 Upon arrival, strands are recommended for dilution at 100 µM using MilliQ water.
 
 To anneal strands into dARTs, each strand is added into a 200 µL PCR tube at concentrations of 1 µL per strand in 1X NEB RNAPol reaction buffer.

 Mixture is heated to 90 °C, incubated for 5 min, then cooled to 20 °C at a rate of 1 °C min−1.

 dARTs whose aptamers are predicted to form a G-quadruplex* are supplemented with 100 mM of KCl during annealing.
 
 *Not all papers mention whether the aptamer forms a G-quadruplex. If unknown, G-quadruplex predictions are available under: https://bioinformatics.ramapo.edu/QGRS/index.php (Kikin and Bagga et al., *Nuc. Acid Res.* **34**, 676-682 (2006).)

 Detailed methods are available in the main text of the paper.

 - Titration grid of salt concentration vs ligand concentration. ARTISTIC picks the highest salt concentration in the input buffer/salt conditions to generate salt titration grid.
  - Column: Monovalent salts (e.g., Na+, K+) span across 0 to 100 mM, and divalent salts (e.g., Mg2+, Ca2+) span across 0 to 5 mM.
    - RNA aptamers should be supplemented with 2 mM of manganese!
    - PBS (Phosphate-Buffered Saline) buffer contains ~ 137 mM of Na+. ARTISTIC will identify PBS as Na+ titration.
  - Columns 2+: RFU (Relative Fluorescence Units) for each ligand concentration
  - Rows: Ligand concentration based on Kd,reported multipliers (0, 1×, 2×, 2.5×, 5×, 10×, 25x, 50x, 100x Kd)

## 2. Kd,apparent Fit
### Purpose:
Normalize RFU vs time kinetic plots and extract dissociation constants from kinetic plot.

### Inputs:
- Excel sheet (.xlsx) containing time-lapse measurements of fluorescence across varying ligand concentrations
  - Column 1: Time points (minutes)
  - Columns 2+: RFU (Relative Fluorescence Units) for each ligand concentration
- Rows: Ligand concentration based on Kd multipliers (0, 1×, 2×, 2.5×, 5×, 10×, etc.)
- First row needs to be 0× Kd, no ligand added.
- At the end of the experiments, 0.5 µL of a DNA strand fully complementary strand at a final concentration of 2.5 µM needs to be added to obtain a maximum of reporter fluorescence intensity.
- Program normalizes fluoresence based on maximum fluorescence intensity.

### Parameters:
- dART concentration (nM): dART template concentration
- Reporter conc. (nM): Total reporter concentration
- Txn rate (nM/min): Slope of reacted reporter kinetics at specific time point (default: 20 min)
- Kd multipliers: Comma-separated list matching Excel columns (e.g., 0, 1, 2, 2.5, 5, 10, 25, 50, 100)

### Outputs:
- Normalized RFU data (RFU converted to [Reacted Reporter] (nM))
- Fits basal transcription rate (txn rate) from 0 nM ligand condition.
- Computes experimental slopes at specified time point.
- Optimizes Kd to minimize residuals between experimental and simulated slopes.
- Fitted Kd: Kd based off Nelder-Mead optimization ± error
- Fitted x: Basal ranscription rate
- Slope Comparison Table: Experimental slopes for each concentration
- Plots normalized reacted reporter kinetics (Reacted Reporter kinetics) and dose-response curves
- Detailed methods are available in the main text of the paper.

## 3. Analog
### Purpose:
Prescribe experimentation for analog response based on desired detection range and simulate dose-response curves.

### Inputs:
- Kd (nM): Kd,reported for prediction; Kd,apparent if fitting is done.
- Transcription rate (k_txn; (nM/min)): Transcription rate (range: 0.0001 – 1.0)
- dART (nM): dART template concentration (range: 10 – 200 nM)

### Modes:
- Single Curve: Simulate one dose-response curve at specified dART concentration (0 to 100 nM slider).
- Sweep: Generate family of curves across multiple dART concentrations.

### Outputs:
- Metrics Table:
  - L10: Ligand concentration at 10% max response
  - L90: Ligand concentration at 90% max response
- Inflection point (EC50)
- Sigmoid slope (Sensitivity)
- Maximum signal (maxY)
- Ligand Concentration Series: Recommended experimental concentrations.
- Dose-Response Plot: Reacted Reporter vs. Ligand concentration.

## 4. Digital
### Purpose:
Prescribe experimentation for digital response based on desired threshold and simulate dose-response curves

### Inputs:
- Kd,apparent (nM): Kd,reported for prediction; Kd,apparent if fitting is done.
- k_txn Ref template (nM/min): Reference template transcription rate (assuming that ligand does not affect txn rate).
- k_txn Inverter dART (nM/min): Inverter dART basal transcription rate.
- Inverter dART (nM): Fixed at 50 nM. User may toggle if desired.

### Modes:
- Single Curve: Simulate at one Reference template concentration (10 to 100 nM slider).
- Sweep: Generate curves across multiple Reference template values.

### Outputs:
- Threshold concentration: Ligand level triggering switch
- Maximum signal (maxY)
- Ligand titration series: 7-point series around threshold (0, 3x below threshold, 3x above threshold)
- Threshold plot: Reporter signal vs. ligand (log scale)

## Troubleshooting
### Common Issues
- "Python is not installed or not on PATH"

Reinstall Python and check "Add Python to PATH". Verify with: python --version

- "ViennaRNA not found"

Install ViennaRNA with Python bindings. Test: python -c "import RNA"

- "draw_rna not available"

RNA structures won't display visually. Install: pip install draw_rna. The program will still function without it.

- Excel file not loading

Ensure .xlsx format (not .xls). First column must contain numeric time values. Check that Kd multipliers match number of data columns.

## Citation
Original ARTIST paper:

Lee, Schaffter, and Schulman et al. "Plug-and-play protein biosensors using aptamer-regulated in vitro transcription." *Nat. Commun.* **15**, 7973 (2024). https://doi.org/10.1038/s41467-024-51907-4

[ARTISTIC Reference TBD]

## Contact & Support
For questions, bug reports, or feature requests:
- Schulman Lab @ Johns Hopkins University : Heonjoon Lee (hlee260@jhmi.edu) and Rebecca Schulman (rschulm3@jhu.edu)
- Cellular Engineering Group @ NIST: Samuel W. Schaffter (samuel.schaffter@nist.gov)
- GitHub: https://github.com/hlee260/ARTIST-Integration-Companion
