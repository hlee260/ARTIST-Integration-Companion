import numpy as np
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import warnings
import RNA

# ── Constants ────────────────────────────────────────────────────────────────
Prom_nt = "TTCTAATACGACTCACTATA"
Prom_t  = "TATAGTGAGTCGTATTAGAA"

O1_nt   = "CTACATCCACATACTAATTAAC"
O1_t    = "GTTAATTAGTATGTGGATGTAG"
O2_nt   = "CTACTTTCACTTCACAACATCA"
O2_t    = "TGATGTTGTGAAGTGAAAGTAG"
O3_nt   = "TACCATCACATTCAATAATCCT"
O3_t    = "AGGATTATTGAATGTGATGGTA"

ALL_INS = [
    "GGGATG","GGGAGT","GGGAGA","GGGAAA","GGGAAC","GGGAAG","GGGAAT","GGGACA","GGGACG","GGGACT",
    "GGGAGC","GGGAGG","GGGATA","GGGATC","GGGATT","GGGCAA","GGGCAC","GGGCAG","GGGCAT","GGGCCA",
    "GGGCCC","GGGCCG","GGGCCT","GGGCGA","GGGCGC","GGGCGG","GGGCGT","GGGCTA","GGGCTC","GGGCTG",
    "GGGCTT","GGGTAA","GGGTAC","GGGTAG","GGGTAT","GGGTCA","GGGTCC","GGGTCG","GGGTCT","GGGTGA",
    "GGGTGC","GGGTGG","GGGTGT","GGGTTA"
]

RC = {"A":"T","T":"A","C":"G","G":"C","U":"A"}

def rev_comp(seq: str) -> str:
    return "".join(RC.get(b, b) for b in reversed(seq))

ALL_COMP = [rev_comp(ins) for ins in ALL_INS]
# One exception: GGGACT paired with TGACCC (not rev_comp AGTCCC) — keep original
ALL_COMP[9] = "TGACCC"


# ── Utilities ────────────────────────────────────────────────────────────────
def sanitize_sequence(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in "ATCGU")

def get_rna_transcript(dna_seq: str) -> str:
    return dna_seq.upper().replace("T", "U")

def build_dart(apt: str, ins: str, ins_comp: str, gate: str):
    out_t = O2_t if gate == "O2" else O3_t if gate == "O3" else O1_t
    encoded = out_t + ins + apt + ins_comp
    template = encoded + Prom_t
    prom_nt = Prom_nt + ins
    output_nt = ins_comp + rev_comp(out_t)
    return {
        "template": template,
        "nonTemplate": prom_nt,
        "rna": get_rna_transcript(rev_comp(encoded)),
        "ins": ins,
        "insComp": ins_comp,
        "encoded": encoded,
        "out_t": out_t,
        "Output_nt": output_nt
    }

DUMMY_SEQ   = "TTCTTCTTCTTCTTCTTCTTCTTC"
DUMMY_INS   = "GGGATG"
DUMMY_COMP  = "CATCCC"
DUMMY_PROM_NT = Prom_nt + DUMMY_INS   # "TTCTAATACGACTCACTATAGGGATG"

def build_reference_dart(gate: str) -> dict:
    """
    Reference (dummy) dART — uses a fixed dummy payload and GGGATG insulation.
    Same for all gate choices (O1/O2/O3) except out_t changes.
    """
    out_t = O2_t if gate == "O2" else O3_t if gate == "O3" else O1_t
    dum_template    = out_t + DUMMY_INS + DUMMY_SEQ + DUMMY_COMP + Prom_t
    dummy_output_nt = DUMMY_COMP + rev_comp(out_t)
    return {
        "gate":             gate,
        "Dummy_Template":   dum_template,
        "Dummy_Prom-nt":    DUMMY_PROM_NT,
        "Dummy_Output-nt":  dummy_output_nt,
    }

def build_inverter_dart(gate: str, ins: str, ins_comp: str) -> dict:
    """
    Inverter dART — encodes the dummy payload but triggered by the aptamer dART output.
    Invert-nt = out_t[:-4]  (output trigger minus last 4 bases).
    Template  = rev_comp(Invert-nt) + ins + dummy + ins_comp + Prom_t
    """
    out_t      = O2_t if gate == "O2" else O3_t if gate == "O3" else O1_t
    invert_nt  = out_t[:-4]
    inv_template   = rev_comp(invert_nt) + ins + DUMMY_SEQ + ins_comp + Prom_t
    inv_prom_nt    = Prom_nt + ins
    inv_output_nt  = ins_comp + invert_nt
    return {
        "gate":              gate,
        "Invert-nt":         invert_nt,
        "Inv_Template":      inv_template,
        "Inv_Prom-nt":       inv_prom_nt,
        "Inv_Output-nt":     inv_output_nt,
    }

def get_digital_strands(gate: str, ins: str, ins_comp: str) -> list:
    """
    Return (label, sequence, notes) rows for the Digital tab order table.
    Includes reference dART strands and inverter dART strands.
    """
    ref  = build_reference_dart(gate)
    inv  = build_inverter_dart(gate, ins, ins_comp)
    return [
        ("Ref dART — Template",
        ref["Dummy_Template"],
        "Reference template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Ref dART — Prom-nt",
        ref["Dummy_Prom-nt"],
        "Reference template promoter non-template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Ref dART — Output-nt",
        ref["Dummy_Output-nt"],
        "Reference template output non-template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Inverter dART — Template",
        inv["Inv_Template"],
        "Inverter dART template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Inverter dART — Prom-nt",
        inv["Inv_Prom-nt"],
        "Inverter dART promoter non-template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Inverter dART — Output-nt",
        inv["Inv_Output-nt"],
        "Inverter dART output non-template strand — order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
    ]

def calculate_unbound_dART(L0, dART0, K_D):
    m = dART0 + K_D + L0
    disc = m * m - 4 * dART0 * L0
    if disc < 0:
        return 0.0
    dART_bound = (m - np.sqrt(disc)) / 2
    return float(dART0 - dART_bound)

# ── Candidate search ─────────────────────────────────────────────────────────
def test_insulations(aptamer_seq: str, gate="O1", bp_threshold=0.5):
    apt = sanitize_sequence(aptamer_seq)
    best = None
    RNA.cvar.temperature = 37

    for ins, ins_comp in zip(ALL_INS, ALL_COMP):
        cand = build_dart(apt, ins, ins_comp, gate)
        rna = cand["rna"]

        struct, mfe = RNA.fold(rna)
        fc = RNA.fold_compound(rna)
        fc.pf()
        bpp = fc.bpp()

        start = range(7)
        target = range(6 + len(apt), 6 + len(apt) + 8)
        count = sum(1 for i in start for j in target if bpp[i][j] > bp_threshold)

        candidate = {
            **cand,
            "structure": struct,
            "mfe": mfe,
            "bpp": bpp,
            "count": count,
        }

        if best is None or mfe < best["mfe"]:
            best = candidate

    return best

# ── Strands to order ─────────────────────────────────────────────────────────
def get_order_strands(best: dict, gate: str) -> list:
    """Return list of (label, sequence, notes) for ordering."""
    out_nt = (
        O2_nt if gate == "O2" else
        O3_nt if gate == "O3" else
        O1_nt
    )
    strands = [
        ("Template strand",
         best["template"],
         "Order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Prom-nt strand",
         best["nonTemplate"],
         "Order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("Output-nt strand",
         best["Output_nt"],
         "Order as DNA oligo, standard desalted, PAGE, or HPLC purified"),
        ("RNA transcript (predicted)",
         best["rna"],
         "Transcribed by dART"),
    ]
    return strands

# ── Salt / titration helpers ──────────────────────────────────────────────────
def parse_salt(salt_str: str):
    """
    Parse a salt string like '150 mM NaCl' or 'NaCl' or '5 mM MgCl2'.
    Also recognises PBS / phosphate-buffered saline → NaCl 150 mM.
    Returns (salt_name, conc_mM_or_None).
    """
    import re
    s = salt_str.strip()

    # PBS / phosphate-buffered saline → treat as 150 mM NaCl
    if re.search(r'\bPBS\b', s, re.I) or re.search(r'phosphate[- ]buffered\s+saline', s, re.I):
        return "NaCl", 150.0

    m = re.match(r"(\d+(?:\.\d+)?)\s*mM\s+(\S+)", s, re.I)
    if m:
        return m.group(2), float(m.group(1))
    m = re.match(r"(\S+)\s+(\d+(?:\.\d+)?)\s*mM", s, re.I)
    if m:
        return m.group(1), float(m.group(2))
    # bare name
    return s if s else "NaCl", None

def salt_rows(salt_name: str = None):
    """
    Fixed titration series — always capped regardless of reported condition.
    Monovalent: 0, 5, 10, 25, 50, 100 mM
    Divalent:   0, 1, 2, 3, 4, 5 mM
    """
    if not salt_name:
        salt_name = "NaCl"
    divalent = salt_name.strip() in ("MgCl2", "CaCl2", "MnCl2", "ZnCl2")
    series = [0, 1, 2, 3, 4, 5] if divalent else [0, 5, 10, 25, 50, 100]
    return [{"salt": salt_name.strip(), "conc": c} for c in series]

def lig_series(kd):
    try:
        k = float(kd)
    except Exception:
        return []
    if not k or np.isnan(k):
        return []
    multipliers = [0, 1, 2, 2.5, 5, 10, 25, 50, 100]
    return [
        {
            "label": "0 (blank)" if m == 0 else f"{m}× Kd",
            "conc": round(m * k, 3),
        }
        for m in multipliers
    ]

def digital_titration_series(threshold):
    if threshold is None or threshold <= 0:
        return []
    below = np.array([threshold / 3.0, threshold / 2.0, threshold * 0.85])
    above = np.array([threshold * 1.15, threshold * 2.0, threshold * 3.0])
    series = [0.0] + below.tolist() + above.tolist()
    return [
        {"label": "0 (blank)" if i == 0 else f"Pt {i}", "conc": float(v)}
        for i, v in enumerate(series)
    ]

# ── ODE models ───────────────────────────────────────────────────────────────
def analog_sensor(y, t, k_txn, k_sd, k_deg, k_hyb, dART_unbound):
    RNAc, FQ, RNAQ, Q, F = y
    return [
        k_txn * dART_unbound - k_sd * FQ * RNAc - k_hyb * RNAc * Q,
        k_hyb * F * Q - k_sd * FQ * RNAc,
        k_sd * FQ * RNAc - k_deg * RNAQ + k_hyb * RNAc * Q,
        k_deg * RNAQ - k_hyb * F * Q - k_hyb * RNAc * Q,
        k_sd * RNAc * FQ - k_hyb * F * Q
    ]

def interpolate(xData, yData, targetY):
    for i in range(len(yData) - 1):
        y1, y2 = yData[i], yData[i + 1]
        if (y1 <= targetY <= y2) or (y2 <= targetY <= y1):
            if abs(y2 - y1) < 1e-10:
                return float(xData[i])
            return float(xData[i] + (targetY - y1) / (y2 - y1) * (xData[i + 1] - xData[i]))
    return None

def fit_sigmoid(xData, yData):
    A = float(np.min(yData))
    B = float(np.max(yData))

    def sig(x, C, D):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return A + (B - A) / (1 + (x / C) ** D)

    best = {"C": None, "D": None, "A": A, "B": B}
    try:
        popt, _ = curve_fit(
            sig, xData, yData,
            p0=[500, 1.0],
            bounds=([1e-9, 1e-9], [np.inf, np.inf]),
            maxfev=20000,
        )
        best["C"] = float(popt[0])
        best["D"] = float(popt[1])
    except Exception:
        pass
    return best

# ── Analog simulation ────────────────────────────────────────────────────────
def simulate_analog_curve(dart0, kd, k_txn_user):
    # k_txn_user is per-second rate; multiply by 60 to get per-minute (matches HTML)
    k_txn = k_txn_user * 60
    k_hyb = 1e6 / 1e9 * 60
    k_deg = 0.001 * 60
    k_sd  = 1e4 / 1e9 * 60
    FQ0   = 250
    # HTML: length=120, t = i*120/999  → t in [0, 119*120/999] ≈ [0, 14.29]
    t = np.array([i * 120 / 999 for i in range(120)])
    L0_list = np.arange(0, 10000)

    y_final = []
    for L0 in L0_list:
        u = calculate_unbound_dART(L0, dart0, kd)
        sol = odeint(analog_sensor, [0, FQ0, 0, 0, 0], t, args=(k_txn, k_sd, k_deg, k_hyb, u))
        y_final.append(sol[-1, 4])

    y_final = np.array(y_final)
    maxY = float(np.max(y_final))
    y10 = maxY * 0.1
    y90 = maxY * 0.9

    L10 = interpolate(L0_list, y_final, y10)
    L90 = interpolate(L0_list, y_final, y90)
    sigmoid = fit_sigmoid(L0_list, y_final)

    metrics = {
        "L10": L10,
        "L90": L90,
        "inflection": sigmoid["C"],
        "slope": sigmoid["D"],
        "maxY": maxY,
    }

    return L0_list, y_final, metrics

# ── Digital simulation ───────────────────────────────────────────────────────
def digital_sensor(y, t, k_txn_ref, k_txn_apt, k_th, k_sd, k_deg,
                   Ref_dART, Apt_dART_unbound, Reporter0):
    """HTML-matched ODE — includes k_deg degradation on both RNA species."""
    Ref_RNA, Apt_RNA, Rep = y
    dRef = (k_txn_ref * Ref_dART
            - k_th  * Ref_RNA * Apt_RNA
            - k_sd  * Ref_RNA * (Reporter0 - Rep)
            - k_deg * Ref_RNA)
    dApt = (k_txn_apt * Apt_dART_unbound
            - k_th  * Ref_RNA * Apt_RNA
            - k_deg * Apt_RNA)
    dRep = k_sd * Ref_RNA * (Reporter0 - Rep)
    return [dRef, dApt, dRep]

def simulate_digital_curve(ref_dart, apt_dart0, kd, k_txn_ref_user, k_txn_apt_user):
    # multiply by 60 to convert per-second → per-minute (matches HTML)
    k_txn_ref = k_txn_ref_user * 60
    k_txn_apt = k_txn_apt_user * 60
    k_th      = 1e6 / 1e9 * 60
    k_sd      = 1e4 / 1e9 * 60
    k_deg     = 0.001 * 60
    Reporter0 = 250
    # HTML: length=240, t = i*240/2399  → t in [0, 239*240/2399] ≈ [0, 23.91]
    t = np.array([i * 240 / 2399 for i in range(240)])
    L0_list = np.logspace(0, 4, 101)

    y_final = []
    for L0 in L0_list:
        apt_unbound = calculate_unbound_dART(L0, apt_dart0, kd)
        sol = odeint(
            digital_sensor,
            [0, 0, 0],
            t,
            args=(k_txn_ref, k_txn_apt, k_th, k_sd, k_deg,
                  ref_dart, apt_unbound, Reporter0)
        )
        y_final.append(max(0.0, sol[-1, 2]))

    y_final = np.array(y_final)

    threshold = None
    maxDeriv = 0
    for i in range(1, len(y_final) - 1):
        d = abs((y_final[i + 1] - y_final[i - 1]) / (np.log10(L0_list[i + 1]) - np.log10(L0_list[i - 1])))
        if d > maxDeriv:
            maxDeriv = d
            threshold = float(L0_list[i])

    return L0_list, y_final, {"threshold": threshold, "maxY": float(np.max(y_final))}

# ── Kd Fitting from experimental data ────────────────────────────────────────
from scipy.optimize import minimize

KD_FIT_MULTIPLIERS = [0, 1, 2, 2.5, 5, 10, 25, 50, 100]  # ×Kd rows expected

def ARTIST_rxn(y, t, k_txn, k_sd, dART, Reporter0):
    RNA, Rep_reacted = y
    dRNA = k_txn * dART - k_sd * RNA * Reporter0 + k_sd * RNA * Rep_reacted
    dRep = k_sd * RNA * Reporter0 - k_sd * RNA * Rep_reacted
    return [dRNA, dRep]

def _slope_at(time_arr, data_arr, idx):
    """Central finite difference slope at index idx."""
    if idx <= 0 or idx >= len(time_arr) - 1:
        return None
    dt = time_arr[idx + 1] - time_arr[idx - 1]
    if dt == 0:
        return None
    return (data_arr[idx + 1] - data_arr[idx - 1]) / dt

def normalize_rfu(rfu_matrix, reporter0):
    """
    Normalise RFU to nM for display purposes only (per-column normalisation).
    Uses the formula from the notebook:
        reporter0 * (RFU - col_min) / (col_max - col_min)
    This is used for plotting traces; NOT for slope-based Kd fitting.
    """
    col_min = rfu_matrix.min(axis=0)
    col_max = rfu_matrix.max(axis=0)
    rng = col_max - col_min
    rng[rng == 0] = 1.0
    return reporter0 * (rfu_matrix - col_min) / rng

def fit_basal_ktxn(time_arr, raw_blank, dart0, reporter0):
    """
    Fit k_txn from the raw (un-normalised, offset-to-zero) blank trace.
    Converts RFU to nM by scaling so max(raw_blank) maps to the simulated
    steady-state reporter signal at t=max(time_arr) with k_txn initial guess.
    Returns fitted k_txn (per minute).
    """
    k_sd = 1e4 / 1e9 * 60

    # RFU → relative units: just normalise to [0,1] range; the ODE shape
    # determines k_txn independently of absolute RFU amplitude
    rng = np.max(raw_blank) - np.min(raw_blank)
    if rng <= 0:
        return 0.001 * 60
    unit_trace = (raw_blank - np.min(raw_blank)) / rng  # 0..1

    def residuals(params):
        k_txn = params[0]
        if k_txn <= 0:
            return 1e12
        sol = odeint(ARTIST_rxn, [0, 0], time_arr,
                     args=(k_txn, k_sd, dart0, reporter0))
        sim = sol[:, 1]
        sim_rng = np.max(sim) - np.min(sim)
        if sim_rng <= 0:
            return 1e12
        unit_sim = (sim - np.min(sim)) / sim_rng  # 0..1
        return float(np.sum((unit_trace - unit_sim) ** 2))

    best_cost, best_ktxn = np.inf, 0.001 * 60
    for x0 in np.logspace(-4, 1, 30):
        r = minimize(residuals, [x0], method='Nelder-Mead',
                     options={'xatol': 1e-10, 'fatol': 1e-10, 'maxiter': 10000})
        if r.fun < best_cost:
            best_cost = r.fun
            best_ktxn = float(r.x[0])
    return best_ktxn

def fit_kd_from_slopes(time_arr, rfu_raw, dart0, reporter0, k_txn,
                       kd_multipliers, slope_time_min=20):
    """
    Fit Kd by matching slopes at `slope_time_min` minutes.
    Uses raw RFU values for experimental slopes (relative shape matters, not units).
    Simulated slopes are computed in nM using the ARTIST_rxn ODE.

    Parameters
    ----------
    time_arr       : 1-D array of time points (minutes), starting at 0
    rfu_raw        : 2-D array (n_timepoints × n_ligand_concs), raw RFU
    dart0          : dART concentration (nM)
    reporter0      : total reporter concentration (nM)
    k_txn          : fitted transcription rate (per minute)
    kd_multipliers : list of Kd multipliers matching columns of rfu_raw
    slope_time_min : time point at which to compute slope (default 20 min)
    """
    k_sd = 1e4 / 1e9 * 60
    t_sim = np.linspace(0, max(float(time_arr[-1]), slope_time_min + 5), 1000)
    slope_idx_exp = int(np.abs(time_arr - slope_time_min).argmin())
    slope_idx_sim = int(np.abs(t_sim - slope_time_min).argmin())

    # Experimental slopes from raw RFU
    exp_slopes_raw = np.array([
        _slope_at(time_arr, rfu_raw[:, col], slope_idx_exp) or 0.0
        for col in range(rfu_raw.shape[1])
    ])

    # Scale experimental slopes to nM/min using the blank (0× Kd) column
    # so that units match simulated slopes
    def _sim_slopes(K_D):
        sim = []
        for mult in kd_multipliers:
            L0 = mult * K_D
            u = calculate_unbound_dART(L0, dart0, K_D) if mult > 0 else float(dart0)
            sol = odeint(ARTIST_rxn, [0, 0], t_sim,
                         args=(k_txn, k_sd, u, reporter0))
            s = _slope_at(t_sim, sol[:, 1], slope_idx_sim)
            sim.append(s if s is not None else 0.0)
        return np.array(sim)

    # Find scale factor: ratio of blank sim slope to blank exp slope
    sim_blank = _sim_slopes(dart0 * 10)[0]   # rough K_D estimate for scale
    exp_blank = exp_slopes_raw[0]
    scale = (sim_blank / exp_blank) if abs(exp_blank) > 1e-12 else 1.0
    exp_slopes = exp_slopes_raw * scale

    def cost(params):
        K_D = params[0]
        if K_D <= 0:
            return 1e12
        return float(np.sum((exp_slopes - _sim_slopes(K_D)) ** 2))

    best_cost, best_kd = np.inf, float(dart0)
    for x0 in np.logspace(-1, 5, 20):
        r = minimize(cost, [x0], method='Nelder-Mead',
                     options={'xatol': 1e-4, 'fatol': 1e-4, 'maxiter': 10000})
        if r.fun < best_cost:
            best_cost = r.fun
            best_kd = float(r.x[0])

    kd_fit = best_kd
    sim_slopes = _sim_slopes(kd_fit)
    residuals_arr = exp_slopes - sim_slopes
    kd_error = float(np.sqrt(np.var(residuals_arr)) / max(np.sqrt(len(residuals_arr)), 1))
    concentrations = [m * kd_fit for m in kd_multipliers]

    return {
        "kd_fit":         kd_fit,
        "kd_error":       kd_error,
        "exp_slopes":     exp_slopes.tolist(),
        "sim_slopes":     sim_slopes.tolist(),
        "concentrations": concentrations,
        "k_txn":          k_txn,
    }
