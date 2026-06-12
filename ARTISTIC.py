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
# GGGATG is the canonical default insulation domain (ALL_INS[0] / ALL_COMP[0]).
# It is ALWAYS tried first. The remaining 43 sequences are only attempted as a
# fallback when GGGATG fails to insulate the aptamer.
#
# "Works" criterion: the 6-bp insulation stem (5' GGGAUG ... CAUCCC) must form in
# the predicted transcript. We measure this with `count` — the number of base
# pairs (from the partition function) between the 5' insulation region and its
# downstream complement. WORKS_MIN_PAIRS is how many of those ~6 pairs must
# register for the default to be accepted. Lower it to make GGGATG "stick" more
# often; raise it to be stricter about insulation before accepting the default.
WORKS_MIN_PAIRS = 4

def _evaluate_insulation(apt: str, ins: str, ins_comp: str, gate: str,
                         bp_threshold: float) -> dict:
    """Fold one insulation candidate and score how well its insulation stem forms."""
    cand = build_dart(apt, ins, ins_comp, gate)
    rna = cand["rna"]

    struct, mfe = RNA.fold(rna)
    fc = RNA.fold_compound(rna)
    fc.pf()
    bpp = fc.bpp()

    start = range(7)
    target = range(6 + len(apt), 6 + len(apt) + 8)
    count = sum(1 for i in start for j in target if bpp[i][j] > bp_threshold)

    return {
        **cand,
        "structure": struct,
        "mfe": mfe,
        "bpp": bpp,
        "count": count,
    }

def test_insulations(aptamer_seq: str, gate="O1", bp_threshold=0.5,
                     works_min_pairs: int = WORKS_MIN_PAIRS):
    """
    Design a dART, defaulting to the GGGATG insulation domain.

    Strategy:
      1. Build the GGGATG (default) design first.
      2. If its insulation stem forms well enough (count >= works_min_pairs),
         return it immediately — no further search is performed.
      3. Otherwise, sweep the remaining insulation sequences and return the
         best fallback, ranked first by how well the insulation stem forms
         (count) and then by minimum free energy (mfe) as a tie-breaker.

    The returned dict carries two extra keys:
      - "used_default": True if GGGATG was accepted, False if a fallback was used.
      - "searched_alternatives": number of fallback sequences evaluated (0 when
        the default worked).
    """
    apt = sanitize_sequence(aptamer_seq)
    RNA.cvar.temperature = 37

    # 1) Always evaluate the canonical default (GGGATG) first.
    default_ins, default_comp = ALL_INS[0], ALL_COMP[0]   # "GGGATG" / "CATCCC"
    default_cand = _evaluate_insulation(apt, default_ins, default_comp, gate, bp_threshold)

    if default_cand["count"] >= works_min_pairs:
        default_cand["used_default"] = True
        default_cand["searched_alternatives"] = 0
        return default_cand

    # 2) Default insulation did not form well — fall back to the other sequences.
    best = default_cand
    n_searched = 0
    for ins, ins_comp in zip(ALL_INS[1:], ALL_COMP[1:]):
        cand = _evaluate_insulation(apt, ins, ins_comp, gate, bp_threshold)
        n_searched += 1
        # Prefer better insulation (higher count); break ties by lower MFE.
        if (cand["count"] > best["count"]
                or (cand["count"] == best["count"] and cand["mfe"] < best["mfe"])):
            best = cand

    best["used_default"] = (best["ins"] == default_ins)
    best["searched_alternatives"] = n_searched
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

def analog_titration_series(detect_low, detect_high, n=7):
    """
    Recommended experimental ligand concentrations spanning a desired analog
    detection range. Returns a blank plus (n-1) log-spaced points from
    detect_low to detect_high inclusive.
    """
    if not detect_low or not detect_high or detect_low <= 0 or detect_high <= 0:
        return []
    lo, hi = sorted((float(detect_low), float(detect_high)))
    pts = np.logspace(np.log10(lo), np.log10(hi), max(int(n) - 1, 1))
    series = [0.0] + pts.tolist()
    return [
        {"label": "0 (blank)" if i == 0 else f"Pt {i}", "conc": float(v)}
        for i, v in enumerate(series)
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
def simulate_analog_curve(dart0, kd, k_txn_user, L0_list=None):
    # k_txn_user is per-second rate; multiply by 60 to get per-minute (matches HTML)
    k_txn = k_txn_user * 60
    k_hyb = 1e6 / 1e9 * 60
    k_deg = 0.001 * 60
    k_sd  = 1e4 / 1e9 * 60
    FQ0   = 250
    t = np.linspace(0, 120, 120)
    # Ligand grid. Default reproduces the original dense linear grid; callers
    # (e.g. the dART sweep) may pass a coarser log-spaced grid for speed.
    if L0_list is None:
        L0_list = np.arange(0, 10000)
    else:
        L0_list = np.asarray(L0_list, dtype=float)

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
    k_deg     = 0.000 * 60
    Reporter0 = 250
    t = np.linspace(0, 120, 120)
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

# ── Design recommendation (Analog): detection range → best dART ───────────────
# Candidate dART concentrations swept when recommending an analog design.
ANALOG_DART_SWEEP = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]

def recommend_analog_dart(detect_low, detect_high, kd, k_txn,
                          dart_candidates=None, sweep_points=160):
    """
    Recommend the dART concentration for a desired analog detection range.

    The user specifies the ligand detection range they want to sense across,
    [detect_low, detect_high] (nM). We sweep dART concentrations, simulate each
    dose-response curve, and read off its detection range (L10–L90), Hill
    coefficient (sigmoid slope) and EC50.

    Objective:
      Among dART concentrations whose L10–L90 span covers the requested range
      (L10 <= detect_low and L90 >= detect_high), pick the one with the LARGEST
      Hill coefficient — i.e. the steepest analog response that still spans the
      whole requested range, giving the best signal change per unit ligand.
      If no candidate fully covers the range, fall back to the one whose
      detection range is closest to the request (smallest endpoint mismatch).

    Returns a dict:
      {
        "sweep":        [ {dart, L10, L90, hill, ec50, maxY}, ... ],
        "recommended":  <the chosen sweep row>,
        "reason":       <str explaining why it was chosen>,
        "covered":      <bool, whether the requested range was fully covered>,
        "detect_low":   detect_low,
        "detect_high":  detect_high,
      }
    """
    if detect_low is None or detect_high is None or detect_low <= 0 or detect_high <= 0:
        raise ValueError("Detection range must be two positive ligand concentrations.")
    if detect_high < detect_low:
        detect_low, detect_high = detect_high, detect_low

    if dart_candidates is None:
        dart_candidates = ANALOG_DART_SWEEP

    # Coarser log-spaced ligand grid keeps the multi-curve sweep fast while still
    # resolving L10/L90 across several decades.
    grid = np.logspace(0, 4, int(sweep_points))

    sweep = []
    for dart in dart_candidates:
        _, _, m = simulate_analog_curve(float(dart), kd, k_txn, L0_list=grid)
        # The analog sensor is an inverse (OFF) response: output falls as ligand
        # rises, so L10 sits at high ligand and L90 at low ligand. Store a sorted
        # detection span [lo, hi] so range logic is direction-agnostic.
        lo = hi = None
        if m["L10"] is not None and m["L90"] is not None:
            lo, hi = sorted((m["L10"], m["L90"]))
        sweep.append({
            "dart": float(dart),
            "L10":  m["L10"],
            "L90":  m["L90"],
            "lo":   lo,
            "hi":   hi,
            "hill": m["slope"],
            "ec50": m["inflection"],
            "maxY": m["maxY"],
        })

    def covers(s):
        return (s["lo"] is not None and s["hi"] is not None
                and s["lo"] <= detect_low and s["hi"] >= detect_high)

    covering = [s for s in sweep if covers(s) and s["hill"] is not None]
    if covering:
        recommended = max(covering, key=lambda s: s["hill"])
        reason = ("Steepest response (max Hill coefficient) whose detection span "
                  "still covers the requested range.")
        covered = True
    else:
        def mismatch(s):
            if s["lo"] is None or s["hi"] is None:
                return float("inf")
            return abs(s["lo"] - detect_low) + abs(s["hi"] - detect_high)
        recommended = min(sweep, key=mismatch)
        reason = ("No dART fully covers the requested range; chose the dART whose "
                  "detection span is closest to the request.")
        covered = False

    return {
        "sweep":       sweep,
        "recommended": recommended,
        "reason":      reason,
        "covered":     covered,
        "detect_low":  float(detect_low),
        "detect_high": float(detect_high),
    }

# ── Design recommendation (Digital): threshold → best reference template ──────
# Candidate reference-template concentrations swept when recommending a digital design.
DIGITAL_REF_SWEEP = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

def recommend_digital_ref(target_threshold, kd, k_txn_ref, k_txn_apt, apt_dart,
                          ref_candidates=None):
    """
    Recommend the reference-template concentration for a desired digital threshold.

    The user specifies the ligand concentration at which the output should switch
    fully ON (the threshold of the digital response). We sweep reference-template
    concentrations, compute each curve's threshold, then interpolate (or, outside
    the swept span, linearly extrapolate) to find the reference-template
    concentration whose threshold matches the requested value.

    Returns a dict:
      {
        "sweep":            [ {ref, threshold, maxY}, ... ],
        "recommended_ref":  <float, interpolated/extrapolated ref template nM>,
        "predicted_threshold": <float, threshold expected at recommended_ref>,
        "nearest":          <the swept row whose threshold is closest to target>,
        "extrapolated":     <bool, True if target lies outside the swept thresholds>,
        "target_threshold": target_threshold,
      }
    """
    if target_threshold is None or target_threshold <= 0:
        raise ValueError("Target threshold must be a positive ligand concentration.")

    if ref_candidates is None:
        ref_candidates = DIGITAL_REF_SWEEP

    sweep = []
    for ref in ref_candidates:
        _, _, m = simulate_digital_curve(float(ref), apt_dart, kd, k_txn_ref, k_txn_apt)
        if m["threshold"] is not None:
            sweep.append({
                "ref":       float(ref),
                "threshold": float(m["threshold"]),
                "maxY":      float(m["maxY"]),
            })

    if not sweep:
        raise ValueError("Could not determine a threshold for any reference-template "
                         "concentration with these parameters.")

    # Nearest swept point (always available as a robust fallback).
    nearest = min(sweep, key=lambda s: abs(s["threshold"] - target_threshold))

    # Build (threshold -> ref) relationship, sorted and de-duplicated by threshold,
    # so we can interpolate the ref that yields the requested threshold.
    pts = sorted(sweep, key=lambda s: s["threshold"])
    th = np.array([p["threshold"] for p in pts])
    rf = np.array([p["ref"] for p in pts])
    # Collapse duplicate thresholds (keep mean ref) so interpolation is well-defined.
    uniq_th, idx = np.unique(th, return_inverse=True)
    uniq_rf = np.array([rf[idx == k].mean() for k in range(len(uniq_th))])

    extrapolated = bool(target_threshold < uniq_th[0] or target_threshold > uniq_th[-1])

    if len(uniq_th) == 1:
        recommended_ref = float(uniq_rf[0])
    elif not extrapolated:
        recommended_ref = float(np.interp(target_threshold, uniq_th, uniq_rf))
    else:
        # Linear extrapolation from the two nearest end points.
        if target_threshold < uniq_th[0]:
            x0, x1, y0, y1 = uniq_th[0], uniq_th[1], uniq_rf[0], uniq_rf[1]
        else:
            x0, x1, y0, y1 = uniq_th[-2], uniq_th[-1], uniq_rf[-2], uniq_rf[-1]
        slope = (y1 - y0) / (x1 - x0) if x1 != x0 else 0.0
        recommended_ref = float(y0 + slope * (target_threshold - x0))

    # Predicted threshold at the recommended ref (linear, extrapolating at ends
    # so it stays consistent with an extrapolated recommended_ref).
    ref_sorted = sorted(sweep, key=lambda s: s["ref"])
    rr = np.array([p["ref"] for p in ref_sorted])
    tt = np.array([p["threshold"] for p in ref_sorted])
    if len(rr) >= 2:
        if recommended_ref <= rr[0]:
            x0, x1, y0, y1 = rr[0], rr[1], tt[0], tt[1]
        elif recommended_ref >= rr[-1]:
            x0, x1, y0, y1 = rr[-2], rr[-1], tt[-2], tt[-1]
        else:
            predicted_threshold = float(np.interp(recommended_ref, rr, tt))
            x0 = None
        if x0 is not None:
            slope = (y1 - y0) / (x1 - x0) if x1 != x0 else 0.0
            predicted_threshold = float(y0 + slope * (recommended_ref - x0))
    else:
        predicted_threshold = float(tt[0])

    return {
        "sweep":               sweep,
        "recommended_ref":     recommended_ref,
        "predicted_threshold": predicted_threshold,
        "nearest":             nearest,
        "extrapolated":        extrapolated,
        "target_threshold":    float(target_threshold),
    }

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
