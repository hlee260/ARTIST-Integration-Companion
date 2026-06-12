import sys
import io
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox,
    QTextEdit, QTableWidget, QTableWidgetItem,
    QHeaderView, QDoubleSpinBox, QSlider,
    QSizePolicy, QGroupBox, QSplitter, QButtonGroup,
    QRadioButton, QScrollArea, QFileDialog, QProgressBar
)
from PySide6.QtCore import Qt, QThread, Signal, QRegularExpression
from PySide6.QtGui import QPixmap, QImage, QRegularExpressionValidator

from ARTISTIC import (
    sanitize_sequence, test_insulations, get_order_strands, get_digital_strands,
    simulate_digital_curve, simulate_analog_curve,
    lig_series, salt_rows, parse_salt,
    digital_titration_series, analog_titration_series,
    recommend_analog_dart, recommend_digital_ref,
)

# ── Worker thread ─────────────────────────────────────────────────────────────
class DesignWorker(QThread):
    done = Signal(object)
    error = Signal(str)
    def __init__(self, seq, gate):
        super().__init__()
        self.seq = seq
        self.gate = gate
    def run(self):
        try:
            self.done.emit(test_insulations(self.seq, gate=self.gate))
        except Exception as e:
            self.error.emit(str(e))

# ── RNA structure → QPixmap ───────────────────────────────────────────────────
def rna_to_pixmap(rna_seq, structure, w=500, h=400):
    try:
        from draw_rna.ipynb_draw import draw_struct
    except ImportError:
        return None
    fig, ax = plt.subplots(figsize=(w/100, h/100), dpi=100)
    fig.patch.set_facecolor("white")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        draw_struct(rna_seq, structure, ax=ax)
    ax.axis("off")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight",
                facecolor="white", pad_inches=0.1)
    plt.close(fig)
    buf.seek(0)
    return QPixmap.fromImage(QImage.fromData(buf.read()))

# ── Plot → QPixmap helper ─────────────────────────────────────────────────────
def fig_to_pixmap(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return QPixmap.fromImage(QImage.fromData(buf.read()))

# ── Sweep colors (matches HTML rainbow) ──────────────────────────────────────
SWEEP_CONCS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
SWEEP_COLORS = ["#ff4444","#ff6644","#ff8844","#ffaa44","#ffcc44",
                "#ccff44","#88ff44","#44ff88","#44ffaa","#44ffcc"]

def make_analog_recommend_plot(rec_curve, sweep_curves, result):
    """
    Plot the recommended analog curve prominently over a faded sweep, with the
    requested detection range shaded and the recommended L10/L90/EC50 marked.

    rec_curve    : (x, y) arrays for the recommended dART (high-res)
    sweep_curves : list of (dart, x, y) for context (faded)
    result       : dict returned by recommend_analog_dart
    """
    rec = result["recommended"]
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")

    # Faded sweep for context
    for dart, x, y in sweep_curves:
        ax.plot(x, y, color="#bbbbbb", linewidth=1.0, alpha=0.55, zorder=1)

    # Shaded requested detection range
    lo, hi = result["detect_low"], result["detect_high"]
    ax.axvspan(lo, hi, color="#ffe08a", alpha=0.35, zorder=0,
               label=f"Requested range {lo:.0f}–{hi:.0f} nM")

    # Recommended curve
    xr, yr = rec_curve
    ax.plot(xr, yr, color="#e6a800", linewidth=2.8, zorder=3,
            label=f"Recommended dART = {rec['dart']:.0f} nM")

    if rec.get("lo") is not None and rec.get("hi") is not None:
        ax.axvline(rec["lo"], color="red", lw=1.4, ls="--", zorder=2,
                   label=f"L90 = {rec['lo']:.1f} nM")
        ax.axvline(rec["hi"], color="red", lw=1.4, ls="--", zorder=2,
                   label=f"L10 = {rec['hi']:.1f} nM")
    if rec.get("ec50"):
        ax.axvline(rec["ec50"], color="green", lw=1.4, ls=":", zorder=2,
                   label=f"EC50 = {rec['ec50']:.1f} nM")

    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    hill = rec.get("hill")
    htxt = f"{hill:.2f}" if hill is not None else "N/A"
    ax.set_title(f"Analog design — recommended dART {rec['dart']:.0f} nM",
                 color="black", fontsize=12)
    ax.tick_params(colors="black"); ax.spines[:].set_color("black")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig_to_pixmap(fig)

def make_digital_recommend_plot(rec_curve, sweep_curves, result):
    """
    Plot the recommended digital curve over a faded reference-template sweep,
    with the requested threshold marked.
    """
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")

    for ref, x, y in sweep_curves:
        ax.plot(x, y, color="#bbbbbb", linewidth=1.0, alpha=0.55, zorder=1)

    target = result["target_threshold"]
    ax.axvline(target, color="#cc0000", lw=2.0, ls="--", zorder=2,
               label=f"Requested threshold = {target:.0f} nM")

    xr, yr = rec_curve
    rref = result["recommended_ref"]
    ax.plot(xr, yr, color="#e6a800", linewidth=2.8, zorder=3,
            label=f"Recommended Ref template = {rref:.0f} nM")

    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    ax.set_title(f"Digital design — recommended Ref template {rref:.0f} nM",
                 color="black", fontsize=12)
    ax.tick_params(colors="black"); ax.spines[:].set_color("black")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig_to_pixmap(fig)


# ── Main window ───────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ARTISTIC")
        self.resize(1400, 900)
        self._last_best = None
        self._last_gate = "O1"

        tabs = QTabWidget()
        tabs.addTab(self._scrollable(self._build_design_tab()), "dART Design")
        tabs.addTab(self._scrollable(self._build_kd_fit_tab()), "Kd,apparent fit")
        tabs.addTab(self._scrollable(self._build_analog_tab()), "Analog")
        tabs.addTab(self._scrollable(self._build_digital_tab()), "Digital")
        self.setCentralWidget(tabs)

    def _scrollable(self, widget):
        """Wrap a widget in a QScrollArea so it scrolls when the window is small."""
        scroll = QScrollArea()
        scroll.setWidget(widget)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        return scroll

    def _sync_kd(self, value, source):
        for spin in (self.design_kd, self.ana_kd, self.dig_kd):
            if spin is not source:
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)

    # ═══════════════════════════════════════════════════════════════════════════
    # DESIGN TAB
    # ═══════════════════════════════════════════════════════════════════════════
    def _build_design_tab(self):
        outer = QWidget()
        layout = QVBoxLayout(outer)

        input_group = QGroupBox("Inputs")
        grid = QGridLayout(input_group)

        grid.addWidget(QLabel("Aptamer sequence (DNA/RNA):"), 0, 0)
        self.apt_seq = QLineEdit()
        self.apt_seq.setPlaceholderText("e.g. GGGTTGGTGTGGTTGG")
        v = QRegularExpressionValidator()
        v.setRegularExpression(QRegularExpression("[ATCGUatcgu]*"))
        self.apt_seq.setValidator(v)
        self.apt_seq.textChanged.connect(self._force_upper)
        grid.addWidget(self.apt_seq, 0, 1)

        grid.addWidget(QLabel("Reported Kd (nM):"), 1, 0)
        self.design_kd = QDoubleSpinBox()
        self.design_kd.setRange(0.001, 1_000_000_000_000); self.design_kd.setDecimals(3)
        self.design_kd.setValue(25.0); self.design_kd.setSuffix(" nM")
        grid.addWidget(self.design_kd, 1, 1)


        grid.addWidget(QLabel("Salt / buffer condition:"), 2, 0)
        self.design_salt = QLineEdit("150 mM NaCl")
        self.design_salt.setPlaceholderText("e.g. 150 mM NaCl  or  5 mM MgCl2")
        grid.addWidget(self.design_salt, 2, 1)

        grid.addWidget(QLabel("Output domain:"), 3, 0)
        self.gate_select = QComboBox()
        self.gate_select.addItems(["O1", "O2", "O3"])
        grid.addWidget(self.gate_select, 3, 1)

        self.generate_btn = QPushButton("Generate Best Design")
        self.generate_btn.setFixedHeight(36)
        self.generate_btn.clicked.connect(self._on_generate)
        grid.addWidget(self.generate_btn, 4, 0, 1, 2)
        layout.addWidget(input_group)

        splitter = QSplitter(Qt.Horizontal)

        # Left
        left = QWidget(); ll = QVBoxLayout(left)

        rg = QGroupBox("Design Result")
        rgl = QVBoxLayout(rg)
        self.design_out = QTextEdit(); self.design_out.setReadOnly(True)
        self.design_out.setFontFamily("Courier New")
        rgl.addWidget(self.design_out)
        ll.addWidget(rg)

        og = QGroupBox("Strands to Order")
        ogl = QVBoxLayout(og)
        self.order_table = QTableWidget()
        self.order_table.setColumnCount(3)
        self.order_table.setHorizontalHeaderLabels(["Strand", "Sequence (5'→3')", "Notes"])
        self.order_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.order_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.order_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.order_table.setWordWrap(True); self.order_table.setMinimumHeight(140)
        ogl.addWidget(self.order_table)
        ll.addWidget(og)

        self.tit_group = QGroupBox("Titration Grid  (Salt × Ligand)")
        tgl = QVBoxLayout(self.tit_group)
        self.tit_grid = QTableWidget()
        self.tit_grid.setMinimumHeight(180)
        # make all cells read-only grey background so user fills them in
        self.tit_grid.setEditTriggers(QTableWidget.NoEditTriggers)
        tgl.addWidget(self.tit_grid)
        ll.addWidget(self.tit_group)
        splitter.addWidget(left)

        # Right: RNA structure
        right = QWidget(); rl = QVBoxLayout(right)
        sg2 = QGroupBox("RNA Transcript Structure")
        sg2l = QVBoxLayout(sg2)
        self.struct_label = QLabel("Structure will appear here after generation.")
        self.struct_label.setAlignment(Qt.AlignCenter)
        self.struct_label.setMinimumSize(400, 350)
        self.struct_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sg2l.addWidget(self.struct_label)
        rl.addWidget(sg2)
        splitter.addWidget(right)

        splitter.setSizes([650, 550])
        layout.addWidget(splitter)
        return outer

    def _force_upper(self, text):
        cur = self.apt_seq.cursorPosition()
        self.apt_seq.blockSignals(True)
        self.apt_seq.setText(text.upper())
        self.apt_seq.setCursorPosition(cur)
        self.apt_seq.blockSignals(False)

    def _on_generate(self):
        seq = sanitize_sequence(self.apt_seq.text())
        gate = self.gate_select.currentText()
        if len(seq) < 10:
            self.design_out.setText("⚠  Aptamer must be at least 10 nucleotides.")
            return
        self.generate_btn.setEnabled(False)
        self.generate_btn.setText("Running folding search…")
        self.design_out.setText("Designing dART — please wait…")
        self._last_gate = gate
        self._worker = DesignWorker(seq, gate)
        self._worker.done.connect(self._on_design_done)
        self._worker.error.connect(self._on_design_error)
        self._worker.start()

    def _on_design_error(self, msg):
        self.design_out.setText(f"Error: {msg}")
        self.generate_btn.setEnabled(True)
        self.generate_btn.setText("Generate Best Design")

    def _on_design_done(self, best):
        self.generate_btn.setEnabled(True)
        self.generate_btn.setText("Generate Best Design")
        if not best:
            self.design_out.setText("No valid design found.")
            return
        self._last_best = best
        gate = self._last_gate
        salt_str = self.design_salt.text().strip()

        if best.get("used_default"):
            ins_status = "GGGATG default — accepted (insulation stem forms)"
        else:
            ins_status = (f"GGGATG default rejected — fallback after searching "
                          f"{best.get('searched_alternatives', 0)} alternatives")

        lines = [
            f"Gate:          {gate}",
            f"Insulation:    {best['ins']}",
            f"Ins Comp:      {best['insComp']}",
            f"Ins status:    {ins_status}",
            f"Template len:  {len(best['template'])} nt",
            f"RNA len:       {len(best['rna'])} nt",
            f"MFE:           {best['mfe']:.2f} kcal/mol",
            f"Pair count:    {best['count']}",
            f"",
            f"Secondary structure (dot-bracket):",
            f"{best['structure']}",
            f"",
            f"RNA transcript:",
            f"{best['rna']}",
        ]
        self.design_out.setText("\n".join(lines))

        strands = get_order_strands(best, gate)
        # Reporter sequences by gate
        reporters = {
            "O1": [
                ("O1 Reporter F", "/5HEX/CTACATCCACATACTA",          "Fluorophore strand - order standard desalted or HPLC purified"),
                ("O1 Reporter Q", "GTTAATTAGTATGTGGATGTAG/3IAbRQSp/", "Quencher strand - order HPLC purified"),
            ],
            "O2": [
                ("O2 Reporter F", "/56-FAM/CTACTTTCACTTCACAA",        "Fluorophore strand  - order standard desalted or HPLC purified"),
                ("O2 Reporter Q", "TGATGTTGTGAAGTGAAAGTAG/3IABkFQ/",  "Quencher strand  - order HPLC purified"),
            ],
            "O3": [
                ("O3 Reporter F", "/56-FAM/TACCATCACATTCAAT",         "Fluorophore strand - order standard desalted or HPLC purified"),
                ("O3 Reporter Q", "AGGATTATTGAATGTGATGGTA/3IABkFQ/",  "Quencher strand - order HPLC purified"),
            ],
        }
        strands = strands + reporters[gate]
        self.order_table.setRowCount(len(strands))
        for i, (label, seq, notes) in enumerate(strands):
            self.order_table.setItem(i, 0, QTableWidgetItem(label))
            self.order_table.setItem(i, 1, QTableWidgetItem(seq))
            self.order_table.setItem(i, 2, QTableWidgetItem(notes))
        self.order_table.resizeRowsToContents()

        salt_name, salt_conc = parse_salt(salt_str)
        s_rows = salt_rows(salt_name)
        # If aptamer is RNA (contains U), add a fixed 2 mM MnCl2
        is_rna = "U" in sanitize_sequence(self.apt_seq.text())
        kd = self.design_kd.value()
        ligs = lig_series(kd)
        divalent = salt_name.strip() in ("MgCl2", "CaCl2", "MnCl2", "ZnCl2")
        kind = "divalent" if divalent else "monovalent"
        rna_note = " RNA aptamer used - add 2 mM MnCl2 to all conditions! " if is_rna else ""
        self.tit_group.setTitle(
            f"Titration Grid — {salt_name} ({kind}) × Ligand (0–100× Kd={kd:.1f} nM){rna_note}"
        )
        # columns = ligand concentrations, rows = salt concentrations
        self.tit_grid.setRowCount(len(s_rows))
        self.tit_grid.setColumnCount(len(ligs))
        self.tit_grid.setHorizontalHeaderLabels([l["label"] for l in ligs])
        self.tit_grid.setVerticalHeaderLabels(
            [f"{r['conc']} mM {r['salt']}" for r in s_rows]
        )
        
        for ri, srow in enumerate(s_rows):
            for ci, lig in enumerate(ligs):
                item = QTableWidgetItem(f"{srow['conc']} mM salt\n{lig['conc']:.3f} nM ligand")
                item.setTextAlignment(Qt.AlignCenter)
                self.tit_grid.setItem(ri, ci, item)
        self.tit_grid.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.tit_grid.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)

        px = rna_to_pixmap(best["rna"], best["structure"])
        if px:
            self.struct_label.setPixmap(
                px.scaled(self.struct_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )
        else:
            self.struct_label.setText("draw_rna not available — pip install draw_rna")

    # ═══════════════════════════════════════════════════════════════════════════
    # ANALOG SENSOR TAB
    # ═══════════════════════════════════════════════════════════════════════════
    def _build_analog_tab(self):
        w = QWidget(); layout = QVBoxLayout(w)

        # Parameters
        pg = QGroupBox("Analog Sensor Parameters")
        grid = QGridLayout(pg)

        grid.addWidget(QLabel("Kd (nM):"), 0, 0)
        self.ana_kd = QDoubleSpinBox()
        self.ana_kd.setRange(0.001, 1_000_000_000_000); self.ana_kd.setDecimals(3)
        self.ana_kd.setValue(25.0); self.ana_kd.setSuffix(" nM")
        grid.addWidget(self.ana_kd, 0, 1)

        grid.addWidget(QLabel("Transcription rate (k_txn):"), 1, 0)
        self.ana_ktxn = QDoubleSpinBox()
        self.ana_ktxn.setRange(0.0001, 1.0); self.ana_ktxn.setDecimals(4)
        self.ana_ktxn.setValue(0.006)
        grid.addWidget(self.ana_ktxn, 1, 1)

        # Desired detection range (the user input that drives the design)
        grid.addWidget(QLabel("Desired detection range (nM):"), 2, 0)
        range_row = QWidget(); rrl = QHBoxLayout(range_row); rrl.setContentsMargins(0,0,0,0)
        self.ana_detect_low = QDoubleSpinBox()
        self.ana_detect_low.setRange(0.001, 1_000_000); self.ana_detect_low.setDecimals(3)
        self.ana_detect_low.setValue(10.0); self.ana_detect_low.setSuffix(" nM")
        self.ana_detect_high = QDoubleSpinBox()
        self.ana_detect_high.setRange(0.001, 1_000_000); self.ana_detect_high.setDecimals(3)
        self.ana_detect_high.setValue(300.0); self.ana_detect_high.setSuffix(" nM")
        rrl.addWidget(QLabel("L90")); rrl.addWidget(self.ana_detect_low)
        rrl.addWidget(QLabel("L10")); rrl.addWidget(self.ana_detect_high)
        grid.addWidget(range_row, 2, 1)

        self.ana_run_btn = QPushButton("Find best dART for this range")
        self.ana_run_btn.setFixedHeight(36)
        self.ana_run_btn.clicked.connect(self._on_run_analog)
        grid.addWidget(self.ana_run_btn, 3, 0, 1, 2)

        self.ana_progress = QProgressBar()
        self.ana_progress.setRange(0, 0); self.ana_progress.setVisible(False)
        grid.addWidget(self.ana_progress, 4, 0, 1, 2)
        layout.addWidget(pg)

        # Results row
        results = QHBoxLayout()

        # Left: tables
        left = QWidget(); ll = QVBoxLayout(left)

        mg = QGroupBox("Recommended dART concentration")
        mgl = QVBoxLayout(mg)
        self.ana_metrics_table = QTableWidget()
        self.ana_metrics_table.setColumnCount(2)
        self.ana_metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.ana_metrics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ana_metrics_table.setMaximumHeight(230)
        mgl.addWidget(self.ana_metrics_table)
        ll.addWidget(mg)

        sg = QGroupBox("dART Sweep (detection range and EC50)")
        sgl = QVBoxLayout(sg)
        self.ana_sweep_table = QTableWidget()
        self.ana_sweep_table.setColumnCount(4)
        self.ana_sweep_table.setHorizontalHeaderLabels(
            ["dART (nM)", "L90 (nM)", "L10 (nM)", "EC50 (nM)"])
        self.ana_sweep_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ana_sweep_table.setMaximumHeight(220)
        sgl.addWidget(self.ana_sweep_table)
        ll.addWidget(sg)

        lg = QGroupBox("Ligand Titration Series")
        lgl = QVBoxLayout(lg)
        self.ana_lig_table = QTableWidget()
        self.ana_lig_table.setColumnCount(2)
        self.ana_lig_table.setHorizontalHeaderLabels(["Series point", "Concentration (nM)"])
        self.ana_lig_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        lgl.addWidget(self.ana_lig_table)
        ll.addWidget(lg)

        results.addWidget(left, 1)

        # Right: plot
        plot_g = QGroupBox("Dose-Response Curve")
        pgl = QVBoxLayout(plot_g)
        self.ana_plot_label = QLabel("Plot will appear after running.")
        self.ana_plot_label.setAlignment(Qt.AlignCenter)
        self.ana_plot_label.setMinimumSize(520, 380)
        self.ana_plot_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pgl.addWidget(self.ana_plot_label)
        results.addWidget(plot_g, 2)

        layout.addLayout(results)
        return w

    def _on_run_analog(self):
        kd    = self.ana_kd.value()
        k_txn = self.ana_ktxn.value()
        lo    = self.ana_detect_low.value()
        hi    = self.ana_detect_high.value()
        if hi <= lo:
            self.ana_plot_label.setText("⚠  Detection range 'L10' must exceed 'L90'.")
            return

        self.ana_run_btn.setEnabled(False)
        self.ana_progress.setVisible(True)
        self.ana_plot_label.setText("Sweeping dART concentrations — please wait…")

        self._ana_worker = _AnalogRecWorker(lo, hi, kd, k_txn)
        self._ana_worker.done.connect(self._on_analog_done)
        self._ana_worker.error.connect(self._on_analog_error)
        self._ana_worker.start()

    def _on_analog_error(self, msg):
        self.ana_run_btn.setEnabled(True)
        self.ana_progress.setVisible(False)
        self.ana_plot_label.setText(f"Error: {msg}")

    def _on_analog_done(self, payload):
        self.ana_run_btn.setEnabled(True)
        self.ana_progress.setVisible(False)

        result = payload["result"]
        rec    = result["recommended"]
        kd     = payload["kd"]
        k_txn  = payload["k_txn"]

        span_txt = (f"{rec['lo']:.2f} – {rec['hi']:.2f} nM"
                    if rec["lo"] is not None and rec["hi"] is not None else "N/A")
        cover_txt = "Yes — fully covers request" if result["covered"] else "No — closest available"
        rows = [
            ("Recommended dART", f"{rec['dart']:.0f} nM"),
            ("Requested range", f"{result['detect_low']:.2f} – {result['detect_high']:.2f} nM"),
            ("Achieved detection span", span_txt),
            ("Covers requested range?", cover_txt),
            ("Hill coefficient", f"{rec['hill']:.3f}" if rec["hill"] is not None else "N/A"),
            ("EC50", f"{rec['ec50']:.2f} nM" if rec["ec50"] else "N/A"),
            ("Max signal", f"{rec['maxY']:.2f} nM"),
            ("Kd", f"{kd:.3f} nM"),
            ("k_txn", f"{k_txn:.4f}"),
            ("RNase H", "6 U/mL recommended"),
        ]
        self.ana_metrics_table.setRowCount(len(rows))
        for i, (k, v) in enumerate(rows):
            self.ana_metrics_table.setItem(i, 0, QTableWidgetItem(k))
            self.ana_metrics_table.setItem(i, 1, QTableWidgetItem(v))

        # Sweep table (highlight recommended row)
        sweep = result["sweep"]
        self.ana_sweep_table.setRowCount(len(sweep))
        for i, s in enumerate(sweep):
            lo_s  = f"{s['lo']:.1f}" if s["lo"] is not None else "—"
            hi_s  = f"{s['hi']:.1f}" if s["hi"] is not None else "—"
            ec50_s = f"{s['ec50']:.1f}" if s["ec50"] is not None else "—"
            vals = [f"{s['dart']:.0f}", lo_s, hi_s, ec50_s]
            for j, v in enumerate(vals):
                item = QTableWidgetItem(v)
                if s["dart"] == rec["dart"]:
                    item.setBackground(Qt.green)
                self.ana_sweep_table.setItem(i, j, item)

        # Ligand series spanning the requested detection range
        ligs = analog_titration_series(result["detect_low"], result["detect_high"])
        self.ana_lig_table.setRowCount(len(ligs))
        for i, row in enumerate(ligs):
            self.ana_lig_table.setItem(i, 0, QTableWidgetItem(row["label"]))
            self.ana_lig_table.setItem(i, 1, QTableWidgetItem(f"{row['conc']:.3f}"))

        px = make_analog_recommend_plot(
            payload["rec_curve"], payload["sweep_curves"], result)
        if px:
            self.ana_plot_label.setPixmap(
                px.scaled(self.ana_plot_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )

    # ═══════════════════════════════════════════════════════════════════════════
    # DIGITAL SENSOR TAB
    # ═══════════════════════════════════════════════════════════════════════════
    def _build_digital_tab(self):
        w = QWidget(); layout = QVBoxLayout(w)

        pg = QGroupBox("Digital Sensor Parameters")
        grid = QGridLayout(pg)

        grid.addWidget(QLabel("Kd (nM):"), 0, 0)
        self.dig_kd = QDoubleSpinBox()
        self.dig_kd.setRange(0.001, 1_000_000_000_000); self.dig_kd.setDecimals(3)
        self.dig_kd.setValue(25.0); self.dig_kd.setSuffix(" nM")
        grid.addWidget(self.dig_kd, 0, 1)

        grid.addWidget(QLabel("k_txn Ref template:"), 1, 0)
        self.dig_ktxn_ref = QDoubleSpinBox()
        self.dig_ktxn_ref.setRange(0.0001, 1.0); self.dig_ktxn_ref.setDecimals(4)
        self.dig_ktxn_ref.setValue(0.01)
        grid.addWidget(self.dig_ktxn_ref, 1, 1)

        grid.addWidget(QLabel("k_txn Inverter dART:"), 2, 0)
        self.dig_ktxn_apt = QDoubleSpinBox()
        self.dig_ktxn_apt.setRange(0.0001, 1.0); self.dig_ktxn_apt.setDecimals(4)
        self.dig_ktxn_apt.setValue(0.01)
        grid.addWidget(self.dig_ktxn_apt, 2, 1)

        grid.addWidget(QLabel("Aptamer dART (nM, fixed):"), 3, 0)
        self.dig_apt = QDoubleSpinBox()
        self.dig_apt.setRange(0.001, 1000); self.dig_apt.setValue(50.0)
        self.dig_apt.setSuffix(" nM")
        grid.addWidget(self.dig_apt, 3, 1)

        # Desired threshold concentration (the user input that drives the design)
        grid.addWidget(QLabel("Desired threshold (nM):"), 4, 0)
        self.dig_threshold = QDoubleSpinBox()
        self.dig_threshold.setRange(0.001, 1_000_000); self.dig_threshold.setDecimals(3)
        self.dig_threshold.setValue(100.0); self.dig_threshold.setSuffix(" nM")
        self.dig_threshold.setToolTip(
            "Ligand concentration at which the output should switch fully ON.")
        grid.addWidget(self.dig_threshold, 4, 1)

        self.dig_run_btn = QPushButton("Find best Ref template for this threshold")
        self.dig_run_btn.setFixedHeight(36)
        self.dig_run_btn.clicked.connect(self._on_run_digital)
        grid.addWidget(self.dig_run_btn, 5, 0, 1, 2)

        self.dig_progress = QProgressBar()
        self.dig_progress.setRange(0, 0); self.dig_progress.setVisible(False)
        grid.addWidget(self.dig_progress, 6, 0, 1, 2)
        layout.addWidget(pg)

        results = QHBoxLayout()

        left = QWidget(); ll = QVBoxLayout(left)

        mg = QGroupBox("Recommended Design")
        mgl = QVBoxLayout(mg)
        self.dig_metrics_table = QTableWidget()
        self.dig_metrics_table.setColumnCount(2)
        self.dig_metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.dig_metrics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.dig_metrics_table.setMaximumHeight(200)
        mgl.addWidget(self.dig_metrics_table)
        ll.addWidget(mg)

        sg = QGroupBox("Reference Template Sweep (threshold per ref template)")
        sgl = QVBoxLayout(sg)
        self.dig_sweep_table = QTableWidget()
        self.dig_sweep_table.setColumnCount(2)
        self.dig_sweep_table.setHorizontalHeaderLabels(
            ["Ref template (nM)", "Threshold (nM)"])
        self.dig_sweep_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.dig_sweep_table.setMaximumHeight(200)
        sgl.addWidget(self.dig_sweep_table)
        ll.addWidget(sg)

        lg = QGroupBox("Ligand Titration Series")
        lgl = QVBoxLayout(lg)
        self.dig_lig_table = QTableWidget()
        self.dig_lig_table.setColumnCount(2)
        self.dig_lig_table.setHorizontalHeaderLabels(["Series point", "Concentration (nM)"])
        self.dig_lig_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        lgl.addWidget(self.dig_lig_table)
        ll.addWidget(lg)

        dsg = QGroupBox("Strands to Order (Reference and Inverter dART)")
        dsgl = QVBoxLayout(dsg)
        self.dig_strands_table = QTableWidget()
        self.dig_strands_table.setColumnCount(3)
        self.dig_strands_table.setHorizontalHeaderLabels(["Strand", "Sequence (5'→3')", "Notes"])
        self.dig_strands_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.dig_strands_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.dig_strands_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.dig_strands_table.setWordWrap(True)
        dsgl.addWidget(self.dig_strands_table)
        ll.addWidget(dsg)

        results.addWidget(left, 1)

        plot_g = QGroupBox("Dose-Response Curve")
        pgl = QVBoxLayout(plot_g)
        self.dig_plot_label = QLabel("Plot will appear after running.")
        self.dig_plot_label.setAlignment(Qt.AlignCenter)
        self.dig_plot_label.setMinimumSize(520, 380)
        self.dig_plot_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pgl.addWidget(self.dig_plot_label)
        results.addWidget(plot_g, 2)

        layout.addLayout(results)
        return w

    def _on_run_digital(self):
        kd        = self.dig_kd.value()
        k_txn_ref = self.dig_ktxn_ref.value()
        k_txn_apt = self.dig_ktxn_apt.value()
        apt_dart  = self.dig_apt.value()
        target    = self.dig_threshold.value()

        self.dig_run_btn.setEnabled(False)
        self.dig_progress.setVisible(True)
        self.dig_plot_label.setText("Sweeping reference-template concentrations — please wait…")

        self._dig_worker = _DigitalRecWorker(target, kd, k_txn_ref, k_txn_apt, apt_dart)
        self._dig_worker.done.connect(self._on_digital_done)
        self._dig_worker.error.connect(self._on_digital_error)
        self._dig_worker.start()

    def _on_digital_error(self, msg):
        self.dig_run_btn.setEnabled(True)
        self.dig_progress.setVisible(False)
        self.dig_plot_label.setText(f"Error: {msg}")

    def _on_digital_done(self, payload):
        self.dig_run_btn.setEnabled(True)
        self.dig_progress.setVisible(False)

        result    = payload["result"]
        kd        = payload["kd"]
        apt_dart  = payload["apt_dart"]
        rref      = result["recommended_ref"]
        target    = result["target_threshold"]
        predicted = result["predicted_threshold"]

        # Honest reporting: warn when the requested threshold can't be matched.
        if result["extrapolated"]:
            match_note = "extrapolated beyond swept range"
        elif abs(predicted - target) > 0.1 * max(target, 1e-9):
            match_note = "approximate"
        else:
            match_note = "matched"

        rows = [
            ("Recommended Ref template", f"{rref:.1f} nM"),
            ("Requested threshold", f"{target:.2f} nM"),
            ("Predicted threshold at ref. template", f"{predicted:.2f} nM ({match_note})"),
        ]
        self.dig_metrics_table.setRowCount(len(rows))
        for i, (k, v) in enumerate(rows):
            self.dig_metrics_table.setItem(i, 0, QTableWidgetItem(k))
            self.dig_metrics_table.setItem(i, 1, QTableWidgetItem(v))

        # Sweep table (highlight the swept row nearest the recommended ref).
        # recommended_ref is interpolated/extrapolated, so it may not equal any
        # swept value exactly — flag the closest swept template to it.
        sweep = result["sweep"]
        best_row = min(range(len(sweep)),
                       key=lambda k: abs(sweep[k]["ref"] - rref)) if sweep else -1
        self.dig_sweep_table.setRowCount(len(sweep))
        for i, s in enumerate(sweep):
            vals = [f"{s['ref']:.0f}", f"{s['threshold']:.1f}", f"{s['maxY']:.1f}"]
            for j, v in enumerate(vals):
                item = QTableWidgetItem(v)
                if i == best_row:
                    item.setBackground(Qt.green)
                    if result["extrapolated"]:
                        item.setToolTip(
                            f"Closest swept template to the recommended "
                            f"{rref:.1f} nM (recommendation is extrapolated "
                            f"beyond this sweep).")
                    else:
                        item.setToolTip(
                            f"Closest swept template to the recommended {rref:.1f} nM.")
                self.dig_sweep_table.setItem(i, j, item)

        # Titration series centered on the requested threshold
        series = digital_titration_series(target)
        self.dig_lig_table.setRowCount(len(series))
        for i, row in enumerate(series):
            self.dig_lig_table.setItem(i, 0, QTableWidgetItem(row["label"]))
            self.dig_lig_table.setItem(i, 1, QTableWidgetItem(f"{row['conc']:.3f}"))

        # Strands to order (uses design tab result if available)
        gate     = self._last_gate
        ins      = self._last_best["ins"]      if self._last_best else "GGGATG"
        ins_comp = self._last_best["insComp"]  if self._last_best else "CATCCC"
        strand_rows = get_digital_strands(gate, ins, ins_comp)
        self.dig_strands_table.setRowCount(len(strand_rows))
        for i, (label, seq, notes) in enumerate(strand_rows):
            self.dig_strands_table.setItem(i, 0, QTableWidgetItem(label))
            self.dig_strands_table.setItem(i, 1, QTableWidgetItem(seq))
            self.dig_strands_table.setItem(i, 2, QTableWidgetItem(notes))
        self.dig_strands_table.resizeRowsToContents()

        px = make_digital_recommend_plot(
            payload["rec_curve"], payload["sweep_curves"], result)
        if px:
            self.dig_plot_label.setPixmap(
                px.scaled(self.dig_plot_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )
    # ═══════════════════════════════════════════════════════════════════════════
    # Kd FITTING TAB
    # ═══════════════════════════════════════════════════════════════════════════
    def _build_kd_fit_tab(self):
        w = QWidget(); layout = QVBoxLayout(w)

        # ── Parameters ───────────────────────────────────────────────────────
        pg = QGroupBox("Experimental Parameters")
        grid = QGridLayout(pg)

        grid.addWidget(QLabel("dART concentration (nM):"), 0, 0)
        self.kdf_dart = QDoubleSpinBox()
        self.kdf_dart.setRange(0.001, 10000); self.kdf_dart.setDecimals(2)
        self.kdf_dart.setValue(25.0); self.kdf_dart.setSuffix(" nM")
        grid.addWidget(self.kdf_dart, 0, 1)

        grid.addWidget(QLabel("Total reporter concentration (nM):"), 1, 0)
        self.kdf_reporter = QDoubleSpinBox()
        self.kdf_reporter.setRange(0.001, 10000); self.kdf_reporter.setDecimals(2)
        self.kdf_reporter.setValue(250.0); self.kdf_reporter.setSuffix(" nM")
        grid.addWidget(self.kdf_reporter, 1, 1)

        grid.addWidget(QLabel("Slope time point (min):"), 2, 0)
        self.kdf_slope_time = QDoubleSpinBox()
        self.kdf_slope_time.setRange(1, 120); self.kdf_slope_time.setDecimals(1)
        self.kdf_slope_time.setValue(20.0); self.kdf_slope_time.setSuffix(" min")
        grid.addWidget(self.kdf_slope_time, 2, 1)

        grid.addWidget(QLabel("Reported Kd multipliers\n(comma-separated, matching Excel rows):"), 3, 0)
        self.kdf_multipliers = QLineEdit("0, 1, 2, 2.5, 5, 10, 25, 50, 100")
        self.kdf_multipliers.setToolTip(
            "One value per row in Excel (Kd multiples).\n"
            "First row must be 0 (blank)."
        )
        grid.addWidget(self.kdf_multipliers, 3, 1)

        # File picker
        grid.addWidget(QLabel("Excel file (.xlsx):"), 4, 0)
        file_row = QWidget(); frl = QHBoxLayout(file_row); frl.setContentsMargins(0,0,0,0)
        self.kdf_file_lbl = QLabel("No file selected")
        self.kdf_file_lbl.setStyleSheet("color: grey;")
        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(90)
        browse_btn.clicked.connect(self._kdf_browse)
        frl.addWidget(self.kdf_file_lbl, 1)
        frl.addWidget(browse_btn)
        grid.addWidget(file_row, 4, 1)

        self.kdf_run_btn = QPushButton("Load & Fit")
        self.kdf_run_btn.setFixedHeight(36)
        self.kdf_run_btn.clicked.connect(self._on_kdf_run)
        grid.addWidget(self.kdf_run_btn, 5, 0, 1, 2)

        self.kdf_progress = QProgressBar()
        self.kdf_progress.setRange(0, 0)   # indeterminate spinner
        self.kdf_progress.setVisible(False)
        grid.addWidget(self.kdf_progress, 6, 0, 1, 2)

        layout.addWidget(pg)

        # ── Results ───────────────────────────────────────────────────────────
        splitter = QSplitter(Qt.Horizontal)

        # Left: result text + slope comparison table
        left = QWidget(); ll = QVBoxLayout(left)

        rg = QGroupBox("Fit Results")
        rgl = QVBoxLayout(rg)
        self.kdf_result_text = QTextEdit()
        self.kdf_result_text.setReadOnly(True)
        self.kdf_result_text.setFontFamily("Courier New")
        self.kdf_result_text.setMaximumHeight(110)
        rgl.addWidget(self.kdf_result_text)
        ll.addWidget(rg)

        sg = QGroupBox("Slope Comparison (experimental vs simulated)")
        sgl = QVBoxLayout(sg)
        self.kdf_slope_table = QTableWidget()
        self.kdf_slope_table.setColumnCount(4)
        self.kdf_slope_table.setHorizontalHeaderLabels(
            ["Kd multiplier", "[Ligand] (nM)", "Exp slope", "Sim slope"]
        )
        self.kdf_slope_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        sgl.addWidget(self.kdf_slope_table)
        ll.addWidget(sg)
        splitter.addWidget(left)

        # Right: plots
        right = QWidget(); rl = QVBoxLayout(right)

        pg2 = QGroupBox("Normalised Traces")
        pg2l = QVBoxLayout(pg2)
        self.kdf_trace_label = QLabel("Traces will appear after fitting.")
        self.kdf_trace_label.setAlignment(Qt.AlignCenter)
        self.kdf_trace_label.setMinimumSize(480, 300)
        self.kdf_trace_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pg2l.addWidget(self.kdf_trace_label)
        rl.addWidget(pg2)

        pg3 = QGroupBox("Slope vs [Ligand]")
        pg3l = QVBoxLayout(pg3)
        self.kdf_slope_label = QLabel("Slope plot will appear after fitting.")
        self.kdf_slope_label.setAlignment(Qt.AlignCenter)
        self.kdf_slope_label.setMinimumSize(480, 260)
        self.kdf_slope_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pg3l.addWidget(self.kdf_slope_label)
        rl.addWidget(pg3)

        splitter.addWidget(right)
        splitter.setSizes([480, 580])
        layout.addWidget(splitter)

        self._kdf_filepath = None
        return w

    def _kdf_browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select Excel file", "", "Excel files (*.xlsx *.xls)"
        )
        if path:
            self._kdf_filepath = path
            self.kdf_file_lbl.setText(path.split("/")[-1].split("\\")[-1])
            self.kdf_file_lbl.setStyleSheet("")

    def _on_kdf_run(self):
        if not self._kdf_filepath:
            self.kdf_result_text.setText("⚠  Please select an Excel file first.")
            return
        try:
            mults_raw = [x.strip() for x in self.kdf_multipliers.text().split(",")]
            multipliers = [float(m) for m in mults_raw if m]
        except ValueError:
            self.kdf_result_text.setText("⚠  Invalid Kd multipliers — use comma-separated numbers.")
            return

        self.kdf_run_btn.setEnabled(False)
        self.kdf_progress.setVisible(True)
        self.kdf_result_text.setText("Loading and fitting — please wait…")

        dart     = self.kdf_dart.value()
        reporter = self.kdf_reporter.value()
        slope_t  = self.kdf_slope_time.value()

        self._kdf_worker = _KdFitWorker(
            self._kdf_filepath, dart, reporter, slope_t, multipliers
        )
        self._kdf_worker.done.connect(self._on_kdf_done)
        self._kdf_worker.error.connect(self._on_kdf_error)
        self._kdf_worker.start()

    def _on_kdf_error(self, msg):
        self.kdf_run_btn.setEnabled(True)
        self.kdf_progress.setVisible(False)
        self.kdf_result_text.setText(f"Error: {msg}")

    def _on_kdf_done(self, result):
        self.kdf_run_btn.setEnabled(True)
        self.kdf_progress.setVisible(False)

        time_arr   = result["time"]
        norm_mat   = result["norm_matrix"]
        fit        = result["fit"]
        multipliers = result["multipliers"]
        k_txn      = fit["k_txn"]
        kd_fit     = fit["kd_fit"]
        kd_error   = fit["kd_error"]

        # ── Result text ───────────────────────────────────────────────────
        self.kdf_result_text.setText(
            f"Fitted k_txn (basal):  {k_txn:.6f} min⁻¹\n"
            f"Fitted Kd:             {kd_fit:.2f} ± {kd_error:.2f} nM\n"
            f"Slope time point:      {result['slope_time']} min\n"
            f"Data dimensions:       {norm_mat.shape[0]} timepoints × {norm_mat.shape[1]} concentrations"
        )

        # ── Slope table ───────────────────────────────────────────────────
        exp_s = fit["exp_slopes"]
        sim_s = fit["sim_slopes"]
        concs = fit["concentrations"]
        self.kdf_slope_table.setRowCount(len(multipliers))
        for i, (mult, conc, es, ss) in enumerate(zip(multipliers, concs, exp_s, sim_s)):
            lbl = "0 (blank)" if mult == 0 else f"{mult}× Kd"
            self.kdf_slope_table.setItem(i, 0, QTableWidgetItem(lbl))
            self.kdf_slope_table.setItem(i, 1, QTableWidgetItem(f"{conc:.2f}"))
            self.kdf_slope_table.setItem(i, 2, QTableWidgetItem(f"{es:.4f}"))
            self.kdf_slope_table.setItem(i, 3, QTableWidgetItem(f"{ss:.4f}"))

        # ── Trace plot ────────────────────────────────────────────────────
        colors = plt.cm.viridis(np.linspace(0, 1, norm_mat.shape[1]))
        fig, ax = plt.subplots(figsize=(5.5, 3.5), dpi=110)
        ax.set_facecolor("white"); fig.patch.set_facecolor("white")
        for ci in range(norm_mat.shape[1]):
            lbl = "0 (blank)" if multipliers[ci] == 0 else f"{multipliers[ci]}× Kd"
            ax.plot(time_arr, norm_mat[:, ci], color=colors[ci], lw=1.8, label=lbl)
        ax.set_xlabel("Time (min)", color="black", fontsize=10)
        ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=10)
        ax.set_title("Normalized Reacted reporter kinetics", color="black", fontsize=11)
        ax.tick_params(colors="black"); ax.spines[:].set_color("black")
        ax.legend(fontsize=6.5, labelcolor="black", ncol=2)
        fig.tight_layout()
        px = fig_to_pixmap(fig)
        if px:
            self.kdf_trace_label.setPixmap(
                px.scaled(self.kdf_trace_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )

        # ── Slope comparison plot ─────────────────────────────────────────
        fig2, ax2 = plt.subplots(figsize=(5.5, 3.0), dpi=110)
        ax2.set_facecolor("white"); fig2.patch.set_facecolor("white")
        ax2.plot(concs, exp_s, "o", color="#ffcc44", lw=0, ms=6, label="Experimental")
        ax2.plot(concs, sim_s, linestyle = "--", color="red", lw=2, ms=6, label=f"Fit (Kd={kd_fit:.1f} nM)")
        ax2.set_xlabel("[Ligand] (nM)", color="black", fontsize=10)
        ax2.set_xscale("log")
        ax2.set_ylabel("Txn rate at t="+f"{result['slope_time']:.0f}"+" min (nM/min)", color="black", fontsize=10)
        ax2.set_title("Dose-response curve", color="black", fontsize=11)
        ax2.tick_params(colors="black"); ax2.spines[:].set_color("black")
        ax2.legend(fontsize=8, labelcolor="black",
                   facecolor="white")
        fig2.tight_layout()
        px2 = fig_to_pixmap(fig2)
        if px2:
            self.kdf_slope_label.setPixmap(
                px2.scaled(self.kdf_slope_label.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            )


# ── Kd Fitting worker thread ──────────────────────────────────────────────────
class _KdFitWorker(QThread):
    done  = Signal(object)
    error = Signal(str)

    def __init__(self, filepath, dart, reporter, slope_time, multipliers):
        super().__init__()
        self.filepath    = filepath
        self.dart        = dart
        self.reporter    = reporter
        self.slope_time  = slope_time
        self.multipliers = multipliers

    def run(self):
        try:
            import pandas as pd
            from ARTISTIC import normalize_rfu, fit_basal_ktxn, fit_kd_from_slopes

            df = pd.read_excel(self.filepath, header=None)

            # First column = time, remaining columns = ligand concentrations
            first_col = pd.to_numeric(df.iloc[:, 0], errors="coerce")
            first_data_row = first_col.first_valid_index()
            if first_data_row is None:
                raise ValueError("No numeric time values found in the first column.")
            if first_data_row > 0:
                df = df.iloc[first_data_row:].reset_index(drop=True)
            time_raw = pd.to_numeric(df.iloc[:, 0], errors="coerce").dropna()
            time_arr = time_raw.values.astype(float)
            time_arr = time_arr - time_arr[0]
            rfu_df = df.iloc[:len(time_arr), 1:].apply(pd.to_numeric, errors="coerce").fillna(0)
            rfu_raw = rfu_df.values.astype(float)

            # Trim to match multipliers count
            n_cols = len(self.multipliers)
            if rfu_raw.shape[1] < n_cols:
                raise ValueError(
                    f"Excel has {rfu_raw.shape[1]} data columns but "
                    f"{n_cols} Kd multipliers were specified."
                )
            rfu_raw = rfu_raw[:, :n_cols]

            norm_mat = normalize_rfu(rfu_raw, self.reporter)

            # Fit k_txn from the raw (un-normalised) blank trace — normalisation
            # removes the amplitude information that k_txn depends on
            raw_blank = rfu_raw[:, 0].astype(float)
            raw_blank = raw_blank - raw_blank[0]   # offset to zero
            k_txn = fit_basal_ktxn(time_arr, raw_blank,
                                   self.dart, self.reporter)

            # Fit Kd from raw RFU slopes (slope shape preserved regardless of RFU scale)
            fit = fit_kd_from_slopes(
                time_arr, rfu_raw, self.dart, self.reporter,
                k_txn, self.multipliers, self.slope_time
            )

            self.done.emit({
                "time":        time_arr,
                "norm_matrix": norm_mat,
                "fit":         fit,
                "multipliers": self.multipliers,
                "slope_time":  self.slope_time,
            })
        except Exception as e:
            import traceback
            self.error.emit(traceback.format_exc())


# ── Analog recommendation worker thread ───────────────────────────────────────
class _AnalogRecWorker(QThread):
    done  = Signal(object)
    error = Signal(str)

    def __init__(self, detect_low, detect_high, kd, k_txn):
        super().__init__()
        self.detect_low  = detect_low
        self.detect_high = detect_high
        self.kd          = kd
        self.k_txn       = k_txn

    def run(self):
        try:
            from ARTISTIC import recommend_analog_dart, simulate_analog_curve
            result = recommend_analog_dart(
                self.detect_low, self.detect_high, self.kd, self.k_txn)
            rec = result["recommended"]

            # High-resolution recommended curve for the foreground plot.
            fine = np.logspace(0, 4, 400)
            xr, yr, _ = simulate_analog_curve(rec["dart"], self.kd, self.k_txn, L0_list=fine)

            # A handful of faded context curves (subset of the sweep).
            coarse = np.logspace(0, 4, 160)
            sweep = result["sweep"]
            pick = sweep[::max(len(sweep) // 6, 1)]
            sweep_curves = []
            for s in pick:
                x, y, _ = simulate_analog_curve(s["dart"], self.kd, self.k_txn, L0_list=coarse)
                sweep_curves.append((s["dart"], x, y))

            self.done.emit({
                "result":       result,
                "kd":           self.kd,
                "k_txn":        self.k_txn,
                "rec_curve":    (xr, yr),
                "sweep_curves": sweep_curves,
            })
        except Exception as e:
            import traceback
            self.error.emit(traceback.format_exc())


# ── Digital recommendation worker thread ──────────────────────────────────────
class _DigitalRecWorker(QThread):
    done  = Signal(object)
    error = Signal(str)

    def __init__(self, target_threshold, kd, k_txn_ref, k_txn_apt, apt_dart):
        super().__init__()
        self.target_threshold = target_threshold
        self.kd               = kd
        self.k_txn_ref        = k_txn_ref
        self.k_txn_apt        = k_txn_apt
        self.apt_dart         = apt_dart

    def run(self):
        try:
            from ARTISTIC import recommend_digital_ref, simulate_digital_curve
            result = recommend_digital_ref(
                self.target_threshold, self.kd,
                self.k_txn_ref, self.k_txn_apt, self.apt_dart)
            rref = result["recommended_ref"]

            xr, yr, _ = simulate_digital_curve(
                rref, self.apt_dart, self.kd, self.k_txn_ref, self.k_txn_apt)

            sweep = result["sweep"]
            pick = sweep[::max(len(sweep) // 6, 1)]
            sweep_curves = []
            for s in pick:
                x, y, _ = simulate_digital_curve(
                    s["ref"], self.apt_dart, self.kd, self.k_txn_ref, self.k_txn_apt)
                sweep_curves.append((s["ref"], x, y))

            self.done.emit({
                "result":       result,
                "kd":           self.kd,
                "k_txn_ref":    self.k_txn_ref,
                "k_txn_apt":    self.k_txn_apt,
                "apt_dart":     self.apt_dart,
                "rec_curve":    (xr, yr),
                "sweep_curves": sweep_curves,
            })
        except Exception as e:
            import traceback
            self.error.emit(traceback.format_exc())


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
