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
    digital_titration_series
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

def make_analog_plot(x, y, metrics, dart_nM, title_suffix=""):
    """Single analog curve with annotation lines."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")
    ax.plot(x, y, color="#ffcc44", linewidth=2.5)
    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    ax.set_title(f"Analog response — {dart_nM} nM dART{title_suffix}", color="black", fontsize=12)
    ax.tick_params(colors="black"); ax.spines[:].set_color("black")
    ax.grid(True, alpha=0.15, color="black")
    if metrics.get("L10") is not None:
        ax.axvline(metrics["L10"], color="#00ffa0", lw=1.8, ls="--",
                   label=f"L10 = {metrics['L10']:.1f} nM")
    if metrics.get("L90") is not None:
        ax.axvline(metrics["L90"], color="#00ffa0", lw=1.8, ls="--",
                   label=f"L90 = {metrics['L90']:.1f} nM")
    if metrics.get("inflection"):
        ax.axvline(metrics["inflection"], color="#ff44ff", lw=1.8, ls=":",
                   label=f"EC50 = {metrics['inflection']:.1f} nM")
    ax.legend(fontsize=15)
    fig.tight_layout()
    return fig_to_pixmap(fig)

def make_analog_sweep_plot(kd, k_txn):
    """Multi-curve analog sweep: dART 10–100 nM."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")
    for dart, color in zip(SWEEP_CONCS, SWEEP_COLORS):
        x, y, _ = simulate_analog_curve(dart, kd, k_txn)
        ax.plot(x, y, color=color, linewidth=1.8, alpha=0.85, label=f"{dart} nM")
    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    ax.set_title("Analog response — Sweep dART 10–100 nM", color="black", fontsize=12)
    ax.tick_params(colors="black"); ax.spines[:].set_color("black")
    ax.grid(True, alpha=0.15, color="black")
    ax.legend(fontsize=8, ncol=2, title="[dART]", title_fontsize=8)
    fig.tight_layout()
    return fig_to_pixmap(fig)

def make_digital_plot(x, y, metrics, ref_nM):
    """Single digital curve."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")
    ax.plot(x, y, color="#ffcc44", linewidth=2.5)
    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    ax.set_title(f"Digital response — Reference template {ref_nM} nM", color="black", fontsize=12)
    ax.tick_params(colors="black"); ax.spines[:].set_color("black")
    ax.grid(True, alpha=0.15, color="black", which="both")
    if metrics.get("threshold"):
        ax.axvline(metrics["threshold"], color="#ff4444", lw=1.8, ls="--",
                   label=f"Threshold ≈ {metrics['threshold']:.1f} nM")
    ax.legend(fontsize=15)
    fig.tight_layout()
    return fig_to_pixmap(fig)

def make_digital_sweep_plot(kd, k_txn_ref, k_txn_apt, apt_dart):
    """Multi-curve digital sweep: Ref dART 10–100 nM."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=110)
    ax.set_facecolor("white"); fig.patch.set_facecolor("white")
    for ref, color in zip(SWEEP_CONCS, SWEEP_COLORS):
        x, y, metrics = simulate_digital_curve(ref, apt_dart, kd, k_txn_ref, k_txn_apt)
        ax.plot(x, y, color=color, linewidth=1.8, alpha=0.85, label=f"Ref {ref} nM")
        if metrics.get("threshold"):
            ax.axvline(metrics["threshold"], color=color, lw=1.2, ls="--", alpha=0.45)
    ax.set_xscale("log")
    ax.set_xlabel("[Ligand] (nM, log scale)", color="black", fontsize=11)
    ax.set_ylabel("[Reacted Reporter] (nM)", color="black", fontsize=11)
    ax.set_title("Digital response — Sweep Ref dART 10–100 nM", color="black", fontsize=12)
    ax.tick_params(colors="#7faaa0"); ax.spines[:].set_color("black")
    ax.grid(True, alpha=0.15, color="black", which="both")
    ax.legend(fontsize=8, labelcolor="#c8ede7", ncol=2,
              title="[Ref template]", title_fontsize=8)
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
        tabs.addTab(self._scrollable(self._build_design_tab()), "Design")
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

        lines = [
            f"Gate:          {gate}",
            f"Insulation:    {best['ins']}",
            f"Ins Comp:      {best['insComp']}",
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
        self.order_table.setRowCount(len(strands))
        for i, (label, seq, notes) in enumerate(strands):
            self.order_table.setItem(i, 0, QTableWidgetItem(label))
            self.order_table.setItem(i, 1, QTableWidgetItem(seq))
            self.order_table.setItem(i, 2, QTableWidgetItem(notes))
        self.order_table.resizeRowsToContents()

        salt_name, salt_conc = parse_salt(salt_str)
        s_rows = salt_rows(salt_name)
        kd = self.design_kd.value()
        ligs = lig_series(kd)
        divalent = salt_name.strip() in ("MgCl2", "CaCl2", "MnCl2", "ZnCl2")
        kind = "divalent" if divalent else "monovalent"
        self.tit_group.setTitle(
            f"Titration Grid — {salt_name} ({kind}) × Ligand (0–100× Kd={kd:.1f} nM)"
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

        # Curve mode: single slider vs sweep
        grid.addWidget(QLabel("Curve mode:"), 3, 0)
        mode_row = QWidget(); mrl = QHBoxLayout(mode_row); mrl.setContentsMargins(0,0,0,0)
        self.ana_single_btn = QRadioButton("Single dART concentration")
        self.ana_sweep_btn  = QRadioButton("Sweep dART 10–100 nM")
        self.ana_single_btn.setChecked(True)
        self.ana_single_btn.toggled.connect(self._ana_mode_changed)
        mrl.addWidget(self.ana_single_btn); mrl.addWidget(self.ana_sweep_btn)
        grid.addWidget(mode_row, 3, 1)

        # Slider row (shown only in single mode)
        self.ana_slider_widget = QWidget()
        slw = QHBoxLayout(self.ana_slider_widget); slw.setContentsMargins(0,0,0,0)
        slw.addWidget(QLabel("dART conc:"))
        self.ana_dart_slider = QSlider(Qt.Horizontal)
        self.ana_dart_slider.setRange(10, 100); self.ana_dart_slider.setSingleStep(5)
        self.ana_dart_slider.setPageStep(10); self.ana_dart_slider.setValue(25)
        self.ana_dart_slider.setTickInterval(10); self.ana_dart_slider.setTickPosition(QSlider.TicksBelow)
        self.ana_dart_val_lbl = QLabel("25 nM")
        self.ana_dart_val_lbl.setMinimumWidth(55)
        self.ana_dart_slider.valueChanged.connect(
            lambda v: self.ana_dart_val_lbl.setText(f"{v} nM"))
        slw.addWidget(self.ana_dart_slider); slw.addWidget(self.ana_dart_val_lbl)
        grid.addWidget(self.ana_slider_widget, 4, 0, 1, 2)

        run_btn = QPushButton("Run Simulation")
        run_btn.setFixedHeight(36); run_btn.clicked.connect(self._on_run_analog)
        grid.addWidget(run_btn, 5, 0, 1, 2)
        layout.addWidget(pg)

        # Results row
        results = QHBoxLayout()

        # Left: tables
        left = QWidget(); ll = QVBoxLayout(left)

        mg = QGroupBox("Metrics")
        mgl = QVBoxLayout(mg)
        self.ana_metrics_table = QTableWidget()
        self.ana_metrics_table.setColumnCount(2)
        self.ana_metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.ana_metrics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ana_metrics_table.setMaximumHeight(190)
        mgl.addWidget(self.ana_metrics_table)
        ll.addWidget(mg)

        lg = QGroupBox("Ligand Concentration Series")
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
        self.ana_plot_label = QLabel("Plot will appear after simulation.")
        self.ana_plot_label.setAlignment(Qt.AlignCenter)
        self.ana_plot_label.setMinimumSize(520, 380)
        self.ana_plot_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pgl.addWidget(self.ana_plot_label)
        results.addWidget(plot_g, 2)

        layout.addLayout(results)
        return w

    def _ana_mode_changed(self):
        self.ana_slider_widget.setVisible(self.ana_single_btn.isChecked())

    def _on_run_analog(self):
        kd     = self.ana_kd.value()
        k_txn  = self.ana_ktxn.value()
        sweep  = self.ana_sweep_btn.isChecked()

        if sweep:
            # Sweep mode — no single-curve metrics to show
            self.ana_metrics_table.setRowCount(0)
            px = make_analog_sweep_plot(kd, k_txn)
        else:
            dart = float(self.ana_dart_slider.value())
            x, y, metrics = simulate_analog_curve(dart, kd, k_txn)

            rows = [
                ("dART concentration", f"{dart:.0f} nM"),
                ("Kd", f"{kd:.3f} nM"),
                ("k_txn", f"{k_txn:.4f}"),
                ("RNase H", "6 U/mL recommended"),
                ("Detection range (L10–L90)",
                 f"{metrics['L10']:.2f} – {metrics['L90']:.2f} nM"
                 if metrics["L10"] is not None and metrics["L90"] is not None else "N/A"),
                ("EC50",
                 f"{metrics['inflection']:.2f} nM" if metrics["inflection"] else "N/A"),
                ("Sensitivity",
                 f"{metrics['slope']:.3f}" if metrics["slope"] else "N/A"),
                ("Max signal", f"{metrics['maxY']:.2f}"),
            ]
            self.ana_metrics_table.setRowCount(len(rows))
            for i, (k, v) in enumerate(rows):
                self.ana_metrics_table.setItem(i, 0, QTableWidgetItem(k))
                self.ana_metrics_table.setItem(i, 1, QTableWidgetItem(v))

            # ligand series
            ligs = lig_series(kd)
            self.ana_lig_table.setRowCount(len(ligs))
            for i, row in enumerate(ligs):
                self.ana_lig_table.setItem(i, 0, QTableWidgetItem(row["label"]))
                self.ana_lig_table.setItem(i, 1, QTableWidgetItem(str(row["conc"])))

            px = make_analog_plot(x, y, metrics, dart)

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

        # Curve mode
        grid.addWidget(QLabel("Curve mode:"), 5, 0)
        mode_row = QWidget(); mrl = QHBoxLayout(mode_row); mrl.setContentsMargins(0,0,0,0)
        self.dig_single_btn = QRadioButton("Single Ref dART concentration")
        self.dig_sweep_btn  = QRadioButton("Sweep Ref dART 10–100 nM")
        self.dig_single_btn.setChecked(True)
        self.dig_single_btn.toggled.connect(self._dig_mode_changed)
        mrl.addWidget(self.dig_single_btn); mrl.addWidget(self.dig_sweep_btn)
        grid.addWidget(mode_row, 5, 1)

        # Slider
        self.dig_slider_widget = QWidget()
        slw = QHBoxLayout(self.dig_slider_widget); slw.setContentsMargins(0,0,0,0)
        slw.addWidget(QLabel("Ref dART:"))
        self.dig_ref_slider = QSlider(Qt.Horizontal)
        self.dig_ref_slider.setRange(10, 100); self.dig_ref_slider.setSingleStep(5)
        self.dig_ref_slider.setPageStep(10); self.dig_ref_slider.setValue(25)
        self.dig_ref_slider.setTickInterval(10); self.dig_ref_slider.setTickPosition(QSlider.TicksBelow)
        self.dig_ref_val_lbl = QLabel("25 nM")
        self.dig_ref_val_lbl.setMinimumWidth(55)
        self.dig_ref_slider.valueChanged.connect(
            lambda v: self.dig_ref_val_lbl.setText(f"{v} nM"))
        slw.addWidget(self.dig_ref_slider); slw.addWidget(self.dig_ref_val_lbl)
        grid.addWidget(self.dig_slider_widget, 6, 0, 1, 2)

        run_btn = QPushButton("Run Simulation")
        run_btn.setFixedHeight(36); run_btn.clicked.connect(self._on_run_digital)
        grid.addWidget(run_btn, 7, 0, 1, 2)
        layout.addWidget(pg)

        results = QHBoxLayout()

        left = QWidget(); ll = QVBoxLayout(left)

        mg = QGroupBox("Metrics")
        mgl = QVBoxLayout(mg)
        self.dig_metrics_table = QTableWidget()
        self.dig_metrics_table.setColumnCount(2)
        self.dig_metrics_table.setHorizontalHeaderLabels(["Metric", "Value"])
        self.dig_metrics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.dig_metrics_table.setMaximumHeight(160)
        mgl.addWidget(self.dig_metrics_table)
        ll.addWidget(mg)

        lg = QGroupBox("Threshold-Centered Titration Series")
        lgl = QVBoxLayout(lg)
        self.dig_lig_table = QTableWidget()
        self.dig_lig_table.setColumnCount(2)
        self.dig_lig_table.setHorizontalHeaderLabels(["Series point", "Concentration (nM)"])
        self.dig_lig_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        lgl.addWidget(self.dig_lig_table)
        ll.addWidget(lg)

        dsg = QGroupBox("Strands to Order (Reference & Inverter dART)")
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
        self.dig_plot_label = QLabel("Plot will appear after simulation.")
        self.dig_plot_label.setAlignment(Qt.AlignCenter)
        self.dig_plot_label.setMinimumSize(520, 380)
        self.dig_plot_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        pgl.addWidget(self.dig_plot_label)
        results.addWidget(plot_g, 2)

        layout.addLayout(results)
        return w

    def _dig_mode_changed(self):
        self.dig_slider_widget.setVisible(self.dig_single_btn.isChecked())

    def _on_run_digital(self):
        kd        = self.dig_kd.value()
        k_txn_ref = self.dig_ktxn_ref.value()
        k_txn_apt = self.dig_ktxn_apt.value()
        apt_dart  = self.dig_apt.value()
        sweep     = self.dig_sweep_btn.isChecked()

        if sweep:
            self.dig_metrics_table.setRowCount(0)
            px = make_digital_sweep_plot(kd, k_txn_ref, k_txn_apt, apt_dart)
        else:
            ref = float(self.dig_ref_slider.value())
            x, y, metrics = simulate_digital_curve(ref, apt_dart, kd, k_txn_ref, k_txn_apt)
            threshold = metrics["threshold"]

            rows = [
                ("Kd", f"{kd:.3f} nM"),
                ("Reference dART", f"{ref:.0f} nM"),
                ("Aptamer dART", f"{apt_dart:.0f} nM"),
                ("k_txn ref", f"{self.dig_ktxn_ref.value():.4f}"),
                ("k_txn inverter", f"{self.dig_ktxn_apt.value():.4f}"),
                ("Threshold", f"{threshold:.2f} nM" if threshold else "N/A"),
                ("Max signal", f"{metrics['maxY']:.2f}"),
            ]
            self.dig_metrics_table.setRowCount(len(rows))
            for i, (k, v) in enumerate(rows):
                self.dig_metrics_table.setItem(i, 0, QTableWidgetItem(k))
                self.dig_metrics_table.setItem(i, 1, QTableWidgetItem(v))

            series = digital_titration_series(threshold)
            self.dig_lig_table.setRowCount(len(series))
            for i, row in enumerate(series):
                self.dig_lig_table.setItem(i, 0, QTableWidgetItem(row["label"]))
                self.dig_lig_table.setItem(i, 1, QTableWidgetItem(f"{row['conc']:.3f}"))

            px = make_digital_plot(x, y, metrics, ref)

        # ── Strands to order (always populated, uses design tab result if available)
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
        ax.set_title("Normalised Reacted reporter kinetics", color="black", fontsize=11)
        ax.tick_params(colors="black"); ax.spines[:].set_color("black")
        ax.grid(True, alpha=0.15, color="#7faaa0")
        ax.legend(fontsize=6.5, labelcolor="#c8ede7", ncol=2)
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
        ax2.plot(concs, sim_s, linestyle = "--", color="#00ffa0", lw=2, ms=6, label=f"Fit (Kd={kd_fit:.1f} nM)")
        ax2.set_xlabel("[Ligand] (nM)", color="black", fontsize=10)
        ax2.set_xscale("log")
        ax2.set_ylabel("Txn rate at t="+f"{result['slope_time']:.0f}"+" min (nM/min)", color="black", fontsize=10)
        ax2.set_title("Dose-response curve", color="black", fontsize=11)
        ax2.tick_params(colors="black"); ax2.spines[:].set_color("black")
        ax2.grid(True, alpha=0.15, color="#7faaa0")
        ax2.legend(fontsize=8, labelcolor="#c8ede7",
                   facecolor="#071a17", edgecolor="#ffcc4433")
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


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
