#!/usr/bin/env python3
"""
Oligo Sequence Deconstruction Tool for qPCR Design
===================================================
A professional tool for analyzing and consolidating oligonucleotide sequences
with support for ambiguity codes and G-A wobble base pairing.

Author: Claude (Anthropic)
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from typing import Optional

# Import the analysis engine
# If oligo_engine.py is not found, fall back to inline definitions
try:
    from oligo_engine import (
        OligoAnalyzer, AnalysisResult, QualityReport,
        format_sequence, format_results_table
    )
except ImportError:
    # Inline fallback - full engine code included below
    from collections import Counter
    from typing import Dict, List, Tuple, Set, Callable
    from dataclasses import dataclass, field

    IUPAC_CODES = {
        'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
        'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'},
        'W': {'A', 'T'}, 'K': {'G', 'T'}, 'M': {'A', 'C'},
        'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'},
    }
    BASES_TO_IUPAC = {frozenset(v): k for k, v in IUPAC_CODES.items()}
    AMBIGUOUS_BASES = set('RYSWKMBDHVN')
    VALID_BASES = set('ACGTRYSWKMBDHVN')
    GAP_CHARS = set('-.')

    @dataclass
    class AnalysisResult:
        variants: List[Tuple[str, int, float]] = field(default_factory=list)
        total_sequences: int = 0
        coverage: float = 0.0
        uncovered_count: int = 0
        uncovered_percentage: float = 0.0
        ambiguity_count: int = 0
        message: str = ""

    @dataclass
    class QualityReport:
        original_count: int = 0
        valid_count: int = 0
        removed_ambiguous: int = 0
        removed_gaps: int = 0
        removed_wrong_length: int = 0
        majority_length: int = 0
        valid_sequences: List[str] = field(default_factory=list)

    def format_sequence(seq: str, group_size: int = 3) -> str:
        return ' '.join(seq[i:i+group_size] for i in range(0, len(seq), group_size))

    class OligoAnalyzer:
        def __init__(self):
            self.sequences: List[str] = []
            self.quality_report: Optional[QualityReport] = None

        def parse_fasta(self, text: str) -> List[str]:
            sequences = []
            current_seq = []
            lines = text.strip().split('\n')
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(''.join(current_seq).upper().replace(' ', ''))
                        current_seq = []
                else:
                    current_seq.append(line.upper().replace(' ', ''))
            if current_seq:
                sequences.append(''.join(current_seq).upper().replace(' ', ''))
            if not sequences and lines:
                for line in lines:
                    line = line.strip().upper().replace(' ', '')
                    if line and not line.startswith('>'):
                        sequences.append(line)
            return sequences

        def quality_filter(self, sequences: List[str]) -> QualityReport:
            report = QualityReport()
            report.original_count = len(sequences)
            if not sequences:
                return report
            length_counts = Counter(len(seq) for seq in sequences)
            report.majority_length = length_counts.most_common(1)[0][0]
            valid = []
            for seq in sequences:
                if any(c in GAP_CHARS for c in seq):
                    report.removed_gaps += 1
                    continue
                if any(c in AMBIGUOUS_BASES for c in seq):
                    report.removed_ambiguous += 1
                    continue
                if not all(c in VALID_BASES for c in seq):
                    report.removed_ambiguous += 1
                    continue
                if len(seq) != report.majority_length:
                    report.removed_wrong_length += 1
                    continue
                valid.append(seq)
            report.valid_sequences = valid
            report.valid_count = len(valid)
            return report

        def load_sequences(self, text: str) -> QualityReport:
            raw_sequences = self.parse_fasta(text)
            self.quality_report = self.quality_filter(raw_sequences)
            self.sequences = self.quality_report.valid_sequences
            return self.quality_report

        def count_variants(self) -> AnalysisResult:
            result = AnalysisResult()
            result.total_sequences = len(self.sequences)
            if not self.sequences:
                result.message = "No sequences to analyze"
                return result
            counts = Counter(self.sequences)
            total = len(self.sequences)
            for seq, count in counts.most_common():
                result.variants.append((seq, count, (count / total) * 100))
            result.coverage = 100.0
            return result

        def _get_ambiguity_code(self, bases: Set[str], treat_g_as_a: bool = False) -> str:
            if treat_g_as_a:
                if bases == {'A', 'G'} or bases == {'A'} or bases == {'G'}:
                    return 'A'
            return BASES_TO_IUPAC.get(frozenset(bases), 'N')

        def _sequence_matches_consensus(self, seq: str, consensus: str, treat_g_as_a: bool = False) -> bool:
            if len(seq) != len(consensus):
                return False
            for s, c in zip(seq, consensus):
                if c in IUPAC_CODES:
                    allowed = IUPAC_CODES[c]
                    if treat_g_as_a:
                        if 'A' in allowed:
                            allowed = allowed | {'G'}
                        if 'G' in allowed:
                            allowed = allowed | {'A'}
                    if s not in allowed:
                        return False
                elif s != c:
                    if not (treat_g_as_a and {s, c} == {'A', 'G'}):
                        return False
            return True

        def _create_consensus(self, sequences: List[str], treat_g_as_a: bool = False) -> Tuple[str, int]:
            if not sequences:
                return "", 0
            seq_len = len(sequences[0])
            consensus = []
            ambiguity_count = 0
            for pos in range(seq_len):
                bases_at_pos = set(seq[pos] for seq in sequences)
                if treat_g_as_a:
                    if bases_at_pos == {'A', 'G'} or bases_at_pos == {'A'} or bases_at_pos == {'G'}:
                        consensus.append('A')
                        continue
                    if 'A' in bases_at_pos and 'G' in bases_at_pos:
                        bases_at_pos = (bases_at_pos - {'G'}) | {'A'}
                if len(bases_at_pos) == 1:
                    consensus.append(bases_at_pos.pop())
                else:
                    code = self._get_ambiguity_code(bases_at_pos, treat_g_as_a)
                    consensus.append(code)
                    ambiguity_count += 1
            return ''.join(consensus), ambiguity_count

        def find_minimum_variants_greedy(self, max_ambiguities: int, treat_g_as_a: bool = False,
                                         progress_callback: Optional[Callable[[str], None]] = None) -> AnalysisResult:
            result = AnalysisResult()
            result.total_sequences = len(self.sequences)
            if not self.sequences:
                result.message = "No sequences to analyze"
                return result
            seq_counts = Counter(self.sequences)
            unique_seqs = list(seq_counts.keys())
            total = len(self.sequences)
            uncovered = set(unique_seqs)
            variants = []
            while uncovered:
                if progress_callback:
                    progress_callback(f"Finding variant {len(variants) + 1}... ({len(uncovered)} remaining)")
                best_consensus = None
                best_coverage = set()
                best_ambiguities = float('inf')
                uncovered_counts = sorted([(seq, seq_counts[seq]) for seq in uncovered], key=lambda x: -x[1])
                for seed_seq, _ in uncovered_counts[:min(50, len(uncovered_counts))]:
                    potential_group = [seed_seq]
                    for other_seq in uncovered:
                        if other_seq == seed_seq:
                            continue
                        combined = potential_group + [other_seq]
                        consensus, amb_count = self._create_consensus(combined, treat_g_as_a)
                        if amb_count <= max_ambiguities:
                            potential_group.append(other_seq)
                    consensus, amb_count = self._create_consensus(potential_group, treat_g_as_a)
                    coverage = set()
                    for seq in uncovered:
                        if self._sequence_matches_consensus(seq, consensus, treat_g_as_a):
                            coverage.add(seq)
                    coverage_score = sum(seq_counts[s] for s in coverage)
                    if coverage_score > sum(seq_counts[s] for s in best_coverage) or \
                       (coverage_score == sum(seq_counts[s] for s in best_coverage) and amb_count < best_ambiguities):
                        best_consensus = consensus
                        best_coverage = coverage
                        best_ambiguities = amb_count
                if best_consensus is None or not best_coverage:
                    most_freq = max(uncovered, key=lambda x: seq_counts[x])
                    best_consensus = most_freq
                    best_coverage = {most_freq}
                    best_ambiguities = 0
                count = sum(seq_counts[s] for s in best_coverage)
                variants.append((best_consensus, count, (count / total) * 100))
                result.ambiguity_count = max(result.ambiguity_count, best_ambiguities)
                uncovered -= best_coverage
            result.variants = variants
            result.coverage = 100.0
            return result

        def find_top_n_variants(self, n: int, treat_g_as_a: bool = False) -> AnalysisResult:
            result = AnalysisResult()
            result.total_sequences = len(self.sequences)
            if not self.sequences:
                result.message = "No sequences to analyze"
                return result
            counts = Counter(self.sequences)
            total = len(self.sequences)
            top_variants = counts.most_common(n)
            covered = 0
            for seq, count in top_variants:
                result.variants.append((seq, count, (count / total) * 100))
                covered += count
            result.coverage = (covered / total) * 100
            result.uncovered_count = total - covered
            result.uncovered_percentage = 100 - result.coverage
            return result


# ============================================================================
# GUI APPLICATION
# ============================================================================
class OligoAnalyzerGUI:
    """Professional Tkinter GUI for the Oligo Analyzer"""
    
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Oligo Sequence Deconstruction Tool for qPCR Design")
        self.root.geometry("1200x800")
        self.root.minsize(900, 600)
        
        # Initialize analyzer
        self.analyzer = OligoAnalyzer()
        
        # Configure styles
        self.setup_styles()
        
        # Create main layout
        self.create_widgets()
        
        # Bind keyboard shortcuts
        self.root.bind('<Control-o>', lambda e: self.load_file())
        self.root.bind('<Control-r>', lambda e: self.run_analysis())
    
    def setup_styles(self):
        """Configure ttk styles for professional appearance"""
        style = ttk.Style()
        style.theme_use('clam')
        
        # Configure colors
        style.configure('TFrame', background='#f5f5f5')
        style.configure('TLabel', background='#f5f5f5', font=('Segoe UI', 10))
        style.configure('TButton', font=('Segoe UI', 10), padding=5)
        style.configure('Header.TLabel', font=('Segoe UI', 12, 'bold'))
        style.configure('Title.TLabel', font=('Segoe UI', 14, 'bold'), foreground='#2c3e50')
        style.configure('Status.TLabel', font=('Segoe UI', 9), foreground='#666666')
        style.configure('Success.TLabel', foreground='#27ae60')
        style.configure('Warning.TLabel', foreground='#e67e22')
        style.configure('Error.TLabel', foreground='#e74c3c')
    
    def create_widgets(self):
        """Create all GUI widgets"""
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky='nsew')
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
        # Title
        title_frame = ttk.Frame(main_frame)
        title_frame.grid(row=0, column=0, sticky='ew', pady=(0, 10))
        ttk.Label(title_frame, text="Oligo Sequence Deconstruction Tool",
                  style='Title.TLabel').pack(side='left')
        
        # Create notebook for tabs
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.grid(row=1, column=0, sticky='nsew')
        
        # Input tab
        self.create_input_tab()
        
        # Analysis tab
        self.create_analysis_tab()
        
        # Results tab
        self.create_results_tab()
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready. Load sequences to begin.")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, style='Status.TLabel')
        status_bar.grid(row=2, column=0, sticky='ew', pady=(10, 0))
    
    def create_input_tab(self):
        """Create the input/data tab"""
        input_frame = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(input_frame, text=" Input Data ")
        
        input_frame.columnconfigure(0, weight=1)
        input_frame.rowconfigure(1, weight=1)
        
        # Control buttons frame
        btn_frame = ttk.Frame(input_frame)
        btn_frame.grid(row=0, column=0, sticky='ew', pady=(0, 10))
        
        ttk.Button(btn_frame, text="📂 Load FASTA File",
                   command=self.load_file).pack(side='left', padx=(0, 5))
        ttk.Button(btn_frame, text="📋 Paste from Clipboard",
                   command=self.paste_clipboard).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="🗑️ Clear",
                   command=self.clear_input).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="📊 Load Example",
                   command=self.load_example).pack(side='left', padx=5)
        
        # Input text area with label
        input_label_frame = ttk.Frame(input_frame)
        input_label_frame.grid(row=1, column=0, sticky='nsew')
        input_label_frame.columnconfigure(0, weight=1)
        input_label_frame.rowconfigure(1, weight=1)
        
        ttk.Label(input_label_frame, text="Sequence Data (FASTA format):",
                  style='Header.TLabel').grid(row=0, column=0, sticky='w', pady=(0, 5))
        
        self.input_text = scrolledtext.ScrolledText(
            input_label_frame, wrap=tk.NONE, font=('Consolas', 10),
            height=20, width=80
        )
        self.input_text.grid(row=1, column=0, sticky='nsew')
        
        # Add horizontal scrollbar
        h_scroll = ttk.Scrollbar(input_label_frame, orient='horizontal',
                                  command=self.input_text.xview)
        h_scroll.grid(row=2, column=0, sticky='ew')
        self.input_text.configure(xscrollcommand=h_scroll.set)
        
        # Quality report frame
        quality_frame = ttk.LabelFrame(input_frame, text="Quality Report", padding="10")
        quality_frame.grid(row=2, column=0, sticky='ew', pady=(10, 0))
        
        self.quality_text = tk.Text(quality_frame, height=6, font=('Consolas', 9),
                                     state='disabled', bg='#f8f8f8')
        self.quality_text.pack(fill='x')
    
    def create_analysis_tab(self):
        """Create the analysis options tab"""
        analysis_frame = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(analysis_frame, text=" Analysis Options ")
        
        analysis_frame.columnconfigure(1, weight=1)
        
        # Analysis type selection
        type_frame = ttk.LabelFrame(analysis_frame, text="Analysis Type", padding="10")
        type_frame.grid(row=0, column=0, columnspan=2, sticky='ew', pady=(0, 10))
        type_frame.columnconfigure(1, weight=1)
        
        self.analysis_type = tk.StringVar(value="all_variants")
        
        types = [
            ("all_variants", "Find all unique sequence variants (no ambiguities)"),
            ("min_variants", "Find minimum variants with up to N ambiguities"),
            ("top_n", "Find top N most frequent variants"),
        ]
        
        for i, (value, text) in enumerate(types):
            rb = ttk.Radiobutton(type_frame, text=text, variable=self.analysis_type,
                                  value=value, command=self.update_options_state)
            rb.grid(row=i, column=0, columnspan=2, sticky='w', pady=2)
        
        # Options frame
        options_frame = ttk.LabelFrame(analysis_frame, text="Options", padding="10")
        options_frame.grid(row=1, column=0, columnspan=2, sticky='ew', pady=(0, 10))
        options_frame.columnconfigure(1, weight=1)
        
        # Max ambiguities
        ttk.Label(options_frame, text="Maximum ambiguities (N):").grid(
            row=0, column=0, sticky='w', pady=5)
        self.max_ambiguities_var = tk.StringVar(value="1")
        self.max_ambiguities_spin = ttk.Spinbox(
            options_frame, from_=0, to=50, width=10,
            textvariable=self.max_ambiguities_var
        )
        self.max_ambiguities_spin.grid(row=0, column=1, sticky='w', pady=5)
        
        # Top N variants
        ttk.Label(options_frame, text="Number of top variants:").grid(
            row=1, column=0, sticky='w', pady=5)
        self.top_n_var = tk.StringVar(value="5")
        self.top_n_spin = ttk.Spinbox(
            options_frame, from_=1, to=100, width=10,
            textvariable=self.top_n_var
        )
        self.top_n_spin.grid(row=1, column=1, sticky='w', pady=5)
        
        # G-A wobble option
        self.treat_g_as_a_var = tk.BooleanVar(value=False)
        self.treat_g_as_a_check = ttk.Checkbutton(
            options_frame,
            text="Treat G as A (G-T wobble base pairing is not considered a mismatch)",
            variable=self.treat_g_as_a_var
        )
        self.treat_g_as_a_check.grid(row=2, column=0, columnspan=2, sticky='w', pady=5)
        
        # Additional options (expandable)
        extra_frame = ttk.LabelFrame(analysis_frame, text="Additional Options", padding="10")
        extra_frame.grid(row=2, column=0, columnspan=2, sticky='ew', pady=(0, 10))
        extra_frame.columnconfigure(1, weight=1)
        
        # Format spacer option
        self.add_spacer_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            extra_frame,
            text="Add spaces every 3 bases in output (visual grouping)",
            variable=self.add_spacer_var
        ).grid(row=0, column=0, columnspan=2, sticky='w', pady=2)
        
        # Show percentages
        self.show_percentages_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            extra_frame,
            text="Show percentages in output",
            variable=self.show_percentages_var
        ).grid(row=1, column=0, columnspan=2, sticky='w', pady=2)
        
        # Run button
        run_frame = ttk.Frame(analysis_frame)
        run_frame.grid(row=3, column=0, columnspan=2, sticky='ew', pady=10)
        
        self.run_btn = ttk.Button(run_frame, text="▶ Run Analysis",
                                   command=self.run_analysis)
        self.run_btn.pack(side='left', padx=(0, 10))
        
        self.progress_var = tk.StringVar(value="")
        ttk.Label(run_frame, textvariable=self.progress_var,
                  style='Status.TLabel').pack(side='left')
        
        # Update initial state
        self.update_options_state()
    
    def create_results_tab(self):
        """Create the results display tab"""
        results_frame = ttk.Frame(self.notebook, padding="10")
        self.notebook.add(results_frame, text=" Results ")
        
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(1, weight=1)
        
        # Control buttons
        btn_frame = ttk.Frame(results_frame)
        btn_frame.grid(row=0, column=0, sticky='ew', pady=(0, 10))
        
        ttk.Button(btn_frame, text="📋 Copy Results",
                   command=self.copy_results).pack(side='left', padx=(0, 5))
        ttk.Button(btn_frame, text="💾 Export to File",
                   command=self.export_results).pack(side='left', padx=5)
        ttk.Button(btn_frame, text="🗑️ Clear Results",
                   command=self.clear_results).pack(side='left', padx=5)
        
        # Results text area
        results_label_frame = ttk.Frame(results_frame)
        results_label_frame.grid(row=1, column=0, sticky='nsew')
        results_label_frame.columnconfigure(0, weight=1)
        results_label_frame.rowconfigure(1, weight=1)
        
        ttk.Label(results_label_frame, text="Analysis Results:",
                  style='Header.TLabel').grid(row=0, column=0, sticky='w', pady=(0, 5))
        
        self.results_text = scrolledtext.ScrolledText(
            results_label_frame, wrap=tk.NONE, font=('Consolas', 10),
            height=25, width=100
        )
        self.results_text.grid(row=1, column=0, sticky='nsew')
        
        # Add horizontal scrollbar
        h_scroll = ttk.Scrollbar(results_label_frame, orient='horizontal',
                                  command=self.results_text.xview)
        h_scroll.grid(row=2, column=0, sticky='ew')
        self.results_text.configure(xscrollcommand=h_scroll.set)
        
        # Summary frame
        summary_frame = ttk.LabelFrame(results_frame, text="Summary", padding="10")
        summary_frame.grid(row=2, column=0, sticky='ew', pady=(10, 0))
        
        self.summary_text = tk.Text(summary_frame, height=4, font=('Consolas', 9),
                                     state='disabled', bg='#f8f8f8')
        self.summary_text.pack(fill='x')
    
    def update_options_state(self):
        """Update the state of options based on analysis type"""
        analysis_type = self.analysis_type.get()
        
        # Enable/disable based on analysis type
        if analysis_type == "all_variants":
            self.max_ambiguities_spin.configure(state='disabled')
            self.top_n_spin.configure(state='disabled')
        elif analysis_type == "min_variants":
            self.max_ambiguities_spin.configure(state='normal')
            self.top_n_spin.configure(state='disabled')
        elif analysis_type == "top_n":
            self.max_ambiguities_spin.configure(state='disabled')
            self.top_n_spin.configure(state='normal')
    
    def load_file(self):
        """Load sequences from a FASTA file"""
        filename = filedialog.askopenfilename(
            title="Open FASTA File",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.fas"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
        )
        
        if filename:
            try:
                with open(filename, 'r') as f:
                    content = f.read()
                self.input_text.delete('1.0', tk.END)
                self.input_text.insert('1.0', content)
                self.process_input()
                self.status_var.set(f"Loaded: {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file:\n{str(e)}")
    
    def paste_clipboard(self):
        """Paste sequences from clipboard"""
        try:
            content = self.root.clipboard_get()
            self.input_text.delete('1.0', tk.END)
            self.input_text.insert('1.0', content)
            self.process_input()
            self.status_var.set("Pasted from clipboard")
        except tk.TclError:
            messagebox.showwarning("Warning", "Clipboard is empty or contains invalid data")
    
    def clear_input(self):
        """Clear input text area"""
        self.input_text.delete('1.0', tk.END)
        self.update_quality_report(None)
        self.status_var.set("Input cleared")
    
    def load_example(self):
        """Load example data for testing"""
        example = """>Sequence1
AAAAAAAAA
>Sequence2
ATAAAAAAA
>Sequence3
ATAAAAAAA
>Sequence4
AAAAAAAAA
>Sequence5
AAAAAAAAA
>Sequence6
AAAAGAAAA"""
        
        self.input_text.delete('1.0', tk.END)
        self.input_text.insert('1.0', example)
        self.process_input()
        self.status_var.set("Example data loaded")
    
    def process_input(self):
        """Process the input text and update quality report"""
        text = self.input_text.get('1.0', tk.END)
        report = self.analyzer.load_sequences(text)
        self.update_quality_report(report)
    
    def update_quality_report(self, report: Optional[QualityReport]):
        """Update the quality report display"""
        self.quality_text.configure(state='normal')
        self.quality_text.delete('1.0', tk.END)
        
        if report is None:
            self.quality_text.insert('1.0', "No data loaded")
        else:
            lines = [
                f"Original sequences:        {report.original_count:,}",
                f"Valid sequences:           {report.valid_count:,}",
                f"Majority sequence length:  {report.majority_length} bp",
                "",
                "Removed due to quality issues:",
                f"  - Ambiguous bases:       {report.removed_ambiguous:,}",
                f"  - Gaps (- or .):         {report.removed_gaps:,}",
                f"  - Wrong length:          {report.removed_wrong_length:,}",
                f"  - Total removed:         {report.removed_ambiguous + report.removed_gaps + report.removed_wrong_length:,}",
            ]
            self.quality_text.insert('1.0', '\n'.join(lines))
        
        self.quality_text.configure(state='disabled')
    
    def format_sequence(self, seq: str) -> str:
        """Format sequence with optional spacing"""
        if self.add_spacer_var.get():
            return ' '.join(seq[i:i+3] for i in range(0, len(seq), 3))
        return seq
    
    def run_analysis(self):
        """Run the selected analysis"""
        if not self.analyzer.sequences:
            self.process_input()
            if not self.analyzer.sequences:
                messagebox.showwarning("Warning", "No valid sequences to analyze.\nPlease load sequence data first.")
                return
        
        analysis_type = self.analysis_type.get()
        treat_g_as_a = self.treat_g_as_a_var.get()
        
        # Disable run button during analysis
        self.run_btn.configure(state='disabled')
        self.progress_var.set("Analyzing...")
        self.root.update()
        
        try:
            if analysis_type == "all_variants":
                result = self.analyzer.count_variants()
            elif analysis_type == "min_variants":
                max_amb = int(self.max_ambiguities_var.get())
                result = self.analyzer.find_minimum_variants_greedy(
                    max_amb, treat_g_as_a,
                    progress_callback=lambda msg: self._update_progress(msg)
                )
            elif analysis_type == "top_n":
                top_n = int(self.top_n_var.get())
                result = self.analyzer.find_top_n_variants(top_n, treat_g_as_a)
            else:
                result = AnalysisResult()
            
            self.display_results(result, analysis_type)
            self.notebook.select(2)  # Switch to results tab
            self.status_var.set("Analysis complete")
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed:\n{str(e)}")
            self.status_var.set("Analysis failed")
        
        finally:
            self.run_btn.configure(state='normal')
            self.progress_var.set("")
    
    def _update_progress(self, message: str):
        """Update progress message"""
        self.progress_var.set(message)
        self.root.update()
    
    def display_results(self, result: AnalysisResult, analysis_type: str):
        """Display analysis results in a clean, copy-friendly format"""
        self.results_text.delete('1.0', tk.END)
        
        # Header
        lines = []
        lines.append("=" * 80)
        lines.append("OLIGO SEQUENCE ANALYSIS RESULTS")
        lines.append("=" * 80)
        lines.append("")
        
        # Analysis info
        type_names = {
            "all_variants": "All Unique Variants (No Ambiguities)",
            "min_variants": f"Minimum Variants (Max {self.max_ambiguities_var.get()} Ambiguities)",
            "top_n": f"Top {self.top_n_var.get()} Most Frequent Variants",
        }
        lines.append(f"Analysis Type:    {type_names.get(analysis_type, analysis_type)}")
        lines.append(f"Total Sequences:  {result.total_sequences:,}")
        lines.append(f"Variants Found:   {len(result.variants):,}")
        
        if self.treat_g_as_a_var.get():
            lines.append("G-A Wobble:       Enabled (G treated as A)")
        
        lines.append("")
        lines.append("-" * 80)
        lines.append("")
        
        # Results table header (tab-separated for easy Excel paste)
        show_pct = self.show_percentages_var.get()
        if show_pct:
            lines.append("Variant\tSequence\tCount\tPercentage")
            lines.append("-" * 80)
        else:
            lines.append("Variant\tSequence\tCount")
            lines.append("-" * 80)
        
        # Results
        for i, (seq, count, pct) in enumerate(result.variants, 1):
            formatted_seq = self.format_sequence(seq)
            if show_pct:
                lines.append(f"Variant {i}\t{formatted_seq}\t{count:,}\t{pct:.1f}%")
            else:
                lines.append(f"Variant {i}\t{formatted_seq}\t{count:,}")
        
        # Coverage info for top_n analysis
        if analysis_type == "top_n" and result.uncovered_count > 0:
            lines.append("")
            lines.append("-" * 80)
            lines.append(f"Sequences with mismatches: {result.uncovered_count:,} ({result.uncovered_percentage:.1f}%)")
        
        lines.append("")
        lines.append("=" * 80)
        
        self.results_text.insert('1.0', '\n'.join(lines))
        
        # Update summary
        self.update_summary(result, analysis_type)
    
    def update_summary(self, result: AnalysisResult, analysis_type: str):
        """Update the summary display"""
        self.summary_text.configure(state='normal')
        self.summary_text.delete('1.0', tk.END)
        
        lines = [
            f"Total sequences analyzed: {result.total_sequences:,}",
            f"Number of variants: {len(result.variants):,}",
            f"Coverage: {result.coverage:.1f}%",
        ]
        
        if result.uncovered_count > 0:
            lines.append(f"Uncovered sequences: {result.uncovered_count:,} ({result.uncovered_percentage:.1f}%)")
        
        self.summary_text.insert('1.0', '\n'.join(lines))
        self.summary_text.configure(state='disabled')
    
    def copy_results(self):
        """Copy results to clipboard"""
        content = self.results_text.get('1.0', tk.END)
        self.root.clipboard_clear()
        self.root.clipboard_append(content)
        self.status_var.set("Results copied to clipboard")
    
    def export_results(self):
        """Export results to a file"""
        filename = filedialog.asksaveasfilename(
            title="Export Results",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("TSV files", "*.tsv"),
                ("All files", "*.*")
            ]
        )
        
        if filename:
            try:
                content = self.results_text.get('1.0', tk.END)
                with open(filename, 'w') as f:
                    f.write(content)
                self.status_var.set(f"Results exported to: {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to export:\n{str(e)}")
    
    def clear_results(self):
        """Clear results display"""
        self.results_text.delete('1.0', tk.END)
        self.summary_text.configure(state='normal')
        self.summary_text.delete('1.0', tk.END)
        self.summary_text.configure(state='disabled')
        self.status_var.set("Results cleared")


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================
def main():
    """Main entry point for the application"""
    root = tk.Tk()
    
    # Set icon if available (optional)
    try:
        # On Windows, you could set an icon
        # root.iconbitmap('icon.ico')
        pass
    except:
        pass
    
    app = OligoAnalyzerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
