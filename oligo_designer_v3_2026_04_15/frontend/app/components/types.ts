// Types shared between the codon optimization and oligo designer panes.

export interface CodonTableEntry {
  amino_acid: string;
  codon: string;
  per_thousand: number;
}

/** Everything the backend's cross-hyb repair pass needs in order to propose
 *  minimally-invasive codon swaps that resolve oligo cross-hybridization.
 *
 *  This is built on the codon optimization pane immediately after a
 *  successful optimization completes and is cleared when the optimization
 *  becomes stale (inputs change). The oligo designer pane forwards it
 *  verbatim to /api/design; if it's missing, the repair pass is skipped.
 */
export interface RepairContext {
  protein: string;
  codon_table: CodonTableEntry[];
  min_codon_frequency: number;
  gc_min: number;         // fraction (0-1)
  gc_max: number;
  gc_window: number;
  avoid_homopolymers_gc: number;
  avoid_homopolymers_at: number;
  avoid_patterns: string[];
}

/** One row of the repair log returned by /api/design.  Kind "applied" means
 *  a codon swap was performed; "no_move_found" / "skipped" describe why a
 *  repair step did NOT change anything.
 */
export interface RepairLogEntry {
  iteration: number;
  kind: "applied" | "no_move_found" | "skipped";
  codon_index?: number | null;
  insert_position?: number | null;
  amino_acid?: string | null;
  old_codon?: string | null;
  new_codon?: string | null;
  old_freq_per_k?: number | null;
  new_freq_per_k?: number | null;
  gc_before?: number | null;
  gc_after?: number | null;
  fixed_pairs: number[][];
  remaining_pairs: number[][];
  remaining_max_tm?: number | null;
  threshold_tm?: number | null;
  notes: string;
}

export interface RepairReport {
  ran: boolean;
  threshold_tm?: number | null;
  issues_before: number;
  issues_after: number;
  log: RepairLogEntry[];
  modified_dna?: string | null;
  skipped_reason?: string | null;
}
