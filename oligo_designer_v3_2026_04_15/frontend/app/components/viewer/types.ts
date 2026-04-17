export interface OligoData {
  index: number;           // position within the owning reaction (0-based)
  seq: string;
  start: number;
  end: number;
  strand: string;
  length: number;
  is_first: boolean;
  is_last: boolean;
  gc: number;
  // In multi-reaction assembly, this is the 0-based reaction index the oligo
  // came from. Undefined when there's only one reaction (legacy view).
  rxnIndex?: number;
}

export interface OverlapData {
  index: number;           // position within the owning reaction (1-based, from the backend)
  seq: string;
  start: number;
  end: number;
  gc: number;
  tm: number | null;
  issues: { kind: string; message: string; severity: string; start?: number; end?: number }[];
  // See OligoData.rxnIndex.
  rxnIndex?: number;
  // True when this is the Gibson handoff between this reaction and the NEXT
  // one (shared last L bp / first L bp). Rendered with a distinct style so
  // the user can spot where one reaction tube ends and the next begins.
  is_junction?: boolean;
}

// Label oligos / overlaps as "<num><letter>" in multi-reaction mode, where
// letter is A for reaction 0, B for reaction 1, etc. Single-reaction mode
// falls back to just the number. Keeps the shorter "1A, 2A, 1B" form the
// user asked for so cross-reaction identity is obvious at a glance.
function rxnLetter(rxnIndex: number | undefined): string {
  if (rxnIndex == null) return "";
  return String.fromCharCode(65 + rxnIndex);
}

export function oligoLabel(oligo: Pick<OligoData, "index" | "rxnIndex">): string {
  return `${oligo.index + 1}${rxnLetter(oligo.rxnIndex)}`;
}

export function overlapLabel(overlap: Pick<OverlapData, "index">): string {
  // Overlaps are numbered globally along the construct (1..N by position)
  // regardless of which reaction they live in — a single unbroken sequence
  // of overlap regions reads more naturally than "5A / 1B". The parent is
  // responsible for rewriting `index` into that global numbering before
  // handing overlaps to the viewer / tables.
  return `${overlap.index}`;
}

export interface SequenceIssue {
  type: string;
  message: string;
  severity: string;
  start?: number;
  end?: number;
}

export interface TrackLayout {
  RULER_Y: number;
  ISSUE_TRACK_Y: number;
  SENSE_Y: number;
  OLIGO_H: number;
  CENTER_Y: number;
  ANTISENSE_Y: number;
  OVERLAP_H: number;
}
