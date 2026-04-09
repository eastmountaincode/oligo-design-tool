export interface OligoData {
  index: number;
  seq: string;
  start: number;
  end: number;
  strand: string;
  length: number;
  is_first: boolean;
  is_last: boolean;
  gc: number;
}

export interface OverlapData {
  index: number;
  seq: string;
  start: number;
  end: number;
  gc: number;
  tm: number | null;
  issues: { kind: string; message: string; severity: string; start?: number; end?: number }[];
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
