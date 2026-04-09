import type { OligoData, OverlapData, SequenceIssue } from "./types";

// Issue type color map — each type gets a distinct hue for quick visual identification
export const ISSUE_COLORS: Record<string, { fill: string; label: string }> = {
  homopolymer:          { fill: "#d93025", label: "Homopolymer run (4+ same base)" },
  g_quadruplex:         { fill: "#7b1fa2", label: "G-quadruplex (4+ G-runs of 3+ in 30bp)" },
  cross_hybridization:  { fill: "#1565c0", label: "Cross-hybridization (off-target Tm too high)" },
  low_tm:               { fill: "#0097a7", label: "Low Tm (on-target annealing may be weak)" },
};

const DEFAULT_ISSUE_COLOR = "#d93025";

export function getIssueColor(issue: SequenceIssue): string {
  // issue.type may be prefixed with "overlap_" from merged overlap issues
  const rawType = issue.type.replace(/^overlap_/, "");
  return ISSUE_COLORS[rawType]?.fill ?? DEFAULT_ISSUE_COLOR;
}

// Google-style palette
// When showBases is true, use lighter fills so base letters are readable
export function getOligoColor(
  o: OligoData,
  isSelected: boolean,
  showBases: boolean = false
): string {
  if (showBases) {
    // Light pastel so dark base-letter colors are legible
    if (isSelected) return o.strand === "sense" ? "#d2e3fc" : "#fce8cd";
    return o.strand === "sense" ? "#e8f0fe" : "#fef3e0";
  }
  if (isSelected) return o.strand === "sense" ? "#1a73e8" : "#e8710a";
  return o.strand === "sense" ? "#4285f4" : "#fa7b17";
}

export function getOverlapColor(
  o: OverlapData,
  isSelected: boolean
): string {
  if (o.issues.length > 0) return isSelected ? "#d93025" : "#ea4335";
  return isSelected ? "#188038" : "#34a853";
}
