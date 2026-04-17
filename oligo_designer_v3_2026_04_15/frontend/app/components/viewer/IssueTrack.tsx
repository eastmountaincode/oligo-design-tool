import type { SequenceIssue } from "./types";
import { getIssueColor } from "./colors";

interface Props {
  issues: SequenceIssue[];
  bp2x: (bp: number) => number;
  y: number;
  selectedIndex: number | null;
  onSelect: (index: number | null) => void;
}

const ROW_HEIGHT = 10;
const ROW_GAP = 2;

/** Assign each issue to a row so overlapping issues stack vertically. */
function assignRows(issues: { start?: number; end?: number }[]): number[] {
  const rows: number[] = [];
  const rowEnds: number[] = []; // tracks the rightmost end in each row
  for (const issue of issues) {
    if (issue.start === undefined || issue.end === undefined) {
      rows.push(0);
      continue;
    }
    let placed = false;
    for (let r = 0; r < rowEnds.length; r++) {
      if (issue.start >= rowEnds[r]) {
        rows.push(r);
        rowEnds[r] = issue.end;
        placed = true;
        break;
      }
    }
    if (!placed) {
      rows.push(rowEnds.length);
      rowEnds.push(issue.end);
    }
  }
  return rows;
}

export function getIssueTrackHeight(issues: SequenceIssue[]): number {
  const regions = issues.filter((i) => i.start !== undefined && i.end !== undefined);
  if (regions.length === 0) return 0;
  const rows = assignRows(regions);
  const numRows = Math.max(...rows) + 1;
  return numRows * (ROW_HEIGHT + ROW_GAP) + 12; // 12 for the label
}

export default function IssueTrack({ issues, bp2x, y, selectedIndex, onSelect }: Props) {
  const regions = issues.filter((i) => i.start !== undefined && i.end !== undefined);
  if (regions.length === 0) return null;

  const rows = assignRows(regions);

  return (
    <g>
      {regions.map((issue, i) => {
        const isSelected = selectedIndex === i;
        const row = rows[i];
        const ry = y + 12 + row * (ROW_HEIGHT + ROW_GAP);
        return (
          <rect
            key={`issue-${i}`}
            x={bp2x(issue.start!)}
            y={ry}
            width={Math.max(bp2x(issue.end!) - bp2x(issue.start!), 3)}
            height={ROW_HEIGHT}
            fill={getIssueColor(issue)}
            opacity={isSelected ? 1 : 0.85}
            stroke={isSelected ? "#202124" : "none"}
            strokeWidth={isSelected ? 1.5 : 0}
            style={{ cursor: "pointer" }}
            onClick={(e) => {
              e.stopPropagation();
              onSelect(isSelected ? null : i);
            }}
          />
        );
      })}
    </g>
  );
}
