import type { OverlapData, TrackLayout } from "./types";
import { getOverlapColor } from "./colors";

interface Props {
  overlap: OverlapData;
  bp2x: (bp: number) => number;
  layout: TrackLayout;
  showNumbers: boolean;
  isSelected: boolean;
  onSelect: (overlap: OverlapData | null) => void;
}

export default function OverlapRegion({
  overlap,
  bp2x,
  layout,
  showNumbers,
  isSelected,
  onSelect,
}: Props) {
  const x = bp2x(overlap.start);
  const w = Math.max(bp2x(overlap.end) - bp2x(overlap.start), 2);
  const color = getOverlapColor(overlap, isSelected);
  const regionH = layout.ANTISENSE_Y + layout.OLIGO_H - layout.SENSE_Y;

  return (
    <g
      onClick={(e) => {
        e.stopPropagation();
        onSelect(isSelected ? null : overlap);
      }}
      style={{ cursor: "pointer" }}
    >
      <rect
        x={x}
        y={layout.SENSE_Y - 2}
        width={w}
        height={regionH + 4}
        fill={color}
        opacity={isSelected ? 0.25 : 0.08}
      />
      <line
        x1={x}
        y1={layout.SENSE_Y}
        x2={x}
        y2={layout.ANTISENSE_Y + layout.OLIGO_H}
        stroke={color}
        strokeWidth={isSelected ? 1.5 : 1}
        opacity={isSelected ? 0.6 : 0.2}
        strokeDasharray="3,3"
      />
      <line
        x1={x + w}
        y1={layout.SENSE_Y}
        x2={x + w}
        y2={layout.ANTISENSE_Y + layout.OLIGO_H}
        stroke={color}
        strokeWidth={isSelected ? 1.5 : 1}
        opacity={isSelected ? 0.6 : 0.2}
        strokeDasharray="3,3"
      />
      {showNumbers && w > 20 && (
        <>
          <text
            x={x + w / 2}
            y={layout.CENTER_Y - 1}
            textAnchor="middle"
            fill={color}
            fontSize={10}
            fontWeight="bold"
            opacity={0.85}
          >
            {overlap.index}
          </text>
          <text
            x={x + w / 2}
            y={layout.CENTER_Y + 10}
            textAnchor="middle"
            fill={color}
            fontSize={9}
            opacity={0.7}
          >
            {overlap.end - overlap.start}bp
          </text>
        </>
      )}
    </g>
  );
}
