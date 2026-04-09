import type { OligoData } from "./types";
import { getOligoColor } from "./colors";
import { BASE_COLORS, complementBase } from "./sequence-utils";

interface Props {
  oligo: OligoData;
  bp2x: (bp: number) => number;
  pxPerBp: number;
  y: number;
  height: number;
  fullSeq: string;
  showBases: boolean;
  showNumbers: boolean;
  isSelected: boolean;
  onSelect: (oligo: OligoData | null) => void;
}

export default function OligoArrow({
  oligo,
  bp2x,
  pxPerBp,
  y,
  height,
  fullSeq,
  showBases,
  showNumbers,
  isSelected,
  onSelect,
}: Props) {
  const x = bp2x(oligo.start);
  const w = Math.max(bp2x(oligo.end) - bp2x(oligo.start), 4);
  const arrowW = Math.min(6, w * 0.15);
  const isSense = oligo.strand === "sense";

  const points = isSense
    ? `${x},${y} ${x + w - arrowW},${y} ${x + w},${y + height / 2} ${x + w - arrowW},${y + height} ${x},${y + height}`
    : `${x},${y + height / 2} ${x + arrowW},${y} ${x + w},${y} ${x + w},${y + height} ${x + arrowW},${y + height}`;

  return (
    <g
      onClick={(e) => {
        e.stopPropagation();
        onSelect(isSelected ? null : oligo);
      }}
      style={{ cursor: "pointer" }}
    >
      <polygon
        points={points}
        fill={getOligoColor(oligo, isSelected, showBases)}
        stroke={isSelected ? "#202124" : showBases ? "#80868b" : "#ffffff"}
        strokeWidth={isSelected ? 1.5 : 1}
        opacity={isSelected ? 1 : 0.85}
      />
      {showBases
        ? renderBases(oligo, bp2x, pxPerBp, y, height, fullSeq)
        : showNumbers &&
          w > 20 && (
            <text
              x={x + w / 2}
              y={y + height / 2 + 4}
              fill="#202124"
              fontSize={11}
              fontWeight="bold"
              textAnchor="middle"
              style={{ pointerEvents: "none" }}
            >
              {oligo.index}
            </text>
          )}
    </g>
  );
}

function renderBases(
  oligo: OligoData,
  bp2x: (bp: number) => number,
  pxPerBp: number,
  y: number,
  height: number,
  fullSeq: string
) {
  const isAntisense = oligo.strand === "antisense";
  const bases = [];
  for (let bp = oligo.start; bp < oligo.end && bp < fullSeq.length; bp++) {
    const senseBase = fullSeq[bp] || "N";
    const displayBase = isAntisense ? complementBase(senseBase) : senseBase;
    bases.push(
      <text
        key={bp}
        x={bp2x(bp) + pxPerBp / 2}
        y={y + height / 2 + 4}
        textAnchor="middle"
        fill={BASE_COLORS[displayBase] || "#202124"}
        fontSize={Math.min(pxPerBp * 0.85, 12)}
        fontFamily="monospace"
        fontWeight="bold"
        style={{ pointerEvents: "none" }}
      >
        {displayBase}
      </text>
    );
  }
  return <>{bases}</>;
}
