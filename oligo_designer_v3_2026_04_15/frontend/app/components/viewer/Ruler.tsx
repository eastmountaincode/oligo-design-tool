interface Props {
  totalLength: number;
  bp2x: (bp: number) => number;
  containerWidth: number;
  pxPerBp: number;
  y: number;
}

export default function Ruler({ totalLength, bp2x, containerWidth, pxPerBp, y }: Props) {
  const viewBp = containerWidth / pxPerBp;
  const tickInterval =
    viewBp <= 50 ? 5
    : viewBp <= 100 ? 10
    : viewBp <= 300 ? 20
    : viewBp <= 600 ? 50
    : viewBp <= 2000 ? 100
    : viewBp <= 5000 ? 500
    : 1000;

  // Ticks are rendered at 1-indexed biology positions (1, 100, 200, ...,
  // totalLength) to match the 1-indexed convention used throughout the
  // rest of the UI. SVG x-coordinates still come from the 0-indexed
  // internal representation: 1-indexed position P sits at bp2x(P - 1).
  const ticks: number[] = [1];
  for (let p = tickInterval; p < totalLength; p += tickInterval) ticks.push(p);
  if (ticks[ticks.length - 1] !== totalLength) ticks.push(totalLength);

  return (
    <g>
      <line
        x1={bp2x(0)}
        y1={y}
        x2={bp2x(totalLength)}
        y2={y}
        stroke="#dadce0"
        strokeWidth={1}
      />
      {ticks.map((p) => (
        <g key={`tick-${p}`}>
          <line
            x1={bp2x(p - 1)}
            y1={y - 4}
            x2={bp2x(p - 1)}
            y2={y + 4}
            stroke="#80868b"
            strokeWidth={1}
          />
          <text
            x={bp2x(p - 1)}
            y={y - 8}
            textAnchor="middle"
            fill="#5f6368"
            fontSize={10}
          >
            {p}
          </text>
        </g>
      ))}
    </g>
  );
}
