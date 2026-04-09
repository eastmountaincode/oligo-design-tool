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

  const ticks: number[] = [];
  for (let t = 0; t <= totalLength; t += tickInterval) ticks.push(t);

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
      {ticks.map((t) => (
        <g key={`tick-${t}`}>
          <line
            x1={bp2x(t)}
            y1={y - 4}
            x2={bp2x(t)}
            y2={y + 4}
            stroke="#80868b"
            strokeWidth={1}
          />
          <text
            x={bp2x(t)}
            y={y - 8}
            textAnchor="middle"
            fill="#5f6368"
            fontSize={10}
          >
            {t}
          </text>
        </g>
      ))}
    </g>
  );
}
