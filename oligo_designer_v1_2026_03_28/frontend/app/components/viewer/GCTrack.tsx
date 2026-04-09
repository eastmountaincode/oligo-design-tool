interface Props {
  sequence: string;
  totalLength: number;
  bp2x: (bp: number) => number;
  y: number;
  height: number;
  windowSize?: number;
}

const GC_TRACK_HEIGHT = 40;

export function getGCTrackHeight() {
  return GC_TRACK_HEIGHT + 12; // 12 for label
}

export default function GCTrack({
  sequence,
  totalLength,
  bp2x,
  y,
  height,
  windowSize = 20,
}: Props) {
  // Compute GC% in sliding windows
  const trackY = y + 12;
  const trackH = height;
  const seq = sequence.toUpperCase();

  // If we don't have the full sequence, use totalLength to draw an empty track
  if (seq.length < windowSize) return null;

  // Compute GC values for each window position
  const gcValues: { pos: number; gc: number }[] = [];
  const halfWin = Math.floor(windowSize / 2);

  // Initialize first window
  let gcCount = 0;
  for (let i = 0; i < windowSize && i < seq.length; i++) {
    if (seq[i] === "G" || seq[i] === "C") gcCount++;
  }
  gcValues.push({ pos: halfWin, gc: gcCount / windowSize });

  // Slide window
  for (let i = 1; i <= seq.length - windowSize; i++) {
    const dropped = seq[i - 1];
    const added = seq[i + windowSize - 1];
    if (dropped === "G" || dropped === "C") gcCount--;
    if (added === "G" || added === "C") gcCount++;
    gcValues.push({ pos: i + halfWin, gc: gcCount / windowSize });
  }

  // Build bar path — each window center gets a bar proportional to GC%
  // Color bars based on GC% — blue for balanced, orange-red for extreme
  const midY = trackY + trackH / 2; // 50% line

  return (
    <g>
      {/* Background */}
      <rect
        x={bp2x(0)}
        y={trackY}
        width={bp2x(totalLength) - bp2x(0)}
        height={trackH}
        fill="#f8f9fa"
        stroke="#e8eaed"
        strokeWidth={0.5}
      />

      {/* 50% reference line */}
      <line
        x1={bp2x(0)}
        y1={midY}
        x2={bp2x(totalLength)}
        y2={midY}
        stroke="#dadce0"
        strokeWidth={1}
        strokeDasharray="4 2"
      />

      {/* 25% and 75% threshold lines */}
      <line
        x1={bp2x(0)}
        y1={trackY + trackH * 0.25}
        x2={bp2x(totalLength)}
        y2={trackY + trackH * 0.25}
        stroke="#e8eaed"
        strokeWidth={0.5}
      />
      <line
        x1={bp2x(0)}
        y1={trackY + trackH * 0.75}
        x2={bp2x(totalLength)}
        y2={trackY + trackH * 0.75}
        stroke="#e8eaed"
        strokeWidth={0.5}
      />

      {/* GC bars */}
      {gcValues.map(({ pos, gc }, i) => {
        const x = bp2x(pos);
        const nextX = i < gcValues.length - 1 ? bp2x(gcValues[i + 1].pos) : x + 1;
        const barW = Math.max(nextX - x, 0.5);

        // Height from the 50% line — bars extend up for >50%, down for <50%
        const deviation = gc - 0.5;
        const barH = Math.abs(deviation) * trackH;
        const barY = deviation >= 0 ? midY - barH : midY;

        // Color: green in ideal range (40-60%), orange outside, red at extremes
        let fill: string;
        if (gc >= 0.40 && gc <= 0.60) {
          fill = "#34a853"; // ideal — green
        } else if (gc >= 0.30 && gc <= 0.70) {
          fill = "#e8710a"; // warning — orange
        } else {
          fill = "#d93025"; // extreme — red
        }

        return (
          <rect
            key={i}
            x={x}
            y={barY}
            width={barW}
            height={Math.max(barH, 0.5)}
            fill={fill}
            opacity={0.7}
          />
        );
      })}

    </g>
  );
}
