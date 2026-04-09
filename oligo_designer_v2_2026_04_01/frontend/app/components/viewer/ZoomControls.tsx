interface Props {
  pxPerBp: number;
  minPxPerBp: number;
  maxPxPerBp: number;
  totalLength: number;
  onZoom: (newPxPerBp: number) => void;
  onFit: () => void;
}

export default function ZoomControls({
  pxPerBp,
  minPxPerBp,
  maxPxPerBp,
  totalLength,
  onZoom,
  onFit,
}: Props) {
  const zoomPercent =
    ((Math.log(pxPerBp) - Math.log(minPxPerBp)) /
      (Math.log(maxPxPerBp) - Math.log(minPxPerBp))) *
    100;

  function handleSlider(val: number) {
    const logMin = Math.log(minPxPerBp);
    const logMax = Math.log(maxPxPerBp);
    onZoom(Math.exp(logMin + (val / 100) * (logMax - logMin)));
  }

  const btnClass =
    "w-7 h-7 flex items-center justify-center bg-[#f1f3f4] border border-[#dadce0] text-[#202124] hover:bg-[#e8eaed] text-sm font-medium cursor-pointer";

  return (
    <div className="flex items-center gap-3">
      <span className="text-xs text-[#5f6368]">Construct: {totalLength} bp</span>
      <div className="flex items-center gap-2">
        <button
          onClick={() => onZoom(Math.max(minPxPerBp, pxPerBp / 1.5))}
          className={btnClass}
          title="Zoom out"
        >
          -
        </button>
        <input
          type="range"
          min={0}
          max={100}
          value={zoomPercent}
          onChange={(e) => handleSlider(Number(e.target.value))}
          className="w-48 h-1 accent-[#1a73e8] cursor-pointer"
          title={`${pxPerBp.toFixed(1)} px/bp`}
        />
        <button
          onClick={() => onZoom(Math.min(maxPxPerBp, pxPerBp * 1.5))}
          className={btnClass}
          title="Zoom in"
        >
          +
        </button>
        <button
          onClick={onFit}
          className="text-xs bg-[#f1f3f4] border border-[#dadce0] text-[#202124] hover:bg-[#e8eaed] px-2 py-1 ml-1 cursor-pointer"
          title="Fit to window"
        >
          Fit
        </button>
      </div>
    </div>
  );
}
