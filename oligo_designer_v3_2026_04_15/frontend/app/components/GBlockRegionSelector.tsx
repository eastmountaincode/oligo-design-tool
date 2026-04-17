"use client";

import { useState, useRef, useCallback } from "react";

export interface GBlockRegion {
  id: string;
  startBp: number;
  endBp: number;
  label: string;
}

const GBLOCK_MIN_BP = 125;
const GBLOCK_MAX_BP = 3000;

interface Props {
  dnaLength: number;  // total DNA length in bp
  regions: GBlockRegion[];
  onChange: (regions: GBlockRegion[]) => void;
}

function nextId(): string {
  return "gb_" + Math.random().toString(36).slice(2, 8);
}

export default function GBlockRegionSelector({
  dnaLength,
  regions,
  onChange,
}: Props) {
  const barRef = useRef<HTMLDivElement>(null);
  const [dragging, setDragging] = useState<{
    regionId: string;
    edge: "left" | "right" | "body";
    startX: number;
    origStart: number;
    origEnd: number;
  } | null>(null);

  const bpToFrac = useCallback(
    (bp: number) => (dnaLength > 0 ? bp / dnaLength : 0),
    [dnaLength]
  );

  // Snap to codon boundary (multiple of 3)
  const snapToCodon = (bp: number) => Math.round(bp / 3) * 3;

  const addRegion = () => {
    const center = Math.floor(dnaLength / 2);
    const halfLen = Math.ceil(GBLOCK_MIN_BP / 2);
    const start = snapToCodon(Math.max(0, center - halfLen));
    const end = snapToCodon(Math.min(dnaLength, start + GBLOCK_MIN_BP));
    onChange([
      ...regions,
      { id: nextId(), startBp: start, endBp: end, label: "" },
    ]);
  };

  const removeRegion = (id: string) => {
    onChange(regions.filter((r) => r.id !== id));
  };

  const updateRegion = (id: string, patch: Partial<GBlockRegion>) => {
    onChange(regions.map((r) => (r.id === id ? { ...r, ...patch } : r)));
  };

  const onPointerDown = useCallback(
    (
      e: React.PointerEvent,
      regionId: string,
      edge: "left" | "right" | "body"
    ) => {
      e.preventDefault();
      e.stopPropagation();
      (e.target as HTMLElement).setPointerCapture(e.pointerId);
      const region = regions.find((r) => r.id === regionId)!;
      setDragging({
        regionId,
        edge,
        startX: e.clientX,
        origStart: region.startBp,
        origEnd: region.endBp,
      });
    },
    [regions]
  );

  const onPointerMove = useCallback(
    (e: React.PointerEvent) => {
      if (!dragging || !barRef.current) return;
      const rect = barRef.current.getBoundingClientRect();
      const pxPerBp = rect.width / dnaLength;
      const dx = e.clientX - dragging.startX;
      const dBp = snapToCodon(Math.round(dx / pxPerBp));

      let newStart = dragging.origStart;
      let newEnd = dragging.origEnd;

      if (dragging.edge === "left") {
        newStart = snapToCodon(Math.max(0, Math.min(dragging.origStart + dBp, newEnd - 3)));
      } else if (dragging.edge === "right") {
        newEnd = snapToCodon(Math.min(dnaLength, Math.max(dragging.origEnd + dBp, newStart + 3)));
      } else {
        const len = dragging.origEnd - dragging.origStart;
        newStart = snapToCodon(Math.max(0, Math.min(dragging.origStart + dBp, dnaLength - len)));
        newEnd = newStart + len;
      }

      updateRegion(dragging.regionId, { startBp: newStart, endBp: newEnd });
    },
    [dragging, dnaLength, updateRegion]
  );

  const onPointerUp = useCallback(() => {
    setDragging(null);
  }, []);

  if (dnaLength <= 0) return null;

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <label className="text-sm font-medium text-[#202124]">
          gBlock Regions
        </label>
        <button
          type="button"
          onClick={addRegion}
          className="text-xs text-[#1a73e8] hover:text-[#1967d2] cursor-pointer"
        >
          + Add gBlock region
        </button>
      </div>

      {regions.length > 0 && (
        <>
          {/* Visual bar */}
          <div
            ref={barRef}
            className="relative h-8 bg-[#f1f3f4] border border-[#dadce0] rounded select-none touch-none"
            onPointerMove={onPointerMove}
            onPointerUp={onPointerUp}
          >
            {/* Tick marks every 300 bp */}
            {Array.from(
              { length: Math.floor(dnaLength / 300) },
              (_, i) => (i + 1) * 300
            ).map((bp) => (
              <div
                key={bp}
                className="absolute top-0 h-full border-l border-[#dadce0]"
                style={{ left: `${bpToFrac(bp) * 100}%` }}
              >
                <span className="absolute -top-3.5 -translate-x-1/2 text-[9px] text-[#80868b]">
                  {bp}
                </span>
              </div>
            ))}

            {/* gBlock region rectangles */}
            {regions.map((r) => {
              const left = bpToFrac(r.startBp) * 100;
              const width = bpToFrac(r.endBp - r.startBp) * 100;
              const lenBp = r.endBp - r.startBp;
              const tooSmall = lenBp < GBLOCK_MIN_BP;
              const tooLarge = lenBp > GBLOCK_MAX_BP;
              const bgColor = tooSmall
                ? "rgba(234, 67, 53, 0.25)"
                : "rgba(66, 133, 244, 0.25)";
              const borderColor = tooSmall ? "#ea4335" : "#4285f4";

              return (
                <div
                  key={r.id}
                  className="absolute top-0 h-full"
                  style={{
                    left: `${left}%`,
                    width: `${width}%`,
                    backgroundColor: bgColor,
                    borderLeft: `2px solid ${borderColor}`,
                    borderRight: `2px solid ${borderColor}`,
                    cursor: dragging ? "grabbing" : "grab",
                    minWidth: 4,
                  }}
                  onPointerDown={(e) => onPointerDown(e, r.id, "body")}
                >
                  <div
                    className="absolute left-0 top-0 h-full w-2 cursor-ew-resize z-10"
                    style={{ transform: "translateX(-50%)" }}
                    onPointerDown={(e) => onPointerDown(e, r.id, "left")}
                  />
                  <div
                    className="absolute right-0 top-0 h-full w-2 cursor-ew-resize z-10"
                    style={{ transform: "translateX(50%)" }}
                    onPointerDown={(e) => onPointerDown(e, r.id, "right")}
                  />
                  {width > 5 && (
                    <span className="absolute inset-0 flex items-center justify-center text-[10px] font-medium text-[#202124] pointer-events-none overflow-hidden whitespace-nowrap">
                      {r.label || `${lenBp} bp`}
                    </span>
                  )}
                </div>
              );
            })}
          </div>

          {/* Per-region detail rows */}
          {regions.map((r) => {
            const lenBp = r.endBp - r.startBp;
            const tooSmall = lenBp < GBLOCK_MIN_BP;
            const tooLarge = lenBp > GBLOCK_MAX_BP;
            return (
              <div
                key={r.id}
                className="flex items-center gap-2 text-xs text-[#5f6368]"
              >
                <input
                  type="text"
                  value={r.label}
                  onChange={(e) =>
                    updateRegion(r.id, { label: e.target.value })
                  }
                  placeholder="Label"
                  className="w-24 px-1 py-0.5 border border-[#dadce0] rounded text-xs"
                />
                <span>nt</span>
                <input
                  type="number"
                  value={r.startBp}
                  min={0}
                  max={r.endBp - 3}
                  step={3}
                  onChange={(e) =>
                    updateRegion(r.id, {
                      startBp: snapToCodon(Math.max(0, Math.min(Number(e.target.value), r.endBp - 3))),
                    })
                  }
                  className="w-16 px-1 py-0.5 border border-[#dadce0] rounded text-xs font-mono"
                />
                <span>–</span>
                <input
                  type="number"
                  value={r.endBp}
                  min={r.startBp + 3}
                  max={dnaLength}
                  step={3}
                  onChange={(e) =>
                    updateRegion(r.id, {
                      endBp: snapToCodon(Math.min(dnaLength, Math.max(Number(e.target.value), r.startBp + 3))),
                    })
                  }
                  className="w-16 px-1 py-0.5 border border-[#dadce0] rounded text-xs font-mono"
                />
                <span
                  className={`font-mono ${tooSmall || tooLarge ? "text-[#ea4335]" : "text-[#5f6368]"}`}
                >
                  {lenBp} bp
                  {tooSmall && ` (min ${GBLOCK_MIN_BP})`}
                  {tooLarge && ` (max ${GBLOCK_MAX_BP})`}
                </span>
                <span className="font-mono">
                  ~${(lenBp * 0.07).toFixed(2)}
                </span>
                <button
                  type="button"
                  onClick={() => removeRegion(r.id)}
                  className="text-[#ea4335] hover:text-[#c5221f] cursor-pointer ml-auto"
                >
                  ×
                </button>
              </div>
            );
          })}
        </>
      )}

      {regions.length === 0 && (
        <p className="text-xs text-[#80868b]">
          No gBlock regions. Add one to mark a section that will be ordered as a
          pre-synthesized IDT gBlock fragment instead of being tiled into
          overlapping oligos. IDT gBlocks: 125–3000 bp, ~$0.07/bp.
        </p>
      )}
    </div>
  );
}
