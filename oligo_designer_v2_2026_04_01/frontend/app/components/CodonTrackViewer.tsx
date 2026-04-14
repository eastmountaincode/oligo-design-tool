"use client";

import { useState, useRef, useMemo, useCallback, useEffect } from "react";
import { createPortal } from "react-dom";

interface CodonDetail {
  amino_acid: string;
  codon: string;
  per_thousand: number;
}

interface CodonTableEntry {
  amino_acid: string;
  codon: string;
  per_thousand: number;
}

interface SequenceWarning {
  kind: string;
  message: string;
  start: number;
  end: number;
  group?: number | null;
}


interface Props {
  codons: CodonDetail[];
  codonTable: CodonTableEntry[];
  warnings: SequenceWarning[];
  originalWarnings?: SequenceWarning[];
  onDnaChange?: (newDna: string) => void;
  onUndo?: () => void;
  canUndo?: boolean;
}

// Layout constants
const CODON_W = 36;
const LEFT_MARGIN = 60;
const RULER_H = 20;
const WARNING_TRACK_H = 16;
const WARNING_GAP = 4;
const AA_ROW_H = 16;
const CODON_ROW_H = 18;
const VALUE_ROW_H = 14;
const ROW_GAP = 1;
const BOTTOM_PAD = 8;

function codonRank(
  perThousand: number,
  alts: { per_thousand: number }[]
): "best" | "worst" | "middle" {
  if (alts.length <= 1) return "best";
  const max = Math.max(...alts.map((a) => a.per_thousand));
  const min = Math.min(...alts.map((a) => a.per_thousand));
  if (max === min) return "best";
  if (perThousand === max) return "best";
  if (perThousand === min) return "worst";
  return "middle";
}

function rankBg(rank: "best" | "worst" | "middle"): string {
  if (rank === "best") return "#e6f4ea";
  if (rank === "middle") return "#fef7e0";
  return "#fce8e6";
}

function rankColor(rank: "best" | "worst" | "middle"): string {
  if (rank === "best") return "#188038";
  if (rank === "middle") return "#e37400";
  return "#d93025";
}

const WARNING_COLORS: Record<string, string> = {
  homopolymer: "#d93025",
  gc_window: "#0097a7",
  repeat_kmer: "#7b1fa2",
  low_codon_frequency: "#e37400",
};

export default function CodonTrackViewer({ codons, codonTable, warnings, originalWarnings, onDnaChange, onUndo, canUndo }: Props) {
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);
  const [zoom, setZoom] = useState(1);
  const [dropdownPos, setDropdownPos] = useState<{ left: number; top: number } | null>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  const aaAlternatives = useMemo(() => {
    const map: Record<string, { codon: string; per_thousand: number }[]> = {};
    for (const e of codonTable) {
      if (!map[e.amino_acid]) map[e.amino_acid] = [];
      map[e.amino_acid].push({ codon: e.codon, per_thousand: e.per_thousand });
    }
    for (const aa of Object.keys(map)) {
      map[aa].sort((a, b) => b.per_thousand - a.per_thousand);
    }
    return map;
  }, [codonTable]);

  // Compute dropdown screen position when selection changes
  useEffect(() => {
    if (selectedIndex === null || !svgRef.current) {
      setDropdownPos(null);
      return;
    }
    const svgRect = svgRef.current.getBoundingClientRect();
    const cW = CODON_W * zoom;
    const x = selectedIndex * cW;
    // Recompute contentStartY locally since hooks run before render variables
    const wRows = assignWarningRows(warnings, cW);
    const nWarnRows = wRows.length > 0 ? Math.max(...wRows) + 1 : 0;
    const wTrackH = nWarnRows > 0 ? nWarnRows * (WARNING_TRACK_H + 2) + WARNING_GAP : 0;
    const cStartY = RULER_H + wTrackH;
    const left = svgRect.left + x + window.scrollX;
    const top = svgRect.top + cStartY + window.scrollY;
    setDropdownPos({ left, top });
  }, [selectedIndex, zoom, warnings]);

  // Close dropdown on click outside
  useEffect(() => {
    if (selectedIndex === null) return;
    const handler = (e: MouseEvent) => {
      const target = e.target as HTMLElement;
      if (target.closest("[data-codon-dropdown]")) return;
      setSelectedIndex(null);
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [selectedIndex]);

  const handleSwap = useCallback((index: number, newCodon: string) => {
    if (!onDnaChange) return;
    const newCodons = codons.map((c, i) => (i === index ? newCodon : c.codon));
    onDnaChange(newCodons.join(""));
    setSelectedIndex(null);
  }, [codons, onDnaChange]);

  // Zoom preserving viewport center
  const handleZoom = useCallback((newZoom: number) => {
    const sc = scrollRef.current;
    if (!sc) {
      setZoom(newZoom);
      return;
    }
    const oldCellW = CODON_W * zoom;
    const viewportCenter = sc.scrollLeft + sc.clientWidth / 2;
    const centerCodonIndex = viewportCenter / oldCellW;
    setZoom(newZoom);
    requestAnimationFrame(() => {
      const newCellW = CODON_W * newZoom;
      sc.scrollLeft = centerCodonIndex * newCellW - sc.clientWidth / 2;
    });
  }, [zoom]);

  if (!codons || codons.length === 0) return null;

  const cellW = CODON_W * zoom;
  const totalW = codons.length * cellW;

  // Use original warnings for layout (stable positions), mark which are still active
  const displayWarnings = useMemo(() => {
    const base = originalWarnings && originalWarnings.length > 0 ? originalWarnings : warnings;
    return base.map((w) => {
      const stillActive = warnings.some(
        (aw) => aw.kind === w.kind && aw.start === w.start && aw.end === w.end
      );
      return { ...w, resolved: !stillActive };
    });
  }, [warnings, originalWarnings]);

  // Also include any NEW warnings that weren't in the original
  const newWarnings = useMemo(() => {
    if (!originalWarnings) return [];
    return warnings.filter(
      (w) => !originalWarnings.some(
        (ow) => ow.kind === w.kind && ow.start === w.start && ow.end === w.end
      )
    ).map((w) => ({ ...w, resolved: false }));
  }, [warnings, originalWarnings]);

  const allTrackWarnings = useMemo(() => [
    ...displayWarnings,
    ...newWarnings,
  ], [displayWarnings, newWarnings]);

  // Warning track: assign rows to avoid overlap
  // Reserve space for at least 2 rows so the layout doesn't jump when warnings appear/disappear
  const MIN_WARNING_ROWS = 2;
  const warningRows = assignWarningRows(allTrackWarnings, cellW);
  const numWarningRows = warningRows.length > 0 ? Math.max(...warningRows) + 1 : 0;
  const displayRows = Math.max(numWarningRows, MIN_WARNING_ROWS);
  const warningTrackH = displayRows * (WARNING_TRACK_H + 2) + WARNING_GAP;

  const contentStartY = RULER_H + warningTrackH;
  const totalH = contentStartY + AA_ROW_H + ROW_GAP + CODON_ROW_H + ROW_GAP + VALUE_ROW_H + BOTTOM_PAD;

  function bp2x(bp: number): number {
    return (bp / 3) * cellW;
  }

  function assignWarningRows(ws: SequenceWarning[], cw: number): number[] {
    const rows: number[] = [];
    const rowEnds: number[] = [];
    for (const w of ws) {
      let placed = false;
      for (let r = 0; r < rowEnds.length; r++) {
        if (w.start >= rowEnds[r]) {
          rows.push(r);
          rowEnds[r] = w.end;
          placed = true;
          break;
        }
      }
      if (!placed) {
        rows.push(rowEnds.length);
        rowEnds.push(w.end);
      }
    }
    return rows;
  }

  return (
    <div className="bg-white border border-[#dadce0] p-4">
      <div className="flex items-center justify-between mb-3">
        <div>
          <h2 className="text-lg font-medium text-[#202124]">Codon Track</h2>
          <p className="text-xs text-[#5f6368]">
            Click a codon to swap. Color = quality relative to alternatives for the same amino acid.
          </p>
        </div>
        <div className="flex items-center gap-3">
          {onUndo && (
            <button
              onClick={onUndo}
              disabled={!canUndo}
              title="Undo last codon swap"
              className={`px-2 h-7 border border-[#dadce0] bg-white text-xs flex items-center gap-1 ${
                canUndo
                  ? "text-[#202124] cursor-pointer hover:bg-[#f1f3f4]"
                  : "text-[#c0c0c0] cursor-default"
              }`}
            >
              <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                <polyline points="1 4 1 10 7 10" />
                <path d="M3.51 15a9 9 0 1 0 2.13-9.36L1 10" />
              </svg>
              Undo
            </button>
          )}
          <div className="flex items-center gap-2">
            <button
              onClick={() => handleZoom(Math.max(0.5, zoom - 0.25))}
              className="w-7 h-7 border border-[#dadce0] bg-white text-[#202124] flex items-center justify-center cursor-pointer hover:bg-[#f1f3f4]"
            >
              -
            </button>
            <span className="text-xs text-[#5f6368] w-10 text-center">{(zoom * 100).toFixed(0)}%</span>
            <button
              onClick={() => handleZoom(Math.min(3, zoom + 0.25))}
              className="w-7 h-7 border border-[#dadce0] bg-white text-[#202124] flex items-center justify-center cursor-pointer hover:bg-[#f1f3f4]"
            >
              +
            </button>
          </div>
        </div>
      </div>

      <div className="flex">
        {/* Fixed labels column */}
        <div style={{ width: LEFT_MARGIN, flexShrink: 0 }}>
          <svg width={LEFT_MARGIN} height={totalH} className="select-none">
            <text x={4} y={RULER_H - 4} fontSize={9} fill="#5f6368">nt</text>
            {warningTrackH > 0 && (
              <text x={4} y={RULER_H + warningTrackH / 2 + 3} fontSize={9} fill="#5f6368">Issues</text>
            )}
            <text x={4} y={contentStartY + AA_ROW_H - 3} fontSize={9} fill="#5f6368">AA</text>
            <text x={4} y={contentStartY + AA_ROW_H + ROW_GAP + CODON_ROW_H - 3} fontSize={9} fill="#5f6368">Codon</text>
            <text x={4} y={contentStartY + AA_ROW_H + ROW_GAP + CODON_ROW_H + ROW_GAP + VALUE_ROW_H - 3} fontSize={9} fill="#5f6368">/1000</text>
          </svg>
        </div>

        {/* Scrollable content */}
        <div ref={scrollRef} className="overflow-x-auto flex-1">
        <svg
          ref={svgRef}
          width={totalW}
          height={totalH}
          style={{ display: "block" }}
          onClick={() => setSelectedIndex(null)}
        >
          {/* Ruler */}
          <g>
            {codons.map((_, i) => {
              const x = i * cellW;
              const bp = i * 3 + 1;
              const showLabel = i % Math.max(1, Math.round(5 / zoom)) === 0;
              return (
                <g key={`ruler-${i}`}>
                  {showLabel && (
                    <text
                      x={x + cellW / 2}
                      y={RULER_H - 4}
                      textAnchor="middle"
                      fontSize={9}
                      fill="#5f6368"
                    >
                      {bp}
                    </text>
                  )}
                  <line
                    x1={x}
                    y1={RULER_H - 2}
                    x2={x}
                    y2={RULER_H}
                    stroke="#dadce0"
                    strokeWidth={0.5}
                  />
                </g>
              );
            })}
            <line
              x1={0}
              y1={RULER_H}
              x2={totalW}
              y2={RULER_H}
              stroke="#dadce0"
              strokeWidth={0.5}
            />
          </g>

          {/* Warning track — active warnings are solid, resolved ones are faded in place */}
          {allTrackWarnings.map((w, i) => {
            const row = warningRows[i] ?? 0;
            const x1 = bp2x(w.start);
            const x2 = bp2x(w.end);
            const y = RULER_H + 2 + row * (WARNING_TRACK_H + 2);
            const color = WARNING_COLORS[w.kind] ?? "#e37400";
            return (
              <g key={`warn-${i}`}>
                <rect
                  x={x1}
                  y={y}
                  width={Math.max(x2 - x1, 3)}
                  height={WARNING_TRACK_H}
                  fill={color}
                  opacity={w.resolved ? 0.15 : 0.7}
                  rx={1}
                />
                {cellW >= 20 && (
                  <text
                    x={x1 + 3}
                    y={y + WARNING_TRACK_H - 4}
                    fontSize={8}
                    fill={w.resolved ? "#c0c0c0" : "white"}
                    fontWeight="bold"
                  >
                    {w.message.split(" at ")[0]}
                  </text>
                )}
              </g>
            );
          })}

          {/* Codon cells */}
          {codons.map((c, i) => {
            const alts = aaAlternatives[c.amino_acid] || [];
            const rank = codonRank(c.per_thousand, alts);
            const x = i * cellW;
            const isSelected = selectedIndex === i;

            return (
              <g
                key={`codon-${i}`}
                onClick={(e) => {
                  e.stopPropagation();
                  setSelectedIndex(isSelected ? null : i);
                }}
                style={{ cursor: "pointer" }}
              >
                {/* Background */}
                <rect
                  x={x}
                  y={contentStartY}
                  width={cellW}
                  height={AA_ROW_H + ROW_GAP + CODON_ROW_H + ROW_GAP + VALUE_ROW_H}
                  fill={rankBg(rank)}
                  stroke={isSelected ? "#1a73e8" : "#dadce0"}
                  strokeWidth={isSelected ? 2 : 0.5}
                />

                {/* Amino acid */}
                <text
                  x={x + cellW / 2}
                  y={contentStartY + AA_ROW_H - 3}
                  textAnchor="middle"
                  fontSize={10}
                  fill="#5f6368"
                >
                  {c.amino_acid}
                </text>

                {/* Codon */}
                <text
                  x={x + cellW / 2}
                  y={contentStartY + AA_ROW_H + ROW_GAP + CODON_ROW_H - 3}
                  textAnchor="middle"
                  fontSize={cellW >= 30 ? 11 : 9}
                  fontFamily="monospace"
                  fontWeight="600"
                  fill={rankColor(rank)}
                >
                  {c.codon}
                </text>

                {/* /1000 value */}
                <text
                  x={x + cellW / 2}
                  y={contentStartY + AA_ROW_H + ROW_GAP + CODON_ROW_H + ROW_GAP + VALUE_ROW_H - 3}
                  textAnchor="middle"
                  fontSize={8}
                  fill="#80868b"
                >
                  {c.per_thousand.toFixed(1)}
                </text>
              </g>
            );
          })}
        </svg>

      </div>
      </div>

      {/* Dropdown rendered via portal so it's never clipped by overflow */}
      {selectedIndex !== null && dropdownPos && (() => {
        const c = codons[selectedIndex];
        const alts = aaAlternatives[c.amino_acid] || [];
        if (alts.length <= 1) return null;

        return createPortal(
          <div
            data-codon-dropdown
            style={{
              position: "absolute",
              left: dropdownPos.left,
              top: dropdownPos.top,
              transform: "translateY(-100%)",
              zIndex: 9999,
            }}
            className="bg-white border border-[#dadce0] shadow-lg min-w-[110px]"
          >
            <div className="px-2 py-1 text-[10px] text-[#5f6368] border-b border-[#dadce0]">
              {c.amino_acid} codons
            </div>
            {alts.map((alt) => {
              const altRank = codonRank(alt.per_thousand, alts);
              const isCurrent = alt.codon === c.codon;
              return (
                <button
                  key={alt.codon}
                  onClick={() => {
                    if (!isCurrent) handleSwap(selectedIndex, alt.codon);
                    else setSelectedIndex(null);
                  }}
                  className={`w-full text-left px-2 py-1 text-xs flex items-center justify-between gap-3 hover:bg-[#f1f3f4] cursor-pointer ${
                    isCurrent ? "bg-[#e8f0fe]" : ""
                  }`}
                >
                  <span className="font-mono" style={{ color: rankColor(altRank) }}>
                    {alt.codon}
                  </span>
                  <span className="text-[#5f6368]">{alt.per_thousand.toFixed(1)}/1000</span>
                </button>
              );
            })}
          </div>,
          document.body
        );
      })()}

      {/* Legend */}
      <div className="mt-3 text-xs text-[#5f6368] space-y-1">
        <div className="flex items-center gap-4">
          <span className="text-[#80868b] w-20 flex-shrink-0">Codon:</span>
          <div className="flex items-center gap-1">
            <span className="inline-block w-3 h-3" style={{ backgroundColor: "#e6f4ea" }} />
            Most frequent
          </div>
          <div className="flex items-center gap-1">
            <span className="inline-block w-3 h-3" style={{ backgroundColor: "#fef7e0" }} />
            Intermediate
          </div>
          <div className="flex items-center gap-1">
            <span className="inline-block w-3 h-3" style={{ backgroundColor: "#fce8e6" }} />
            Least frequent
          </div>
        </div>
        {warnings.length > 0 && (
          <div className="flex items-center gap-4 flex-wrap">
            <span className="text-[#80868b] w-20 flex-shrink-0">Warnings:</span>
            <div className="flex items-center gap-1">
              <span className="inline-block w-3 h-3" style={{ backgroundColor: "#d93025" }} />
              Homopolymer run
            </div>
            <div className="flex items-center gap-1">
              <span className="inline-block w-3 h-3" style={{ backgroundColor: "#0097a7" }} />
              GC content
            </div>
            <div className="flex items-center gap-1">
              <span className="inline-block w-3 h-3" style={{ backgroundColor: "#7b1fa2" }} />
              Repeated k-mer
            </div>
            <div className="flex items-center gap-1">
              <span className="inline-block w-3 h-3" style={{ backgroundColor: "#e37400" }} />
              Low codon frequency
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
