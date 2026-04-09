"use client";

import { useState, useRef, useEffect, useCallback, useMemo } from "react";
import type { OligoData, OverlapData, SequenceIssue, TrackLayout } from "./viewer/types";
import { reconstructSequence } from "./viewer/sequence-utils";
import ZoomControls from "./viewer/ZoomControls";
import Ruler from "./viewer/Ruler";
import IssueTrack, { getIssueTrackHeight } from "./viewer/IssueTrack";
import GCTrack, { getGCTrackHeight } from "./viewer/GCTrack";
import OligoArrow from "./viewer/OligoArrow";
import OverlapRegion from "./viewer/OverlapRegion";
import OverlapTooltip from "./viewer/OverlapTooltip";
import OligoDetailPanel from "./viewer/OligoDetailPanel";
import { ISSUE_COLORS } from "./viewer/colors";

interface Props {
  oligos: OligoData[];
  overlaps: OverlapData[];
  totalLength: number;
  sequenceIssues: SequenceIssue[];
}

const MIN_PX_PER_BP = 0.3;
const MAX_PX_PER_BP = 14;
const LABEL_WIDTH = 55;

export default function OligoViewer({ oligos, overlaps, totalLength, sequenceIssues }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const initializedRef = useRef(false);

  const [containerWidth, setContainerWidth] = useState(900);
  const [pxPerBp, setPxPerBp] = useState(0.5);
  const [selectedOligo, setSelectedOligo] = useState<OligoData | null>(null);
  const [selectedOverlap, setSelectedOverlap] = useState<OverlapData | null>(null);
  const [selectedIssueIndex, setSelectedIssueIndex] = useState<number | null>(null);
  const [isPanning, setIsPanning] = useState(false);
  const [panStart, setPanStart] = useState({ x: 0, scrollLeft: 0 });
  const [showIssues, setShowIssues] = useState(true);
  const [showGC, setShowGC] = useState(true);

  // Clear other selections when one is made
  function selectOligo(oligo: OligoData | null) {
    setSelectedOligo(oligo);
    if (oligo) { setSelectedOverlap(null); setSelectedIssueIndex(null); }
  }
  function selectOverlap(overlap: OverlapData | null) {
    setSelectedOverlap(overlap);
    if (overlap) { setSelectedOligo(null); setSelectedIssueIndex(null); }
  }
  function selectIssue(index: number | null) {
    setSelectedIssueIndex(index);
    if (index !== null) { setSelectedOligo(null); setSelectedOverlap(null); }
  }

  // --- Resize: set initial zoom once ---
  useEffect(() => {
    if (!containerRef.current) return;
    const obs = new ResizeObserver((entries) => {
      const w = entries[0].contentRect.width;
      setContainerWidth(w);
      if (!initializedRef.current) {
        initializedRef.current = true;
        const fitPx = (w - LABEL_WIDTH - 20) / totalLength;
        setPxPerBp(Math.max(MIN_PX_PER_BP, Math.min(fitPx, MAX_PX_PER_BP)));
      }
    });
    obs.observe(containerRef.current);
    return () => obs.disconnect();
  }, [totalLength]);

  // --- Pan handlers ---
  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const sc = scrollRef.current;
    if (!sc || sc.scrollWidth <= sc.clientWidth) return;
    setIsPanning(true);
    setPanStart({ x: e.clientX, scrollLeft: sc.scrollLeft });
  }, []);

  const handleMouseMove = useCallback(
    (e: React.MouseEvent) => {
      if (!isPanning) return;
      const sc = scrollRef.current;
      if (sc) sc.scrollLeft = panStart.scrollLeft - (e.clientX - panStart.x);
    },
    [isPanning, panStart]
  );

  const handleMouseUp = useCallback(() => setIsPanning(false), []);

  // --- Derived layout ---
  // For the scrollable content area, bp2x starts at 0 (no label offset)
  const bp2x = useCallback((bp: number) => bp * pxPerBp, [pxPerBp]);
  const contentWidth = Math.max(totalLength * pxPerBp + 40, containerWidth - LABEL_WIDTH - 32);
  const showBases = pxPerBp >= 6;
  const showNumbers = pxPerBp >= 1.5;
  const fullSeq = useMemo(
    () => (showBases ? reconstructSequence(oligos, totalLength) : reconstructSequence(oligos, totalLength)),
    [oligos, totalLength]
  );

  const senseOligos = oligos.filter((o) => o.strand === "sense");
  const antisenseOligos = oligos.filter((o) => o.strand === "antisense");

  // Merge overlap issues into the issue track
  const overlapIssuesAsSequenceIssues: SequenceIssue[] = [];
  for (const ovl of overlaps) {
    for (const iss of ovl.issues) {
      const issStart = iss.start ?? ovl.start;
      const issEnd = iss.end ?? ovl.end;
      const isDuplicate = sequenceIssues.some((si) =>
        si.start !== undefined && si.end !== undefined &&
        si.type === iss.kind &&
        si.start <= issStart && si.end >= issEnd
      );
      if (!isDuplicate) {
        overlapIssuesAsSequenceIssues.push({
          type: `overlap_${iss.kind}`,
          message: `Overlap ${ovl.index}: ${iss.message}`,
          severity: iss.severity,
          start: issStart,
          end: issEnd,
        });
      }
    }
  }
  const allIssues = [...sequenceIssues, ...overlapIssuesAsSequenceIssues];

  // --- Compute track Y positions ---
  const RULER_Y = 20;
  const RULER_H = 30;
  let currentY = RULER_Y + RULER_H;

  const issueTrackY = currentY;
  const issueTrackH = showIssues ? getIssueTrackHeight(allIssues) : 0;
  if (showIssues && issueTrackH > 0) currentY += issueTrackH + 6;

  const gcTrackY = currentY;
  const gcTrackH = showGC ? getGCTrackHeight() : 0;
  if (showGC) currentY += gcTrackH + 6;

  const senseY = currentY;
  const oligoH = showBases ? 22 : 26;
  currentY += oligoH;

  const antisenseY = currentY + 40;
  currentY = antisenseY + oligoH + 15;

  const layout: TrackLayout = {
    RULER_Y,
    ISSUE_TRACK_Y: issueTrackY,
    SENSE_Y: senseY,
    OLIGO_H: oligoH,
    CENTER_Y: senseY + oligoH + 20,
    ANTISENSE_Y: antisenseY,
    OVERLAP_H: antisenseY + oligoH - senseY,
  };
  const svgHeight = currentY;

  // Zoom centered on viewport middle
  const handleZoom = useCallback((newPx: number) => {
    const sc = scrollRef.current;
    if (!sc) {
      setPxPerBp(newPx);
      return;
    }
    const viewportCenter = sc.scrollLeft + sc.clientWidth / 2;
    const centerBp = viewportCenter / pxPerBp;
    setPxPerBp(newPx);
    requestAnimationFrame(() => {
      sc.scrollLeft = centerBp * newPx - sc.clientWidth / 2;
    });
  }, [pxPerBp]);

  function handleFit() {
    const fitPx = (containerWidth - LABEL_WIDTH - 40) / totalLength;
    setPxPerBp(Math.max(MIN_PX_PER_BP, Math.min(fitPx, MAX_PX_PER_BP)));
    if (scrollRef.current) scrollRef.current.scrollLeft = 0;
  }

  const issueRegions = allIssues.filter((i) => i.start !== undefined && i.end !== undefined);
  const selectedIssue = selectedIssueIndex !== null ? issueRegions[selectedIssueIndex] : null;

  return (
    <div ref={containerRef} className="bg-white border border-[#dadce0] p-4">
      {/* Header */}
      <div className="flex items-center justify-between mb-3">
        <h2 className="text-lg font-medium text-[#202124]">Assembly View</h2>
        <ZoomControls
          pxPerBp={pxPerBp}
          minPxPerBp={MIN_PX_PER_BP}
          maxPxPerBp={MAX_PX_PER_BP}
          totalLength={totalLength}
          onZoom={handleZoom}
          onFit={handleFit}
        />
      </div>

      {/* Track toggles */}
      <div className="text-xs text-[#5f6368] mb-2">
        <span className="text-[#202124] font-medium">Tracks:</span>
        <div className="mt-1 space-y-1">
          <label className="flex items-center gap-1.5 cursor-pointer">
            <input
              type="checkbox"
              checked={showIssues}
              onChange={(e) => setShowIssues(e.target.checked)}
              className="accent-[#1a73e8] cursor-pointer"
            />
            Issues
          </label>
          <label className="flex items-center gap-1.5 cursor-pointer">
            <input
              type="checkbox"
              checked={showGC}
              onChange={(e) => setShowGC(e.target.checked)}
              className="accent-[#1a73e8] cursor-pointer"
            />
            GC% (20bp window)
          </label>
        </div>
      </div>

      {/* Main viewer: fixed labels + scrollable content + fixed right scale */}
      <div className="flex">
        {/* Fixed labels column (left) */}
        <div style={{ width: LABEL_WIDTH, flexShrink: 0 }}>
          <svg width={LABEL_WIDTH} height={svgHeight} className="select-none">
            <text x={4} y={RULER_Y + RULER_H - 12} fill="#5f6368" fontSize={9}>
              nt
            </text>
            {showIssues && issueTrackH > 0 && (
              <text x={4} y={issueTrackY + 9} fill="#5f6368" fontSize={9}>
                Issues
              </text>
            )}

            {showGC && (
              <>
                <text x={4} y={gcTrackY + 9} fill="#5f6368" fontSize={9}>
                  GC%
                </text>
                <text x={4} y={gcTrackY + 19} fill="#5f6368" fontSize={7}>
                  20bp win
                </text>
              </>
            )}

            <text x={4} y={senseY + oligoH / 2 + 4} fill="#1a73e8" fontSize={10}>
              5&apos;-3&apos;
            </text>
            <text x={4} y={antisenseY + oligoH / 2 + 4} fill="#e8710a" fontSize={10}>
              3&apos;-5&apos;
            </text>
          </svg>
        </div>

        {/* Scrollable content (middle) */}
        <div
          ref={scrollRef}
          className="overflow-x-auto overflow-y-hidden flex-1"
          style={{
            cursor: isPanning ? "grabbing" : contentWidth > containerWidth - LABEL_WIDTH - 32 ? "grab" : "default",
          }}
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseUp}
        >
          <svg width={contentWidth} height={svgHeight} className="select-none">
            <Ruler
              totalLength={totalLength}
              bp2x={bp2x}
              containerWidth={containerWidth - LABEL_WIDTH}
              pxPerBp={pxPerBp}
              y={layout.RULER_Y}
            />

            {showIssues && (
              <IssueTrack
                issues={allIssues}
                bp2x={bp2x}
                y={layout.ISSUE_TRACK_Y}
                selectedIndex={selectedIssueIndex}
                onSelect={selectIssue}
              />
            )}

            {showGC && (
              <GCTrack
                sequence={fullSeq}
                totalLength={totalLength}
                bp2x={bp2x}
                y={gcTrackY}
                height={40}
              />
            )}

            {/* Overlaps (behind oligos) */}
            {overlaps.map((ovl) => (
              <OverlapRegion
                key={ovl.index}
                overlap={ovl}
                bp2x={bp2x}
                layout={layout}
                showNumbers={showNumbers}
                isSelected={selectedOverlap?.index === ovl.index}
                onSelect={selectOverlap}
              />
            ))}

            {/* Sense oligos */}
            {senseOligos.map((o) => (
              <OligoArrow
                key={o.index}
                oligo={o}
                bp2x={bp2x}
                pxPerBp={pxPerBp}
                y={layout.SENSE_Y}
                height={layout.OLIGO_H}
                fullSeq={fullSeq}
                showBases={showBases}
                showNumbers={showNumbers}
                isSelected={selectedOligo?.index === o.index}
                onSelect={selectOligo}
              />
            ))}

            {/* Antisense oligos */}
            {antisenseOligos.map((o) => (
              <OligoArrow
                key={o.index}
                oligo={o}
                bp2x={bp2x}
                pxPerBp={pxPerBp}
                y={layout.ANTISENSE_Y}
                height={layout.OLIGO_H}
                fullSeq={fullSeq}
                showBases={showBases}
                showNumbers={showNumbers}
                isSelected={selectedOligo?.index === o.index}
                onSelect={selectOligo}
              />
            ))}
          </svg>
        </div>

        {/* Fixed scale column (right) */}
        {showGC && (
          <div style={{ width: 32, flexShrink: 0 }}>
            <svg width={32} height={svgHeight} className="select-none">
              <text x={4} y={gcTrackY + 12 + 7} fill="#5f6368" fontSize={7}>
                100%
              </text>
              <text x={4} y={gcTrackY + 12 + 23} fill="#5f6368" fontSize={7}>
                50%
              </text>
              <text x={4} y={gcTrackY + 12 + 40} fill="#5f6368" fontSize={7}>
                0%
              </text>
            </svg>
          </div>
        )}
      </div>

      {/* Selected overlap detail */}
      {selectedOverlap && (
        <OverlapTooltip overlap={selectedOverlap} onClose={() => setSelectedOverlap(null)} />
      )}

      {/* Selected oligo detail */}
      {selectedOligo && (
        <OligoDetailPanel oligo={selectedOligo} onClose={() => setSelectedOligo(null)} />
      )}

      {/* Selected issue detail */}
      {selectedIssue && (() => {
        const rawType = selectedIssue.type.replace(/^overlap_/, "");
        const colorInfo = ISSUE_COLORS[rawType];
        const color = colorInfo?.fill ?? "#d93025";
        const typeLabel = colorInfo?.label?.split("(")[0]?.trim() ?? rawType;
        return (
          <div className="mt-3 bg-white border border-[#dadce0] p-4">
            <div className="flex items-center justify-between mb-3">
              <div className="flex items-center gap-2">
                <span
                  className="inline-block w-3 h-3 flex-shrink-0"
                  style={{ backgroundColor: color }}
                />
                <span className="font-medium text-[#202124]">{typeLabel}</span>
                <span className="text-xs text-[#5f6368]">
                  {selectedIssue.severity === "error" ? "Error" : "Warning"}
                </span>
              </div>
              <button
                onClick={() => setSelectedIssueIndex(null)}
                className="text-[#5f6368] hover:text-[#202124] text-sm cursor-pointer"
              >
                close
              </button>
            </div>
            <table className="text-sm">
              <tbody>
                {selectedIssue.start !== undefined && (
                  <tr>
                    <td className="text-[#5f6368] pr-6 py-0.5 align-top">Position</td>
                    <td className="py-0.5">{selectedIssue.start}-{selectedIssue.end}</td>
                  </tr>
                )}
                <tr>
                  <td className="text-[#5f6368] pr-6 py-0.5 align-top">Detail</td>
                  <td className="py-0.5">{selectedIssue.message}</td>
                </tr>
              </tbody>
            </table>
          </div>
        );
      })()}
    </div>
  );
}
