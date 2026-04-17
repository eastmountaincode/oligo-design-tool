"use client";

import { useState, useRef, useEffect, useLayoutEffect, useCallback, useMemo } from "react";
import type { OligoData, OverlapData, SequenceIssue, TrackLayout } from "./viewer/types";
import { reconstructSequence, BASE_COLORS, complementBase } from "./viewer/sequence-utils";
import ZoomControls from "./viewer/ZoomControls";
import Ruler from "./viewer/Ruler";
import IssueTrack, { getIssueTrackHeight } from "./viewer/IssueTrack";
import GCTrack, { getGCTrackHeight } from "./viewer/GCTrack";
import OligoArrow from "./viewer/OligoArrow";
import OverlapRegion from "./viewer/OverlapRegion";
import OverlapTooltip from "./viewer/OverlapTooltip";
import OligoDetailPanel from "./viewer/OligoDetailPanel";
import GBlockDetailPanel from "./viewer/GBlockDetailPanel";
import { ISSUE_COLORS } from "./viewer/colors";

interface GBlockRegionView {
  start_bp: number;
  end_bp: number;
  label: string;
  len_bp: number;
  // Optional richer fields — if provided, clicking the gBlock opens a
  // detail panel (same UX as oligos/overlaps). When missing, the gBlock
  // renders but is not selectable.
  index?: number;
  seq?: string;
  gc?: number;
}

interface Props {
  oligos: OligoData[];
  overlaps: OverlapData[];
  totalLength: number;
  sequenceIssues: SequenceIssue[];
  gblockRegions?: GBlockRegionView[];
  fullSequence?: string;
  // Insert boundaries within the full construct (full = upstream + insert + downstream).
  // [0, insertStart) and [insertEnd, totalLength) are plasmid flanks.
  // Defaults to the entire construct being insert (no flanks) for back-compat.
  insertStart?: number;
  insertEnd?: number;
}

const MIN_PX_PER_BP = 0.3;
const MAX_PX_PER_BP = 14;
const LABEL_WIDTH = 55;

const GBLOCK_COLOR = "#00897b";        // teal — distinct from blue (sense), orange (antisense), green (overlaps)
const GBLOCK_FILL = "#e0f2f1";        // light teal for base-visible mode
const GBLOCK_FILL_SOLID = "#00897b";  // solid teal when zoomed out

// Plasmid flank styling — neutral gray so it reads as "fixed, not designed by us"
const FLANK_FILL = "rgba(128, 134, 139, 0.12)";
const FLANK_STROKE = "#80868b";
const FLANK_LABEL_COLOR = "#5f6368";

export default function OligoViewer({
  oligos,
  overlaps,
  totalLength,
  sequenceIssues,
  gblockRegions,
  fullSequence,
  insertStart = 0,
  insertEnd,
}: Props) {
  const effectiveInsertEnd = insertEnd ?? totalLength;
  const hasUpstreamFlank = insertStart > 0;
  const hasDownstreamFlank = effectiveInsertEnd < totalLength;
  const containerRef = useRef<HTMLDivElement>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const initializedRef = useRef(false);

  const [containerWidth, setContainerWidth] = useState(900);
  const [pxPerBp, setPxPerBp] = useState(0.5);
  const [selectedOligo, setSelectedOligo] = useState<OligoData | null>(null);
  const [selectedOverlap, setSelectedOverlap] = useState<OverlapData | null>(null);
  const [selectedIssueIndex, setSelectedIssueIndex] = useState<number | null>(null);
  const [selectedGblockKey, setSelectedGblockKey] = useState<string | null>(null);
  const [isPanning, setIsPanning] = useState(false);
  const [panStart, setPanStart] = useState({ x: 0, scrollLeft: 0 });
  const [showIssues, setShowIssues] = useState(true);
  const [showGC, setShowGC] = useState(true);

  // Clear other selections when one is made
  function selectOligo(oligo: OligoData | null) {
    setSelectedOligo(oligo);
    if (oligo) { setSelectedOverlap(null); setSelectedIssueIndex(null); setSelectedGblockKey(null); }
  }
  function selectOverlap(overlap: OverlapData | null) {
    setSelectedOverlap(overlap);
    if (overlap) { setSelectedOligo(null); setSelectedIssueIndex(null); setSelectedGblockKey(null); }
  }
  function selectIssue(index: number | null) {
    setSelectedIssueIndex(index);
    if (index !== null) { setSelectedOligo(null); setSelectedOverlap(null); setSelectedGblockKey(null); }
  }
  function selectGblock(key: string | null) {
    setSelectedGblockKey(key);
    if (key !== null) { setSelectedOligo(null); setSelectedOverlap(null); setSelectedIssueIndex(null); }
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

  const pxPerBpRef = useRef(pxPerBp);
  useEffect(() => { pxPerBpRef.current = pxPerBp; }, [pxPerBp]);

  // Pending scroll value written by the wheel handler, consumed by the
  // layout effect after React commits.  Prevents stale-DOM reads during
  // rapid pinch bursts.
  const pendingScrollRef = useRef<number | null>(null);

  useLayoutEffect(() => {
    const el = scrollRef.current;
    const pending = pendingScrollRef.current;
    if (el == null || pending == null) return;
    pendingScrollRef.current = null;
    el.scrollLeft = pending;
  }, [pxPerBp]);

  useEffect(() => {
    const el = scrollRef.current;
    if (!el) return;

    // Pinch gesture state.  We track:
    //   - anchorMouseX / anchorBp: the pointer-anchored zoom point,
    //     captured on the first event of a gesture.  The bp under the
    //     cursor should stay visually pinned to that screen position.
    //   - initialPx: the pxPerBp at gesture start.  All zoom levels in
    //     the gesture are computed from this base, not from the running
    //     value — so accumulated floating-point and rAF timing drift
    //     can't stack up.
    //   - cumulativeDelta: the running sum of e.deltaY for the gesture.
    //     Zoom target = initialPx * exp(-cumulativeDelta * k).  Because
    //     we sum raw deltas (not multiply factors), trackpad sensor
    //     noise that oscillates in sign cancels in the sum instead of
    //     producing visible zoom-back.
    let anchorMouseX: number | null = null;
    let anchorBp: number | null = null;
    let initialPx = 0;
    let cumulativeDelta = 0;
    let idleTimer: ReturnType<typeof setTimeout> | null = null;
    let rafId: number | null = null;

    const endGesture = () => {
      anchorMouseX = null;
      anchorBp = null;
      cumulativeDelta = 0;
      idleTimer = null;
    };

    const flush = () => {
      rafId = null;
      if (anchorMouseX === null || anchorBp === null) return;
      const target = initialPx * Math.exp(-cumulativeDelta * 0.01);
      const newPx = Math.max(MIN_PX_PER_BP, Math.min(MAX_PX_PER_BP, target));
      if (newPx === pxPerBpRef.current) return;
      pendingScrollRef.current = anchorBp * newPx - anchorMouseX;
      pxPerBpRef.current = newPx;
      setPxPerBp(newPx);
    };

    const onWheel = (e: WheelEvent) => {
      if (!e.ctrlKey && !e.metaKey) return;
      e.preventDefault();

      if (anchorMouseX === null || anchorBp === null) {
        const rect = el.getBoundingClientRect();
        const mouseX = e.clientX - rect.left;
        const scrollLeft = pendingScrollRef.current ?? el.scrollLeft;
        initialPx = pxPerBpRef.current;
        anchorMouseX = mouseX;
        anchorBp = (scrollLeft + mouseX) / initialPx;
        cumulativeDelta = 0;
      }

      cumulativeDelta += e.deltaY;

      if (rafId === null) rafId = requestAnimationFrame(flush);

      if (idleTimer) clearTimeout(idleTimer);
      idleTimer = setTimeout(endGesture, 150);
    };

    el.addEventListener("wheel", onWheel, { passive: false });
    return () => {
      el.removeEventListener("wheel", onWheel);
      if (idleTimer) clearTimeout(idleTimer);
      if (rafId !== null) cancelAnimationFrame(rafId);
    };
  }, []);

  // --- Derived layout ---
  // For the scrollable content area, bp2x starts at 0 (no label offset)
  const bp2x = useCallback((bp: number) => bp * pxPerBp, [pxPerBp]);
  const contentWidth = Math.max(totalLength * pxPerBp + 40, containerWidth - LABEL_WIDTH - 32);
  const showBases = pxPerBp >= 6;
  const showNumbers = pxPerBp >= 1.5;
  const fullSeq = useMemo(
    () => fullSequence || reconstructSequence(oligos, totalLength),
    [oligos, totalLength, fullSequence]
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

  // Amino acid translation track
  const AA_TRACK_H = 16;
  const showAA = pxPerBp >= 1.0;
  const aaTrackY = currentY;
  if (showAA) currentY += AA_TRACK_H + 8;
  currentY += 20; // extra padding so the scrollbar doesn't cover the last track

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
            {showAA && (
              <text x={4} y={aaTrackY + AA_TRACK_H / 2 + 4} fill="#5f6368" fontSize={9}>
                AA
              </text>
            )}
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

            {/* Plasmid flank shading — drawn behind everything else */}
            {(hasUpstreamFlank || hasDownstreamFlank) && (() => {
              const sY = layout.SENSE_Y;
              const aY = layout.ANTISENSE_Y;
              const oH = layout.OLIGO_H;
              const top = sY - 4;
              const bot = aY + oH + 4;
              const h = bot - top;
              const showLabel = pxPerBp >= 0.6;
              const flanks: { start: number; end: number; label: string; anchor: "start" | "end" }[] = [];
              if (hasUpstreamFlank) {
                flanks.push({ start: 0, end: insertStart, label: "upstream plasmid flank", anchor: "start" });
              }
              if (hasDownstreamFlank) {
                flanks.push({ start: effectiveInsertEnd, end: totalLength, label: "downstream plasmid flank", anchor: "end" });
              }
              return flanks.map((f, i) => {
                const x = bp2x(f.start);
                const w = Math.max(bp2x(f.end) - x, 1);
                const labelX = f.anchor === "start" ? x + 4 : x + w - 4;
                const textAnchor = f.anchor === "start" ? "start" : "end";
                return (
                  <g key={`flank-${i}`}>
                    <rect
                      x={x}
                      y={top}
                      width={w}
                      height={h}
                      fill={FLANK_FILL}
                      stroke={FLANK_STROKE}
                      strokeWidth={0.5}
                      strokeDasharray="3 2"
                      style={{ pointerEvents: "none" }}
                    />
                    {showLabel && w >= 40 && (
                      <text
                        x={labelX}
                        y={top - 2}
                        textAnchor={textAnchor}
                        fontSize={9}
                        fill={FLANK_LABEL_COLOR}
                        fontStyle="italic"
                        style={{ pointerEvents: "none" }}
                      >
                        {f.label} ({f.end - f.start} bp)
                      </text>
                    )}
                  </g>
                );
              });
            })()}

            {/* gBlock regions — double-stranded, styled like oligos */}
            {(gblockRegions ?? []).map((gb, i) => {
              const x = bp2x(gb.start_bp);
              const w = Math.max(bp2x(gb.end_bp) - x, 2);
              const sY = layout.SENSE_Y;
              const aY = layout.ANTISENSE_Y;
              const oH = layout.OLIGO_H;
              const gbKey = `${gb.start_bp}-${gb.end_bp}`;
              const isSelected = selectedGblockKey === gbKey;
              const selectable = gb.seq !== undefined && gb.gc !== undefined && gb.index !== undefined;
              const fill = showBases ? GBLOCK_FILL : GBLOCK_FILL_SOLID;
              const stroke = isSelected ? "#202124" : showBases ? GBLOCK_COLOR : "#ffffff";
              const strokeW = isSelected ? 1.5 : 1;

              // Render bases for one strand
              const renderStrandBases = (yPos: number, isComplement: boolean) => {
                const bases = [];
                for (let bp = gb.start_bp; bp < gb.end_bp && bp < fullSeq.length; bp++) {
                  const senseBase = fullSeq[bp] || "N";
                  const displayBase = isComplement ? complementBase(senseBase) : senseBase;
                  bases.push(
                    <text
                      key={bp}
                      x={bp2x(bp) + pxPerBp / 2}
                      y={yPos + oH / 2 + 4}
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
                return bases;
              };

              return (
                <g
                  key={`gblock-asm-${i}`}
                  onClick={selectable ? (e) => {
                    e.stopPropagation();
                    selectGblock(isSelected ? null : gbKey);
                  } : undefined}
                  style={selectable ? { cursor: "pointer" } : undefined}
                >
                  {/* Sense strand bar */}
                  <rect
                    x={x} y={sY} width={w} height={oH}
                    fill={fill}
                    stroke={stroke}
                    strokeWidth={strokeW}
                    opacity={showBases ? 1 : 0.85}
                  />
                  {/* Antisense strand bar */}
                  <rect
                    x={x} y={aY} width={w} height={oH}
                    fill={fill}
                    stroke={stroke}
                    strokeWidth={strokeW}
                    opacity={showBases ? 1 : 0.85}
                  />
                  {/* Connecting side bars */}
                  <line x1={x} y1={sY + oH} x2={x} y2={aY} stroke={isSelected ? "#202124" : showBases ? GBLOCK_COLOR : "#ffffff"} strokeWidth={isSelected ? 1.5 : 1} />
                  <line x1={x + w} y1={sY + oH} x2={x + w} y2={aY} stroke={isSelected ? "#202124" : showBases ? GBLOCK_COLOR : "#ffffff"} strokeWidth={isSelected ? 1.5 : 1} />

                  {showBases ? (
                    <>
                      {renderStrandBases(sY, false)}
                      {renderStrandBases(aY, true)}
                    </>
                  ) : showNumbers && w > 30 ? (
                    <text
                      x={x + w / 2} y={sY + oH + (aY - sY - oH) / 2 + 4}
                      textAnchor="middle"
                      fontSize={11} fill="white" fontWeight="bold"
                      style={{ pointerEvents: "none" }}
                    >
                      {gb.label || "gBlock"}
                    </text>
                  ) : null}
                </g>
              );
            })}

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

            {/* Amino acid translation track — frame anchored to the insert,
                not to position 0, so the upstream flank doesn't shift the frame. */}
            {showAA && fullSeq.length - insertStart >= 3 && (() => {
              const codonTable: Record<string, string> = {
                TTT:"F",TTC:"F",TTA:"L",TTG:"L",CTT:"L",CTC:"L",CTA:"L",CTG:"L",
                ATT:"I",ATC:"I",ATA:"I",ATG:"M",GTT:"V",GTC:"V",GTA:"V",GTG:"V",
                TCT:"S",TCC:"S",TCA:"S",TCG:"S",CCT:"P",CCC:"P",CCA:"P",CCG:"P",
                ACT:"T",ACC:"T",ACA:"T",ACG:"T",GCT:"A",GCC:"A",GCA:"A",GCG:"A",
                TAT:"Y",TAC:"Y",TAA:"*",TAG:"*",CAT:"H",CAC:"H",CAA:"Q",CAG:"Q",
                AAT:"N",AAC:"N",AAA:"K",AAG:"K",GAT:"D",GAC:"D",GAA:"E",GAG:"E",
                TGT:"C",TGC:"C",TGA:"*",TGG:"W",CGT:"R",CGC:"R",CGA:"R",CGG:"R",
                AGT:"S",AGC:"S",AGA:"R",AGG:"R",GGT:"G",GGC:"G",GGA:"G",GGG:"G",
              };
              const aas = [];
              for (let i = insertStart; i + 2 < fullSeq.length && i + 2 < effectiveInsertEnd; i += 3) {
                const codon = fullSeq.slice(i, i + 3).toUpperCase();
                const aa = codonTable[codon] || "?";
                const x = bp2x(i) + bp2x(3) / 2;
                const w = bp2x(3);
                if (showBases) {
                  aas.push(
                    <text
                      key={`aa-${i}`}
                      x={x}
                      y={aaTrackY + AA_TRACK_H / 2 + 4}
                      textAnchor="middle"
                      fontSize={Math.min(pxPerBp * 2, 12)}
                      fontFamily="monospace"
                      fontWeight="bold"
                      fill="#5f6368"
                      style={{ pointerEvents: "none" }}
                    >
                      {aa}
                    </text>
                  );
                } else if (w >= 8) {
                  aas.push(
                    <text
                      key={`aa-${i}`}
                      x={x}
                      y={aaTrackY + AA_TRACK_H / 2 + 4}
                      textAnchor="middle"
                      fontSize={9}
                      fontFamily="monospace"
                      fill="#80868b"
                      style={{ pointerEvents: "none" }}
                    >
                      {aa}
                    </text>
                  );
                }
              }
              return (
                <g>
                  <line x1={0} y1={aaTrackY - 2} x2={contentWidth} y2={aaTrackY - 2} stroke="#dadce0" strokeWidth={0.5} />
                  {aas}
                </g>
              );
            })()}
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

      {/* Selected gBlock detail */}
      {selectedGblockKey && (() => {
        const gb = (gblockRegions ?? []).find(
          (g) => `${g.start_bp}-${g.end_bp}` === selectedGblockKey
        );
        if (!gb || gb.seq === undefined || gb.gc === undefined || gb.index === undefined) {
          return null;
        }
        return (
          <GBlockDetailPanel
            gblock={{
              index: gb.index,
              label: gb.label,
              start: gb.start_bp,
              end: gb.end_bp,
              length: gb.len_bp,
              gc: gb.gc,
              seq: gb.seq,
            }}
            onClose={() => setSelectedGblockKey(null)}
          />
        );
      })()}

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
                    <td className="text-[#5f6368] pr-6 py-0.5 align-top">
                      Position <span className="text-xs text-[#9aa0a6]">(nt, 1-indexed)</span>
                    </td>
                    <td className="py-0.5">{selectedIssue.start + 1}-{selectedIssue.end ?? selectedIssue.start + 1}</td>
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
