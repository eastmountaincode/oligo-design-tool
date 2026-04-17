"use client";

import { Fragment, useState, useEffect, useRef, useCallback } from "react";
import OligoViewer from "./OligoViewer";
import { ISSUE_COLORS } from "./viewer/colors";
import { oligoLabel, overlapLabel } from "./viewer/types";
import type { RepairContext, RepairReport } from "./types";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

interface OligoResult {
  index: number;
  seq: string;
  start: number;
  end: number;
  strand: string;
  length: number;
  is_first: boolean;
  is_last: boolean;
  gc: number;
}

interface OverlapResult {
  index: number;
  seq: string;
  start: number;
  end: number;
  gc: number;
  tm: number | null;
  issues: { kind: string; message: string; severity: string; start?: number; end?: number }[];
}

interface GBlockResult {
  index: number;
  seq: string;
  start: number;
  end: number;
  length: number;
  label: string;
  gc: number;
  fragment_type: string;
}

interface ReactionResult {
  index: number;                // 0-based
  insert_start: number;         // insert-bp coords
  insert_end: number;
  left_flank: string;
  right_flank: string;
  left_is_junction: boolean;
  right_is_junction: boolean;
  num_oligos: number;
  num_overlaps: number;
  avg_overlap_tm: number | null;
  oligos: OligoResult[];
  overlaps: OverlapResult[];
  gblocks: GBlockResult[];
  repair?: RepairReport;
}

interface DesignResult {
  insert_length: number;
  total_length: number;
  num_oligos: number;
  num_overlaps: number;
  avg_overlap_tm: number | null;
  sequence_issues: { type: string; message: string; severity: string }[];
  report_text: string;
  reactions: ReactionResult[];
}

interface GBlockInfo {
  start_bp: number;
  end_bp: number;
  label: string;
  len_bp: number;
}

interface OligoDesignerPageProps {
  liveDna?: string;
  liveGblocks?: GBlockInfo[];
  // Plasmid flanks are owned by the parent so codon optimization and oligo
  // design stay in sync. The oligo designer no longer renders an input for
  // them — editing lives on the codon opt pane, where flanks actually
  // influence the optimization.
  plasmidUpstream: string;
  plasmidDownstream: string;
  // Optional: when the codon-opt pane has produced a result, it emits a
  // RepairContext that the backend needs in order to run its cross-hyb
  // repair pass. We forward it verbatim.
  repairContext?: RepairContext | null;
}

/** Compute evenly-distributed default split positions for N reactions.
 *
 * Splits are returned in insert-bp coordinates; each is a multiple of 3
 * (codon-aligned) and is nudged to the nearest valid codon outside of any
 * gBlock. Also respects a minimum distance `overlapLength` from each insert
 * edge (otherwise the junction flank has nothing to cut from).
 *
 * If n <= 1 or the insert is too short to be split, returns [].
 */
function defaultSplits(
  n: number,
  insertLen: number,
  gblocks: GBlockInfo[],
  overlapLength: number,
): number[] {
  if (n <= 1 || insertLen <= 0) return [];

  // Snap a position to the nearest codon boundary that is (a) inside the
  // admissible band [overlapLength, insertLen - overlapLength], and (b) not
  // strictly inside any gBlock. We scan outward from the ideal position
  // alternating +3/-3 until we find a valid slot.
  const lo = Math.ceil(overlapLength / 3) * 3;
  const hi = Math.floor((insertLen - overlapLength) / 3) * 3;
  const inBand = (p: number) => p >= lo && p <= hi;
  const inGblock = (p: number) =>
    gblocks.some((g) => g.start_bp < p && p < g.end_bp);
  const isValid = (p: number) => inBand(p) && !inGblock(p);

  function nearestValid(ideal: number): number | null {
    const start = Math.round(ideal / 3) * 3;
    if (isValid(start)) return start;
    for (let d = 3; d <= insertLen; d += 3) {
      const up = start + d;
      const down = start - d;
      if (isValid(up)) return up;
      if (isValid(down)) return down;
    }
    return null;
  }

  const results: number[] = [];
  for (let k = 1; k <= n - 1; k++) {
    const ideal = Math.round((insertLen * k) / n);
    const snapped = nearestValid(ideal);
    if (snapped != null && !results.includes(snapped)) {
      results.push(snapped);
    }
  }
  results.sort((a, b) => a - b);
  return results;
}

export default function OligoDesignerPage({
  liveDna,
  liveGblocks,
  plasmidUpstream,
  plasmidDownstream,
  repairContext,
}: OligoDesignerPageProps) {
  const [sequence, setSequence] = useState("");
  const [maxOligoLength, setMaxOligoLength] = useState(60);
  const [overlapLength, setOverlapLength] = useState(20);
  // Local aliases used throughout the rest of the component. These are
  // effectively read-only — the source of truth is the parent's state.
  const upstream = plasmidUpstream;
  const downstream = plasmidDownstream;
  const [mvConc, setMvConc] = useState(50.0);
  const [dvConc, setDvConc] = useState(2.0);
  const [dntpConc, setDntpConc] = useState(0.8);
  const [dnaConc, setDnaConc] = useState(250.0);
  const [annealingTemp, setAnnealingTemp] = useState(50.0);

  // Multi-reaction assembly: user-chosen split count and the resulting split
  // positions (in INSERT bp coords). When numReactions=1, splits=[] and the
  // backend tiles the whole insert as one reaction (legacy behavior). When
  // numReactions>1, cross-hybridization constraints are only enforced WITHIN
  // each reaction — separating a problematic overlap across reactions
  // relaxes it to zero. See the user's note: "split into n reactions" is
  // explicit, user-chosen, and super cheap to iterate on.
  const [numReactions, setNumReactions] = useState(1);
  const [splits, setSplits] = useState<number[]>([]);
  // Draft string for the "Split into N reactions" input so the user can
  // transiently clear the field while retyping (otherwise Number("") → 0
  // → clamped back to 1 mid-edit, which traps the caret).
  const [numReactionsDraft, setNumReactionsDraft] = useState("1");

  const [result, setResult] = useState<DesignResult | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Ref-mirror of splits so `runDesign` can read the current splits without
  // itself depending on the `splits` state. If `runDesign` took `splits` as
  // a useCallback dep, recreating it would re-trigger the liveDna effect,
  // which sets splits, which recreates runDesign, ... infinite loop.
  const splitsRef = useRef<number[]>(splits);
  useEffect(() => {
    splitsRef.current = splits;
  }, [splits]);

  const runDesign = useCallback(async (seq: string, overrideSplits?: number[]) => {
    if (!seq.trim()) return;
    setError("");
    setLoading(true);
    try {
      const res = await fetch(`${API_URL}/api/design`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence: seq,
          max_oligo_length: maxOligoLength,
          overlap_length: overlapLength,
          plasmid_upstream: upstream,
          plasmid_downstream: downstream,
          mv_conc: mvConc,
          dv_conc: dvConc,
          dntp_conc: dntpConc,
          dna_conc: dnaConc,
          annealing_temp: annealingTemp,
          gblock_regions: liveGblocks ?? [],
          repair_context: repairContext ?? null,
          splits: overrideSplits ?? splitsRef.current,
        }),
      });
      if (!res.ok) {
        const data = await res.json();
        throw new Error(data.detail || "Design failed");
      }
      setResult(await res.json());
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : "Unknown error");
    } finally {
      setLoading(false);
    }
  }, [maxOligoLength, overlapLength, upstream, downstream, mvConc, dvConc, dntpConc, dnaConc, annealingTemp, liveGblocks, repairContext]);

  // Auto-run design when liveDna or liveGblocks change (debounced). Also
  // recomputes default splits whenever the DNA or gBlock layout changes,
  // because the user's previous split positions may no longer be valid
  // (e.g. a new protein is shorter, or a gBlock now covers the old split).
  useEffect(() => {
    if (!liveDna) return;
    setSequence(liveDna);
    const freshSplits = defaultSplits(numReactions, liveDna.length, liveGblocks ?? [], overlapLength);
    setSplits(freshSplits);
    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => {
      runDesign(liveDna, freshSplits);
    }, 300);
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, [liveDna, liveGblocks, numReactions, overlapLength, runDesign]);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    runDesign(sequence);
  }

  // Change N: recompute default splits and re-run design.
  const handleNumReactionsChange = useCallback((n: number) => {
    const safeN = Math.max(1, Math.min(10, Math.round(n)));
    setNumReactions(safeN);
    const insertLen = liveDna?.length ?? sequence.length;
    const fresh = defaultSplits(safeN, insertLen, liveGblocks ?? [], overlapLength);
    setSplits(fresh);
    if (sequence.trim()) runDesign(sequence, fresh);
  }, [liveDna, sequence, liveGblocks, overlapLength, runDesign]);

  // Called by OligoViewer when the user drags a split handle and releases.
  // We receive raw INSERT-bp positions (one per split); we snap each to the
  // nearest codon boundary outside any gBlock, dedupe, sort, and re-run
  // tiling with the cleaned-up split list.
  const handleSplitsChanged = useCallback((rawSplits: number[]) => {
    const insertLen = liveDna?.length ?? sequence.length;
    if (insertLen <= 0) return;
    const lo = Math.ceil(overlapLength / 3) * 3;
    const hi = Math.floor((insertLen - overlapLength) / 3) * 3;
    const gbs = liveGblocks ?? [];
    const inBand = (p: number) => p >= lo && p <= hi;
    const inGblock = (p: number) => gbs.some((g) => g.start_bp < p && p < g.end_bp);
    const isValid = (p: number) => inBand(p) && !inGblock(p);
    function snap(bp: number): number | null {
      const start = Math.max(0, Math.min(insertLen, Math.round(bp / 3) * 3));
      if (isValid(start)) return start;
      for (let d = 3; d <= insertLen; d += 3) {
        if (isValid(start + d)) return start + d;
        if (isValid(start - d)) return start - d;
      }
      return null;
    }
    const snapped = rawSplits.map((s) => snap(s)).filter((s): s is number => s != null);
    const deduped = Array.from(new Set(snapped)).sort((a, b) => a - b);
    // Also enforce min spacing between splits — drop any split that falls
    // within overlapLength of its neighbor (the backend would reject it
    // anyway; better to show the user a sane result here).
    const enforced: number[] = [];
    for (const s of deduped) {
      if (enforced.length === 0 || s - enforced[enforced.length - 1] >= overlapLength * 2) {
        enforced.push(s);
      }
    }
    setSplits(enforced);
    if (sequence.trim()) runDesign(sequence, enforced);
  }, [liveDna, sequence, liveGblocks, overlapLength, runDesign]);

  return (
    <div>
      <h1 className="text-xl font-semibold mb-2 text-[#202124]">Oligo Designer</h1>
      <p className="text-[#5f6368] mb-8">
        Design overlapping oligos for gene synthesis
      </p>

      <form onSubmit={handleSubmit} className="space-y-4 mb-8">
        <div>
          <label className="block text-sm font-medium mb-1 text-[#202124]">
            DNA Sequence
          </label>
          <textarea
            value={sequence}
            onChange={(e) => setSequence(e.target.value)}
            placeholder="Paste DNA sequence (ATCGs only)..."
            rows={6}
            className="w-full bg-white border border-[#dadce0] p-3 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
          />
        </div>

        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div>
            <label className="block text-sm font-medium mb-1 text-[#202124]">
              Max oligo length (bp)
            </label>
            <input
              type="number"
              value={maxOligoLength}
              onChange={(e) => setMaxOligoLength(Number(e.target.value))}
              min={30}
              max={200}
              className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
            />
          </div>
          <div>
            <label className="block text-sm font-medium mb-1 text-[#202124]">
              Overlap length (bp)
            </label>
            <input
              type="number"
              value={overlapLength}
              onChange={(e) => setOverlapLength(Number(e.target.value))}
              min={15}
              max={40}
              className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
            />
          </div>
          <div>
            <label className="block text-sm font-medium mb-1 text-[#202124]">
              Split into N reactions
              <span className="ml-1 text-xs text-[#9aa0a6] font-normal">(1 = don&apos;t split)</span>
            </label>
            <input
              type="number"
              value={numReactionsDraft}
              onChange={(e) => {
                const raw = e.target.value;
                setNumReactionsDraft(raw);
                // Only commit when the draft parses to a valid integer in
                // range; otherwise let the user keep typing.
                if (raw === "") return;
                const parsed = Number(raw);
                if (!Number.isFinite(parsed)) return;
                const n = Math.round(parsed);
                if (n < 1 || n > 10) return;
                if (n !== numReactions) handleNumReactionsChange(n);
              }}
              onBlur={() => {
                // On blur, snap the draft back to the committed value so
                // empty / out-of-range drafts don't linger.
                setNumReactionsDraft(String(numReactions));
              }}
              min={1}
              max={10}
              className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
            />
            {numReactions > 1 && splits.length > 0 && (
              <div className="text-xs text-[#5f6368] mt-1 font-mono">
                splits @ {splits.map((s) => s).join(", ")}
              </div>
            )}
          </div>
        </div>

        {(upstream || downstream) && (
          <div className="text-xs text-[#5f6368] border-l-2 border-[#dadce0] pl-3">
            <div className="mb-1">
              Plasmid flanks in use (edit on the codon optimization pane):
            </div>
            <div className="font-mono text-[#202124]">
              5&prime; {upstream || <span className="text-[#9aa0a6]">—</span>}
              {" · "}
              3&prime; {downstream || <span className="text-[#9aa0a6]">—</span>}
            </div>
          </div>
        )}

        <details className="text-sm">
          <summary className="cursor-pointer text-[#5f6368] hover:text-[#202124]">
            Reaction conditions
          </summary>
          <div className="grid grid-cols-2 md:grid-cols-5 gap-4 mt-2">
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Na⁺/K⁺ (mM)
              </label>
              <input
                type="number"
                value={mvConc}
                onChange={(e) => setMvConc(Number(e.target.value))}
                step={1}
                min={0}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Mg²⁺ (mM)
              </label>
              <input
                type="number"
                value={dvConc}
                onChange={(e) => setDvConc(Number(e.target.value))}
                step={0.1}
                min={0}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                dNTPs (mM)
              </label>
              <input
                type="number"
                value={dntpConc}
                onChange={(e) => setDntpConc(Number(e.target.value))}
                step={0.1}
                min={0}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Oligo conc (nM)
              </label>
              <input
                type="number"
                value={dnaConc}
                onChange={(e) => setDnaConc(Number(e.target.value))}
                step={10}
                min={0}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Anneal temp (°C)
              </label>
              <input
                type="number"
                value={annealingTemp}
                onChange={(e) => setAnnealingTemp(Number(e.target.value))}
                step={1}
                min={0}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
          </div>
        </details>

        <div className="flex items-center gap-3">
          <button
            type="submit"
            disabled={loading || !sequence.trim()}
            className="bg-[#1a73e8] hover:bg-[#1967d2] text-white disabled:bg-[#dadce0] disabled:text-[#80868b] px-6 py-2 font-medium transition-colors cursor-pointer"
          >
            {loading ? "Designing..." : "Design Oligos"}
          </button>
        </div>
      </form>

      {error && (
        <div className="bg-[#fce8e6] border border-[#d93025] p-4 mb-6">
          <p className="text-[#d93025]">{error}</p>
        </div>
      )}

      {result && (() => {
        // Flatten per-reaction arrays for components that still render a
        // single merged view (assembly viewer, issue legend, summary stats).
        // Each flattened item carries the rxnIndex it came from so the UI
        // can show R1/R2/... badges.
        const reactions = result.reactions ?? [];
        const flatOligos = reactions.flatMap((r) => r.oligos.map((o) => ({ ...o, rxnIndex: r.index })));
        const flatOverlaps = reactions.flatMap((r) => r.overlaps.map((o) => ({ ...o, rxnIndex: r.index })));
        const flatGblocks = reactions.flatMap((r) => r.gblocks.map((g) => ({ ...g, rxnIndex: r.index })));
        const multiReaction = reactions.length > 1;
        const rxnBadge = (idx: number) => `R${idx + 1}`;

        return (
        <div className="space-y-6">
          {/* Summary */}
          <div className="bg-white p-4 border border-[#dadce0]">
            <h2 className="text-lg font-medium mb-3 text-[#202124]">Summary</h2>
            <table className="text-sm">
              <tbody>
                <tr>
                  <td className="text-[#5f6368] pr-6 py-0.5">Insert length</td>
                  <td className="py-0.5">{result.insert_length} bp</td>
                </tr>
                <tr>
                  <td className="text-[#5f6368] pr-6 py-0.5">Total construct length</td>
                  <td className="py-0.5">{result.total_length} bp</td>
                </tr>
                {multiReaction && (
                  <tr>
                    <td className="text-[#5f6368] pr-6 py-0.5">Reactions</td>
                    <td className="py-0.5">
                      {reactions.length}
                      <span className="text-[#5f6368] text-xs ml-2 font-mono">
                        {reactions.map((r) => `R${r.index + 1}: insert ${r.insert_start + 1}-${r.insert_end} (${r.num_oligos}o+${r.gblocks.length}g)`).join("  ·  ")}
                      </span>
                    </td>
                  </tr>
                )}
                <tr>
                  <td className="text-[#5f6368] pr-6 py-0.5">Number of oligos</td>
                  <td className="py-0.5">{result.num_oligos}</td>
                </tr>
                <tr>
                  <td className="text-[#5f6368] pr-6 py-0.5">Number of overlaps</td>
                  <td className="py-0.5">{result.num_overlaps}</td>
                </tr>
                {result.avg_overlap_tm != null && (
                  <>
                    <tr>
                      <td className="text-[#5f6368] pr-6 py-0.5">Overlap Tm range</td>
                      <td className="py-0.5">
                        {(() => {
                          const tms = flatOverlaps.map((o) => o.tm).filter((t): t is number => t != null);
                          if (tms.length === 0) return "—";
                          const min = Math.min(...tms);
                          const max = Math.max(...tms);
                          return `${min.toFixed(1)}°C – ${max.toFixed(1)}°C`;
                        })()}
                      </td>
                    </tr>
                    <tr>
                      <td className="text-[#5f6368] pr-6 py-0.5">Avg overlap Tm</td>
                      <td className="py-0.5">
                        {result.avg_overlap_tm.toFixed(1)}°C
                        <span className="text-[#5f6368] text-xs ml-2">
                          (cross-hyb threshold: {(result.avg_overlap_tm - 20).toFixed(0)}°C)
                        </span>
                      </td>
                    </tr>
                  </>
                )}
              </tbody>
            </table>
          </div>

          {/* Cross-hybridization repair log — one block per reaction that
              actually ran the repair pass. When N=1 the rendering is
              identical to before; when N>1 each reaction's repair log gets
              its own "R{n}: Cross-hybridization repair" panel so you can
              see which reaction needed which swaps. */}
          {reactions.map((r) =>
            r.repair && r.repair.ran ? (
              <RepairLogPanel
                key={`repair-${r.index}`}
                repair={r.repair}
                reactionLabel={multiReaction ? rxnBadge(r.index) : undefined}
              />
            ) : null,
          )}

          {/* IGV-style assembly viewer */}
          <OligoViewer
            oligos={flatOligos}
            overlaps={flatOverlaps}
            totalLength={result.total_length}
            sequenceIssues={result.sequence_issues}
            gblockRegions={flatGblocks.map((g) => ({
              start_bp: g.start,
              end_bp: g.end,
              label: g.label,
              len_bp: g.length,
              index: g.index,
              seq: g.seq,
              gc: g.gc,
            }))}
            fullSequence={(upstream || "").toUpperCase() + sequence + (downstream || "").toUpperCase()}
            insertStart={(upstream || "").length}
            insertEnd={(upstream || "").length + sequence.length}
            splits={splits}
            onSplitsChanged={handleSplitsChanged}
          />

          {/* Issue legend */}
          {(() => {
            const allIssueTypes = [
              ...result.sequence_issues.map((i) => i.type),
              ...flatOverlaps.flatMap((o) => o.issues.map((i) => i.kind)),
            ];
            const seenLabels = new Set<string>();
            const entries: { color: string; label: string }[] = [];
            for (const t of allIssueTypes) {
              const info = ISSUE_COLORS[t];
              if (!info || seenLabels.has(info.label)) continue;
              seenLabels.add(info.label);
              entries.push({ color: info.fill, label: info.label });
            }
            if (entries.length === 0) return null;
            return (
              <div className="bg-white border border-[#dadce0] px-4 py-3 text-xs text-[#5f6368]">
                <span className="font-medium text-[#202124]">Issue Color Legend</span>
                <div className="mt-2 space-y-1.5">
                  {entries.map((e) => (
                    <div key={e.label} className="flex items-center gap-1.5">
                      <span
                        className="inline-block w-3 h-3"
                        style={{ backgroundColor: e.color }}
                      />
                      <span>{e.label}</span>
                    </div>
                  ))}
                </div>
              </div>
            );
          })()}

          {/* Sequence warnings */}
          {result.sequence_issues.length > 0 && (
            <div className="bg-[#fef7e0] border border-[#f9ab00] p-4">
              <h2 className="text-lg font-medium mb-2 text-[#e37400]">
                Sequence Warnings
              </h2>
              <ul className="space-y-1">
                {result.sequence_issues.map((issue, i) => {
                  const color =
                    ISSUE_COLORS[issue.type.replace(/^overlap_/, "")]?.fill ??
                    "#5f6368";
                  return (
                    <li
                      key={i}
                      className="text-sm text-[#202124] flex items-start gap-2"
                    >
                      <span
                        className="inline-block w-2.5 h-2.5 mt-1.5 rounded-sm flex-shrink-0"
                        style={{ backgroundColor: color }}
                        aria-hidden
                      />
                      <span className="flex-1">
                        {issue.message}
                        {issue.severity === "error" && (
                          <span className="ml-2 text-xs uppercase tracking-wide text-[#c5221f] font-medium">
                            severe
                          </span>
                        )}
                      </span>
                    </li>
                  );
                })}
              </ul>
            </div>
          )}

          {/* Fragments table (oligos + gBlocks). When there are multiple
              reactions, rows are grouped under a "Reaction N" header row
              so it's clear which oligos go into which tube. gBlock labels
              stay globally stable (G1, G2, ...) via the index emitted by
              the backend. */}
          <div className="bg-white border border-[#dadce0] overflow-x-auto">
            <h2 className="text-lg font-medium text-[#202124] p-4 pb-2">
              Fragments
              <span className="text-sm font-normal text-[#5f6368] ml-2">
                {result.num_oligos} oligo{result.num_oligos !== 1 ? "s" : ""}
                {flatGblocks.length > 0 && ` + ${flatGblocks.length} gBlock${flatGblocks.length !== 1 ? "s" : ""}`}
              </span>
            </h2>
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-[#5f6368] border-b border-[#dadce0]">
                  <th className="px-4 py-2">#</th>
                  <th className="px-4 py-2">Type</th>
                  <th className="px-4 py-2">
                    Position
                    <span className="ml-1 text-xs text-[#9aa0a6] font-normal">(nt, 1-indexed)</span>
                  </th>
                  <th className="px-4 py-2">Length</th>
                  <th className="px-4 py-2">GC%</th>
                  <th className="px-4 py-2">Sequence</th>
                </tr>
              </thead>
              <tbody>
                {reactions.map((r) => {
                  const items = [
                    ...r.oligos.map((o) => ({ kind: "oligo" as const, pos: o.start, data: o })),
                    ...r.gblocks.map((g) => ({ kind: "gblock" as const, pos: g.start, data: g })),
                  ].sort((a, b) => a.pos - b.pos);
                  return (
                    <Fragment key={`rxn-frag-${r.index}`}>
                      {multiReaction && (
                        <tr className="bg-[#f8f9fa]">
                          <td colSpan={6} className="px-4 py-1.5 text-xs font-medium text-[#5f6368] uppercase tracking-wide">
                            Reaction {r.index + 1}
                            <span className="ml-2 text-[#9aa0a6] normal-case tracking-normal font-normal">
                              insert {r.insert_start + 1}–{r.insert_end}
                              {" · "}
                              {r.left_is_junction ? "junction" : "plasmid"} 5&prime; flank / {r.right_is_junction ? "junction" : "plasmid"} 3&prime; flank
                            </span>
                          </td>
                        </tr>
                      )}
                      {items.map((item) => {
                        if (item.kind === "gblock") {
                          const g = item.data as GBlockResult;
                          return (
                            <tr
                              key={`gb-${g.index}`}
                              className="border-b border-[#dadce0]/50 hover:bg-[#e8f5e9]"
                              style={{ backgroundColor: "rgba(52, 168, 83, 0.06)" }}
                            >
                              <td className="px-4 py-2">G{g.index + 1}</td>
                              <td className="px-4 py-2">
                                <span className="text-[#34a853] font-medium">gBlock</span>
                                {multiReaction && (
                                  <span className="ml-2 text-[10px] bg-[#e8f0fe] text-[#1967d2] px-1.5 py-0.5 rounded">
                                    {rxnBadge(r.index)}
                                  </span>
                                )}
                              </td>
                              <td className="px-4 py-2 text-[#5f6368]">
                                {g.start + 1}-{g.end}
                              </td>
                              <td className="px-4 py-2">{g.length}</td>
                              <td className="px-4 py-2">{(g.gc * 100).toFixed(0)}%</td>
                              <td className="px-4 py-2 font-mono text-xs break-all">{g.seq}</td>
                            </tr>
                          );
                        }
                        const o = item.data as OligoResult;
                        return (
                          <tr
                            key={`ol-${r.index}-${o.index}`}
                            className="border-b border-[#dadce0]/50 hover:bg-[#f1f3f4]"
                          >
                            <td className="px-4 py-2">
                              {oligoLabel({ index: o.index, rxnIndex: multiReaction ? r.index : undefined })}
                            </td>
                            <td className="px-4 py-2">
                              <span
                                className={
                                  o.strand === "sense"
                                    ? "text-[#1a73e8]"
                                    : "text-[#e8710a]"
                                }
                              >
                                {o.strand === "sense" ? "5\u2032 \u2192 3\u2032" : "3\u2032 \u2190 5\u2032"}
                              </span>
                              {multiReaction && (
                                <span className="ml-2 text-[10px] bg-[#e8f0fe] text-[#1967d2] px-1.5 py-0.5 rounded">
                                  {rxnBadge(r.index)}
                                </span>
                              )}
                            </td>
                            <td className="px-4 py-2 text-[#5f6368]">
                              {o.start + 1}-{o.end}
                            </td>
                            <td className="px-4 py-2">{o.length}</td>
                            <td className="px-4 py-2">{(o.gc * 100).toFixed(0)}%</td>
                            <td className="px-4 py-2 font-mono text-xs break-all">{o.seq}</td>
                          </tr>
                        );
                      })}
                    </Fragment>
                  );
                })}
              </tbody>
            </table>
          </div>

          {/* Overlaps */}
          <div className="bg-white border border-[#dadce0] overflow-x-auto">
            <h2 className="text-lg font-medium p-4 pb-2 text-[#202124]">Overlaps</h2>
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-[#5f6368] border-b border-[#dadce0]">
                  <th className="px-4 py-2">#</th>
                  <th className="px-4 py-2">
                    Position
                    <span className="ml-1 text-xs text-[#9aa0a6] font-normal">(nt, 1-indexed)</span>
                  </th>
                  <th className="px-4 py-2">GC%</th>
                  <th className="px-4 py-2">Tm</th>
                  <th className="px-4 py-2">Sequence</th>
                  <th className="px-4 py-2">Issues</th>
                </tr>
              </thead>
              <tbody>
                {reactions.map((r) => (
                  <Fragment key={`ovl-frag-${r.index}`}>
                    {multiReaction && (
                      <tr className="bg-[#f8f9fa]">
                        <td colSpan={6} className="px-4 py-1.5 text-xs font-medium text-[#5f6368] uppercase tracking-wide">
                          Reaction {r.index + 1}
                          {r.avg_overlap_tm != null && (
                            <span className="ml-2 text-[#9aa0a6] normal-case tracking-normal font-normal">
                              avg Tm {r.avg_overlap_tm.toFixed(1)}°C
                            </span>
                          )}
                        </td>
                      </tr>
                    )}
                    {r.overlaps.map((o) => (
                      <tr
                        key={`ovl-${r.index}-${o.index}`}
                        className={`border-b border-[#dadce0]/50 ${
                          o.issues.length > 0
                            ? "bg-[#fce8e6]"
                            : "hover:bg-[#f1f3f4]"
                        }`}
                      >
                        <td className="px-4 py-2">
                          {overlapLabel({ index: o.index, rxnIndex: multiReaction ? r.index : undefined })}
                        </td>
                        <td className="px-4 py-2 text-[#5f6368]">
                          {o.start + 1}-{o.end}
                        </td>
                        <td className="px-4 py-2">{(o.gc * 100).toFixed(0)}%</td>
                        <td className="px-4 py-2">
                          {o.tm != null ? `${o.tm.toFixed(1)}°C` : "—"}
                        </td>
                        <td className="px-4 py-2 font-mono text-xs">{o.seq}</td>
                        <td className="px-4 py-2">
                          {o.issues.length === 0 ? (
                            <span className="text-[#188038]">OK</span>
                          ) : (
                            o.issues.map((iss, i) => (
                              <p
                                key={i}
                                className={
                                  iss.severity === "error"
                                    ? "text-[#d93025]"
                                    : "text-[#e37400]"
                                }
                              >
                                {iss.message}
                              </p>
                            ))
                          )}
                        </td>
                      </tr>
                    ))}
                  </Fragment>
                ))}
              </tbody>
            </table>
          </div>

          {/* Tm equation reference */}
          <details className="bg-white border border-[#dadce0]">
            <summary className="p-4 cursor-pointer text-sm text-[#5f6368] hover:text-[#202124]">
              Tm calculation method
            </summary>
            <div className="px-4 pb-4 text-sm text-[#202124]">
              <p className="mb-3">
                Overlap melting temperatures are computed using the <strong>SantaLucia (1998) nearest-neighbor model</strong> for
                base-stacking thermodynamics, with the <strong>Owczarzy (2008) salt correction</strong> applied for Mg²⁺,
                via the primer3 library (<code className="text-xs bg-[#f8f9fa] px-1">calc_tm</code>).
              </p>

              <div className="bg-[#f8f9fa] border border-[#dadce0] p-3 font-mono text-xs mb-3">
                <div className="mb-2 text-[#202124]">
                  Tm = <span className="text-[#1a73e8]">ΔH</span> / (<span className="text-[#1a73e8]">ΔS</span> + R × ln(C<sub>t</sub>/4)) − 273.15
                </div>
                <div className="text-[#5f6368] space-y-1">
                  <div><span className="text-[#1a73e8]">ΔH</span> = sum of nearest-neighbor enthalpies + initiation terms (kcal/mol)</div>
                  <div><span className="text-[#1a73e8]">ΔS</span> = sum of nearest-neighbor entropies + initiation + salt correction (cal/mol·K)</div>
                  <div>R = 1.987 cal/mol·K (gas constant)</div>
                  <div>C<sub>t</sub> = total strand concentration</div>
                </div>
              </div>

              <p className="mb-2 text-[#5f6368] text-xs">
                The nearest-neighbor model breaks the sequence into overlapping dinucleotide pairs
                (e.g., AT, TG, GC, ...) and looks up the stacking enthalpy (ΔH) and entropy (ΔS) for each pair
                from experimentally determined tables. This captures how adjacent bases interact —
                e.g., a GC pair next to another GC pair is more stable than a GC pair next to AT.
              </p>

              <table className="text-xs mb-3">
                <tbody>
                  <tr>
                    <td className="text-[#5f6368] pr-4 py-0.5">Na⁺/K⁺ concentration</td>
                    <td className="font-mono text-[#1a73e8]">{mvConc} mM</td>
                  </tr>
                  <tr>
                    <td className="text-[#5f6368] pr-4 py-0.5">Mg²⁺ concentration</td>
                    <td className="font-mono text-[#1a73e8]">{dvConc} mM</td>
                  </tr>
                  <tr>
                    <td className="text-[#5f6368] pr-4 py-0.5">dNTP concentration</td>
                    <td className="font-mono text-[#1a73e8]">{dntpConc} mM</td>
                  </tr>
                  <tr>
                    <td className="text-[#5f6368] pr-4 py-0.5">Oligo concentration</td>
                    <td className="font-mono text-[#1a73e8]">{dnaConc} nM</td>
                  </tr>
                  <tr>
                    <td className="text-[#5f6368] pr-4 py-0.5">Annealing temperature</td>
                    <td className="font-mono text-[#1a73e8]">{annealingTemp}°C</td>
                  </tr>
                </tbody>
              </table>

              <p className="text-[#5f6368] text-xs mb-2">
                Salt correction: ΔS is adjusted for ionic strength using the <strong>Owczarzy 2008</strong> empirical
                model (Owczarzy et al., <em>Biochemistry</em> 47:5336), which captures Mg²⁺ effects more accurately
                than the SantaLucia Na⁺-only correction. This matters because assembly buffers contain Mg²⁺
                and dNTPs that bind Mg²⁺, both of which shift the effective ionic strength.
              </p>

              <p className="text-[#5f6368] text-xs">
                Cross-hybridization Tm is computed via <code className="bg-[#f8f9fa] px-1">calc_heterodimer</code>,
                which finds the best possible alignment between two sequences (including mismatches and bulges)
                and reports the Tm of that interaction. Threshold: avg overlap Tm − 20°C
                (Rouillard et al. 2004, Nucleic Acids Res. 32:W176).
              </p>
            </div>
          </details>

          {/* Raw report toggle */}
          <details className="bg-white border border-[#dadce0]">
            <summary className="p-4 cursor-pointer text-sm text-[#5f6368] hover:text-[#202124]">
              Raw text report
            </summary>
            <div className="px-4 pb-4">
              <pre className="text-xs font-mono whitespace-pre-wrap text-[#202124] max-h-96 overflow-y-auto bg-[#f8f9fa] border border-[#dadce0] p-3">
                {result.report_text}
              </pre>
            </div>
          </details>
        </div>
        );
      })()}
    </div>
  );
}


/** Repair log — shown when the backend ran its cross-hybridization repair
 *  pass. Each "applied" row is a synonymous codon swap the backend performed
 *  to push one or more overlap pairs below the heterodimer-Tm threshold;
 *  "no_move_found" rows describe why a repair step gave up before resolving
 *  every flagged pair. Row order reflects the iteration order of the
 *  repair loop so users can follow the chain of edits.
 */
function RepairLogPanel({
  repair,
  reactionLabel,
}: {
  repair: RepairReport;
  reactionLabel?: string;
}) {
  const applied = repair.log.filter((e) => e.kind === "applied");
  const failed = repair.log.filter((e) => e.kind !== "applied");
  const allResolved = repair.issues_after === 0;
  return (
    <div
      className={
        "border p-4 " +
        (allResolved
          ? "bg-[#e6f4ea] border-[#34a853]"
          : "bg-[#fef7e0] border-[#f9ab00]")
      }
    >
      <h2
        className={
          "text-lg font-medium mb-2 " +
          (allResolved ? "text-[#188038]" : "text-[#e37400]")
        }
      >
        {reactionLabel ? `${reactionLabel}: ` : ""}Cross-hybridization repair
      </h2>
      <p className="text-sm text-[#202124] mb-3">
        {allResolved ? (
          <>
            Resolved <strong>{repair.issues_before}</strong> flagged pair
            {repair.issues_before !== 1 ? "s" : ""} using{" "}
            <strong>{applied.length}</strong> synonymous codon swap
            {applied.length !== 1 ? "s" : ""}.
          </>
        ) : (
          <>
            Started with <strong>{repair.issues_before}</strong> flagged pair
            {repair.issues_before !== 1 ? "s" : ""}; resolved{" "}
            <strong>{repair.issues_before - repair.issues_after}</strong>,{" "}
            <strong>{repair.issues_after}</strong> remain.
          </>
        )}
        {repair.threshold_tm != null && (
          <span className="text-[#5f6368] ml-2">
            (threshold: {repair.threshold_tm.toFixed(1)}°C)
          </span>
        )}
      </p>

      {applied.length > 0 && (
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="text-left text-[#5f6368] border-b border-[#dadce0]">
                <th className="px-2 py-1.5">Step</th>
                <th className="px-2 py-1.5">
                  Codon
                  <span className="ml-1 text-xs text-[#9aa0a6] font-normal">(1-indexed)</span>
                </th>
                <th className="px-2 py-1.5">AA</th>
                <th className="px-2 py-1.5">Swap</th>
                <th className="px-2 py-1.5">Freq (/1000)</th>
                <th className="px-2 py-1.5">GC</th>
                <th className="px-2 py-1.5">Fixed pairs</th>
              </tr>
            </thead>
            <tbody>
              {applied.map((e) => (
                <tr
                  key={`step-${e.iteration}`}
                  className="border-b border-[#dadce0]/50"
                >
                  <td className="px-2 py-1.5 font-medium">{e.iteration}</td>
                  <td className="px-2 py-1.5">
                    #{(e.codon_index ?? 0) + 1}
                    <span className="text-[#9aa0a6] text-xs ml-1">
                      (nt {(e.insert_position ?? 0) + 1}–
                      {(e.insert_position ?? 0) + 3})
                    </span>
                  </td>
                  <td className="px-2 py-1.5">{e.amino_acid ?? "—"}</td>
                  <td className="px-2 py-1.5 font-mono">
                    {e.old_codon} → {e.new_codon}
                  </td>
                  <td className="px-2 py-1.5 text-xs">
                    {e.old_freq_per_k != null ? e.old_freq_per_k.toFixed(0) : "—"}
                    {" → "}
                    {e.new_freq_per_k != null ? e.new_freq_per_k.toFixed(0) : "—"}
                  </td>
                  <td className="px-2 py-1.5 text-xs">
                    {e.gc_before != null ? (e.gc_before * 100).toFixed(1) + "%" : "—"}
                    {" → "}
                    {e.gc_after != null ? (e.gc_after * 100).toFixed(1) + "%" : "—"}
                  </td>
                  <td className="px-2 py-1.5 text-xs font-mono">
                    {e.fixed_pairs.length === 0
                      ? "—"
                      : e.fixed_pairs.map((p) => `#${p[0]}↔#${p[1]}`).join(", ")}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {failed.length > 0 && (
        <div className="mt-3 text-sm text-[#5f6368]">
          {failed.map((e, i) => (
            <div key={`fail-${i}`} className="mb-1">
              <strong className="text-[#d93025]">Step {e.iteration}:</strong>{" "}
              {e.notes || "No eligible swap found."}
              {e.remaining_pairs.length > 0 && (
                <span className="ml-1 font-mono">
                  Remaining: {e.remaining_pairs.map((p) => `#${p[0]}↔#${p[1]}`).join(", ")}
                  {e.remaining_max_tm != null && ` (max Tm ${e.remaining_max_tm.toFixed(1)}°C)`}
                </span>
              )}
            </div>
          ))}
        </div>
      )}

      <p className="text-xs text-[#5f6368] mt-3">
        The repair pass runs after codon optimization and oligo tiling. It
        searches for synonymous codon swaps inside flagged overlap regions
        that drop the heterodimer Tm below the threshold without breaking any
        other constraint (GC window, homopolymer limits, avoid patterns, or
        the minimum codon frequency). Ranking prefers swaps that change GC
        content the least and preserve high-frequency codons.
      </p>
    </div>
  );
}
