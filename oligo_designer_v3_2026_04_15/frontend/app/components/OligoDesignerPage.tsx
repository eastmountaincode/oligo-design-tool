"use client";

import { useState, useEffect, useRef, useCallback } from "react";
import OligoViewer from "./OligoViewer";
import { ISSUE_COLORS } from "./viewer/colors";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

// Example plasmid flanks for the pcDNA3.1(+) MCS, between HindIII and XhoI
// (the most common insertion window for mammalian expression). The first
// oligo will overlap into HindIII; the last oligo into XhoI — the Gibson
// overlap case Slim described in the 2026-04-07 meeting.
const PCDNA3_UPSTREAM = "AAGCTTGGTACCGAGCTCGG";
const PCDNA3_DOWNSTREAM = "CTCGAGTCTAGAGGGCCCGT";

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

interface DesignResult {
  insert_length: number;
  total_length: number;
  num_oligos: number;
  num_overlaps: number;
  avg_overlap_tm: number | null;
  oligos: OligoResult[];
  overlaps: OverlapResult[];
  gblocks: GBlockResult[];
  sequence_issues: { type: string; message: string; severity: string }[];
  report_text: string;
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
}

export default function OligoDesignerPage({ liveDna, liveGblocks }: OligoDesignerPageProps) {
  const [sequence, setSequence] = useState("");
  const [oligoLength, setOligoLength] = useState(60);
  const [overlapLength, setOverlapLength] = useState(20);
  const [upstream, setUpstream] = useState("");
  const [downstream, setDownstream] = useState("");
  const [mvConc, setMvConc] = useState(50.0);
  const [dvConc, setDvConc] = useState(2.0);
  const [dntpConc, setDntpConc] = useState(0.8);
  const [dnaConc, setDnaConc] = useState(250.0);
  const [annealingTemp, setAnnealingTemp] = useState(50.0);
  const [result, setResult] = useState<DesignResult | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [showRaw, setShowRaw] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const runDesign = useCallback(async (seq: string) => {
    if (!seq.trim()) return;
    setError("");
    setLoading(true);
    try {
      const res = await fetch(`${API_URL}/api/design`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence: seq,
          oligo_length: oligoLength,
          overlap_length: overlapLength,
          plasmid_upstream: upstream,
          plasmid_downstream: downstream,
          mv_conc: mvConc,
          dv_conc: dvConc,
          dntp_conc: dntpConc,
          dna_conc: dnaConc,
          annealing_temp: annealingTemp,
          gblock_regions: liveGblocks ?? [],
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
  }, [oligoLength, overlapLength, upstream, downstream, mvConc, dvConc, dntpConc, dnaConc, annealingTemp, liveGblocks]);

  // Auto-run design when liveDna or liveGblocks change (debounced)
  useEffect(() => {
    if (!liveDna) return;
    setSequence(liveDna);
    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => {
      runDesign(liveDna);
    }, 300);
    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, [liveDna, liveGblocks, runDesign]);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    runDesign(sequence);
  }

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
              Oligo length (bp)
            </label>
            <input
              type="number"
              value={oligoLength}
              onChange={(e) => setOligoLength(Number(e.target.value))}
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
        </div>

        <details className="text-sm">
          <summary className="cursor-pointer text-[#5f6368] hover:text-[#202124]">
            Plasmid flanking sequences (optional)
          </summary>
          <div className="flex items-center gap-3 text-xs mt-2">
            <span className="text-[#5f6368]">Examples:</span>
            <button
              type="button"
              onClick={() => {
                setUpstream(PCDNA3_UPSTREAM);
                setDownstream(PCDNA3_DOWNSTREAM);
              }}
              className="text-[#1a73e8] hover:text-[#1967d2] cursor-pointer"
            >
              pcDNA3.1 (HindIII / XhoI)
            </button>
            <button
              type="button"
              onClick={() => {
                setUpstream("");
                setDownstream("");
              }}
              className="text-[#5f6368] hover:text-[#202124] cursor-pointer"
            >
              clear
            </button>
          </div>
          <div className="grid grid-cols-2 gap-4 mt-2">
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Upstream flank (~20bp)
              </label>
              <input
                type="text"
                value={upstream}
                onChange={(e) => setUpstream(e.target.value)}
                placeholder="e.g. GCTAGCTAGCTAGCTAGCTA"
                className="w-full bg-white border border-[#dadce0] p-2 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Downstream flank (~20bp)
              </label>
              <input
                type="text"
                value={downstream}
                onChange={(e) => setDownstream(e.target.value)}
                placeholder="e.g. TTAAGCTTGCATGCCTGCAG"
                className="w-full bg-white border border-[#dadce0] p-2 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
          </div>
        </details>

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

      {result && (
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
                          const tms = result.overlaps.map((o) => o.tm).filter((t): t is number => t != null);
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

          {/* IGV-style assembly viewer */}
          <OligoViewer
            oligos={result.oligos}
            overlaps={result.overlaps}
            totalLength={result.total_length}
            sequenceIssues={result.sequence_issues}
            gblockRegions={(result.gblocks ?? []).map(g => ({
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
          />

          {/* Issue legend */}
          {(() => {
            const allIssueTypes = [
              ...result.sequence_issues.map((i) => i.type),
              ...result.overlaps.flatMap((o) => o.issues.map((i) => i.kind)),
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

          {/* Fragments table (oligos + gBlocks) */}
          <div className="bg-white border border-[#dadce0] overflow-x-auto">
            <h2 className="text-lg font-medium text-[#202124] p-4 pb-2">
              Fragments
              <span className="text-sm font-normal text-[#5f6368] ml-2">
                {result.oligos.length} oligo{result.oligos.length !== 1 ? "s" : ""}
                {result.gblocks?.length > 0 && ` + ${result.gblocks.length} gBlock${result.gblocks.length !== 1 ? "s" : ""}`}
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
                {/* Merge oligos and gBlocks into one sorted-by-position list */}
                {[
                  ...result.oligos.map((o) => ({ kind: "oligo" as const, pos: o.start, data: o })),
                  ...(result.gblocks ?? []).map((g) => ({ kind: "gblock" as const, pos: g.start, data: g })),
                ]
                  .sort((a, b) => a.pos - b.pos)
                  .map((item, i) => {
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
                          </td>
                          <td className="px-4 py-2 text-[#5f6368]">
                            {g.start + 1}-{g.end}
                          </td>
                          <td className="px-4 py-2">{g.length}</td>
                          <td className="px-4 py-2">{(g.gc * 100).toFixed(0)}%</td>
                          <td className="px-4 py-2 font-mono text-xs break-all">
                            {g.seq}
                          </td>
                        </tr>
                      );
                    }
                    const o = item.data as OligoResult;
                    return (
                      <tr
                        key={`ol-${o.index}`}
                        className="border-b border-[#dadce0]/50 hover:bg-[#f1f3f4]"
                      >
                        <td className="px-4 py-2">{o.index + 1}</td>
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
                        </td>
                        <td className="px-4 py-2 text-[#5f6368]">
                          {o.start + 1}-{o.end}
                        </td>
                        <td className="px-4 py-2">{o.length}</td>
                        <td className="px-4 py-2">{(o.gc * 100).toFixed(0)}%</td>
                        <td className="px-4 py-2 font-mono text-xs break-all">
                          {o.seq}
                        </td>
                      </tr>
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
                {result.overlaps.map((o) => (
                  <tr
                    key={o.index}
                    className={`border-b border-[#dadce0]/50 ${
                      o.issues.length > 0
                        ? "bg-[#fce8e6]"
                        : "hover:bg-[#f1f3f4]"
                    }`}
                  >
                    <td className="px-4 py-2">{o.index}</td>
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
      )}
    </div>
  );
}
