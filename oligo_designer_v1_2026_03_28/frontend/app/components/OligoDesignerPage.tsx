"use client";

import { useState, useEffect } from "react";
import OligoViewer from "./OligoViewer";
import { ISSUE_COLORS } from "./viewer/colors";

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

interface DesignResult {
  insert_length: number;
  total_length: number;
  num_oligos: number;
  num_overlaps: number;
  avg_overlap_tm: number | null;
  oligos: OligoResult[];
  overlaps: OverlapResult[];
  sequence_issues: { type: string; message: string; severity: string }[];
  report_text: string;
}

interface OligoDesignerPageProps {
  prefillSequence?: string;
  onPrefillConsumed?: () => void;
}

export default function OligoDesignerPage({ prefillSequence, onPrefillConsumed }: OligoDesignerPageProps) {
  const [sequence, setSequence] = useState("");

  useEffect(() => {
    if (prefillSequence) {
      setSequence(prefillSequence);
      onPrefillConsumed?.();
    }
  }, [prefillSequence, onPrefillConsumed]);
  const [oligoLength, setOligoLength] = useState(45);
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

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    setError("");
    setResult(null);
    setLoading(true);

    try {
      const res = await fetch(`${API_URL}/api/design`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          sequence,
          oligo_length: oligoLength,
          overlap_length: overlapLength,
          plasmid_upstream: upstream,
          plasmid_downstream: downstream,
          mv_conc: mvConc,
          dv_conc: dvConc,
          dntp_conc: dntpConc,
          dna_conc: dnaConc,
          annealing_temp: annealingTemp,
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
  }

  return (
    <div>
      <h1 className="text-3xl font-semibold mb-2 text-[#202124]">Oligo Designer</h1>
      <p className="text-[#5f6368] mb-8">
        Design overlapping oligos for gene synthesis via overlap extension assembly
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

        <button
          type="submit"
          disabled={loading || !sequence.trim()}
          className="bg-[#1a73e8] hover:bg-[#1967d2] text-white disabled:bg-[#dadce0] disabled:text-[#80868b] px-6 py-2 font-medium transition-colors cursor-pointer"
        >
          {loading ? "Designing..." : "Design Oligos"}
        </button>
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
                  <tr>
                    <td className="text-[#5f6368] pr-6 py-0.5">Avg overlap Tm</td>
                    <td className="py-0.5">
                      {result.avg_overlap_tm.toFixed(1)}°C
                      <span className="text-[#5f6368] text-xs ml-2">
                        (cross-hyb threshold: {(result.avg_overlap_tm - 20).toFixed(0)}°C)
                      </span>
                    </td>
                  </tr>
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
              {result.sequence_issues.map((issue, i) => (
                <p key={i} className="text-sm text-[#202124]">
                  {issue.severity === "error" ? "!!" : "*"} {issue.message}
                </p>
              ))}
            </div>
          )}

          {/* Oligos table */}
          <div className="bg-white border border-[#dadce0] overflow-x-auto">
            <h2 className="text-lg font-medium text-[#202124] p-4 pb-2">Oligos</h2>
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-[#5f6368] border-b border-[#dadce0]">
                  <th className="px-4 py-2">#</th>
                  <th className="px-4 py-2">Strand</th>
                  <th className="px-4 py-2">Position</th>
                  <th className="px-4 py-2">Length</th>
                  <th className="px-4 py-2">GC%</th>
                  <th className="px-4 py-2">Sequence</th>
                </tr>
              </thead>
              <tbody>
                {result.oligos.map((o) => (
                  <tr
                    key={o.index}
                    className="border-b border-[#dadce0]/50 hover:bg-[#f1f3f4]"
                  >
                    <td className="px-4 py-2">{o.index}</td>
                    <td className="px-4 py-2">
                      <span
                        className={
                          o.strand === "sense"
                            ? "text-[#1a73e8]"
                            : "text-[#e8710a]"
                        }
                      >
                        {o.strand === "sense" ? "\u2192" : "\u2190"}
                      </span>
                    </td>
                    <td className="px-4 py-2 text-[#5f6368]">
                      {o.start}-{o.end}
                    </td>
                    <td className="px-4 py-2">{o.length}</td>
                    <td className="px-4 py-2">{(o.gc * 100).toFixed(0)}%</td>
                    <td className="px-4 py-2 font-mono text-xs break-all">
                      {o.seq}
                    </td>
                  </tr>
                ))}
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
                  <th className="px-4 py-2">Position</th>
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
                      {o.start}-{o.end}
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

          {/* Raw report toggle */}
          <div className="bg-white border border-[#dadce0] p-4">
            <button
              onClick={() => setShowRaw(!showRaw)}
              className="text-sm text-[#5f6368] hover:text-[#202124] cursor-pointer"
            >
              {showRaw ? "Hide" : "Show"} raw text report
            </button>
            {showRaw && (
              <pre className="mt-3 text-xs font-mono whitespace-pre-wrap text-[#202124] max-h-96 overflow-y-auto bg-[#f8f9fa] border border-[#dadce0] p-3">
                {result.report_text}
              </pre>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
