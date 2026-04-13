"use client";

import { useState, useEffect, useCallback } from "react";
import InfoTip from "./InfoTip";
import CodonTableView from "./CodonTableView";
import CodonTrackViewer from "./CodonTrackViewer";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

// Number input that allows a transient empty string while editing — so the
// user can clear the field and type a new value without it snapping to 0.
function NumericInput({
  value,
  onChange,
  className,
  ...rest
}: {
  value: number;
  onChange: (v: number) => void;
  className?: string;
} & Omit<React.InputHTMLAttributes<HTMLInputElement>, "value" | "onChange" | "type">) {
  const [draft, setDraft] = useState<string | null>(null);
  return (
    <input
      type="number"
      value={draft ?? value}
      className={className}
      onChange={(e) => {
        const v = e.target.value;
        setDraft(v);
        if (v !== "" && !Number.isNaN(Number(v))) {
          onChange(Number(v));
        }
      }}
      onBlur={() => {
        if (draft === "" || draft === null) {
          setDraft(null);
          return;
        }
        const parsed = Number(draft);
        if (!Number.isNaN(parsed)) onChange(parsed);
        setDraft(null);
      }}
      {...rest}
    />
  );
}

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

interface CodonOptResult {
  input_protein: string;
  optimized_dna: string;
  back_translated_protein: string;
  length_dna: number;
  length_protein: number;
  gc_content: number;
  gc_window_min?: number | null;
  gc_window_max?: number | null;
  gc_window_size?: number;
  constraints_pass: boolean;
  constraints_summary: { constraint: string; passing: boolean; message: string; enforced?: boolean }[];
  warnings: SequenceWarning[];
  codons: CodonDetail[];
  codon_table: CodonTableEntry[];
  codon_table_name: string;
  optimization_report?: {
    ideal_score: number;
    chosen_score: number;
    total_loss: number;
    codons_changed: number;
    source?: string;
    chisel_score?: number;
    beam_score?: number;
    chisel_loss?: number;
    beam_loss?: number;
    beam_k?: number;
    min_frequency_achieved?: number;
  };
}

interface CodonOptPageProps {
  onDnaChanged?: (dna: string) => void;
  prefillDna?: string;
  onPrefillConsumed?: () => void;
}

export default function CodonOptPage({ onDnaChanged, prefillDna, onPrefillConsumed }: CodonOptPageProps) {
  const [protein, setProtein] = useState("");
  const [species, setSpecies] = useState("h_sapiens_9606");
  const [customTable, setCustomTable] = useState("");
  const [useCustomTable, setUseCustomTable] = useState(false);
  const [avoidHomopolymersGC, setAvoidHomopolymersGC] = useState(4);
  const [avoidHomopolymersAT, setAvoidHomopolymersAT] = useState(6);
  const [gcMin, setGcMin] = useState(25);
  const [gcMax, setGcMax] = useState(75);
  const [gcWindow, setGcWindow] = useState(50);
  const [uniquifyKmers, setUniquifyKmers] = useState(10);
  const [minCodonFrequency, setMinCodonFrequency] = useState(10);
  const [beamK, setBeamK] = useState(250);
  const [avoidPatterns, setAvoidPatterns] = useState("");
  const [progress, setProgress] = useState<{ stage: string; current: number; total: number } | null>(null);
  const [tables, setTables] = useState<string[]>([]);
  const [customTables, setCustomTables] = useState<string[]>([]);
  const [result, setResult] = useState<CodonOptResult | null>(null);
  const [editedDna, setEditedDna] = useState<string | null>(null);
  const [editedCodons, setEditedCodons] = useState<CodonDetail[] | null>(null);
  const [codonHistory, setCodonHistory] = useState<{ dna: string; codons: CodonDetail[] }[]>([]);
  const [editedConstraints, setEditedConstraints] = useState<{
    constraints_pass: boolean;
    constraints_summary: { constraint: string; passing: boolean; message: string; enforced?: boolean }[];
    gc_content: number;
    warnings: SequenceWarning[];
  } | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [copied, setCopied] = useState(false);

  // Handle DNA sent back from oligo designer
  useEffect(() => {
    if (prefillDna && result) {
      // Set as edited DNA so user can continue tweaking
      handleCodonDnaChange(prefillDna);
      onPrefillConsumed?.();
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [prefillDna]);

  useEffect(() => {
    fetch(`${API_URL}/api/codon-tables`)
      .then((r) => r.json())
      .then((d) => {
        setTables(d.tables);
        setCustomTables(d.custom_tables || []);
        // Default to the first bundled custom table if available
        if (d.custom_tables?.length > 0) {
          setSpecies(d.custom_tables[0]);
        }
      })
      .catch(() => {});
  }, []);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    setError("");
    setResult(null);
    setEditedDna(null);
    setEditedCodons(null);
    setEditedConstraints(null);
    setProgress(null);
    setLoading(true);

    try {
      const patterns = avoidPatterns
        .split(",")
        .map((s) => s.trim())
        .filter(Boolean);

      const body: Record<string, unknown> = {
        protein_sequence: protein,
        avoid_homopolymers_gc: avoidHomopolymersGC,
        avoid_homopolymers_at: avoidHomopolymersAT,
        gc_min: gcMin / 100,
        gc_max: gcMax / 100,
        gc_window: gcWindow,
        uniquify_kmers: uniquifyKmers,
        avoid_patterns: patterns,
        min_codon_frequency: minCodonFrequency,
        beam_k: beamK,
      };

      if (useCustomTable && customTable.trim()) {
        body.custom_codon_table = customTable;
      } else {
        body.species = species;
      }

      const res = await fetch(`${API_URL}/api/codon-optimize-stream`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      if (!res.ok || !res.body) {
        throw new Error(`Optimization failed (${res.status})`);
      }

      const reader = res.body.getReader();
      const decoder = new TextDecoder();
      let buffer = "";
      let finalResult: CodonOptResult | null = null;

      while (true) {
        const { done, value } = await reader.read();
        if (done) break;
        buffer += decoder.decode(value, { stream: true });

        const events = buffer.split("\n\n");
        buffer = events.pop() ?? "";
        for (const block of events) {
          const line = block.split("\n").find((l) => l.startsWith("data: "));
          if (!line) continue;
          const payload = JSON.parse(line.slice(6));
          if (payload.type === "progress") {
            setProgress({ stage: payload.stage, current: payload.current, total: payload.total });
            // Yield to React so it can paint this update before the next
            // setProgress in the same chunk gets batched on top of it.
            await new Promise((r) => setTimeout(r, 0));
          } else if (payload.type === "result") {
            // Snap the bar to 100% so the user visibly sees it complete
            // before the result hides it.
            if (payload.data?.length_protein) {
              setProgress({ stage: "beam", current: payload.data.length_protein, total: payload.data.length_protein });
              await new Promise((r) => setTimeout(r, 0));
            }
            finalResult = payload.data;
          } else if (payload.type === "error") {
            throw new Error(payload.message || "Optimization failed");
          }
        }
      }

      if (!finalResult) {
        throw new Error("No result received from server");
      }
      setResult(finalResult);
      setCodonHistory([]);
      onDnaChanged?.(finalResult.optimized_dna);
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : "Unknown error");
    } finally {
      setLoading(false);
      setProgress(null);
    }
  }

  // Build a lookup from the codon table: { codon -> per_thousand }
  const codonLookup = useCallback((table: CodonTableEntry[]) => {
    const map: Record<string, number> = {};
    for (const e of table) {
      map[e.codon] = e.per_thousand;
    }
    return map;
  }, []);

  function handleCodonDnaChange(newDna: string) {
    if (!result) return;

    // Push current state onto undo history before applying the change
    const currentDna = editedDna ?? result.optimized_dna;
    const currentCodons = editedCodons ?? result.codons;
    setCodonHistory((prev) => [...prev, { dna: currentDna, codons: currentCodons }].slice(-50));

    const lookup = codonLookup(result.codon_table);
    const newCodons: CodonDetail[] = [];
    for (let i = 0; i < newDna.length - 2; i += 3) {
      const codon = newDna.substring(i, i + 3);
      const aa = result.codons[i / 3]?.amino_acid ?? "?";
      newCodons.push({
        amino_acid: aa,
        codon,
        per_thousand: lookup[codon] ?? 0,
      });
    }
    setEditedDna(newDna);
    setEditedCodons(newCodons);
    onDnaChanged?.(newDna);
    recheckConstraints(newDna);
  }

  function recheckConstraints(dnaToCheck: string) {
    const patterns = avoidPatterns.split(",").map((s) => s.trim()).filter(Boolean);
    fetch(`${API_URL}/api/check-constraints`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        dna_sequence: dnaToCheck,
        avoid_homopolymers_gc: avoidHomopolymersGC,
        avoid_homopolymers_at: avoidHomopolymersAT,
        gc_min: gcMin / 100,
        gc_max: gcMax / 100,
        gc_window: gcWindow,
        uniquify_kmers: uniquifyKmers,
        avoid_patterns: patterns,
        min_codon_frequency: minCodonFrequency,
      }),
    })
      .then((r) => r.json())
      .then((d) => setEditedConstraints(d))
      .catch(() => {});
  }

  function handleCodonUndo() {
    if (codonHistory.length === 0) return;
    const prev = codonHistory[codonHistory.length - 1];
    setCodonHistory((h) => h.slice(0, -1));
    setEditedDna(prev.dna);
    setEditedCodons(prev.codons);
    onDnaChanged?.(prev.dna);
    recheckConstraints(prev.dna);
  }

  const currentDna = editedDna ?? result?.optimized_dna ?? "";
  const currentCodons = editedCodons ?? result?.codons ?? [];

  function handleCopyDna() {
    if (currentDna) {
      navigator.clipboard.writeText(currentDna);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    }
  }

  return (
    <div>
      <h1 className="text-xl font-semibold mb-2 text-[#202124]">Codon Optimization <span className="text-sm font-normal text-[#5f6368]">DNA Chisel</span></h1>
      <p className="text-[#5f6368] mb-8">
        Optimize a protein sequence for expression using DNA Chisel
      </p>

      <form onSubmit={handleSubmit} className="space-y-4 mb-8">
        <div>
          <div className="flex items-center justify-between mb-1">
            <label className="text-sm font-medium text-[#202124]">
              Protein Sequence
            </label>
            <button
              type="button"
              onClick={() => setProtein("MKGLVTDDEGQPIPGATIKVNRNKKPVTSSARGEYWRLLLPGNYTLTASAQEYISLTISI\nHIPEKTDPDDWSKPAQRQDFRLVRKSTGVTKPLTNSNNPTIITNEKNVPTELSFPIEPS\nASSINSPPFITTKSNLKPNFFPTELTAQTEEPRHFENKGKQNNKPFQNNSNKDIKNNEN\nPRIEATTRRSNQSRGRGGRRCKCANGKRRCKKCSKRSKTRNRNLRNAQKKQK*")}
              className="text-xs text-[#1a73e8] hover:text-[#1967d2] cursor-pointer"
            >
              Use example sequence
            </button>
          </div>
          <textarea
            value={protein}
            onChange={(e) => setProtein(e.target.value)}
            placeholder="Paste amino acid sequence (single-letter codes)..."
            rows={6}
            className="w-full bg-white border border-[#dadce0] p-3 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
          />
        </div>

        <div className="space-y-3">
          <div className="flex items-center gap-4">
            <label className="flex items-center gap-2 text-sm text-[#202124] cursor-pointer">
              <input
                type="radio"
                checked={!useCustomTable}
                onChange={() => setUseCustomTable(false)}
                className="accent-[#1a73e8]"
              />
              Built-in codon table
            </label>
            <label className="flex items-center gap-2 text-sm text-[#202124] cursor-pointer">
              <input
                type="radio"
                checked={useCustomTable}
                onChange={() => setUseCustomTable(true)}
                className="accent-[#1a73e8]"
              />
              Custom codon table
            </label>
          </div>

          {!useCustomTable ? (
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Species
              </label>
              <select
                value={species}
                onChange={(e) => setSpecies(e.target.value)}
                className="bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              >
                {customTables.length > 0 && (
                  <optgroup label="Lab tables">
                    {customTables.map((t) => (
                      <option key={t} value={t}>
                        {t.replace(/_/g, " ")}
                      </option>
                    ))}
                  </optgroup>
                )}
                <optgroup label="Species">
                  {tables.map((t) => (
                    <option key={t} value={t}>
                      {t.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase())}
                    </option>
                  ))}
                </optgroup>
              </select>
            </div>
          ) : (
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Custom codon table (Kazusa/GCG format)
              </label>
              <textarea
                value={customTable}
                onChange={(e) => setCustomTable(e.target.value)}
                placeholder={"AmAcid  Codon     Number    /1000     Fraction\nGly     GGG      905.00     18.70      0.00\nGly     GGA      527.00     10.89      0.00\n..."}
                rows={6}
                className="w-full bg-white border border-[#dadce0] p-3 font-mono text-xs text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
          )}
        </div>

        <details className="text-sm">
          <summary className="cursor-pointer text-[#5f6368] hover:text-[#202124]">
            Synthesis constraints
          </summary>
          <table className="mt-2 text-sm">
            <tbody>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    No G/C runs &ge;
                    <InfoTip text="Avoid runs of this many or more consecutive G or C bases. G/C runs are especially problematic for synthesis." />
                  </span>
                </td>
                <td className="py-1.5 pr-6 align-middle">
                  <NumericInput value={avoidHomopolymersGC} onChange={setAvoidHomopolymersGC} min={3} max={8} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    No A/T runs &ge;
                    <InfoTip text="Avoid runs of this many or more consecutive A or T bases. Less problematic than G/C but long runs still cause issues." />
                  </span>
                </td>
                <td className="py-1.5 align-middle">
                  <NumericInput value={avoidHomopolymersAT} onChange={setAvoidHomopolymersAT} min={3} max={8} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    GC% min
                    <InfoTip text="Minimum GC content allowed in any sliding window. Regions below this are hard to synthesize." />
                  </span>
                </td>
                <td className="py-1.5 pr-6 align-middle">
                  <NumericInput value={gcMin} onChange={setGcMin} step={5} min={0} max={100} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    GC% max
                    <InfoTip text="Maximum GC content allowed in any sliding window. High GC regions form strong secondary structures." />
                  </span>
                </td>
                <td className="py-1.5 pr-6 align-middle">
                  <NumericInput value={gcMax} onChange={setGcMax} step={5} min={0} max={100} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    Window (bp)
                    <InfoTip text="Size of the sliding window (in base pairs) used to check GC content." />
                  </span>
                </td>
                <td className="py-1.5 align-middle">
                  <NumericInput value={gcWindow} onChange={setGcWindow} min={20} max={200} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    Uniquify k-mers (bp)
                    <InfoTip text="Ensure no two identical k-mers of this length appear in the sequence. Prevents mispriming during assembly." />
                  </span>
                </td>
                <td className="py-1.5 align-middle">
                  <NumericInput value={uniquifyKmers} onChange={setUniquifyKmers} min={6} max={20} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    Min codon frequency (/1000)
                    <InfoTip text="Minimum allowed codon usage frequency. Codons below this threshold are excluded during optimization to avoid rare codons that may reduce expression." />
                  </span>
                </td>
                <td className="py-1.5 align-middle">
                  <NumericInput value={minCodonFrequency} onChange={setMinCodonFrequency} min={0} max={100} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    Beam width (K)
                    <InfoTip text="Beam search width per bucket in the refinement step. Higher values explore more candidate sequences and find lower-loss solutions, but take longer. If it's slow, lower this." />
                  </span>
                </td>
                <td className="py-1.5 align-middle">
                  <NumericInput value={beamK} onChange={setBeamK} min={1} max={500} className="w-16 bg-white border border-[#dadce0] px-2 py-1 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
              <tr>
                <td className="py-1.5 pr-4 text-[#5f6368] whitespace-nowrap align-middle">
                  <span className="flex items-center gap-1">
                    Avoid patterns
                    <InfoTip text="DNA motifs to exclude, e.g. restriction sites. Comma-separated." />
                  </span>
                </td>
                <td colSpan={5} className="py-1.5 align-middle">
                  <input type="text" value={avoidPatterns} onChange={(e) => setAvoidPatterns(e.target.value)} placeholder="e.g. GAATTC, GGATCC" className="w-full bg-white border border-[#dadce0] px-2 py-1 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none" />
                </td>
              </tr>
            </tbody>
          </table>
        </details>

        <button
          type="submit"
          disabled={loading || !protein.trim()}
          className="bg-[#1a73e8] hover:bg-[#1967d2] text-white disabled:bg-[#dadce0] disabled:text-[#80868b] px-6 py-2 font-medium transition-colors cursor-pointer"
        >
          {loading ? "Optimizing..." : "Optimize Codons"}
        </button>

        {loading && (
          <div className="mt-3">
            <div className="flex items-center justify-between mb-1 text-xs text-[#5f6368]">
              <span>
                Beam search
                {progress && progress.total > 0 && ` — codon ${progress.current} of ${progress.total}`}
              </span>
              <span>{progress && progress.total > 0 ? `${Math.round((progress.current / progress.total) * 100)}%` : ""}</span>
            </div>
            <div className="w-full h-1.5 bg-[#e8eaed] overflow-hidden">
              <div
                className="h-full bg-[#1a73e8] transition-all"
                style={{ width: progress && progress.total > 0 ? `${(progress.current / progress.total) * 100}%` : "0%" }}
              />
            </div>
          </div>
        )}
      </form>

      {error && (
        <div className="bg-[#fce8e6] border border-[#d93025] p-4 mb-6">
          <p className="text-[#d93025]">{error}</p>
        </div>
      )}

      {result && (
        <div className="space-y-6">
          <div className="bg-white p-4 border border-[#dadce0]">
            <h2 className="text-lg font-medium mb-3 text-[#202124]">Summary</h2>
            <table className="text-sm">
              <tbody>
                <tr>
                  <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Protein length</td>
                  <td className="py-0.5">{result.length_protein} aa</td>
                </tr>
                <tr>
                  <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">DNA length</td>
                  <td className="py-0.5">{result.length_dna} bp</td>
                </tr>
                <tr>
                  <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">GC content (overall)</td>
                  <td className="py-0.5">{((editedConstraints?.gc_content ?? result.gc_content) * 100).toFixed(1)}%</td>
                </tr>
                {result.gc_window_min != null && result.gc_window_max != null && (
                  <tr>
                    <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">GC range ({result.gc_window_size}bp window)</td>
                    <td className="py-0.5">
                      {(result.gc_window_min * 100).toFixed(1)}% – {(result.gc_window_max * 100).toFixed(1)}%
                      <span className="text-[#80868b] ml-2">
                        (allowed {(gcMin).toFixed(0)}–{(gcMax).toFixed(0)}%)
                      </span>
                    </td>
                  </tr>
                )}
                <tr>
                  <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">All constraints pass</td>
                  <td className="py-0.5">
                    {(() => {
                      const pass = editedConstraints?.constraints_pass ?? result.constraints_pass;
                      return (
                        <span className={pass ? "text-[#188038]" : "text-[#d93025]"}>
                          {pass ? "Yes" : "No"}
                        </span>
                      );
                    })()}
                  </td>
                </tr>
                <tr>
                  <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Translation matches input</td>
                  <td className="py-0.5">
                    <span className={
                      result.input_protein === result.back_translated_protein.replace(/\*$/, "")
                        ? "text-[#188038]"
                        : "text-[#d93025]"
                    }>
                      {result.input_protein === result.back_translated_protein.replace(/\*$/, "")
                        ? "Yes"
                        : "No"}
                    </span>
                  </td>
                </tr>
              </tbody>
            </table>
          </div>

          {/* Constraints + detection checks */}
          {(() => {
            const allChecks = [...(editedConstraints?.constraints_summary ?? result.constraints_summary)];
            // Add client-side min codon frequency check for live edits
            if (editedConstraints && minCodonFrequency > 0) {
              const belowFloor = currentCodons.filter((c) => c.per_thousand < minCodonFrequency);
              allChecks.push({
                constraint: `MinCodonFrequency(min:${minCodonFrequency}/1000)`,
                passing: belowFloor.length === 0,
                message: belowFloor.length === 0 ? "" : `${belowFloor.length} codon(s) below ${minCodonFrequency}/1000`,
                enforced: true,
              });
            }
            const enforced = allChecks.filter((c) => c.enforced !== false);
            const detection = allChecks.filter((c) => c.enforced === false);
            return (
              <div className="bg-white border border-[#dadce0] p-4">
                <div className="flex items-center justify-between mb-3">
                  <h2 className="text-lg font-medium text-[#202124]">Constraints</h2>
                  {editedConstraints && (
                    <span className="text-xs text-[#5f6368]">Live check</span>
                  )}
                </div>
                <p className="text-xs text-[#5f6368] mb-2">Enforced during optimization</p>
                <div className="space-y-1 text-sm">
                  {enforced.map((c, i) => (
                    <div key={i} className="flex items-center gap-2">
                      <span className={c.passing ? "text-[#188038]" : "text-[#d93025]"}>
                        {c.passing ? "PASS" : "FAIL"}
                      </span>
                      <span className="text-[#5f6368]">{c.constraint}</span>
                    </div>
                  ))}
                </div>
                {detection.length > 0 && (
                  <>
                    <div className="border-t border-[#dadce0] my-3" />
                    <p className="text-xs text-[#5f6368] mb-2">Detection only (not enforced)</p>
                    <div className="space-y-1 text-sm">
                      {detection.map((c, i) => (
                        <div key={i} className="flex items-center gap-2">
                          <span className={c.passing ? "text-[#188038]" : "text-[#e37400]"}>
                            {c.passing ? "PASS" : "FAIL"}
                          </span>
                          <span className="text-[#5f6368]">{c.constraint}</span>
                        </div>
                      ))}
                    </div>
                  </>
                )}
              </div>
            );
          })()}

          {/* Warnings */}
          {(() => {
            const originalWarnings = result.warnings ?? [];
            const currentWarnings = editedConstraints?.warnings ?? originalWarnings;

            // Build display list from originals, marking which are resolved
            const displayWarnings = originalWarnings.map((ow) => {
              const stillActive = currentWarnings.some(
                (cw) => cw.kind === ow.kind && cw.start === ow.start && cw.end === ow.end
              );
              return { ...ow, resolved: !stillActive };
            });

            // Add any NEW warnings not in the original set
            // Offset their group numbers so they don't collide with original groups
            const maxOriginalGroup = Math.max(0, ...originalWarnings.map((w) => w.group ?? 0));
            const newWarnings = currentWarnings
              .filter((cw) => !originalWarnings.some(
                (ow) => ow.kind === cw.kind && ow.start === cw.start && ow.end === cw.end
              ))
              .map((w) => ({
                ...w,
                resolved: false,
                group: w.group != null ? w.group + maxOriginalGroup : w.group,
              }));

            const allWarnings = [...displayWarnings, ...newWarnings];
            const activeCount = allWarnings.filter((w) => !w.resolved).length;

            const GROUP_LABELS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"];

            return (
              <div className="bg-white border border-[#dadce0] p-4 h-[480px] flex flex-col">
                <div className="flex items-center justify-between mb-3 flex-shrink-0">
                  <h2 className="text-lg font-medium text-[#202124]">
                    Warnings
                    <span className="text-sm font-normal text-[#5f6368] ml-2">
                      ({activeCount})
                    </span>
                  </h2>
                  {editedConstraints && (
                    <span className="text-xs text-[#5f6368]">Live check</span>
                  )}
                </div>
                <div className="overflow-y-auto flex-1">
                  {allWarnings.length > 0 ? (
                    <div className="space-y-1 text-sm">
                      {allWarnings.map((w, i) => (
                        <div
                          key={i}
                          className="flex items-center gap-2"
                          style={{ opacity: w.resolved ? 0.25 : 1 }}
                        >
                          <span className="text-[#80868b] flex-shrink-0 w-5 text-right text-xs">{i + 1}</span>
                          <span className="text-[#e37400] flex-shrink-0">
                            {w.kind === "homopolymer" ? "HOMOPOLYMER" : w.kind === "gc_window" ? "GC" : "10-MER"}
                          </span>
                          <span className="text-[#202124]">{w.message}</span>
                          {w.group != null && (
                            <span
                              className="inline-flex items-center justify-center w-4 h-4 text-[9px] font-bold text-white flex-shrink-0"
                              style={{ backgroundColor: "#202124", borderRadius: "50%" }}
                              title={`Group ${GROUP_LABELS[(w.group - 1) % GROUP_LABELS.length]} — matching repeats share this label`}
                            >
                              {GROUP_LABELS[(w.group - 1) % GROUP_LABELS.length]}
                            </span>
                          )}
                        </div>
                      ))}
                    </div>
                  ) : (
                    <p className="text-sm text-[#188038]">No warnings</p>
                  )}
                </div>
              </div>
            );
          })()}

          {/* Optimized DNA with copy + send to tiler */}
          <div className="bg-white border border-[#dadce0] p-4">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-lg font-medium text-[#202124]">
                {editedDna ? "Edited DNA" : "Optimized DNA"}
              </h2>
              <div className="flex items-center gap-3">
                {editedDna && (
                  <button
                    onClick={() => {
                      setEditedDna(null);
                      setEditedCodons(null);
                      setEditedConstraints(null);
                      onDnaChanged?.(result.optimized_dna);
                    }}
                    className="text-sm text-[#5f6368] hover:text-[#202124] cursor-pointer"
                  >
                    Reset to optimized
                  </button>
                )}
                <button
                  onClick={handleCopyDna}
                  className="text-sm text-[#1a73e8] hover:text-[#1967d2] cursor-pointer min-w-[120px] text-right"
                >
                  {copied ? "Copied" : "Copy to clipboard"}
                </button>
              </div>
            </div>
            <pre className="font-mono text-xs text-[#202124] whitespace-pre-wrap break-all bg-[#f8f9fa] border border-[#dadce0] p-3 max-h-64 overflow-y-auto">
              {currentDna}
            </pre>
          </div>

          {/* Optimization report */}
          {result.optimization_report && result.optimization_report.ideal_score != null && (
            <div className="bg-white p-4 border border-[#dadce0]">
              <h2 className="text-lg font-medium mb-3 text-[#202124]">Optimization Report</h2>
              <table className="text-sm">
                <tbody>
                  <tr>
                    <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Ideal score (all best codons)</td>
                    <td className="py-0.5">{result.optimization_report.ideal_score.toFixed(1)}</td>
                  </tr>
                  {result.optimization_report.chisel_score != null && result.optimization_report.chisel_loss != null && (
                    <tr>
                      <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">DNA Chisel score</td>
                      <td className="py-0.5">
                        {result.optimization_report.chisel_score.toFixed(1)}
                        <span className="text-[#80868b] ml-2">(loss {result.optimization_report.chisel_loss.toFixed(1)})</span>
                      </td>
                    </tr>
                  )}
                  {result.optimization_report.beam_score != null && result.optimization_report.beam_loss != null && (
                    <tr>
                      <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Beam search score</td>
                      <td className="py-0.5">
                        {result.optimization_report.beam_score.toFixed(1)}
                        <span className="text-[#80868b] ml-2">(loss {result.optimization_report.beam_loss.toFixed(1)})</span>
                      </td>
                    </tr>
                  )}
                  <tr>
                    <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Chosen score</td>
                    <td className="py-0.5">
                      {result.optimization_report.chosen_score.toFixed(1)}
                      <span className="text-[#80868b] ml-2">
                        ({result.optimization_report.ideal_score > 0
                          ? `${(100 * result.optimization_report.chosen_score / result.optimization_report.ideal_score).toFixed(1)}% of ideal, `
                          : ""}
                        {result.optimization_report.source === "beam" ? "beam" :
                         result.optimization_report.source === "dna_chisel_fallback" ? "DNA Chisel fallback" :
                         result.optimization_report.source === "beam_pruned_empty" ? "beam failed (constraints too tight) — using DNA Chisel" :
                         result.optimization_report.source === "beam_no_solution" ? "beam found no solution — using DNA Chisel" :
                         result.optimization_report.source}
                        {result.optimization_report.source === "beam" && result.optimization_report.beam_k != null && ` K=${result.optimization_report.beam_k}`})
                      </span>
                    </td>
                  </tr>
                  <tr>
                    <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Codons adjusted by optimizer</td>
                    <td className="py-0.5">{result.optimization_report.codons_changed} of {result.length_protein}</td>
                  </tr>
                  {result.optimization_report.min_frequency_achieved != null && (
                    <tr>
                      <td className="py-0.5 pr-4 text-[#5f6368] whitespace-nowrap">Min codon frequency in output</td>
                      <td className="py-0.5">
                        {result.optimization_report.min_frequency_achieved.toFixed(1)}/1000
                        <span className="text-[#80868b] ml-2">
                          (requested floor {minCodonFrequency}/1000)
                        </span>
                      </td>
                    </tr>
                  )}
                </tbody>
              </table>
            </div>
          )}

          {/* Interactive codon track */}
          <CodonTrackViewer
            codons={currentCodons}
            codonTable={result.codon_table}
            warnings={editedConstraints?.warnings ?? result.warnings ?? []}
            originalWarnings={result.warnings ?? []}
            onDnaChange={handleCodonDnaChange}
            onUndo={handleCodonUndo}
            canUndo={codonHistory.length > 0}
          />

          {/* Codon table used */}
          <details className="bg-white border border-[#dadce0]">
            <summary className="p-4 cursor-pointer text-lg font-medium text-[#202124] hover:bg-[#f1f3f4]">
              Codon Table
              <span className="text-sm font-normal text-[#5f6368] ml-2">{result.codon_table_name}</span>
            </summary>
            <div className="px-4 pb-4">
              <CodonTableView
                entries={result.codon_table}
                tableName={result.codon_table_name}
                hideHeader
              />
            </div>
          </details>
        </div>
      )}
    </div>
  );
}
