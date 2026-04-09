"use client";

import { useState, useEffect, useCallback } from "react";
import CodonTableView from "./CodonTableView";
import CodonTrackViewer from "./CodonTrackViewer";

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";

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
  constraints_pass: boolean;
  constraints_summary: { constraint: string; passing: boolean; message: string }[];
  warnings: SequenceWarning[];
  codons: CodonDetail[];
  codon_table: CodonTableEntry[];
  codon_table_name: string;
}

interface CodonOptPageProps {
  onSendToTiler?: (dna: string) => void;
}

export default function CodonOptPage({ onSendToTiler }: CodonOptPageProps) {
  const [protein, setProtein] = useState("");
  const [species, setSpecies] = useState("h_sapiens_9606");
  const [customTable, setCustomTable] = useState("");
  const [useCustomTable, setUseCustomTable] = useState(false);
  const [avoidHomopolymers, setAvoidHomopolymers] = useState(5);
  const [gcMin, setGcMin] = useState(0.25);
  const [gcMax, setGcMax] = useState(0.75);
  const [gcWindow, setGcWindow] = useState(50);
  const [uniquifyKmers, setUniquifyKmers] = useState(10);
  const [avoidPatterns, setAvoidPatterns] = useState("");
  const [tables, setTables] = useState<string[]>([]);
  const [customTables, setCustomTables] = useState<string[]>([]);
  const [result, setResult] = useState<CodonOptResult | null>(null);
  const [editedDna, setEditedDna] = useState<string | null>(null);
  const [editedCodons, setEditedCodons] = useState<CodonDetail[] | null>(null);
  const [editedConstraints, setEditedConstraints] = useState<{
    constraints_pass: boolean;
    constraints_summary: { constraint: string; passing: boolean; message: string }[];
    gc_content: number;
    warnings: SequenceWarning[];
  } | null>(null);
  const [error, setError] = useState("");
  const [loading, setLoading] = useState(false);
  const [copied, setCopied] = useState(false);

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
    setLoading(true);

    try {
      const patterns = avoidPatterns
        .split(",")
        .map((s) => s.trim())
        .filter(Boolean);

      const body: Record<string, unknown> = {
        protein_sequence: protein,
        avoid_homopolymers: avoidHomopolymers,
        gc_min: gcMin,
        gc_max: gcMax,
        gc_window: gcWindow,
        uniquify_kmers: uniquifyKmers,
        avoid_patterns: patterns,
      };

      if (useCustomTable && customTable.trim()) {
        body.custom_codon_table = customTable;
      } else {
        body.species = species;
      }

      const res = await fetch(`${API_URL}/api/codon-optimize`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(body),
      });
      if (!res.ok) {
        const data = await res.json();
        throw new Error(data.detail || "Optimization failed");
      }
      setResult(await res.json());
    } catch (err: unknown) {
      setError(err instanceof Error ? err.message : "Unknown error");
    } finally {
      setLoading(false);
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

    // Re-check constraints with current settings
    const patterns = avoidPatterns.split(",").map((s) => s.trim()).filter(Boolean);
    fetch(`${API_URL}/api/check-constraints`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        dna_sequence: newDna,
        avoid_homopolymers: avoidHomopolymers,
        gc_min: gcMin,
        gc_max: gcMax,
        gc_window: gcWindow,
        uniquify_kmers: uniquifyKmers,
        avoid_patterns: patterns,
      }),
    })
      .then((r) => r.json())
      .then((d) => setEditedConstraints(d))
      .catch(() => {});
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
      <h1 className="text-3xl font-semibold mb-2 text-[#202124]">Codon Optimization</h1>
      <p className="text-[#5f6368] mb-8">
        Optimize a protein sequence for expression using DNA Chisel
      </p>

      <form onSubmit={handleSubmit} className="space-y-4 mb-8">
        <div>
          <label className="block text-sm font-medium mb-1 text-[#202124]">
            Protein Sequence
          </label>
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
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mt-2">
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Max homopolymer
              </label>
              <input
                type="number"
                value={avoidHomopolymers}
                onChange={(e) => setAvoidHomopolymers(Number(e.target.value))}
                min={3}
                max={8}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                GC% min
              </label>
              <input
                type="number"
                value={gcMin}
                onChange={(e) => setGcMin(Number(e.target.value))}
                step={0.05}
                min={0}
                max={1}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                GC% max
              </label>
              <input
                type="number"
                value={gcMax}
                onChange={(e) => setGcMax(Number(e.target.value))}
                step={0.05}
                min={0}
                max={1}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                GC window (bp)
              </label>
              <input
                type="number"
                value={gcWindow}
                onChange={(e) => setGcWindow(Number(e.target.value))}
                min={20}
                max={200}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
          </div>
          <div className="grid grid-cols-2 gap-4 mt-2">
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Uniquify k-mers (bp)
              </label>
              <input
                type="number"
                value={uniquifyKmers}
                onChange={(e) => setUniquifyKmers(Number(e.target.value))}
                min={6}
                max={20}
                className="w-full bg-white border border-[#dadce0] p-2 text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
            <div>
              <label className="block text-sm font-medium mb-1 text-[#202124]">
                Avoid patterns (comma-separated)
              </label>
              <input
                type="text"
                value={avoidPatterns}
                onChange={(e) => setAvoidPatterns(e.target.value)}
                placeholder="e.g. GAATTC, GGATCC"
                className="w-full bg-white border border-[#dadce0] p-2 font-mono text-sm text-[#202124] focus:border-[#1a73e8] focus:outline-none"
              />
            </div>
          </div>
        </details>

        <button
          type="submit"
          disabled={loading || !protein.trim()}
          className="bg-[#1a73e8] hover:bg-[#1967d2] text-white disabled:bg-[#dadce0] disabled:text-[#80868b] px-6 py-2 font-medium transition-colors cursor-pointer"
        >
          {loading ? "Optimizing..." : "Optimize Codons"}
        </button>
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
            <div className="space-y-1 text-sm">
              <div>
                <span className="text-[#5f6368]">Protein length:</span>{" "}
                {result.length_protein} aa
              </div>
              <div>
                <span className="text-[#5f6368]">DNA length:</span>{" "}
                {result.length_dna} bp
              </div>
              <div>
                <span className="text-[#5f6368]">GC content:</span>{" "}
                {((editedConstraints?.gc_content ?? result.gc_content) * 100).toFixed(1)}%
              </div>
              <div>
                <span className="text-[#5f6368]">All constraints pass:</span>{" "}
                {(() => {
                  const pass = editedConstraints?.constraints_pass ?? result.constraints_pass;
                  return (
                    <span className={pass ? "text-[#188038]" : "text-[#d93025]"}>
                      {pass ? "Yes" : "No"}
                    </span>
                  );
                })()}
              </div>
              <div>
                <span className="text-[#5f6368]">Translation matches input:</span>{" "}
                <span className={
                  result.input_protein === result.back_translated_protein.replace(/\*$/, "")
                    ? "text-[#188038]"
                    : "text-[#d93025]"
                }>
                  {result.input_protein === result.back_translated_protein.replace(/\*$/, "")
                    ? "Yes"
                    : "No"}
                </span>
              </div>
            </div>
          </div>

          {/* Constraints detail */}
          <div className="bg-white border border-[#dadce0] p-4">
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-lg font-medium text-[#202124]">Constraints</h2>
              {editedConstraints && (
                <span className="text-xs text-[#5f6368]">Live check</span>
              )}
            </div>
            <div className="space-y-1 text-sm">
              {(editedConstraints?.constraints_summary ?? result.constraints_summary).map((c, i) => (
                <div key={i} className="flex items-center gap-2">
                  <span className={c.passing ? "text-[#188038]" : "text-[#d93025]"}>
                    {c.passing ? "PASS" : "FAIL"}
                  </span>
                  <span className="text-[#5f6368]">{c.constraint}</span>
                </div>
              ))}
            </div>
          </div>

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
                    onClick={() => { setEditedDna(null); setEditedCodons(null); }}
                    className="text-sm text-[#5f6368] hover:text-[#202124] cursor-pointer"
                  >
                    Reset to optimized
                  </button>
                )}
                {onSendToTiler && (
                  <button
                    onClick={() => onSendToTiler(currentDna)}
                    className="text-sm text-[#1a73e8] hover:text-[#1967d2] cursor-pointer"
                  >
                    Send to Oligo Designer
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

          {/* Interactive codon track */}
          <CodonTrackViewer
            codons={currentCodons}
            codonTable={result.codon_table}
            warnings={editedConstraints?.warnings ?? result.warnings ?? []}
            originalWarnings={result.warnings ?? []}
            onDnaChange={handleCodonDnaChange}
          />

          {/* Codon table used */}
          <CodonTableView
            entries={result.codon_table}
            tableName={result.codon_table_name}
          />
        </div>
      )}
    </div>
  );
}
