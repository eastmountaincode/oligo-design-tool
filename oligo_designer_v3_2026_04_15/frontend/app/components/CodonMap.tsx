"use client";

import { useState, useMemo } from "react";

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

interface Props {
  codons: CodonDetail[];
  codonTable: CodonTableEntry[];
  onDnaChange?: (newDna: string) => void;
}

/**
 * Color a codon relative to the best/worst choice for its amino acid.
 * Returns a 0-1 ratio where 1 = best codon for this aa, 0 = worst.
 */
function codonRank(
  perThousand: number,
  aaAlternatives: { codon: string; per_thousand: number }[]
): number {
  if (aaAlternatives.length <= 1) return 1;
  const max = Math.max(...aaAlternatives.map((a) => a.per_thousand));
  const min = Math.min(...aaAlternatives.map((a) => a.per_thousand));
  if (max === min) return 1;
  return (perThousand - min) / (max - min);
}

function rankColor(rank: number): string {
  if (rank >= 0.8) return "#188038";
  if (rank >= 0.4) return "#e37400";
  return "#d93025";
}

function rankBg(rank: number): string {
  if (rank >= 0.8) return "#e6f4ea";
  if (rank >= 0.4) return "#fef7e0";
  return "#fce8e6";
}

export default function CodonMap({ codons, codonTable, onDnaChange }: Props) {
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);

  // Build lookup: amino acid -> list of alternative codons with /1000 values
  const aaAlternatives = useMemo(() => {
    const map: Record<string, { codon: string; per_thousand: number }[]> = {};
    for (const e of codonTable) {
      if (!map[e.amino_acid]) map[e.amino_acid] = [];
      map[e.amino_acid].push({ codon: e.codon, per_thousand: e.per_thousand });
    }
    // Sort each by per_thousand descending
    for (const aa of Object.keys(map)) {
      map[aa].sort((a, b) => b.per_thousand - a.per_thousand);
    }
    return map;
  }, [codonTable]);

  function handleSwap(index: number, newCodon: string) {
    if (!onDnaChange) return;
    const newCodons = codons.map((c, i) => (i === index ? newCodon : c.codon));
    onDnaChange(newCodons.join(""));
    setSelectedIndex(null);
  }

  if (!codons || codons.length === 0) return null;

  return (
    <div className="bg-white border border-[#dadce0] p-4">
      <h2 className="text-lg font-medium mb-1 text-[#202124]">Codon Map</h2>
      <p className="text-xs text-[#5f6368] mb-3">
        Color indicates quality of codon choice relative to alternatives for the same amino acid. Click to swap.
      </p>
      <div>
        <div className="flex flex-wrap gap-px font-mono text-xs">
          {codons.map((c, i) => {
            const alts = aaAlternatives[c.amino_acid] || [];
            const rank = codonRank(c.per_thousand, alts);
            const isSelected = selectedIndex === i;

            return (
              <div key={i} className="relative" style={{ width: 32 }}>
                <div
                  className="flex flex-col items-center py-1 cursor-pointer"
                  style={{
                    backgroundColor: rankBg(rank),
                    outline: isSelected ? "2px solid #1a73e8" : "none",
                    outlineOffset: "-2px",
                    width: 32,
                  }}
                  title={`${c.amino_acid}: ${c.codon} (${c.per_thousand}/1000)`}
                  onClick={() => setSelectedIndex(isSelected ? null : i)}
                >
                  <span className="text-[#5f6368] text-[10px]">{c.amino_acid}</span>
                  <span style={{ color: rankColor(rank) }} className="font-medium">
                    {c.codon}
                  </span>
                  <span className="text-[9px] text-[#80868b]">
                    {c.per_thousand.toFixed(1)}
                  </span>
                </div>

                {/* Dropdown for alternatives */}
                {isSelected && alts.length > 1 && (
                  <div className="absolute bottom-full left-0 z-10 mb-1 bg-white border border-[#dadce0] shadow-lg min-w-[100px]">
                    <div className="px-2 py-1 text-[10px] text-[#5f6368] border-b border-[#dadce0]">
                      {c.amino_acid} codons
                    </div>
                    {alts.map((alt) => {
                      const altRank = codonRank(alt.per_thousand, alts);
                      const isCurrent = alt.codon === c.codon;
                      return (
                        <button
                          key={alt.codon}
                          onClick={(e) => {
                            e.stopPropagation();
                            if (!isCurrent) handleSwap(i, alt.codon);
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
                  </div>
                )}
              </div>
            );
          })}
        </div>
      </div>
      <div className="flex items-center gap-4 mt-3 text-xs text-[#5f6368]">
        <div className="flex items-center gap-1">
          <span className="inline-block w-3 h-3" style={{ backgroundColor: "#e6f4ea" }} />
          Best choice for this amino acid
        </div>
        <div className="flex items-center gap-1">
          <span className="inline-block w-3 h-3" style={{ backgroundColor: "#fce8e6" }} />
          Worst choice for this amino acid
        </div>
      </div>
    </div>
  );
}
