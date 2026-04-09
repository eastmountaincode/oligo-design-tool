"use client";

interface CodonTableEntry {
  amino_acid: string;
  codon: string;
  per_thousand: number;
}

interface Props {
  entries: CodonTableEntry[];
  tableName: string;
}

const AA_NAMES: Record<string, string> = {
  A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys",
  E: "Glu", Q: "Gln", G: "Gly", H: "His", I: "Ile",
  L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro",
  S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val",
};

export default function CodonTableView({ entries, tableName }: Props) {
  if (!entries || entries.length === 0 || entries[0].per_thousand === undefined) return null;

  // Group by amino acid
  const grouped: Record<string, CodonTableEntry[]> = {};
  for (const e of entries) {
    if (!grouped[e.amino_acid]) grouped[e.amino_acid] = [];
    grouped[e.amino_acid].push(e);
  }

  // Find max per_thousand across all entries for scaling bars
  const maxVal = Math.max(...entries.map((e) => e.per_thousand));

  const aminoAcids = Object.keys(grouped).sort();

  return (
    <div className="bg-white border border-[#dadce0] p-4">
      <h2 className="text-lg font-medium mb-1 text-[#202124]">
        Codon Table
      </h2>
      <p className="text-xs text-[#5f6368] mb-3">{tableName} — values are per 1000 codons</p>
      <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-4 gap-3">
        {aminoAcids.map((aa) => (
          <div key={aa} className="border border-[#dadce0] p-2">
            <div className="text-xs font-medium text-[#202124] mb-1">
              {aa} <span className="text-[#5f6368]">{AA_NAMES[aa] ?? aa}</span>
            </div>
            <table className="w-full text-xs">
              <tbody>
                {grouped[aa].map((e) => (
                  <tr key={e.codon}>
                    <td className="font-mono pr-2 py-px">{e.codon}</td>
                    <td className="w-full py-px">
                      <div className="flex items-center gap-1">
                        <div
                          className="h-2 bg-[#1a73e8]"
                          style={{
                            width: `${Math.max((e.per_thousand / maxVal) * 100, 2)}%`,
                          }}
                        />
                        <span className="text-[#5f6368] text-[10px] whitespace-nowrap">
                          {e.per_thousand.toFixed(1)}
                        </span>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        ))}
      </div>
    </div>
  );
}
