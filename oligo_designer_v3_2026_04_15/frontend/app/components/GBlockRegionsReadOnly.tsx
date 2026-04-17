"use client";

// Read-only visualization + detail list for gBlock regions returned by the
// beam search / auto-detector. Deliberately has no drag handles, no inputs,
// and no add/remove buttons — gBlocks are output of codon optimization now,
// not user-editable state. If you need to change them, rerun optimization
// with different constraints.
//
// The shape mirrors the `gblock_segments` / `suggested_gblocks` payload the
// backend emits on `CodonOptResult`, so this component can be fed directly
// from either field.

export interface ReadOnlyGBlock {
  // Stable gBlock number — matches the G{index+1} label emitted by the
  // oligo designer's Fragments table so both panels refer to the same
  // physical fragment. Assigned by the backend in start_bp order.
  index: number;
  start_bp: number;
  end_bp: number;
  label: string;
  len_bp: number;
  valid_size: boolean;
  gc?: number;
  max_g_run?: number;
  max_c_run?: number;
  start_aa?: number;
  end_aa?: number;
}

const GBLOCK_MIN_BP = 125;
const GBLOCK_MAX_BP = 3000;

interface Props {
  dnaLength: number;
  regions: ReadOnlyGBlock[];
}

export default function GBlockRegionsReadOnly({ dnaLength, regions }: Props) {
  if (dnaLength <= 0 || regions.length === 0) return null;

  const bpToFrac = (bp: number) => (dnaLength > 0 ? bp / dnaLength : 0);
  const totalBp = regions.reduce((s, r) => s + r.len_bp, 0);
  const totalCost = totalBp * 0.07;

  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <label className="text-sm font-medium text-[#202124]">
          gBlock Regions
          <span className="text-xs font-normal text-[#5f6368] ml-2">
            ({regions.length} region{regions.length === 1 ? "" : "s"}, {totalBp} bp, ~${totalCost.toFixed(2)})
          </span>
        </label>
        <span className="text-xs text-[#80868b]">read-only</span>
      </div>

      {/* Visual bar */}
      <div className="relative h-8 bg-[#f1f3f4] border border-[#dadce0] rounded select-none">
        {/* Tick marks every 300 bp */}
        {Array.from(
          { length: Math.floor(dnaLength / 300) },
          (_, i) => (i + 1) * 300
        ).map((bp) => (
          <div
            key={bp}
            className="absolute top-0 h-full border-l border-[#dadce0]"
            style={{ left: `${bpToFrac(bp) * 100}%` }}
          >
            <span className="absolute -top-3.5 -translate-x-1/2 text-[9px] text-[#80868b]">
              {bp}
            </span>
          </div>
        ))}

        {/* gBlock region rectangles */}
        {regions.map((r, i) => {
          const left = bpToFrac(r.start_bp) * 100;
          const width = bpToFrac(r.end_bp - r.start_bp) * 100;
          const tooSmall = r.len_bp < GBLOCK_MIN_BP;
          const tooLarge = r.len_bp > GBLOCK_MAX_BP;
          const invalid = tooSmall || tooLarge || !r.valid_size;
          const bgColor = invalid
            ? "rgba(234, 67, 53, 0.25)"
            : "rgba(66, 133, 244, 0.25)";
          const borderColor = invalid ? "#ea4335" : "#4285f4";

          return (
            <div
              key={`gb-${i}`}
              className="absolute top-0 h-full"
              style={{
                left: `${left}%`,
                width: `${width}%`,
                backgroundColor: bgColor,
                borderLeft: `2px solid ${borderColor}`,
                borderRight: `2px solid ${borderColor}`,
                minWidth: 4,
              }}
              title={`G${r.index + 1} · nt ${r.start_bp + 1}–${r.end_bp} · ${r.len_bp} bp`}
            >
              {width > 5 && (
                <span className="absolute inset-0 flex items-center justify-center text-[10px] font-medium text-[#202124] pointer-events-none overflow-hidden whitespace-nowrap">
                  G{r.index + 1}
                </span>
              )}
            </div>
          );
        })}
      </div>

      {/* Per-region detail table — columns align so values are easy to scan. */}
      <div className="pt-1 overflow-x-auto">
        <table className="text-xs font-mono">
          <thead>
            <tr className="text-left text-[#80868b] font-sans">
              <th className="pr-4 pb-1 font-normal"></th>
              <th className="pr-4 pb-1 font-normal">Label</th>
              <th className="pr-4 pb-1 font-normal">Amino acids</th>
              <th className="pr-4 pb-1 font-normal">Nucleotides</th>
              <th className="pr-4 pb-1 font-normal text-right">Length</th>
              <th className="pr-4 pb-1 font-normal text-right">GC</th>
              <th className="pb-1 font-normal text-right">Cost</th>
            </tr>
          </thead>
          <tbody>
            {regions.map((r, i) => {
              const tooSmall = r.len_bp < GBLOCK_MIN_BP;
              const tooLarge = r.len_bp > GBLOCK_MAX_BP;
              const invalid = tooSmall || tooLarge || !r.valid_size;
              return (
                <tr key={`row-${i}`} className="align-middle">
                  <td className="pr-2 py-0.5">
                    <span className="inline-block w-2 h-2 rounded-sm"
                          style={{ backgroundColor: invalid ? "#ea4335" : "#4285f4" }} />
                  </td>
                  <td className="pr-4 py-0.5 text-[#202124]">
                    G{r.index + 1}
                  </td>
                  <td className="pr-4 py-0.5 text-[#5f6368]">
                    {r.start_aa != null && r.end_aa != null
                      ? `${r.start_aa}–${r.end_aa}`
                      : "—"}
                  </td>
                  <td className="pr-4 py-0.5 text-[#5f6368]">
                    {r.start_bp + 1}–{r.end_bp}
                  </td>
                  <td className={`pr-4 py-0.5 text-right ${invalid ? "text-[#ea4335]" : "text-[#5f6368]"}`}>
                    {r.len_bp} bp
                    {tooSmall && ` (min ${GBLOCK_MIN_BP})`}
                    {tooLarge && ` (max ${GBLOCK_MAX_BP})`}
                  </td>
                  <td className="pr-4 py-0.5 text-right text-[#5f6368]">
                    {r.gc != null ? `${(r.gc * 100).toFixed(0)}%` : "—"}
                  </td>
                  <td className="py-0.5 text-right text-[#5f6368]">
                    ~${(r.len_bp * 0.07).toFixed(2)}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
