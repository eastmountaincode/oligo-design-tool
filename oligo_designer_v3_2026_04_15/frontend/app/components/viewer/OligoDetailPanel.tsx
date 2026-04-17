import { oligoLabel, type OligoData } from "./types";

interface Props {
  oligo: OligoData;
  onClose: () => void;
}

export default function OligoDetailPanel({ oligo, onClose }: Props) {
  return (
    <div className="mt-3 bg-white border border-[#dadce0] p-4">
      <div className="flex items-center justify-between mb-3">
        <h3 className="font-medium text-[#202124]">
          Oligo {oligoLabel(oligo)}
          <span
            className={`ml-2 text-sm ${
              oligo.strand === "sense" ? "text-[#1a73e8]" : "text-[#e8710a]"
            }`}
          >
            {oligo.strand === "sense" ? "sense (5'-3')" : "antisense (3'-5')"}
          </span>
        </h3>
        <button
          onClick={onClose}
          className="text-[#5f6368] hover:text-[#202124] text-sm cursor-pointer"
        >
          close
        </button>
      </div>

      <div className="grid grid-cols-2 md:grid-cols-5 gap-3 text-sm mb-3 text-[#202124]">
        <div>
          <span className="text-[#5f6368] block">
            Position <span className="text-xs text-[#9aa0a6]">(nt, 1-indexed)</span>
          </span>
          {oligo.start + 1} - {oligo.end}
        </div>
        <div>
          <span className="text-[#5f6368] block">Length</span>
          {oligo.length} bp
        </div>
        <div>
          <span className="text-[#5f6368] block">GC Content</span>
          {(oligo.gc * 100).toFixed(1)}%
        </div>
        <div>
          <span className="text-[#5f6368] block">Role</span>
          {oligo.is_first
            ? "First (has upstream flank)"
            : oligo.is_last
            ? "Last (has downstream flank)"
            : "Internal"}
        </div>
        <div>
          <span className="text-[#5f6368] block">Order as</span>
          {oligo.length <= 45
            ? "Standard oligo"
            : oligo.length <= 60
            ? "Extended oligo"
            : "200-mer / gBlock"}
        </div>
      </div>

      <div>
        <span className="text-[#5f6368] text-sm block mb-1">
          Sequence (5&apos;-3&apos; as ordered)
        </span>
        <div className="font-mono text-xs bg-[#f8f9fa] border border-[#dadce0] p-3 break-all select-all text-[#202124]">
          {oligo.seq}
        </div>
      </div>
    </div>
  );
}
