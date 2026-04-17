export interface GBlockDetail {
  index: number;
  label: string;
  start: number;
  end: number;
  length: number;
  gc: number;
  seq: string;
}

interface Props {
  gblock: GBlockDetail;
  onClose: () => void;
}

export default function GBlockDetailPanel({ gblock, onClose }: Props) {
  return (
    <div className="mt-3 bg-white border border-[#dadce0] p-4">
      <div className="flex items-center justify-between mb-3">
        <h3 className="font-medium text-[#202124]">
          gBlock G{gblock.index + 1}
          {gblock.label && (
            <span className="ml-2 text-sm text-[#5f6368]">({gblock.label})</span>
          )}
          <span className="ml-2 text-sm text-[#34a853]">double-stranded</span>
        </h3>
        <button
          onClick={onClose}
          className="text-[#5f6368] hover:text-[#202124] text-sm cursor-pointer"
        >
          close
        </button>
      </div>

      <div className="grid grid-cols-2 md:grid-cols-4 gap-3 text-sm mb-3 text-[#202124]">
        <div>
          <span className="text-[#5f6368] block">
            Position <span className="text-xs text-[#9aa0a6]">(nt, 1-indexed)</span>
          </span>
          {gblock.start + 1} - {gblock.end}
        </div>
        <div>
          <span className="text-[#5f6368] block">Length</span>
          {gblock.length} bp
        </div>
        <div>
          <span className="text-[#5f6368] block">GC Content</span>
          {(gblock.gc * 100).toFixed(1)}%
        </div>
        <div>
          <span className="text-[#5f6368] block">Order as</span>
          gBlock / synthetic fragment
        </div>
      </div>

      <div>
        <span className="text-[#5f6368] text-sm block mb-1">
          Sequence (5&apos;-3&apos;)
        </span>
        <div className="font-mono text-xs bg-[#f8f9fa] border border-[#dadce0] p-3 break-all select-all text-[#202124]">
          {gblock.seq}
        </div>
      </div>
    </div>
  );
}
