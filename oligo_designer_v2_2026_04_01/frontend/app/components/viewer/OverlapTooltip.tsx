import type { OverlapData } from "./types";
import { ISSUE_COLORS } from "./colors";

interface Props {
  overlap: OverlapData;
  onClose: () => void;
}

export default function OverlapTooltip({ overlap, onClose }: Props) {
  return (
    <div className="mt-3 bg-white border border-[#dadce0] p-4">
      <div className="flex items-center justify-between mb-3">
        <span className="font-medium text-[#202124]">Overlap {overlap.index}</span>
        <button
          onClick={onClose}
          className="text-[#5f6368] hover:text-[#202124] text-sm cursor-pointer"
        >
          close
        </button>
      </div>

      <table className="text-sm mb-3">
        <tbody>
          <tr>
            <td className="text-[#5f6368] pr-6 py-0.5">Position</td>
            <td className="py-0.5">{overlap.start}-{overlap.end}</td>
          </tr>
          <tr>
            <td className="text-[#5f6368] pr-6 py-0.5">GC</td>
            <td className="py-0.5">{(overlap.gc * 100).toFixed(0)}%</td>
          </tr>
          {overlap.tm != null && (
            <tr>
              <td className="text-[#5f6368] pr-6 py-0.5">Tm</td>
              <td className="py-0.5">{overlap.tm.toFixed(1)}°C</td>
            </tr>
          )}
          <tr>
            <td className="text-[#5f6368] pr-6 py-0.5 align-top">Sequence <span className="text-xs">(sense, 5&apos;-3&apos;)</span></td>
            <td className="py-0.5">
              <span className="font-mono text-xs break-all select-all">{overlap.seq}</span>
            </td>
          </tr>
        </tbody>
      </table>

      {overlap.issues.length > 0 ? (
        <div>
          <span className="text-[#5f6368] text-sm block mb-1">Issues</span>
          <div className="space-y-2">
            {overlap.issues.map((iss, i) => {
              const colorInfo = ISSUE_COLORS[iss.kind];
              const color = colorInfo?.fill ?? "#d93025";
              const typeLabel = colorInfo?.label?.split("(")[0]?.trim() ?? iss.kind;
              return (
                <div
                  key={i}
                  className="text-sm p-2 bg-[#f8f9fa] border-l-2"
                  style={{ borderLeftColor: color }}
                >
                  <div className="flex items-center gap-2 mb-0.5">
                    <span
                      className="inline-block w-2 h-2 flex-shrink-0"
                      style={{ backgroundColor: color }}
                    />
                    <span className="text-xs font-medium" style={{ color }}>
                      {typeLabel}
                    </span>
                    <span className="text-xs text-[#5f6368]">
                      {iss.severity === "error" ? "Error" : "Warning"}
                    </span>
                  </div>
                  <div className="text-[#202124]">{iss.message}</div>
                </div>
              );
            })}
          </div>
        </div>
      ) : (
        <div className="text-sm text-[#188038]">No issues</div>
      )}
    </div>
  );
}
