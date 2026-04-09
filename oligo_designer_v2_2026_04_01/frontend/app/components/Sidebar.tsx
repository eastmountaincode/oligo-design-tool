"use client";

export type Tab = "oligo-designer" | "codon-opt";

const TABS: { id: Tab; label: string; sublabel: string }[] = [
  { id: "codon-opt", label: "DNA Chisel", sublabel: "Codon Optimization" },
  { id: "oligo-designer", label: "Oligo Designer", sublabel: "Tiling" },
];

interface SidebarProps {
  activeTab: Tab;
  onTabChange: (tab: Tab) => void;
}

export default function Sidebar({ activeTab, onTabChange }: SidebarProps) {
  return (
    <nav className="w-56 bg-white border-r border-[#dadce0] flex-shrink-0 h-screen sticky top-0 overflow-hidden">
      <div className="px-4 py-5 border-b border-[#dadce0]">
        <h1 className="text-sm tracking-wide uppercase text-[#5f6368]">Oligo Design Tool</h1>
      </div>
      <ul className="py-2">
        {TABS.map((tab) => (
          <li key={tab.id}>
            <button
              onClick={() => onTabChange(tab.id)}
              className={`w-full text-left px-4 py-2.5 text-sm transition-colors cursor-pointer ${
                activeTab === tab.id
                  ? "bg-[#e8f0fe] text-[#1a73e8] border-r-2 border-[#1a73e8]"
                  : "text-[#202124] hover:bg-[#f1f3f4]"
              }`}
            >
              <div className="font-medium">{tab.label}</div>
              <div className={`text-xs mt-0.5 ${
                activeTab === tab.id ? "text-[#1a73e8]/70" : "text-[#5f6368]"
              }`}>
                {tab.sublabel}
              </div>
            </button>
          </li>
        ))}
      </ul>
    </nav>
  );
}
