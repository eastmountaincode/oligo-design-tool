"use client";

import { useState, useCallback } from "react";
import OligoDesignerPage from "./components/OligoDesignerPage";
import CodonOptPage from "./components/CodonOptPage";

export default function Home() {
  // Live DNA flows from DNA Chisel → Oligo Designer automatically
  const [liveDna, setLiveDna] = useState("");

  const handleDnaChanged = useCallback((dna: string) => {
    setLiveDna(dna);
  }, []);

  return (
    <div className="min-h-screen" style={{ overscrollBehavior: "none" }}>
      {/* Top bar */}
      <header className="bg-white border-b border-[#dadce0] px-6 py-3">
        <h1 className="text-sm tracking-wide uppercase text-[#5f6368]">Oligo Design Tool</h1>
      </header>

      {/* Side-by-side panels — fixed height, independent scroll */}
      <div className="flex h-[calc(100vh-49px)]">
        {/* Left: DNA Chisel */}
        <div className="flex-1 border-r border-[#dadce0] overflow-y-auto">
          <div className="px-6 py-6">
            <CodonOptPage onDnaChanged={handleDnaChanged} />
          </div>
        </div>

        {/* Right: Oligo Designer */}
        <div className="flex-1 overflow-y-auto">
          <div className="px-6 py-6">
            <OligoDesignerPage liveDna={liveDna} />
          </div>
        </div>
      </div>
    </div>
  );
}
