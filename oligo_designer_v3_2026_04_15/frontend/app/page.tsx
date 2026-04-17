"use client";

import { useState, useCallback } from "react";
import OligoDesignerPage from "./components/OligoDesignerPage";
import CodonOptPage from "./components/CodonOptPage";
import type { RepairContext } from "./components/types";

interface GBlockInfo {
  start_bp: number;
  end_bp: number;
  label: string;
  len_bp: number;
}

export default function Home() {
  const [liveDna, setLiveDna] = useState("");
  const [liveGblocks, setLiveGblocks] = useState<GBlockInfo[]>([]);

  // Plasmid flanking sequences are lifted to the app level because they now
  // participate in BOTH panels: the codon optimizer uses them during beam
  // search / k-mer uniqueness checks (so problematic flank↔insert clashes
  // are avoided up front, not surfaced later as oligo cross-hybridization),
  // and the oligo designer uses them to define Gibson overlaps at the ends
  // of the insert. Keeping one source of truth avoids the two panels
  // drifting out of sync.
  const [plasmidUpstream, setPlasmidUpstream] = useState("");
  const [plasmidDownstream, setPlasmidDownstream] = useState("");

  // Repair context: codon-opt inputs the oligo designer's repair pass needs
  // (protein, codon table, constraint settings) to propose minimally-invasive
  // codon swaps that resolve oligo cross-hybridization. Populated by the
  // codon opt pane whenever an optimization completes; cleared when the
  // optimization becomes stale (inputs changed).
  const [repairContext, setRepairContext] = useState<RepairContext | null>(null);

  const handleDnaChanged = useCallback((dna: string) => {
    setLiveDna(dna);
  }, []);

  const handleGblocksChanged = useCallback((gblocks: GBlockInfo[]) => {
    setLiveGblocks(gblocks);
  }, []);

  const handleRepairContextChanged = useCallback((ctx: RepairContext | null) => {
    setRepairContext(ctx);
  }, []);

  return (
    <div className="min-h-screen" style={{ overscrollBehavior: "none" }}>
      {/* Top bar */}
      <header className="bg-white border-b border-[#dadce0] px-6 py-3">
        <h1 className="text-sm tracking-wide uppercase text-[#5f6368]">Oligo Design Tool</h1>
      </header>

      {/* Side-by-side panels — fixed height, independent scroll */}
      <div className="flex h-[calc(100vh-49px)]">
        {/* Left: Codon Optimization */}
        <div className="flex-1 border-r border-[#dadce0] overflow-y-auto">
          <div className="px-6 py-6">
            <CodonOptPage
              onDnaChanged={handleDnaChanged}
              onGblocksChanged={handleGblocksChanged}
              onRepairContextChanged={handleRepairContextChanged}
              plasmidUpstream={plasmidUpstream}
              plasmidDownstream={plasmidDownstream}
              setPlasmidUpstream={setPlasmidUpstream}
              setPlasmidDownstream={setPlasmidDownstream}
            />
          </div>
        </div>

        {/* Right: Oligo Designer */}
        <div className="flex-1 overflow-y-auto">
          <div className="px-6 py-6">
            <OligoDesignerPage
              liveDna={liveDna}
              liveGblocks={liveGblocks}
              plasmidUpstream={plasmidUpstream}
              plasmidDownstream={plasmidDownstream}
              repairContext={repairContext}
            />
          </div>
        </div>
      </div>
    </div>
  );
}
