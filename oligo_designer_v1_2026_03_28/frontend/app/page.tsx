"use client";

import { useState } from "react";
import Sidebar, { type Tab } from "./components/Sidebar";
import OligoDesignerPage from "./components/OligoDesignerPage";
import CodonOptPage from "./components/CodonOptPage";

export default function Home() {
  const [activeTab, setActiveTab] = useState<Tab>("codon-opt");
  const [prefillDna, setPrefillDna] = useState("");

  function handleSendToTiler(dna: string) {
    setPrefillDna(dna);
    setActiveTab("oligo-designer");
  }

  return (
    <div className="flex min-h-screen" style={{ overscrollBehavior: "none" }}>
      <Sidebar activeTab={activeTab} onTabChange={setActiveTab} />
      <main className="flex-1 max-w-5xl px-8 py-8">
        {activeTab === "oligo-designer" && (
          <OligoDesignerPage prefillSequence={prefillDna} onPrefillConsumed={() => setPrefillDna("")} />
        )}
        {activeTab === "codon-opt" && (
          <CodonOptPage onSendToTiler={handleSendToTiler} />
        )}
      </main>
    </div>
  );
}
