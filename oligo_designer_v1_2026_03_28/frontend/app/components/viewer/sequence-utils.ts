import type { OligoData } from "./types";

const COMPLEMENT: Record<string, string> = {
  A: "T", T: "A", C: "G", G: "C", N: "N",
  a: "t", t: "a", c: "g", g: "c",
};

// High-contrast base colors — must be legible on light pastel backgrounds
export const BASE_COLORS: Record<string, string> = {
  A: "#188038", // dark green
  T: "#d93025", // dark red
  C: "#1a73e8", // dark blue
  G: "#e37400", // dark amber
  N: "#5f6368", // gray
};

export function reverseComplement(seq: string): string {
  return seq
    .split("")
    .reverse()
    .map((b) => COMPLEMENT[b] || b)
    .join("");
}

export function complementBase(base: string): string {
  return COMPLEMENT[base] || base;
}

export function reconstructSequence(
  oligos: OligoData[],
  totalLength: number
): string {
  const arr = new Array(totalLength).fill("N");
  for (const o of oligos) {
    const seq = o.strand === "antisense" ? reverseComplement(o.seq) : o.seq;
    for (let i = 0; i < seq.length; i++) {
      const pos = o.start + i;
      if (pos < totalLength) arr[pos] = seq[i];
    }
  }
  return arr.join("");
}
