"use client";

import { useState, useRef, useCallback, useLayoutEffect } from "react";
import { createPortal } from "react-dom";
import { Info } from "lucide-react";

interface InfoTipProps {
  text: string;
}

export default function InfoTip({ text }: InfoTipProps) {
  const [visible, setVisible] = useState(false);
  const [anchor, setAnchor] = useState({ x: 0, y: 0 });
  const [nudge, setNudge] = useState(0);
  const iconRef = useRef<HTMLSpanElement>(null);
  const tipRef = useRef<HTMLDivElement>(null);

  const show = useCallback(() => {
    if (iconRef.current) {
      const rect = iconRef.current.getBoundingClientRect();
      setAnchor({ x: rect.left + rect.width / 2, y: rect.top });
      setNudge(0);
    }
    setVisible(true);
  }, []);

  const hide = useCallback(() => setVisible(false), []);

  useLayoutEffect(() => {
    if (visible && tipRef.current) {
      const tipRect = tipRef.current.getBoundingClientRect();
      let shift = 0;
      if (tipRect.left < 8) {
        shift = 8 - tipRect.left;
      } else if (tipRect.right > window.innerWidth - 8) {
        shift = window.innerWidth - 8 - tipRect.right;
      }
      if (shift !== nudge) setNudge(shift);
    }
  }, [visible, anchor, nudge]);

  return (
    <>
      <span
        ref={iconRef}
        onMouseEnter={show}
        onMouseLeave={hide}
        className="inline-flex cursor-help"
      >
        <Info size={14} className="text-[#5f6368]" />
      </span>
      {visible &&
        createPortal(
          <div
            ref={tipRef}
            style={{
              position: "fixed",
              left: anchor.x + nudge,
              top: anchor.y,
              transform: "translate(-50%, -100%)",
              marginTop: -4,
              zIndex: 9999,
            }}
            className="px-2 py-1 text-xs text-white bg-[#3c4043] rounded shadow-lg max-w-xs"
          >
            {text}
          </div>,
          document.body
        )}
    </>
  );
}
