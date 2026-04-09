import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Oligo Designer",
  description: "Design overlapping oligos for gene synthesis via overlap extension assembly",
  icons: {
    icon: "/favicon.svg",
  },
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <head>
        <link rel="stylesheet" href="https://use.typekit.net/fos1ece.css" />
      </head>
      <body className="bg-[#f8f9fa] text-[#202124] min-h-screen" style={{ fontFamily: '"futura-100", "Futura", sans-serif' }}>
        {children}
      </body>
    </html>
  );
}
