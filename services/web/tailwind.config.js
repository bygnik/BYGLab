/** @type {import('tailwindcss').Config} */
export default {
  content: ["./index.html", "./src/**/*.{js,jsx}"],
  theme: {
    extend: {
      fontFamily: {
        sans: ['"IBM Plex Sans"', "system-ui", "sans-serif"],
        mono: ['"IBM Plex Mono"', "ui-monospace", "monospace"],
        display: ['"Barlow Condensed"', "system-ui", "sans-serif"],
      },
      colors: {
        void: "#070A12",
        panel: "#101623",
        panel2: "#152138",
        line: "#2C3A5C",
        hair: "#273149",
        cyan: "#3A7CF6",
        amber: "#F2A93B",
        danger: "#E5644B",
        good: "#7FD858",
        hi: "#EBF1FF",
        lo: "#98A6C0",
      },
    },
  },
  plugins: [],
};
