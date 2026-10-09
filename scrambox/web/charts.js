// Minimal SVG charts (no dependencies): line charts with linear/log axes, stacked areas,
// grouped bars; legend toggling and a hover read-out.
const Charts = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const PALETTE = ["#2f6fdf", "#e0712c", "#1d9a6c", "#c8394a", "#7b5cd6", "#c2a00e", "#16a3b8",
    "#d6589f", "#6b7a8f", "#8c5a2b", "#3fa34d", "#b04ad1", "#e05a5a", "#4a6fa5", "#a3a3a3"];
  const color = i => PALETTE[i % PALETTE.length];

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function fmt(v) {
    if (v === null || v === undefined || !isFinite(v)) return "–";
    const a = Math.abs(v);
    if (a === 0) return "0";
    if (a >= 1e5 || a < 1e-3) return v.toExponential(2);
    if (a >= 100) return v.toFixed(0);
    if (a >= 10) return v.toFixed(1);
    if (a >= 1) return v.toFixed(2);
    return v.toPrecision(3);
  }

  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) { hi = lo + 1; }
    const span = hi - lo, step0 = span / Math.max(n, 1);
    const mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const r = step0 / mag;
    const step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(+v.toPrecision(12));
    return out;
  }

  function logTicks(lo, hi) {
    const out = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
      const v = Math.pow(10, e);
      if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v);
    }
    if (out.length < 2) return [lo, hi];
    return out;
  }

  // series: [{name, x:[], y:[], color?, dash?, width?, area?}]
  function line(container, opt) {
    container.innerHTML = "";
    const W = opt.width || 560, H = opt.height || 300;
    const m = { l: 62, r: 16, t: opt.title ? 28 : 12, b: 42 };
    const wrap = document.createElement("div");
    wrap.className = "chart";
    container.appendChild(wrap);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}` }, wrap);
    const tip = document.createElement("div");
    tip.className = "tip";
    wrap.appendChild(tip);
    const hidden = new Set(opt.hidden || []);
    const series = (opt.series || []).map((s, i) => ({ ...s, color: s.color || color(i) }));

    function draw() {
      svg.innerHTML = "";
      const vis = series.filter(s => !hidden.has(s.name));
      let xs = [], ys = [];
      vis.forEach(s => s.x.forEach((x, i) => {
        const y = s.y[i];
        if (x !== null && isFinite(x) && (!opt.logx || x > 0)) xs.push(x);
        if (y !== null && isFinite(y) && (!opt.logy || y > 0)) ys.push(y);
      }));
      if (opt.stacked) {
        const n = vis.length ? vis[0].y.length : 0;
        for (let i = 0; i < n; i++) ys.push(vis.reduce((a, s) => a + (isFinite(s.y[i]) ? s.y[i] : 0), 0));
      }
      if (!xs.length) { el("text", { x: W / 2, y: H / 2, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.empty || "no data"; return; }
      let x0 = Math.min(...xs), x1 = Math.max(...xs);
      let y0 = opt.ymin !== undefined ? opt.ymin : Math.min(...ys), y1 = opt.ymax !== undefined ? opt.ymax : Math.max(...ys);
      if (!opt.logy && opt.ymin === undefined) y0 = Math.min(0, y0);
      if (opt.logy) { y0 = Math.max(y0, y1 * 1e-6); if (!(y1 > 0)) y1 = 1; }
      if (y1 === y0) { y1 = y0 === 0 ? 1 : y0 * 1.1; }
      if (x1 === x0) { x1 = x0 + 1; }
      const pad = opt.logy ? 1 : 0;
      const X = v => m.l + (opt.logx ? Math.log(v / x0) / Math.log(x1 / x0) : (v - x0) / (x1 - x0)) * (W - m.l - m.r);
      const Y = v => {
        const f = opt.logy ? Math.log(Math.max(v, y0) / y0) / Math.log(y1 / y0) : (v - y0) / (y1 - y0 || 1);
        return H - m.b - f * (H - m.t - m.b);
      };
      if (!opt.logy && !opt.ymax) { y1 = y1 + (y1 - y0) * 0.05; }
      // grid & axes
      const gx = el("g", { class: "axis" }, svg);
      const xt = opt.logx ? logTicks(x0, x1) : niceTicks(x0, x1, 6);
      const yt = opt.logy ? logTicks(y0, y1) : niceTicks(y0, y1, 5);
      yt.forEach(v => {
        el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: "gridl" }, gx);
        el("text", { x: m.l - 6, y: Y(v) + 4, "text-anchor": "end" }, gx).textContent = fmt(v);
      });
      xt.forEach(v => {
        el("line", { x1: X(v), x2: X(v), y1: H - m.b, y2: H - m.b + 4 }, gx);
        el("text", { x: X(v), y: H - m.b + 16, "text-anchor": "middle" }, gx).textContent = fmt(v);
      });
      el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b }, gx);
      if (opt.xlabel) el("text", { x: (m.l + W - m.r) / 2, y: H - 6, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.xlabel;
      if (opt.ylabel) el("text", { x: 14, y: (m.t + H - m.b) / 2, transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.ylabel;
      if (opt.title) el("text", { x: m.l, y: 17, class: "title" }, svg).textContent = opt.title;
      // series
      let base = null;
      vis.forEach(s => {
        const pts = [];
        s.x.forEach((x, i) => {
          let y = s.y[i];
          if (!(x !== null && isFinite(x)) || !(y !== null && isFinite(y))) { pts.push(null); return; }
          if (opt.stacked) y += base ? base[i] : 0;
          pts.push([X(x), Y(y)]);
        });
        if (opt.stacked) {
          const prev = base ? base.slice() : s.y.map(() => 0);
          base = s.y.map((v, i) => prev[i] + (isFinite(v) ? v : 0));
          const top = pts.filter(Boolean), bot = s.x.map((x, i) => [X(x), Y(Math.max(prev[i], opt.logy ? y0 : -Infinity))]).reverse();
          el("path", { d: "M" + top.concat(bot).map(p => p.join(",")).join("L") + "Z", fill: s.color, "fill-opacity": 0.75, stroke: s.color, "stroke-width": 0.8 }, svg);
          return;
        }
        let d = "", pen = false;
        pts.forEach(p => { if (!p) { pen = false; return; } d += (pen ? "L" : "M") + p[0].toFixed(1) + "," + p[1].toFixed(1); pen = true; });
        if (s.area) el("path", { d: d + `L${pts.filter(Boolean).slice(-1)[0][0]},${H - m.b}L${pts.filter(Boolean)[0][0]},${H - m.b}Z`, fill: s.color, "fill-opacity": 0.12 }, svg);
        el("path", { d, fill: "none", stroke: s.color, "stroke-width": s.width || 2, "stroke-dasharray": s.dash || "", "stroke-linejoin": "round" }, svg);
        if (s.markers || pts.length < 3) pts.filter(Boolean).forEach(p => el("circle", { cx: p[0], cy: p[1], r: 2.6, fill: s.color }, svg));
      });
      // hover
      const cross = el("line", { y1: m.t, y2: H - m.b, stroke: "currentColor", "stroke-opacity": 0.25, visibility: "hidden" }, svg);
      const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", ev => {
        const r = svg.getBoundingClientRect();
        const px = (ev.clientX - r.left) * W / r.width;
        const ref = vis[0];
        let bi = 0, bd = Infinity;
        ref.x.forEach((x, i) => { const d = Math.abs(X(x) - px); if (d < bd) { bd = d; bi = i; } });
        cross.setAttribute("x1", X(ref.x[bi])); cross.setAttribute("x2", X(ref.x[bi])); cross.setAttribute("visibility", "visible");
        tip.innerHTML = `<b>${opt.xlabel || "x"}: ${fmt(ref.x[bi])}</b><br>` + vis.map(s =>
          `<span style="color:${s.color}">■</span> ${s.name}: ${fmt(s.y[bi])}`).join("<br>");
        tip.style.display = "block";
        const left = (ev.clientX - r.left) + 14;
        tip.style.left = Math.min(left, r.width - tip.offsetWidth - 4) + "px";
        tip.style.top = (ev.clientY - r.top + 10) + "px";
      });
      hit.addEventListener("mouseleave", () => { tip.style.display = "none"; cross.setAttribute("visibility", "hidden"); });
    }
    draw();
    if (opt.legend !== false && series.length > 1) {
      const lg = document.createElement("div");
      lg.className = "legend";
      series.forEach(s => {
        const sp = document.createElement("span");
        sp.innerHTML = `<i style="background:${s.color}"></i>${s.name}`;
        if (hidden.has(s.name)) sp.classList.add("off");
        sp.onclick = () => { hidden.has(s.name) ? hidden.delete(s.name) : hidden.add(s.name); sp.classList.toggle("off"); draw(); };
        lg.appendChild(sp);
      });
      wrap.appendChild(lg);
    }
  }

  // groups: category labels; series: [{name, values:[]}] -> grouped (or stacked) bars
  function bars(container, opt) {
    container.innerHTML = "";
    const W = opt.width || 560, H = opt.height || 280;
    const m = { l: 62, r: 12, t: opt.title ? 28 : 12, b: opt.rotate ? 78 : 44 };
    const wrap = document.createElement("div");
    wrap.className = "chart";
    container.appendChild(wrap);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}` }, wrap);
    const tip = document.createElement("div"); tip.className = "tip"; wrap.appendChild(tip);
    const series = opt.series.map((s, i) => ({ ...s, color: s.color || color(i) }));
    const ng = opt.groups.length;
    const totals = opt.groups.map((_, g) => opt.stacked ? series.reduce((a, s) => a + Math.max(s.values[g] || 0, 0), 0)
      : Math.max(...series.map(s => s.values[g] || 0)));
    let y1 = Math.max(...totals, 0) * 1.05 || 1, y0 = 0;
    const Y = v => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
    const bw = (W - m.l - m.r) / Math.max(ng, 1);
    const gx = el("g", { class: "axis" }, svg);
    niceTicks(y0, y1, 5).forEach(v => {
      el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: "gridl" }, gx);
      el("text", { x: m.l - 6, y: Y(v) + 4, "text-anchor": "end" }, gx).textContent = fmt(v);
    });
    opt.groups.forEach((gname, g) => {
      const cx = m.l + bw * (g + 0.5);
      const t = el("text", opt.rotate ? { x: cx, y: H - m.b + 12, "text-anchor": "end", transform: `rotate(-35 ${cx} ${H - m.b + 12})` }
        : { x: cx, y: H - m.b + 16, "text-anchor": "middle" }, gx);
      t.textContent = gname;
      let acc = 0;
      series.forEach((s, i) => {
        const v = Math.max(s.values[g] || 0, 0);
        let x, w, ytop, ybot;
        if (opt.stacked) { w = bw * 0.62; x = cx - w / 2; ybot = Y(acc); acc += v; ytop = Y(acc); }
        else { w = bw * 0.8 / series.length; x = cx - bw * 0.4 + i * w; ybot = Y(0); ytop = Y(v); }
        const r = el("rect", { x, y: ytop, width: Math.max(w - 1, 1), height: Math.max(ybot - ytop, 0), fill: s.color, rx: 2 }, svg);
        r.addEventListener("mousemove", ev => {
          const b = svg.getBoundingClientRect();
          tip.innerHTML = `<b>${gname}</b><br>${s.name}: ${fmt(s.values[g])}`;
          tip.style.display = "block"; tip.style.left = (ev.clientX - b.left + 12) + "px"; tip.style.top = (ev.clientY - b.top + 8) + "px";
        });
        r.addEventListener("mouseleave", () => tip.style.display = "none");
      });
    });
    el("line", { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0) }, gx);
    if (opt.ylabel) el("text", { x: 14, y: (m.t + H - m.b) / 2, transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.ylabel;
    if (opt.title) el("text", { x: m.l, y: 17, class: "title" }, svg).textContent = opt.title;
    if (series.length > 1) {
      const lg = document.createElement("div"); lg.className = "legend";
      series.forEach(s => { const sp = document.createElement("span"); sp.innerHTML = `<i style="background:${s.color}"></i>${s.name}`; lg.appendChild(sp); });
      wrap.appendChild(lg);
    }
  }

  // 2-D histogram: points (x, y, weight) on log-x / linear-y bins, colour = log10(weight)
  function heat(container, opt) {
    container.innerHTML = "";
    const W = opt.width || 560, H = opt.height || 300, m = { l: 62, r: 70, t: 12, b: 42 };
    const nx = opt.nx || 40, ny = opt.ny || 25;
    const wrap = document.createElement("div"); wrap.className = "chart"; container.appendChild(wrap);
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}` }, wrap);
    const xs = opt.x, ys = opt.y, ws = opt.w;
    if (!xs.length) { el("text", { x: W / 2, y: H / 2, "text-anchor": "middle", class: "lbl" }, svg).textContent = "no data"; return; }
    const x0 = opt.xmin || Math.min(...xs.filter(v => v > 0)), x1 = opt.xmax || Math.max(...xs);
    const y0 = opt.ymin !== undefined ? opt.ymin : Math.min(...ys), y1 = opt.ymax !== undefined ? opt.ymax : Math.max(...ys);
    const lx0 = Math.log10(x0), lx1 = Math.log10(x1 * 1.0001);
    const grid = Array.from({ length: ny }, () => new Float64Array(nx));
    xs.forEach((x, i) => {
      if (!(x > 0)) return;
      const ix = Math.min(nx - 1, Math.max(0, Math.floor((Math.log10(x) - lx0) / (lx1 - lx0) * nx)));
      const iy = Math.min(ny - 1, Math.max(0, Math.floor((ys[i] - y0) / (y1 - y0 || 1) * ny)));
      grid[iy][ix] += ws[i];
    });
    let vmax = 0; grid.forEach(r => r.forEach(v => { if (v > vmax) vmax = v; }));
    const vmin = vmax * 1e-3;
    const cmap = f => {   // perceptual-ish ramp: dark blue -> teal -> yellow
      const stops = [[49, 54, 149], [69, 117, 180], [116, 173, 209], [171, 217, 233], [254, 224, 144], [253, 174, 97], [244, 109, 67], [215, 48, 39]];
      const p = Math.min(Math.max(f, 0), 1) * (stops.length - 1), i = Math.floor(p), w = p - i, a = stops[i], b = stops[Math.min(i + 1, stops.length - 1)];
      return `rgb(${a.map((v, k) => Math.round(v + w * (b[k] - v))).join(",")})`;
    };
    const pw = (W - m.l - m.r) / nx, ph = (H - m.t - m.b) / ny;
    const tip = document.createElement("div"); tip.className = "tip"; wrap.appendChild(tip);
    grid.forEach((row, iy) => row.forEach((v, ix) => {
      if (v <= vmin) return;
      const f = Math.log10(v / vmin) / Math.log10(vmax / vmin);
      const r = el("rect", { x: m.l + ix * pw, y: H - m.b - (iy + 1) * ph, width: pw + 0.5, height: ph + 0.5, fill: cmap(f) }, svg);
      r.addEventListener("mousemove", ev => {
        const b = svg.getBoundingClientRect();
        tip.innerHTML = `${opt.xlabel}: ${fmt(Math.pow(10, lx0 + (ix + 0.5) * (lx1 - lx0) / nx))}<br>${opt.ylabel}: ${fmt(y0 + (iy + 0.5) * (y1 - y0) / ny)}<br>${opt.wlabel || "N"}: ${fmt(v)}`;
        tip.style.display = "block"; tip.style.left = (ev.clientX - b.left + 12) + "px"; tip.style.top = (ev.clientY - b.top + 8) + "px";
      });
      r.addEventListener("mouseleave", () => tip.style.display = "none");
    }));
    const gx = el("g", { class: "axis" }, svg);
    for (let e = Math.ceil(lx0); e <= Math.floor(lx1); e++) {
      const X = m.l + (e - lx0) / (lx1 - lx0) * (W - m.l - m.r);
      el("line", { x1: X, x2: X, y1: H - m.b, y2: H - m.b + 4 }, gx);
      el("text", { x: X, y: H - m.b + 16, "text-anchor": "middle" }, gx).textContent = fmt(Math.pow(10, e));
    }
    [0, 0.25, 0.5, 0.75, 1].forEach(f => {
      const v = y0 + f * (y1 - y0), Y = H - m.b - f * (H - m.t - m.b);
      el("text", { x: m.l - 6, y: Y + 4, "text-anchor": "end" }, gx).textContent = fmt(v);
    });
    el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "none", stroke: "currentColor", "stroke-opacity": 0.2 }, svg);
    el("text", { x: (m.l + W - m.r) / 2, y: H - 6, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.xlabel;
    el("text", { x: 14, y: (m.t + H - m.b) / 2, transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})`, "text-anchor": "middle", class: "lbl" }, svg).textContent = opt.ylabel;
    for (let k = 0; k < 20; k++) el("rect", { x: W - m.r + 18, y: m.t + (19 - k) * (H - m.t - m.b) / 20, width: 12, height: (H - m.t - m.b) / 20 + 0.5, fill: cmap(k / 19) }, svg);
    el("text", { x: W - m.r + 34, y: m.t + 10, class: "lbl" }, svg).textContent = fmt(vmax);
    el("text", { x: W - m.r + 34, y: H - m.b, class: "lbl" }, svg).textContent = fmt(vmin);
  }

  return { line, bars, heat, color, fmt };
})();
