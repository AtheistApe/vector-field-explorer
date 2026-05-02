"use strict";

// ============================================================
// Vector Field Explorer — Math 253
// ============================================================

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
const readoutEl = document.getElementById("readout");
const cursorInfo = document.getElementById("cursor-info");
const parseErrorEl = document.getElementById("parse-error");

// ------------------------------------------------------------
// State
// ------------------------------------------------------------
const state = {
  // Field expressions
  Pexpr: "-y",
  Qexpr: "x",
  Pcompiled: null,
  Qcompiled: null,
  // Cached symbolic derivatives for analytic curl/div (when possible)
  dPdy: null,
  dQdx: null,
  dPdx: null,
  dQdy: null,
  hasAnalyticDerivs: false,

  // Plot range (world coordinates)
  xmin: -5, xmax: 5, ymin: -5, ymax: 5,

  // Display options
  density: 25,
  length: 1.0,
  showFlow: false,
  colorByMag: true,
  showGrid: true,
  normalize: false,

  // Animation
  animating: false,
  particles: [],
  numParticles: 800,
  speed: 1.0,
  time: 0,
  lastFrame: 0,

  // Tool: trace, paddle, divergence, probe
  tool: "trace",

  // Probes (placed objects)
  curves: [],          // {points: [...], color}
  paddles: [],         // {x, y, x0, y0, angle, omega, fixed}
  divRings: [],        // {x, y, r, flux, div}
  probes: [],          // {x, y}
  loops: [],           // {vertices: [[x,y]], closed, circulation, doubleIntegral}
  loopInProgress: null,// {vertices: [[x,y]]}
  potentialAt: null,   // {x, y, value} - clicked target for potential

  // Overlay
  overlay: "none",     // "none" | "curl" | "div" | "speed"
  conservativeFlag: false,

  // Hover
  hover: null,         // {x, y} during mouse move

  // View interaction
  isPanning: false,
  panStart: null,
};

// ------------------------------------------------------------
// Math expression parsing
// ------------------------------------------------------------
function compileField() {
  parseErrorEl.style.display = "none";
  document.getElementById("fp").classList.remove("error");
  document.getElementById("fq").classList.remove("error");

  const knownSyms = new Set([
    "x", "y", "t",
    "pi", "e", "tau", "phi",
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
    "exp", "log", "log2", "log10", "ln", "sqrt", "cbrt",
    "abs", "sign", "floor", "ceil", "round", "min", "max",
    "pow", "mod", "gcd", "lcm",
  ]);

  function checkUndefined(node, label) {
    let bad = null;
    node.traverse((n) => {
      if (bad) return;
      if (n.type === "SymbolNode" && !knownSyms.has(n.name)) {
        bad = n.name;
      }
      if (n.type === "FunctionNode" && !knownSyms.has(n.name) && !knownSyms.has(n.fn?.name)) {
        bad = n.name;
      }
    });
    if (bad) throw new Error(`${label}: unknown symbol "${bad}". Allowed: x, y, t, sin, cos, exp, sqrt, …`);
  }

  // Keep last-good compilation so a typo mid-edit doesn't blank the field
  let pNode, qNode, pCompiled, qCompiled;
  try {
    pNode = math.parse(state.Pexpr);
    qNode = math.parse(state.Qexpr);
    checkUndefined(pNode, "P");
    checkUndefined(qNode, "Q");
    pCompiled = pNode.compile();
    qCompiled = qNode.compile();
    pCompiled.evaluate({ x: 1, y: 1, t: 0 });
    qCompiled.evaluate({ x: 1, y: 1, t: 0 });
  } catch (err) {
    parseErrorEl.textContent = err.message;
    parseErrorEl.style.display = "block";
    if (/^P:/.test(err.message)) document.getElementById("fp").classList.add("error");
    else if (/^Q:/.test(err.message)) document.getElementById("fq").classList.add("error");
    else { document.getElementById("fp").classList.add("error"); document.getElementById("fq").classList.add("error"); }
    return false;
  }

  state.Pcompiled = pCompiled;
  state.Qcompiled = qCompiled;
  try {
    state.dPdx = math.derivative(pNode, "x").compile();
    state.dPdy = math.derivative(pNode, "y").compile();
    state.dQdx = math.derivative(qNode, "x").compile();
    state.dQdy = math.derivative(qNode, "y").compile();
    state.hasAnalyticDerivs = true;
  } catch (e) {
    state.hasAnalyticDerivs = false;
  }
  // Check conservativeness whenever the field changes
  if (typeof checkConservative === "function") checkConservative();
  return true;
}

function evalField(x, y, t) {
  try {
    const scope = { x, y, t: t || 0 };
    return [
      state.Pcompiled.evaluate(scope),
      state.Qcompiled.evaluate(scope),
    ];
  } catch (e) {
    return [0, 0];
  }
}

// Numerical partials (central difference) — fallback if symbolic fails
function numericalPartials(x, y) {
  const h = 1e-4;
  const [Pxp, Qxp] = evalField(x + h, y);
  const [Pxm, Qxm] = evalField(x - h, y);
  const [Pyp, Qyp] = evalField(x, y + h);
  const [Pym, Qym] = evalField(x, y - h);
  return {
    Px: (Pxp - Pxm) / (2 * h),
    Py: (Pyp - Pym) / (2 * h),
    Qx: (Qxp - Qxm) / (2 * h),
    Qy: (Qyp - Qym) / (2 * h),
  };
}

function partials(x, y) {
  if (state.hasAnalyticDerivs) {
    try {
      const scope = { x, y, t: state.time };
      return {
        Px: state.dPdx.evaluate(scope),
        Py: state.dPdy.evaluate(scope),
        Qx: state.dQdx.evaluate(scope),
        Qy: state.dQdy.evaluate(scope),
      };
    } catch (e) {
      return numericalPartials(x, y);
    }
  }
  return numericalPartials(x, y);
}

function curlAt(x, y) {
  const p = partials(x, y);
  return p.Qx - p.Py;
}
function divAt(x, y) {
  const p = partials(x, y);
  return p.Px + p.Qy;
}

// ------------------------------------------------------------
// Coordinate transforms
// ------------------------------------------------------------
function worldToScreen(x, y) {
  const w = cssWidth, h = cssHeight;
  const sx = (x - state.xmin) / (state.xmax - state.xmin) * w;
  const sy = h - (y - state.ymin) / (state.ymax - state.ymin) * h;
  return [sx, sy];
}
function screenToWorld(sx, sy) {
  const w = cssWidth, h = cssHeight;
  const x = state.xmin + sx / w * (state.xmax - state.xmin);
  const y = state.ymin + (h - sy) / h * (state.ymax - state.ymin);
  return [x, y];
}
function scaleX() { return cssWidth / (state.xmax - state.xmin); }
function scaleY() { return cssHeight / (state.ymax - state.ymin); }

// ------------------------------------------------------------
// Color mapping
// ------------------------------------------------------------
function magToColor(mag, refMag) {
  // Viridis-ish: dark blue → cyan → yellow
  const t = Math.min(1, mag / refMag);
  const r = Math.floor(68 + (253 - 68) * Math.pow(t, 1.2));
  const g = Math.floor(1 + (231 - 1) * Math.pow(t, 0.8));
  const b = Math.floor(84 + (37 - 84) * t);
  return `rgb(${r},${g},${b})`;
}

// Compute a reference magnitude (~95th percentile) for color scaling
function computeReferenceMag() {
  const samples = [];
  const n = 12;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const x = state.xmin + (i + 0.5) / n * (state.xmax - state.xmin);
      const y = state.ymin + (j + 0.5) / n * (state.ymax - state.ymin);
      const [u, v] = evalField(x, y, state.time);
      const m = Math.hypot(u, v);
      if (isFinite(m)) samples.push(m);
    }
  }
  if (samples.length === 0) return 1;
  samples.sort((a, b) => a - b);
  const p95 = samples[Math.floor(samples.length * 0.95)];
  return Math.max(p95, 1e-6);
}

// ------------------------------------------------------------
// Drawing — grid and axes
// ------------------------------------------------------------
function drawGrid() {
  const w = cssWidth, h = cssHeight;
  ctx.save();

  // Choose grid spacing based on view extent
  const xspan = state.xmax - state.xmin;
  const yspan = state.ymax - state.ymin;
  const span = Math.min(xspan, yspan);
  const targetLines = 12;
  const raw = span / targetLines;
  const pow10 = Math.pow(10, Math.floor(Math.log10(raw)));
  const candidates = [1, 2, 5, 10].map(c => c * pow10);
  let step = candidates[0];
  for (const c of candidates) if (Math.abs(c - raw) < Math.abs(step - raw)) step = c;

  // Minor grid
  ctx.strokeStyle = "rgba(255,255,255,0.04)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  const xstart = Math.ceil(state.xmin / step) * step;
  for (let x = xstart; x <= state.xmax; x += step) {
    const [sx] = worldToScreen(x, 0);
    ctx.moveTo(sx, 0); ctx.lineTo(sx, h);
  }
  const ystart = Math.ceil(state.ymin / step) * step;
  for (let y = ystart; y <= state.ymax; y += step) {
    const [, sy] = worldToScreen(0, y);
    ctx.moveTo(0, sy); ctx.lineTo(w, sy);
  }
  ctx.stroke();

  // Axes (if in view)
  ctx.strokeStyle = "rgba(255,255,255,0.25)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  if (0 >= state.xmin && 0 <= state.xmax) {
    const [sx] = worldToScreen(0, 0);
    ctx.moveTo(sx, 0); ctx.lineTo(sx, h);
  }
  if (0 >= state.ymin && 0 <= state.ymax) {
    const [, sy] = worldToScreen(0, 0);
    ctx.moveTo(0, sy); ctx.lineTo(w, sy);
  }
  ctx.stroke();

  // Tick labels
  ctx.fillStyle = "rgba(230,237,243,0.5)";
  ctx.font = "11px 'JetBrains Mono', monospace";
  ctx.textAlign = "center"; ctx.textBaseline = "top";
  const yAxisInView = (0 >= state.ymin && 0 <= state.ymax);
  const xAxisInView = (0 >= state.xmin && 0 <= state.xmax);
  const labelY = yAxisInView ? worldToScreen(0, 0)[1] + 3 : h - 14;
  const labelX = xAxisInView ? worldToScreen(0, 0)[0] + 3 : 4;
  for (let x = xstart; x <= state.xmax; x += step) {
    if (Math.abs(x) < 1e-9) continue;
    const [sx] = worldToScreen(x, 0);
    ctx.fillText(formatTick(x, step), sx, labelY);
  }
  ctx.textAlign = "left"; ctx.textBaseline = "middle";
  for (let y = ystart; y <= state.ymax; y += step) {
    if (Math.abs(y) < 1e-9) continue;
    const [, sy] = worldToScreen(0, y);
    ctx.fillText(formatTick(y, step), labelX, sy);
  }
  ctx.restore();
}

function formatTick(v, step) {
  const decimals = step < 1 ? Math.max(0, -Math.floor(Math.log10(step))) : 0;
  return v.toFixed(decimals);
}

// ------------------------------------------------------------
// Drawing — vector field
// ------------------------------------------------------------
function drawVectorField(refMag) {
  const w = cssWidth, h = cssHeight;
  const n = state.density;
  // Cell sizes in world coordinates
  const dx = (state.xmax - state.xmin) / n;
  const dy = (state.ymax - state.ymin) / n;
  // Maximum arrow length in pixels (cell size based, scaled by user)
  const cellPx = Math.min(scaleX() * dx, scaleY() * dy);
  const maxLenPx = cellPx * 0.85 * state.length;

  ctx.lineCap = "round";

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const x = state.xmin + (i + 0.5) * dx;
      const y = state.ymin + (j + 0.5) * dy;
      const [u, v] = evalField(x, y, state.time);
      if (!isFinite(u) || !isFinite(v)) continue;
      const mag = Math.hypot(u, v);
      if (mag < 1e-9) continue;

      let lenPx;
      if (state.normalize) {
        lenPx = maxLenPx;
      } else {
        // Scale by magnitude relative to refMag
        lenPx = maxLenPx * Math.min(1, mag / refMag);
      }

      const [sx, sy] = worldToScreen(x, y);
      // Direction in screen coords: dy is flipped
      const angle = Math.atan2(-v, u);
      const ex = sx + lenPx * Math.cos(angle);
      const ey = sy + lenPx * Math.sin(angle);

      const color = state.colorByMag ? magToColor(mag, refMag) : "#58a6ff";
      ctx.strokeStyle = color;
      ctx.fillStyle = color;
      ctx.lineWidth = Math.max(0.8, lenPx * 0.06);

      // Shaft
      ctx.beginPath();
      ctx.moveTo(sx, sy);
      ctx.lineTo(ex, ey);
      ctx.stroke();

      // Arrowhead — scales with arrow length
      const headLen = lenPx * 0.32;
      const headW = lenPx * 0.18;
      const cosA = Math.cos(angle), sinA = Math.sin(angle);
      const bx = ex - headLen * cosA;
      const by = ey - headLen * sinA;
      ctx.beginPath();
      ctx.moveTo(ex, ey);
      ctx.lineTo(bx + headW * sinA, by - headW * cosA);
      ctx.lineTo(bx - headW * sinA, by + headW * cosA);
      ctx.closePath();
      ctx.fill();
    }
  }
}

// ------------------------------------------------------------
// Streamlines (when "Flow lines" toggle is on)
// ------------------------------------------------------------
function rk4Step(x, y, dt, t) {
  const [k1u, k1v] = evalField(x, y, t);
  const [k2u, k2v] = evalField(x + 0.5 * dt * k1u, y + 0.5 * dt * k1v, t);
  const [k3u, k3v] = evalField(x + 0.5 * dt * k2u, y + 0.5 * dt * k2v, t);
  const [k4u, k4v] = evalField(x + dt * k3u, y + dt * k3v, t);
  return [
    x + dt * (k1u + 2 * k2u + 2 * k3u + k4u) / 6,
    y + dt * (k1v + 2 * k2v + 2 * k3v + k4v) / 6,
  ];
}

function integrateCurve(x0, y0, dir, maxSteps, refMag) {
  const span = Math.max(state.xmax - state.xmin, state.ymax - state.ymin);
  const dt = dir * span / 600 / Math.max(refMag, 0.5);
  const points = [[x0, y0]];
  let x = x0, y = y0;
  const margin = (state.xmax - state.xmin) * 0.2;
  for (let i = 0; i < maxSteps; i++) {
    const [u, v] = evalField(x, y, state.time);
    if (!isFinite(u) || !isFinite(v)) break;
    const mag = Math.hypot(u, v);
    if (mag < 1e-6) break;
    [x, y] = rk4Step(x, y, dt, state.time);
    if (!isFinite(x) || !isFinite(y)) break;
    if (x < state.xmin - margin || x > state.xmax + margin) break;
    if (y < state.ymin - margin || y > state.ymax + margin) break;
    points.push([x, y]);
  }
  return points;
}

function drawStreamlines(refMag) {
  // Seed a roughly uniform grid of streamlines, but skip cells already near another line
  const n = Math.max(8, Math.floor(state.density * 0.7));
  const occupancy = new Uint8Array(n * n);
  const cellW = cssWidth / n, cellH = cssHeight / n;

  function markCells(points) {
    for (const [x, y] of points) {
      const [sx, sy] = worldToScreen(x, y);
      const ci = Math.floor(sx / cellW);
      const cj = Math.floor(sy / cellH);
      if (ci >= 0 && ci < n && cj >= 0 && cj < n) occupancy[cj * n + ci] = 1;
    }
  }

  ctx.lineWidth = 1.3;
  ctx.lineCap = "round";

  // Visit cells in shuffled order so seeding doesn't bias top-left
  const order = [];
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) order.push([i, j]);
  for (let i = order.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [order[i], order[j]] = [order[j], order[i]];
  }

  for (const [i, j] of order) {
    if (occupancy[j * n + i]) continue;
    const sx = (i + 0.5) * cellW, sy = (j + 0.5) * cellH;
    const [x0, y0] = screenToWorld(sx, sy);
    const fwd = integrateCurve(x0, y0, +1, 500, refMag);
    const bwd = integrateCurve(x0, y0, -1, 500, refMag);
    const all = bwd.slice().reverse().concat(fwd.slice(1));
    if (all.length < 3) continue;
    markCells(all);

    // Draw with magnitude-based color along the line
    ctx.beginPath();
    let prev = null;
    for (let k = 0; k < all.length; k++) {
      const [x, y] = all[k];
      const [px, py] = worldToScreen(x, y);
      if (k === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
      prev = [px, py];
    }
    if (state.colorByMag) {
      const mid = all[Math.floor(all.length / 2)];
      const [u, v] = evalField(mid[0], mid[1], state.time);
      const m = Math.hypot(u, v);
      ctx.strokeStyle = magToColor(m, refMag);
    } else {
      ctx.strokeStyle = "rgba(126, 231, 135, 0.7)";
    }
    ctx.stroke();
  }

  // Add small arrows along streamlines for direction
  ctx.fillStyle = "rgba(255,255,255,0.7)";
  for (const [i, j] of order) {
    if (!occupancy[j * n + i]) continue;
    if ((i + j) % 3 !== 0) continue;
    const sx = (i + 0.5) * cellW, sy = (j + 0.5) * cellH;
    const [x, y] = screenToWorld(sx, sy);
    const [u, v] = evalField(x, y, state.time);
    const m = Math.hypot(u, v);
    if (m < 1e-9 || !isFinite(m)) continue;
    const angle = Math.atan2(-v, u);
    const sz = 4;
    ctx.save();
    ctx.translate(sx, sy);
    ctx.rotate(angle);
    ctx.beginPath();
    ctx.moveTo(sz, 0);
    ctx.lineTo(-sz * 0.6, sz * 0.5);
    ctx.lineTo(-sz * 0.6, -sz * 0.5);
    ctx.closePath();
    ctx.fill();
    ctx.restore();
  }
}

// ------------------------------------------------------------
// Particle animation (flow as water)
// ------------------------------------------------------------
function spawnParticle(p) {
  p.x = state.xmin + Math.random() * (state.xmax - state.xmin);
  p.y = state.ymin + Math.random() * (state.ymax - state.ymin);
  p.life = Math.random() * 3 + 1;
  p.maxLife = p.life;
  // Trail
  p.trail = [[p.x, p.y]];
}

function ensureParticles() {
  while (state.particles.length < state.numParticles) {
    const p = {};
    spawnParticle(p);
    p.life = Math.random() * 4;
    state.particles.push(p);
  }
  while (state.particles.length > state.numParticles) state.particles.pop();
}

function updateParticles(dt, refMag) {
  const margin = (state.xmax - state.xmin) * 0.05;
  const stepScale = state.speed / Math.max(refMag, 0.5);
  for (const p of state.particles) {
    const [u, v] = evalField(p.x, p.y, state.time);
    if (!isFinite(u) || !isFinite(v)) { spawnParticle(p); continue; }
    p.x += u * dt * stepScale;
    p.y += v * dt * stepScale;
    p.life -= dt;
    if (p.x < state.xmin - margin || p.x > state.xmax + margin ||
        p.y < state.ymin - margin || p.y > state.ymax + margin ||
        p.life <= 0) {
      spawnParticle(p);
      continue;
    }
    p.trail.push([p.x, p.y]);
    if (p.trail.length > 12) p.trail.shift();
  }
}

function drawParticles(refMag) {
  ctx.lineWidth = 1.4;
  ctx.lineCap = "round";
  for (const p of state.particles) {
    if (p.trail.length < 2) continue;
    const [u, v] = evalField(p.x, p.y, state.time);
    const m = Math.hypot(u, v);
    const col = state.colorByMag ? magToColor(m, refMag) : "#7ee787";
    // Fade in then fade out
    const lifeFrac = p.life / p.maxLife;
    const alpha = Math.min(1, 1 - Math.abs(lifeFrac - 0.5) * 2) * 0.85;
    ctx.strokeStyle = col;
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    for (let i = 0; i < p.trail.length; i++) {
      const [px, py] = worldToScreen(p.trail[i][0], p.trail[i][1]);
      if (i === 0) ctx.moveTo(px, py);
      else ctx.lineTo(px, py);
    }
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
}

// ------------------------------------------------------------
// Probes — curves, paddle wheels, divergence rings
// ------------------------------------------------------------
function drawCurves() {
  for (const c of state.curves) {
    ctx.strokeStyle = c.color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (let i = 0; i < c.points.length; i++) {
      const [x, y] = c.points[i];
      const [sx, sy] = worldToScreen(x, y);
      if (i === 0) ctx.moveTo(sx, sy);
      else ctx.lineTo(sx, sy);
    }
    ctx.stroke();
    // Mark start point
    const [x0, y0] = c.start;
    const [sx0, sy0] = worldToScreen(x0, y0);
    ctx.fillStyle = "white";
    ctx.beginPath(); ctx.arc(sx0, sy0, 3.5, 0, Math.PI * 2); ctx.fill();
    ctx.strokeStyle = "black"; ctx.lineWidth = 1; ctx.stroke();
  }
}

function drawPaddles() {
  for (const p of state.paddles) {
    const [sx, sy] = worldToScreen(p.x, p.y);
    const r = 22;

    // Outer ring
    ctx.strokeStyle = "rgba(255,216,102,0.5)";
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(sx, sy, r, 0, Math.PI * 2); ctx.stroke();

    // Paddle blades — 4 blades rotated by p.angle
    ctx.save();
    ctx.translate(sx, sy);
    ctx.rotate(-p.angle); // screen rotation is opposite of math
    ctx.strokeStyle = "#ffd866";
    ctx.fillStyle = "rgba(255,216,102,0.18)";
    ctx.lineWidth = 2.5;
    for (let k = 0; k < 4; k++) {
      ctx.save();
      ctx.rotate(k * Math.PI / 2);
      ctx.beginPath();
      ctx.moveTo(0, 0);
      ctx.lineTo(r - 3, 0);
      ctx.stroke();
      ctx.fillRect(r - 8, -3.5, 6, 7);
      ctx.strokeRect(r - 8, -3.5, 6, 7);
      ctx.restore();
    }
    // Center hub
    ctx.beginPath(); ctx.arc(0, 0, 3, 0, Math.PI * 2);
    ctx.fillStyle = "#ffd866"; ctx.fill();
    ctx.restore();

    // Readout: curl, optional translational speed
    const curl = curlAt(p.x, p.y);
    const [u, v] = evalField(p.x, p.y, state.time);
    const speed = Math.hypot(u, v);
    const lines = [`ω = ${curl.toFixed(3)}`];
    if (!p.fixed) lines.push(`|v| = ${speed.toFixed(3)}`);

    ctx.font = "11px 'JetBrains Mono', monospace";
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
    const pinSize = 18;
    const boxX = sx + r + 4, boxY = sy - r - 2;
    const boxW = maxW + pinSize + 14, boxH = lines.length * 14 + 6;

    // Highlight box on hover, so user can see they're about to interact
    p._hovered = state.hover && (
      // Mouse in the readout box
      (Math.abs((cssMouseX ?? -1e9) - (boxX + boxW/2)) < boxW/2 + 4 &&
       Math.abs((cssMouseY ?? -1e9) - (boxY + boxH/2)) < boxH/2 + 4)
      ||
      // Mouse over the paddle disc
      (Math.hypot((cssMouseX ?? -1e9) - sx, (cssMouseY ?? -1e9) - sy) < r + 8)
    );

    ctx.fillStyle = p._hovered ? "rgba(40,55,75,0.95)" : "rgba(15,20,25,0.85)";
    ctx.fillRect(boxX, boxY, boxW, boxH);
    if (p._hovered) {
      ctx.strokeStyle = "rgba(255,216,102,0.6)";
      ctx.lineWidth = 1;
      ctx.strokeRect(boxX, boxY, boxW, boxH);
    }
    ctx.fillStyle = "#ffd866";
    lines.forEach((l, i) => ctx.fillText(l, boxX + 5, boxY + 3 + i * 14));

    // Pin/lock toggle button — top-right of the readout box
    const pinX = boxX + boxW - pinSize - 3;
    const pinY = boxY + (boxH - pinSize) / 2;
    p._pinHitbox = [pinX, pinY, pinSize, pinSize];

    ctx.save();
    // Subtle background circle to make it feel like a button
    ctx.translate(pinX + pinSize / 2, pinY + pinSize / 2);
    if (p.fixed) {
      ctx.fillStyle = "rgba(255,123,114,0.25)";
      ctx.beginPath(); ctx.arc(0, 0, pinSize/2 - 1, 0, 2*Math.PI); ctx.fill();
      ctx.fillStyle = "#ff7b72";
      ctx.strokeStyle = "#ff7b72";
    } else {
      ctx.fillStyle = p._hovered ? "rgba(255,216,102,0.3)" : "rgba(255,216,102,0.15)";
      ctx.beginPath(); ctx.arc(0, 0, pinSize/2 - 1, 0, 2*Math.PI); ctx.fill();
      ctx.fillStyle = "rgba(255,216,102,0.3)";
      ctx.strokeStyle = p._hovered ? "#ffd866" : "rgba(255,216,102,0.8)";
    }
    ctx.lineWidth = 1.4;
    // Thumbtack: triangle head on top, line below
    ctx.beginPath();
    ctx.moveTo(0, -5);
    ctx.lineTo(4, -1);
    ctx.lineTo(-4, -1);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(0, -1);
    ctx.lineTo(0, 5);
    ctx.lineWidth = 1.6;
    ctx.stroke();
    ctx.restore();

    // Trail showing where this paddle has been (only when free-flowing)
    if (!p.fixed && p._trail && p._trail.length > 2) {
      ctx.strokeStyle = "rgba(255,216,102,0.4)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let i = 0; i < p._trail.length; i++) {
        const [tx, ty] = p._trail[i];
        const [tsx, tsy] = worldToScreen(tx, ty);
        if (i === 0) ctx.moveTo(tsx, tsy);
        else ctx.lineTo(tsx, tsy);
      }
      ctx.stroke();
    }
  }
}

function drawDivRings() {
  for (const ring of state.divRings) {
    const [sx, sy] = worldToScreen(ring.x, ring.y);
    const rPx = ring.r * scaleX();

    // Compute flux through the ring and divergence (recompute each frame; cheap)
    const samples = 64;
    let flux = 0;
    for (let i = 0; i < samples; i++) {
      const theta = (i / samples) * 2 * Math.PI;
      const px = ring.x + ring.r * Math.cos(theta);
      const py = ring.y + ring.r * Math.sin(theta);
      const [u, v] = evalField(px, py, state.time);
      // Outward normal is (cos theta, sin theta)
      flux += (u * Math.cos(theta) + v * Math.sin(theta));
    }
    flux *= (2 * Math.PI * ring.r) / samples;
    const div = divAt(ring.x, ring.y);
    ring.flux = flux; ring.div = div;

    // Color: green for outflow (positive), red for inflow (negative)
    const intensity = Math.min(1, Math.abs(div) / 3);
    const col = div >= 0 ? `rgba(126,231,135,${0.25 + 0.5*intensity})`
                         : `rgba(255,123,114,${0.25 + 0.5*intensity})`;
    const stroke = div >= 0 ? "#7ee787" : "#ff7b72";

    ctx.fillStyle = col;
    ctx.strokeStyle = stroke;
    ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.arc(sx, sy, rPx, 0, Math.PI * 2); ctx.fill(); ctx.stroke();

    // Outward normal arrows scaled by F·n
    const arrows = 16;
    for (let i = 0; i < arrows; i++) {
      const theta = (i / arrows) * 2 * Math.PI;
      const px = ring.x + ring.r * Math.cos(theta);
      const py = ring.y + ring.r * Math.sin(theta);
      const [u, v] = evalField(px, py, state.time);
      const fdotn = u * Math.cos(theta) + v * Math.sin(theta);
      const arrowLen = Math.max(-15, Math.min(15, fdotn * 6));
      const startX = sx + rPx * Math.cos(theta);
      const startY = sy - rPx * Math.sin(theta);
      const endX = startX + arrowLen * Math.cos(theta);
      const endY = startY - arrowLen * Math.sin(theta);
      ctx.strokeStyle = fdotn >= 0 ? "#7ee787" : "#ff7b72";
      ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(startX, startY); ctx.lineTo(endX, endY); ctx.stroke();
    }

    // Center label
    ctx.font = "11px 'JetBrains Mono', monospace";
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    const txt = `∇·F = ${div.toFixed(3)}`;
    const tw = ctx.measureText(txt).width;
    ctx.fillStyle = "rgba(15,20,25,0.85)";
    ctx.fillRect(sx + rPx + 4, sy - rPx - 2, tw + 8, 16);
    ctx.fillStyle = stroke;
    ctx.fillText(txt, sx + rPx + 8, sy - rPx);
  }
}

function drawProbes() {
  for (const probe of state.probes) {
    const [sx, sy] = worldToScreen(probe.x, probe.y);
    const [u, v] = evalField(probe.x, probe.y, state.time);
    const mag = Math.hypot(u, v);
    const angle = Math.atan2(-v, u);
    const lenPx = 50;
    ctx.strokeStyle = "#58a6ff";
    ctx.fillStyle = "#58a6ff";
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(sx, sy);
    ctx.lineTo(sx + lenPx * Math.cos(angle), sy + lenPx * Math.sin(angle));
    ctx.stroke();
    // Arrowhead
    const ex = sx + lenPx * Math.cos(angle), ey = sy + lenPx * Math.sin(angle);
    const headLen = 10, headW = 6;
    const cosA = Math.cos(angle), sinA = Math.sin(angle);
    ctx.beginPath();
    ctx.moveTo(ex, ey);
    ctx.lineTo(ex - headLen * cosA + headW * sinA, ey - headLen * sinA - headW * cosA);
    ctx.lineTo(ex - headLen * cosA - headW * sinA, ey - headLen * sinA + headW * cosA);
    ctx.closePath();
    ctx.fill();

    // Crosshair
    ctx.fillStyle = "white";
    ctx.beginPath(); ctx.arc(sx, sy, 3, 0, Math.PI * 2); ctx.fill();
    ctx.strokeStyle = "black"; ctx.lineWidth = 1; ctx.stroke();

    // Label
    const lines = [
      `(${probe.x.toFixed(2)}, ${probe.y.toFixed(2)})`,
      `F = (${u.toFixed(2)}, ${v.toFixed(2)})`,
      `|F| = ${mag.toFixed(3)}`,
    ];
    ctx.font = "11px 'JetBrains Mono', monospace";
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
    ctx.fillStyle = "rgba(15,20,25,0.9)";
    ctx.fillRect(sx + 8, sy + 8, maxW + 12, lines.length * 14 + 6);
    ctx.fillStyle = "#e6edf3";
    lines.forEach((l, i) => ctx.fillText(l, sx + 14, sy + 12 + i * 14));
  }
}

// ------------------------------------------------------------
// Background overlay (curl / div / |F| heatmap)
// ------------------------------------------------------------
const overlayCanvas = document.createElement("canvas");
const overlayCtx = overlayCanvas.getContext("2d");

function drawOverlay() {
  if (state.overlay === "none") return;
  const N = 80;
  overlayCanvas.width = N;
  overlayCanvas.height = N;
  const img = overlayCtx.createImageData(N, N);
  const dx = (state.xmax - state.xmin) / N;
  const dy = (state.ymax - state.ymin) / N;

  const vals = new Float32Array(N * N);
  for (let j = 0; j < N; j++) {
    const y = state.ymin + (j + 0.5) * dy;
    for (let i = 0; i < N; i++) {
      const x = state.xmin + (i + 0.5) * dx;
      let v;
      if (state.overlay === "curl") v = curlAt(x, y);
      else if (state.overlay === "div") v = divAt(x, y);
      else { const [u, w] = evalField(x, y, state.time); v = Math.hypot(u, w); }
      if (!isFinite(v)) v = 0;
      vals[j * N + i] = v;
    }
  }
  // Clamp at 95th percentile so a single huge spike (singularity) doesn't wash out everything else
  const sorted = Array.from(vals).map(Math.abs).sort((a, b) => a - b);
  const scale = Math.max(sorted[Math.floor(sorted.length * 0.95)] || 1, 1e-6);

  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const v = vals[j * N + i];
      const idx = ((N - 1 - j) * N + i) * 4;
      let r, g, b;
      if (state.overlay === "speed") {
        const t = Math.min(1, v / scale);
        r = Math.floor(68 + (253 - 68) * Math.pow(t, 1.2));
        g = Math.floor(1 + (231 - 1) * Math.pow(t, 0.8));
        b = Math.floor(84 + (37 - 84) * t);
      } else {
        const t = Math.max(-1, Math.min(1, v / scale));
        if (t >= 0) { r = 255; g = Math.floor(255 - 175 * t); b = Math.floor(255 - 200 * t); }
        else { r = Math.floor(255 - 200 * (-t)); g = Math.floor(255 - 175 * (-t)); b = 255; }
      }
      img.data[idx] = r;
      img.data[idx + 1] = g;
      img.data[idx + 2] = b;
      img.data[idx + 3] = 200;
    }
  }
  overlayCtx.putImageData(img, 0, 0);

  ctx.save();
  ctx.imageSmoothingEnabled = true;
  ctx.imageSmoothingQuality = "high";
  ctx.globalAlpha = 0.55;
  ctx.drawImage(overlayCanvas, 0, 0, cssWidth, cssHeight);
  ctx.restore();

  drawColorbar(scale);
}

function drawColorbar(scale) {
  const w = 14, h = 140;
  const x = cssWidth - w - 18, y = 18;
  ctx.save();
  for (let i = 0; i < h; i++) {
    const t = 1 - 2 * (i / (h - 1));
    let r, g, b;
    if (state.overlay === "speed") {
      const tt = (t + 1) / 2;
      r = Math.floor(68 + (253 - 68) * Math.pow(tt, 1.2));
      g = Math.floor(1 + (231 - 1) * Math.pow(tt, 0.8));
      b = Math.floor(84 + (37 - 84) * tt);
    } else if (t >= 0) {
      r = 255; g = Math.floor(255 - 175 * t); b = Math.floor(255 - 200 * t);
    } else {
      r = Math.floor(255 - 200 * (-t)); g = Math.floor(255 - 175 * (-t)); b = 255;
    }
    ctx.fillStyle = `rgb(${r},${g},${b})`;
    ctx.fillRect(x, y + i, w, 1);
  }
  ctx.strokeStyle = "rgba(255,255,255,0.3)";
  ctx.lineWidth = 0.5;
  ctx.strokeRect(x, y, w, h);

  ctx.fillStyle = "#e6edf3";
  ctx.font = "10px 'JetBrains Mono', monospace";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  if (state.overlay === "speed") {
    ctx.fillText(scale.toFixed(2), x - 4, y);
    ctx.fillText("0", x - 4, y + h);
  } else {
    ctx.fillText("+" + scale.toFixed(2), x - 4, y);
    ctx.fillText("0", x - 4, y + h / 2);
    ctx.fillText("-" + scale.toFixed(2), x - 4, y + h);
  }
  ctx.textAlign = "center";
  ctx.fillStyle = "rgba(230,237,243,0.7)";
  const labels = { curl: "curl", div: "div F", speed: "|F|" };
  ctx.fillText(labels[state.overlay] || "", x + w / 2, y + h + 14);
  ctx.restore();
}

// ------------------------------------------------------------
// Conservative-field detector
// ------------------------------------------------------------
function checkConservative() {
  const N = 12;
  let maxAbs = 0, totalAbs = 0, count = 0;
  for (let j = 0; j < N; j++) {
    for (let i = 0; i < N; i++) {
      const x = state.xmin + (i + 0.5) / N * (state.xmax - state.xmin);
      const y = state.ymin + (j + 0.5) / N * (state.ymax - state.ymin);
      const c = curlAt(x, y);
      if (!isFinite(c)) continue;
      maxAbs = Math.max(maxAbs, Math.abs(c));
      totalAbs += Math.abs(c);
      count++;
    }
  }
  const meanAbs = count > 0 ? totalAbs / count : 0;
  const refMag = computeReferenceMag();
  const span = Math.max(state.xmax - state.xmin, state.ymax - state.ymin);
  // Tolerance scales with characteristic 1/length, so curl ~ refMag/span is "order 1" relative
  const relTol = 1e-4;
  const isCons = maxAbs < relTol * refMag / span * 1000 && meanAbs < relTol * refMag / span * 100;

  state.conservativeFlag = isCons;
  const badge = document.getElementById("conservative-badge");
  if (!badge) return;
  if (isCons) {
    badge.style.display = "block";
    badge.style.background = "rgba(126,231,135,0.15)";
    badge.style.color = "#7ee787";
    badge.style.border = "1px solid rgba(126,231,135,0.4)";
    badge.textContent = "✓ curl ≈ 0 — likely conservative. Use Potential f tool to find f.";
  } else {
    badge.style.display = "block";
    badge.style.background = "rgba(255,123,114,0.1)";
    badge.style.color = "#ff7b72";
    badge.style.border = "1px solid rgba(255,123,114,0.3)";
    badge.textContent = `✗ not conservative (max |curl| ≈ ${maxAbs.toFixed(3)})`;
  }
}

// ------------------------------------------------------------
// Potential f(x,y) such that ∇f = F (when conservative)
// f(P) = ∫₀^P F·dr along L-shape: (0,0) → (x,0) → (x,y)
// ------------------------------------------------------------
function potentialAt(targetX, targetY) {
  function simpson(f, a, b, n) {
    if (a === b) return 0;
    if (n % 2 === 1) n++;
    if (n < 2) n = 2;
    const h = (b - a) / n;
    let s = f(a) + f(b);
    for (let i = 1; i < n; i++) s += (i % 2 === 0 ? 2 : 4) * f(a + i * h);
    return s * h / 3;
  }
  const I1 = simpson((t) => evalField(t, 0, state.time)[0], 0, targetX, 200);
  const I2 = simpson((s) => evalField(targetX, s, state.time)[1], 0, targetY, 200);
  return I1 + I2;
}

// ------------------------------------------------------------
// Circulation loops (Green's theorem)
// ------------------------------------------------------------
function circulationOfLoop(vertices) {
  const segPerEdge = 30;
  let total = 0;
  for (let i = 0; i < vertices.length; i++) {
    const [x1, y1] = vertices[i];
    const [x2, y2] = vertices[(i + 1) % vertices.length];
    const ex = x2 - x1, ey = y2 - y1;
    const dx = ex / segPerEdge, dy = ey / segPerEdge;
    for (let k = 0; k < segPerEdge; k++) {
      const tm = (k + 0.5) / segPerEdge;
      const xm = x1 + ex * tm, ym = y1 + ey * tm;
      const [u, v] = evalField(xm, ym, state.time);
      total += u * dx + v * dy;
    }
  }
  return total;
}

function pointInPolygon(x, y, verts) {
  let inside = false;
  for (let i = 0, j = verts.length - 1; i < verts.length; j = i++) {
    const [xi, yi] = verts[i];
    const [xj, yj] = verts[j];
    const intersect = ((yi > y) !== (yj > y)) &&
      (x < (xj - xi) * (y - yi) / (yj - yi + 1e-15) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

function curlIntegralOverPolygon(vertices) {
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const [x, y] of vertices) {
    if (x < xmin) xmin = x; if (x > xmax) xmax = x;
    if (y < ymin) ymin = y; if (y > ymax) ymax = y;
  }
  const N = 60;
  const dx = (xmax - xmin) / N;
  const dy = (ymax - ymin) / N;
  let sum = 0;
  for (let j = 0; j < N; j++) {
    const y = ymin + (j + 0.5) * dy;
    for (let i = 0; i < N; i++) {
      const x = xmin + (i + 0.5) * dx;
      if (pointInPolygon(x, y, vertices)) sum += curlAt(x, y) * dx * dy;
    }
  }
  return sum;
}

function drawLoops() {
  if (state.loopInProgress) {
    const v = state.loopInProgress.vertices;
    if (v.length >= 1) {
      ctx.strokeStyle = "#7ee787";
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 4]);
      ctx.beginPath();
      for (let i = 0; i < v.length; i++) {
        const [sx, sy] = worldToScreen(v[i][0], v[i][1]);
        if (i === 0) ctx.moveTo(sx, sy);
        else ctx.lineTo(sx, sy);
      }
      if (state.hover) {
        const [hx, hy] = worldToScreen(state.hover.x, state.hover.y);
        ctx.lineTo(hx, hy);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = "#7ee787";
      for (const [x, y] of v) {
        const [sx, sy] = worldToScreen(x, y);
        ctx.beginPath(); ctx.arc(sx, sy, 4, 0, 2 * Math.PI); ctx.fill();
      }
      if (v.length >= 3 && state.hover) {
        const [hx, hy] = worldToScreen(state.hover.x, state.hover.y);
        ctx.fillStyle = "rgba(126,231,135,0.9)";
        ctx.font = "11px 'JetBrains Mono', monospace";
        ctx.fillText("dbl-click to close", hx + 10, hy - 8);
      }
    }
  }

  for (const loop of state.loops) {
    const v = loop.vertices;
    const circ = circulationOfLoop(v);
    const dbl = curlIntegralOverPolygon(v);
    loop.circulation = circ;
    loop.doubleIntegral = dbl;

    ctx.beginPath();
    for (let i = 0; i < v.length; i++) {
      const [sx, sy] = worldToScreen(v[i][0], v[i][1]);
      if (i === 0) ctx.moveTo(sx, sy);
      else ctx.lineTo(sx, sy);
    }
    ctx.closePath();
    const intensity = Math.min(0.4, Math.abs(circ) * 0.05);
    ctx.fillStyle = circ >= 0 ? `rgba(126,231,135,${0.1 + intensity})`
                              : `rgba(255,123,114,${0.1 + intensity})`;
    ctx.fill();
    ctx.strokeStyle = "#7ee787";
    ctx.lineWidth = 2;
    ctx.stroke();

    ctx.fillStyle = "#7ee787";
    for (let i = 0; i < v.length; i++) {
      const [x1, y1] = v[i];
      const [x2, y2] = v[(i + 1) % v.length];
      const [sx1, sy1] = worldToScreen(x1, y1);
      const [sx2, sy2] = worldToScreen(x2, y2);
      const mx = (sx1 + sx2) / 2, my = (sy1 + sy2) / 2;
      const angle = Math.atan2(sy2 - sy1, sx2 - sx1);
      ctx.save();
      ctx.translate(mx, my);
      ctx.rotate(angle);
      ctx.beginPath();
      ctx.moveTo(6, 0);
      ctx.lineTo(-3, 4);
      ctx.lineTo(-3, -4);
      ctx.closePath();
      ctx.fill();
      ctx.restore();
    }

    let cx = 0, cy = 0;
    for (const [x, y] of v) { cx += x; cy += y; }
    cx /= v.length; cy /= v.length;
    const [scx, scy] = worldToScreen(cx, cy);
    const lines = [
      `∮ F·dr  = ${circ.toFixed(3)}`,
      `∬ curl dA = ${dbl.toFixed(3)}`,
      `Δ = ${Math.abs(circ - dbl).toExponential(1)}`,
    ];
    ctx.font = "11px 'JetBrains Mono', monospace";
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
    const boxX = scx - maxW / 2 - 6;
    const boxY = scy - 24;
    ctx.fillStyle = "rgba(15,20,25,0.92)";
    ctx.fillRect(boxX, boxY, maxW + 12, lines.length * 14 + 6);
    ctx.strokeStyle = "rgba(126,231,135,0.4)";
    ctx.lineWidth = 0.5;
    ctx.strokeRect(boxX, boxY, maxW + 12, lines.length * 14 + 6);
    lines.forEach((l, i) => {
      ctx.fillStyle = i === 2 ? "rgba(230,237,243,0.6)" : "#7ee787";
      ctx.fillText(l, boxX + 6, boxY + 4 + i * 14);
    });
  }
}

function drawPotential() {
  if (!state.potentialAt) return;
  const { x, y, value } = state.potentialAt;
  const [sx, sy] = worldToScreen(x, y);

  const [ox, oy] = worldToScreen(0, 0);
  const [mx, my] = worldToScreen(x, 0);
  ctx.strokeStyle = "rgba(210,168,255,0.9)";
  ctx.lineWidth = 2;
  ctx.setLineDash([4, 3]);
  ctx.beginPath();
  ctx.moveTo(ox, oy);
  ctx.lineTo(mx, my);
  ctx.lineTo(sx, sy);
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "#d2a8ff";
  ctx.beginPath(); ctx.arc(sx, sy, 5, 0, 2 * Math.PI); ctx.fill();
  ctx.strokeStyle = "white"; ctx.lineWidth = 1; ctx.stroke();
  ctx.fillStyle = "white";
  ctx.beginPath(); ctx.arc(ox, oy, 3, 0, 2 * Math.PI); ctx.fill();

  const lines = [
    `f(${x.toFixed(2)}, ${y.toFixed(2)})`,
    `= ${value.toFixed(4)}`,
    `(line integral from origin)`,
  ];
  ctx.font = "11px 'JetBrains Mono', monospace";
  ctx.textAlign = "left"; ctx.textBaseline = "top";
  const maxW = Math.max(...lines.map(l => ctx.measureText(l).width));
  ctx.fillStyle = "rgba(15,20,25,0.92)";
  ctx.fillRect(sx + 10, sy + 10, maxW + 12, lines.length * 14 + 6);
  ctx.strokeStyle = "rgba(210,168,255,0.4)";
  ctx.strokeRect(sx + 10, sy + 10, maxW + 12, lines.length * 14 + 6);
  lines.forEach((l, i) => {
    ctx.fillStyle = i === 2 ? "rgba(230,237,243,0.55)" : "#d2a8ff";
    ctx.fillText(l, sx + 16, sy + 14 + i * 14);
  });
}

// ------------------------------------------------------------
// Main render
// ------------------------------------------------------------
function render() {
  ctx.fillStyle = "#0f1419";
  ctx.fillRect(0, 0, cssWidth, cssHeight);

  drawOverlay();

  if (state.showGrid) drawGrid();

  const refMag = computeReferenceMag();

  if (state.animating) {
    drawParticles(refMag);
  } else if (state.showFlow) {
    drawStreamlines(refMag);
  } else {
    drawVectorField(refMag);
  }

  drawCurves();
  drawDivRings();
  drawPaddles();
  drawLoops();
  drawPotential();
  drawProbes();
}

// ------------------------------------------------------------
// Animation loop
// ------------------------------------------------------------
function animationFrame(now) {
  const dt = Math.min(0.05, (now - state.lastFrame) / 1000 || 0.016);
  state.lastFrame = now;

  if (state.animating) {
    state.time += dt * state.speed;
    ensureParticles();
    const refMag = computeReferenceMag();
    updateParticles(dt, refMag);
  }

  // Update paddle wheel state from local field
  // Always free to scale by user's animation speed slider so paddles feel synced with particles.
  const paddleSpeed = state.speed;
  const paddleRefMag = computeReferenceMag();
  const margin = (state.xmax - state.xmin) * 0.05;
  for (const p of state.paddles) {
    // Spin from local curl at the paddle's CURRENT location
    const c = curlAt(p.x, p.y);
    p.omega = isFinite(c) ? c / 2 : 0;
    p.angle += p.omega * dt * paddleSpeed;

    // Translate along the flow when not pinned and not being hovered (hover gives the
    // user a chance to aim at the pin without chasing a moving target). The paddle
    // continues spinning in either case so the curl readout stays live.
    if (!p.fixed && !p._hovered) {
      const [u, v] = evalField(p.x, p.y, state.time);
      if (isFinite(u) && isFinite(v)) {
        const stepDt = dt * paddleSpeed;
        const [k1u, k1v] = [u, v];
        const [k2u, k2v] = evalField(p.x + 0.5 * stepDt * k1u, p.y + 0.5 * stepDt * k1v, state.time);
        const [k3u, k3v] = evalField(p.x + 0.5 * stepDt * k2u, p.y + 0.5 * stepDt * k2v, state.time);
        const [k4u, k4v] = evalField(p.x + stepDt * k3u, p.y + stepDt * k3v, state.time);
        p.x += stepDt * (k1u + 2*k2u + 2*k3u + k4u) / 6;
        p.y += stepDt * (k1v + 2*k2v + 2*k3v + k4v) / 6;
      }
      // Accumulate trail
      if (!p._trail) p._trail = [];
      p._trail.push([p.x, p.y]);
      if (p._trail.length > 60) p._trail.shift();
      // If paddle drifts off screen, respawn at its original click location
      if (!isFinite(p.x) || !isFinite(p.y) ||
          p.x < state.xmin - margin || p.x > state.xmax + margin ||
          p.y < state.ymin - margin || p.y > state.ymax + margin) {
        p.x = p.x0; p.y = p.y0;
        p._trail = [];
      }
    } else if (p.fixed) {
      // Only clear trail when explicitly pinned, not on transient hover
      p._trail = null;
    }
  }

  render();
  requestAnimationFrame(animationFrame);
}

// ------------------------------------------------------------
// Resize handling
// ------------------------------------------------------------
// We work in CSS pixels everywhere. Canvas backing store is dpr-scaled
// via setTransform so HiDPI screens render crisply without any other math
// changing. cssWidth/cssHeight are the dimensions our coord transforms use.
let cssWidth = 0, cssHeight = 0;

function resize() {
  const rect = canvas.parentElement.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  cssWidth = rect.width;
  cssHeight = rect.height;
  canvas.width = Math.floor(cssWidth * dpr);
  canvas.height = Math.floor(cssHeight * dpr);
  canvas.style.width = cssWidth + "px";
  canvas.style.height = cssHeight + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

// ------------------------------------------------------------
// Mouse / interaction
// ------------------------------------------------------------
// Live mouse position in CSS pixels relative to canvas — used by the renderer
// for hover-state visuals on paddle wheels and other interactive elements.
let cssMouseX = null, cssMouseY = null;

function mousePos(evt) {
  const rect = canvas.getBoundingClientRect();
  return [
    evt.clientX - rect.left,
    evt.clientY - rect.top,
  ];
}

canvas.addEventListener("mousemove", (e) => {
  const [sx, sy] = mousePos(e);
  cssMouseX = sx; cssMouseY = sy;
  const [x, y] = screenToWorld(sx, sy);
  state.hover = { x, y };
  cursorInfo.textContent = `x: ${x.toFixed(2)}   y: ${y.toFixed(2)}`;

  // Show pointer cursor when over a paddle pin
  let overPin = false;
  for (const p of state.paddles) {
    if (!p._pinHitbox) continue;
    const [hx, hy, hw, hh] = p._pinHitbox;
    if (sx >= hx && sx <= hx + hw && sy >= hy && sy <= hy + hh) { overPin = true; break; }
  }
  if (!state.isPanning) {
    canvas.style.cursor = overPin ? "pointer" : "crosshair";
  }

  // Update right-side readout based on tool
  showHoverReadout(x, y);

  // Pan
  if (state.isPanning) {
    const dx = (e.clientX - state.panStart.cx) * (state.xmax - state.xmin) / canvas.getBoundingClientRect().width;
    const dy = (e.clientY - state.panStart.cy) * (state.ymax - state.ymin) / canvas.getBoundingClientRect().height;
    state.xmin = state.panStart.xmin - dx;
    state.xmax = state.panStart.xmax - dx;
    state.ymin = state.panStart.ymin + dy;
    state.ymax = state.panStart.ymax + dy;
    updateRangeInputs();
  }
});

canvas.addEventListener("mouseleave", () => {
  state.hover = null;
  cssMouseX = null; cssMouseY = null;
  readoutEl.style.display = "none";
});

canvas.addEventListener("mousedown", (e) => {
  if (e.button === 1 || e.shiftKey) {
    state.isPanning = true;
    state.panStart = {
      cx: e.clientX, cy: e.clientY,
      xmin: state.xmin, xmax: state.xmax,
      ymin: state.ymin, ymax: state.ymax,
    };
    canvas.style.cursor = "grabbing";
    e.preventDefault();
  }
});

window.addEventListener("mouseup", () => {
  if (state.isPanning) {
    state.isPanning = false;
    canvas.style.cursor = "crosshair";
  }
});

// Drag-to-pan with primary button if not over an interaction
let mouseDownPos = null;
canvas.addEventListener("mousedown", (e) => {
  if (e.button === 0 && !e.shiftKey) {
    mouseDownPos = { cx: e.clientX, cy: e.clientY, t: Date.now() };
  }
});
canvas.addEventListener("mouseup", (e) => {
  if (e.button !== 0 || !mouseDownPos) return;
  const dx = e.clientX - mouseDownPos.cx;
  const dy = e.clientY - mouseDownPos.cy;
  const dt = Date.now() - mouseDownPos.t;
  mouseDownPos = null;
  // Treat as click only if mouse barely moved
  if (Math.hypot(dx, dy) < 4 && dt < 400) {
    handleClick(e);
  }
});

// Drag-with-left-button pans when distance threshold passed
canvas.addEventListener("mousemove", (e) => {
  if (mouseDownPos && !state.isPanning) {
    const dx = e.clientX - mouseDownPos.cx;
    const dy = e.clientY - mouseDownPos.cy;
    if (Math.hypot(dx, dy) > 4) {
      state.isPanning = true;
      state.panStart = {
        cx: mouseDownPos.cx, cy: mouseDownPos.cy,
        xmin: state.xmin, xmax: state.xmax,
        ymin: state.ymin, ymax: state.ymax,
      };
      canvas.style.cursor = "grabbing";
    }
  }
});

canvas.addEventListener("wheel", (e) => {
  e.preventDefault();
  const [sx, sy] = mousePos(e);
  const [wx, wy] = screenToWorld(sx, sy);
  const factor = e.deltaY > 0 ? 1.15 : 1 / 1.15;
  state.xmin = wx + (state.xmin - wx) * factor;
  state.xmax = wx + (state.xmax - wx) * factor;
  state.ymin = wy + (state.ymin - wy) * factor;
  state.ymax = wy + (state.ymax - wy) * factor;
  updateRangeInputs();
}, { passive: false });

function handleClick(e) {
  const [sx, sy] = mousePos(e);
  const [x, y] = screenToWorld(sx, sy);
  const refMag = computeReferenceMag();

  // Check if click hits an existing paddle's pin button — toggles fixed state
  for (const p of state.paddles) {
    if (!p._pinHitbox) continue;
    const [hx, hy, hw, hh] = p._pinHitbox;
    if (sx >= hx && sx <= hx + hw && sy >= hy && sy <= hy + hh) {
      p.fixed = !p.fixed;
      // When unpinning, reset origin to current location so respawn lands here
      if (!p.fixed) { p.x0 = p.x; p.y0 = p.y; p._trail = []; }
      return;
    }
  }

  if (state.tool === "trace") {
    const fwd = integrateCurve(x, y, +1, 800, refMag);
    const bwd = integrateCurve(x, y, -1, 800, refMag);
    state.curves.push({ start: [x, y], points: fwd, color: "#ffa657" });
    state.curves.push({ start: [x, y], points: bwd, color: "#d2a8ff" });
  } else if (state.tool === "paddle") {
    state.paddles.push({ x, y, x0: x, y0: y, angle: 0, omega: 0, fixed: false });
  } else if (state.tool === "divergence") {
    const r = (state.xmax - state.xmin) * 0.04;
    state.divRings.push({ x, y, r, flux: 0, div: 0 });
  } else if (state.tool === "probe") {
    state.probes.push({ x, y });
  } else if (state.tool === "circulation") {
    if (!state.loopInProgress) {
      state.loopInProgress = { vertices: [[x, y]] };
    } else {
      // If clicking near the first vertex (and we have at least 3 points), close
      const v0 = state.loopInProgress.vertices[0];
      const span = Math.max(state.xmax - state.xmin, state.ymax - state.ymin);
      const closeThresh = span * 0.02;
      if (state.loopInProgress.vertices.length >= 3 &&
          Math.hypot(x - v0[0], y - v0[1]) < closeThresh) {
        closeLoopInProgress();
      } else {
        state.loopInProgress.vertices.push([x, y]);
      }
    }
  } else if (state.tool === "potential") {
    const value = potentialAt(x, y);
    state.potentialAt = { x, y, value };
  }
}

function closeLoopInProgress() {
  if (state.loopInProgress && state.loopInProgress.vertices.length >= 3) {
    // Strip a duplicate-trailing vertex if dblclick added one nearly identical to the previous
    const v = state.loopInProgress.vertices;
    const span = Math.max(state.xmax - state.xmin, state.ymax - state.ymin);
    const eps = span * 0.01;
    while (v.length >= 4) {
      const last = v[v.length - 1];
      const prev = v[v.length - 2];
      if (Math.hypot(last[0] - prev[0], last[1] - prev[1]) < eps) v.pop();
      else break;
    }
    if (v.length >= 3) {
      state.loops.push({ vertices: v.slice(), circulation: 0, doubleIntegral: 0 });
    }
    state.loopInProgress = null;
  }
}

canvas.addEventListener("dblclick", (e) => {
  if (state.tool === "circulation") {
    e.preventDefault();
    closeLoopInProgress();
  }
});

window.addEventListener("keydown", (e) => {
  if (e.key === "Enter" && state.tool === "circulation" && state.loopInProgress) {
    closeLoopInProgress();
  }
  if (e.key === "Escape") {
    state.loopInProgress = null;
  }
});

function showHoverReadout(x, y) {
  const [u, v] = evalField(x, y, state.time);
  const mag = Math.hypot(u, v);
  const curl = curlAt(x, y);
  const div = divAt(x, y);
  const sign = (n) => n >= 0 ? "pos" : "neg";
  const fmt = (n) => (isFinite(n) ? n.toFixed(3) : "—");
  const html = `
    <div class="title">at (${x.toFixed(2)}, ${y.toFixed(2)})</div>
    <div class="row"><span class="label">P (=u)</span><span class="val">${fmt(u)}</span></div>
    <div class="row"><span class="label">Q (=v)</span><span class="val">${fmt(v)}</span></div>
    <div class="row"><span class="label">|F|</span><span class="val">${fmt(mag)}</span></div>
    <div class="row"><span class="label">curl·k̂</span><span class="val ${sign(curl)}">${fmt(curl)}</span></div>
    <div class="row"><span class="label">div F</span><span class="val ${sign(div)}">${fmt(div)}</span></div>
  `;
  readoutEl.innerHTML = html;
  readoutEl.style.display = "block";
}

// ------------------------------------------------------------
// UI wiring
// ------------------------------------------------------------
function bindToggle(id, key) {
  const el = document.getElementById(id);
  el.addEventListener("click", () => {
    state[key] = !state[key];
    el.classList.toggle("on", state[key]);
  });
}

function updateRangeInputs() {
  document.getElementById("xmin").value = state.xmin.toFixed(2);
  document.getElementById("xmax").value = state.xmax.toFixed(2);
  document.getElementById("ymin").value = state.ymin.toFixed(2);
  document.getElementById("ymax").value = state.ymax.toFixed(2);
}

function setupUI() {
  // Field expression inputs
  document.getElementById("fp").addEventListener("input", (e) => {
    state.Pexpr = e.target.value;
    if (!compileField()) {
      document.getElementById("fp").classList.add("error");
    }
  });
  document.getElementById("fq").addEventListener("input", (e) => {
    state.Qexpr = e.target.value;
    if (!compileField()) {
      document.getElementById("fq").classList.add("error");
    }
  });

  // Range inputs
  ["xmin", "xmax", "ymin", "ymax"].forEach(id => {
    document.getElementById(id).addEventListener("change", (e) => {
      const v = parseFloat(e.target.value);
      if (isFinite(v)) state[id] = v;
    });
  });

  document.getElementById("reset-view").addEventListener("click", () => {
    state.xmin = -5; state.xmax = 5; state.ymin = -5; state.ymax = 5;
    updateRangeInputs();
  });

  document.getElementById("square-view").addEventListener("click", () => {
    const cx = (state.xmin + state.xmax) / 2;
    const cy = (state.ymin + state.ymax) / 2;
    const xspan = state.xmax - state.xmin;
    const aspect = cssWidth / cssHeight;
    const yspan = xspan / aspect;
    state.ymin = cy - yspan / 2;
    state.ymax = cy + yspan / 2;
    updateRangeInputs();
  });

  // Sliders
  const slider = (id, key, fmt) => {
    const el = document.getElementById(id);
    const out = document.getElementById(id + "-val");
    el.addEventListener("input", (e) => {
      state[key] = parseFloat(e.target.value);
      out.textContent = fmt ? fmt(state[key]) : state[key];
    });
  };
  slider("density", "density");
  slider("length", "length", v => v.toFixed(2));
  slider("particles", "numParticles");
  slider("speed", "speed", v => v.toFixed(2));

  // Toggles
  bindToggle("toggle-flow", "showFlow");
  bindToggle("toggle-color", "colorByMag");
  bindToggle("toggle-grid", "showGrid");
  bindToggle("toggle-normalize", "normalize");

  // Animation
  document.getElementById("anim-toggle").addEventListener("click", (e) => {
    state.animating = !state.animating;
    e.target.textContent = state.animating ? "⏸ Pause" : "▶ Animate";
    if (state.animating) ensureParticles();
  });
  document.getElementById("anim-clear").addEventListener("click", () => {
    state.particles = [];
  });

  // Tool selection (sidebar)
  document.querySelectorAll(".tool-btn-side").forEach(btn => {
    btn.addEventListener("click", () => {
      state.tool = btn.dataset.tool;
      document.querySelectorAll(".tool-btn-side").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".tool-btn[data-toolbar-tool]").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      document.querySelector(`.tool-btn[data-toolbar-tool="${state.tool}"]`).classList.add("active");
      updateToolHint();
    });
  });

  // Toolbar tool selection (top of canvas)
  document.querySelectorAll(".tool-btn[data-toolbar-tool]").forEach(btn => {
    btn.addEventListener("click", () => {
      state.tool = btn.dataset.toolbarTool;
      document.querySelectorAll(".tool-btn[data-toolbar-tool]").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".tool-btn-side").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      document.querySelector(`.tool-btn-side[data-tool="${state.tool}"]`).classList.add("active");
      updateToolHint();
    });
  });

  document.getElementById("clear-tools").addEventListener("click", () => {
    state.curves = [];
    state.paddles = [];
    state.divRings = [];
    state.probes = [];
    state.loops = [];
    state.loopInProgress = null;
    state.potentialAt = null;
  });

  // Overlay dropdown
  document.getElementById("overlay-select").addEventListener("change", (e) => {
    state.overlay = e.target.value;
  });

  // Presets
  const presets = [
    { label: "Rotation", P: "-y", Q: "x" },
    { label: "Source", P: "x", Q: "y" },
    { label: "Sink", P: "-x", Q: "-y" },
    { label: "Saddle", P: "x", Q: "-y" },
    { label: "Shear", P: "y", Q: "0" },
    { label: "Conservative", P: "2x", Q: "2y" },
    { label: "Vortex (1/r)", P: "-y/(x^2+y^2)", Q: "x/(x^2+y^2)" },
    { label: "Wave", P: "sin(y)", Q: "sin(x)" },
    { label: "Spiral in", P: "-x-y", Q: "x-y" },
    { label: "Spiral out", P: "x-y", Q: "x+y" },
    { label: "Time wave", P: "sin(x+t)", Q: "cos(y+t)" },
    { label: "Dipole", P: "(x^2-y^2)/(x^2+y^2)^2", Q: "2*x*y/(x^2+y^2)^2" },
  ];
  const presetGrid = document.getElementById("presets");
  presets.forEach(p => {
    const btn = document.createElement("button");
    btn.className = "btn";
    btn.textContent = p.label;
    btn.addEventListener("click", () => {
      document.getElementById("fp").value = p.P;
      document.getElementById("fq").value = p.Q;
      state.Pexpr = p.P;
      state.Qexpr = p.Q;
      compileField();
    });
    presetGrid.appendChild(btn);
  });
}

function updateToolHint() {
  const hints = {
    trace: "Click in field: trace forward (orange) and backward (purple) integral curves.",
    paddle: "Click to drop a paddle wheel. It rides the flow, spinning at the local curl/2. Click the pin icon to fix it in place.",
    divergence: "Click to place a small circle. Shows ∇·F via outward arrows and flux integral.",
    probe: "Click to place a probe showing F(x,y) and its components at that point.",
    circulation: "Click vertices of a closed loop. Double-click (or press Enter) to close. Shows ∮F·dr vs ∬curl dA — Green's theorem.",
    potential: "Click anywhere to compute f(x,y) = ∫₀ᴾ F·dr along an L-path from origin. (Only meaningful when curl ≈ 0.)",
  };
  document.getElementById("tool-hint").textContent = hints[state.tool];
}

// ------------------------------------------------------------
// Init
// ------------------------------------------------------------
function init() {
  setupUI();
  compileField();
  resize();
  window.addEventListener("resize", resize);
  requestAnimationFrame((t) => {
    state.lastFrame = t;
    animationFrame(t);
  });
}

init();
