// GWAS browser: stacked Manhattan and line tracks on one genome axis.
// Data layout is described in README.md ("Prepared data contract").
import * as core from './core.js';

const params = new URLSearchParams(location.search);
const DATA = (params.get('data') || '/output_data/gwas_browser/').replace(/\/?$/, '/');

const FULL_SPAN = 100e6;       // below this span, fetch full-resolution SNPs
const GENE_MODEL_SPAN = 10e6;  // below this span, draw gene models instead of density
const GENE_LABEL_SPAN = 2e6;   // below this span, label genes
const MIN_SPAN = 200;
const AXIS = 46;               // css px reserved for the y axis on every canvas
const PAD_R = 12;
const HIT_R = 6;
const LRU_BYTES = 600e6;
const GENE_ROW_H = 100;

const $ = sel => document.querySelector(sel);
function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text != null) e.textContent = text;
  return e;
}

const state = {
  manifest: null, scale: 1000, binSize: 10000,
  contigMap: new Map(), dataContigs: new Set(), layout: null,
  gs: 0, ge: 1, hidden: new Set(), yMode: 'track', yMax: null,
  showUps: false, showUnfiltered: false, height: 90,
  sel: null, band: null, drag: null, pendingLoc: null, pendingSel: null,
};
const tracks = [];
let overview, ruler, geneRow;
let geom = null;
let palette = {};
const objByCanvas = new WeakMap();

let geneEntries = [];
const geneIndexByContig = new Map();
const geneCache = new Map();
const genePending = new Map();
const fullCache = new Map();
const fullPending = new Map();
const fullFailed = new Set();
let fullBytes = 0;

// ---------- utilities ----------

// fetch with one retry on network errors (not on HTTP errors)
async function fetchRetry(url) {
  try {
    return await fetch(url);
  } catch (err) {
    await new Promise(r => setTimeout(r, 300));
    return fetch(url);
  }
}

async function fetchBuf(path) {
  const r = await fetchRetry(DATA + path);
  if (!r.ok) throw new Error(`${r.status} ${DATA}${path}`);
  return r.arrayBuffer();
}

async function fetchJSON(path, missingOk = false) {
  const r = await fetchRetry(DATA + path);
  if (!r.ok) {
    if (missingOk && r.status === 404) return null;
    throw new Error(`${r.status} ${DATA}${path}`);
  }
  return r.json();
}

const probe = document.createElement('canvas').getContext('2d', { willReadFrequently: true });
function rgba32(color, alpha = 1) {
  probe.clearRect(0, 0, 1, 1);
  probe.fillStyle = color;
  probe.fillRect(0, 0, 1, 1);
  const [r, g, b, a] = probe.getImageData(0, 0, 1, 1).data;
  return ((Math.round(a * alpha) << 24) | (b << 16) | (g << 8) | r) >>> 0;
}

function readPalette() {
  const cs = getComputedStyle(document.documentElement);
  const v = n => cs.getPropertyValue(n).trim();
  palette = {
    fg: v('--fg'), muted: v('--muted'), line: v('--line'), grid: v('--grid'), band: v('--band'),
    panel: v('--panel'), accent: v('--accent'), hit: v('--hit'), sig: v('--sig'), sug: v('--sug'),
    dotACss: v('--dot-a'), dotBCss: v('--dot-b'),
    dotA: rgba32(v('--dot-a')), dotB: rgba32(v('--dot-b')),
    massA: rgba32(v('--dot-a'), 0.3), massB: rgba32(v('--dot-b'), 0.3),
    peak: v('--peak'), peakFaint: v('--peak-faint'), peakEdge: v('--peak-edge'),
    gene: v('--gene'), geneNc: v('--gene-nc'), lineMean: v('--line-mean'), lineEnv: v('--line-env'),
    ovA: v('--ov-a'), ovB: v('--ov-b'),
  };
}

function fmtNum(v) {
  if (v == null || v !== v) return '—';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a < 1e-3 || a >= 1e5) return v.toExponential(2);
  return String(+v.toPrecision(4));
}

function fmtP(logp) {
  const p = Math.pow(10, -logp);
  return p < 1e-3 ? p.toExponential(2) : String(+p.toPrecision(3));
}

function fmtTick(v, step) {
  if (step >= 1e6) return (v / 1e6).toFixed(0) + ' Mb';
  if (step >= 1e5) return (v / 1e6).toFixed(1) + ' Mb';
  if (step >= 1e4) return (v / 1e6).toFixed(2) + ' Mb';
  if (step >= 1e3) return (v / 1e3).toFixed(0) + ' kb';
  if (step >= 100) return (v / 1e3).toFixed(1) + ' kb';
  return core.formatInt(v);
}

function fmtAxis(v, step) {
  if (v === 0) return '0';
  const a = Math.abs(step);
  if (a >= 1) return String(Math.round(v));
  if (a >= 1e-4) return v.toFixed(Math.max(0, -Math.floor(Math.log10(a))));
  return v.toExponential(0);
}

// Track names are long underscore-joined identifiers: allow line breaks after each underscore.
const displayName = m => (m.label || m.id).replace(/_/g, '_\u200b');

function setStatus(text) { $('#status').textContent = text; }

function updateStatus() {
  const nG = tracks.filter(t => t.kind === 'gwas').length, nL = tracks.length - nG;
  let s = `${nG} GWAS · ${nL} line · ${core.formatInt(geneEntries.length)} genes`;
  if (fullPending.size) s += ` · loading ${fullPending.size} full-res`;
  setStatus(s);
}

function fatal(html) {
  const m = $('#message');
  m.innerHTML = html;
  m.hidden = false;
  setStatus('');
}

// ---------- layout and view ----------

function buildLayout() {
  const isMain = n => /^chr/i.test(n);
  const hasMain = state.manifest.contigs.some(c => isMain(c.name) && state.dataContigs.has(c.name));
  const keep = state.manifest.contigs.filter(c => {
    if (hasMain && isMain(c.name)) return true;
    if (!hasMain) return state.dataContigs.has(c.name);
    return state.showUps && state.dataContigs.has(c.name);
  });
  state.layout = core.buildLayout(keep);
}

function setView(gs, ge) {
  const total = state.layout.total;
  if (!(ge > gs)) ge = gs + MIN_SPAN;
  if (ge - gs < MIN_SPAN) {
    const c = (gs + ge) / 2;
    gs = c - MIN_SPAN / 2; ge = c + MIN_SPAN / 2;
  }
  if (ge - gs >= total) { gs = 0; ge = total; } else {
    if (gs < 0) { ge -= gs; gs = 0; }
    if (ge > total) { gs -= ge - total; ge = total; }
  }
  state.gs = gs; state.ge = ge;
  requestRedraw();
}

function zoomBy(f, anchor) {
  const a = anchor ?? (state.gs + state.ge) / 2;
  setView(a - (a - state.gs) * f, a + (state.ge - a) * f);
}

function panBy(frac) {
  const d = (state.ge - state.gs) * frac;
  setView(state.gs + d, state.ge + d);
}

// Make sure a contig is on the axis (turning on ups scaffolds if needed).
function ensureContig(contig) {
  if (state.layout.byName.has(contig)) return true;
  if (!state.dataContigs.has(contig) && !state.contigMap.has(contig)) return false;
  setUps(true);
  return state.layout.byName.has(contig);
}

function goLocus(contig, start, end) {
  if (!ensureContig(contig)) { setStatus(`${contig} has no data`); return false; }
  setView(core.toGenome(state.layout, contig, start), core.toGenome(state.layout, contig, end) + 1);
  return true;
}

function goGene(e) {
  const pad = Math.max(20000, Math.round((e.end - e.start) * 0.3));
  const len = state.contigMap.get(e.contig)?.length ?? e.end + pad;
  goLocus(e.contig, Math.max(1, e.start - pad), Math.min(len, e.end + pad));
}

function goText(text) {
  const loc = core.parseLocus(text, state.contigMap);
  if (loc) return goLocus(loc.contig, loc.start, loc.end);
  const hit = core.searchGenes(geneEntries, text, 1)[0];
  if (hit) { goGene(hit); return true; }
  setStatus(`Not found: ${text}`);
  return false;
}

function setUps(on) {
  if (state.showUps === on) return;
  const keep = core.encodeView(state.layout, state.gs, state.ge);
  state.showUps = on;
  $('#ups').checked = on;
  buildLayout();
  const v = core.decodeView(state.layout, keep);
  if (v) { state.gs = v[0]; state.ge = v[1]; } else { state.gs = 0; state.ge = state.layout.total; }
  requestRedraw();
}

// ---------- URL hash ----------

function readHash() {
  const h = new URLSearchParams(location.hash.slice(1));
  state.pendingLoc = h.get('loc');
  state.pendingSel = h.get('sel');
  state.hidden = new Set((h.get('hide') || '').split(',').filter(Boolean));
  if (h.get('y')) state.yMode = h.get('y');
  state.yMax = h.get('ymax') ? Number(h.get('ymax')) || null : null;
  state.showUps = h.get('ups') === '1' || (state.pendingLoc != null && !/^chr/i.test(state.pendingLoc));
  state.showUnfiltered = h.get('unf') === '1';
  if (h.get('h')) state.height = Number(h.get('h')) || 90;
}

let hashTimer = 0;
function scheduleHash() {
  clearTimeout(hashTimer);
  hashTimer = setTimeout(writeHash, 250);
}

function writeHash() {
  const parts = [['loc', core.encodeView(state.layout, state.gs, state.ge)]];
  if (state.sel) parts.push(['sel', `${state.sel.contig}:${state.sel.pos}`]);
  if (state.hidden.size) parts.push(['hide', [...state.hidden].join(',')]);
  if (state.yMode !== 'track') parts.push(['y', state.yMode]);
  if (state.yMax) parts.push(['ymax', state.yMax]);
  if (state.showUps) parts.push(['ups', '1']);
  if (state.showUnfiltered) parts.push(['unf', '1']);
  if (state.height !== 90) parts.push(['h', state.height]);
  const enc = v => encodeURIComponent(v).replace(/%3A/g, ':').replace(/%2C/g, ',');
  history.replaceState(null, '', '#' + parts.map(([k, v]) => `${k}=${enc(v)}`).join('&'));
}

// ---------- data loading ----------

async function loadTrack(t) {
  const m = t.meta;
  if (t.kind === 'gwas') {
    const [ov, bins] = await Promise.all([fetchBuf(`gwas/${m.id}/overview.bin`), fetchBuf(`gwas/${m.id}/bins.bin`)]);
    const n = ov.byteLength / 6;
    t.ovPos = new Uint32Array(ov, 0, n);
    t.ovVal = new Uint16Array(ov, 4 * n, n);
    t.bins = new Uint16Array(bins);
  } else {
    const buf = await fetchBuf(`line/${m.id}.bin`);
    const n = m.n;
    t.starts = new Uint32Array(buf, 0, n);
    t.ends = new Uint32Array(buf, 4 * n, n);
    t.vals = new Float32Array(buf, 8 * n, n);
    let maxLen = 0;
    for (let i = 0; i < n; i++) maxLen = Math.max(maxLen, t.ends[i] - t.starts[i]);
    t.maxLen = maxLen;
  }
  t.loaded = true;
}

async function loadGeneIndex() {
  const g = state.manifest.genes;
  if (!g) return;
  const r = await fetchRetry(DATA + g.index);
  if (!r.ok) throw new Error(`${r.status} ${DATA}${g.index}`);
  geneEntries = core.parseGeneIndex(await r.text());
  for (const e of geneEntries) {
    if (!geneIndexByContig.has(e.contig)) geneIndexByContig.set(e.contig, []);
    geneIndexByContig.get(e.contig).push(e);
  }
  for (const [c, list] of geneIndexByContig) {
    list.sort((a, b) => a.start - b.start);
    geneIndexByContig.set(c, { list, starts: Uint32Array.from(list, e => e.start) });
  }
}

function loadFull(t, contig) {
  const key = `${t.meta.id}/${contig}`;
  const hit = fullCache.get(key);
  if (hit) { fullCache.delete(key); fullCache.set(key, hit); return Promise.resolve(hit); }
  if (fullPending.has(key)) return fullPending.get(key);
  if (!t.meta.contigs[contig] || fullFailed.has(key)) return Promise.resolve(null);
  const p = fetchBuf(`gwas/${t.meta.id}/full/${contig}.bin`).then(buf => {
    const n = buf.byteLength / 6;
    const v = { pos: new Uint32Array(buf, 0, n), val: new Uint16Array(buf, 4 * n, n), bytes: buf.byteLength };
    fullCache.set(key, v);
    fullBytes += v.bytes;
    for (const [k, old] of fullCache) {
      if (fullBytes <= LRU_BYTES || fullCache.size <= 1) break;
      fullCache.delete(k);
      fullBytes -= old.bytes;
    }
    return v;
  }).catch(err => {
    console.error(err);
    fullFailed.add(key);
    return null;
  }).finally(() => {
    fullPending.delete(key);
    updateStatus();
    requestRedraw();
  });
  fullPending.set(key, p);
  updateStatus();
  return p;
}

// Full-resolution data if cached; otherwise start loading and return null.
function getFull(t, contig) {
  const key = `${t.meta.id}/${contig}`;
  const hit = fullCache.get(key);
  if (hit) { fullCache.delete(key); fullCache.set(key, hit); return hit; }
  loadFull(t, contig);
  return null;
}

function loadDetail(t, contig) {
  if (!t.detail) t.detail = new Map();
  if (!t.detail.has(contig)) {
    const p = t.meta.contigs[contig]
      ? fetchJSON(`gwas/${t.meta.id}/detail/${contig}.json`, true).catch(err => { console.error(err); return null; })
      : Promise.resolve(null);
    t.detail.set(contig, p);
  }
  return t.detail.get(contig);
}

function loadGenes(contig) {
  if (geneCache.has(contig)) return Promise.resolve(geneCache.get(contig));
  if (genePending.has(contig)) return genePending.get(contig);
  const g = state.manifest.genes;
  if (!g || !g.contigs.includes(contig)) { geneCache.set(contig, []); return Promise.resolve([]); }
  const p = fetchJSON(`genes/${contig}.json`, true)
    .then(list => list || [])
    .catch(err => { console.error(err); return []; })
    .then(list => { geneCache.set(contig, list); genePending.delete(contig); requestRedraw(); return list; });
  genePending.set(contig, p);
  return p;
}

function getGenes(contig) {
  if (geneCache.has(contig)) return geneCache.get(contig);
  loadGenes(contig);
  return null;
}

// SNP at an exact position in a track (highest logp if several), from full data.
function snpAt(full, pos) {
  if (!full) return null;
  let i = core.lowerBound(full.pos, pos), best = -1;
  for (; i < full.pos.length && full.pos[i] === pos; i++) {
    if (best < 0 || full.val[i] > full.val[best]) best = i;
  }
  return best < 0 ? null : { pos, logp: full.val[best] / state.scale };
}

function detailAt(detail, pos) {
  if (!detail || !detail.pos) return null;
  const i = core.lowerBound(detail.pos, pos);
  if (detail.pos[i] !== pos) return null;
  const out = {};
  for (const k of Object.keys(detail)) out[k] = detail[k][i];
  return out;
}

function lineValueAt(t, contig, pos) {
  const c = t.meta.contigs[contig];
  if (!c) return null;
  const lo = core.lowerBound(t.starts, pos - t.maxLen, c.off, c.off + c.n);
  const hi = core.upperBound(t.starts, pos, c.off, c.off + c.n);
  if (t.meta.kind === 'point') {
    // nearest point
    let best = -1, bd = Infinity;
    const a = Math.max(c.off, hi - 1), b = Math.min(c.off + c.n - 1, hi);
    for (let i = a; i <= b; i++) {
      const d = Math.abs(t.starts[i] - pos);
      if (d < bd) { bd = d; best = i; }
    }
    return best < 0 ? null : { value: t.vals[best], start: t.starts[best], end: t.ends[best], dist: bd };
  }
  let sum = 0, n = 0, s = 0, e = 0;
  for (let i = lo; i < hi; i++) {
    if (t.ends[i] >= pos && t.vals[i] === t.vals[i]) {
      sum += t.vals[i]; n++;
      if (n === 1) { s = t.starts[i]; e = t.ends[i]; } else { s = Math.min(s, t.starts[i]); e = Math.max(e, t.ends[i]); }
    }
  }
  return n ? { value: sum / n, start: s, end: e, n } : null;
}

// ---------- DOM ----------

function makeRow(parent, cls, cssH, labelNodes) {
  const row = el('div', 'row ' + cls);
  const label = el('div', 'label');
  label.append(...labelNodes);
  const canvas = el('canvas');
  row.append(label, canvas);
  parent.append(row);
  const obj = { row, label, canvas, ctx: canvas.getContext('2d'), cssH, w: 0, h: 0, dpr: 1 };
  objByCanvas.set(canvas, obj);
  return obj;
}

function sizeCanvas(o) {
  const dpr = window.devicePixelRatio || 1;
  const w = o.canvas.clientWidth, h = o.cssH;
  if (o.w !== w || o.h !== h || o.dpr !== dpr) {
    o.canvas.width = Math.max(1, Math.round(w * dpr));
    o.canvas.height = Math.max(1, Math.round(h * dpr));
    o.canvas.style.height = h + 'px';
    o.w = w; o.h = h; o.dpr = dpr;
  }
  o.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function gwasSub(m) {
  const parts = [m.test || ''];
  if (m.sig != null) parts.push(`sig ${m.sig.toFixed(2)} · sug ${m.sug.toFixed(2)}`);
  else parts.push('no thresholds');
  if (m.peaks && m.peaks.length) parts.push(`${m.peaks.length} peak${m.peaks.length > 1 ? 's' : ''}`);
  return parts.filter(Boolean).join(' · ');
}

function gwasTitle(m) {
  return [
    m.source,
    `${core.formatInt(m.n_snps)} SNPs (${m.p_col})` + (m.min_af ? `, af ≥ ${m.min_af}` : ''),
    m.thresh_source ? `thresholds: ${m.thresh_source}` : 'no permutation thresholds',
    `max −log10 p ${m.max_logp}`,
  ].join('\n');
}

function buildRows() {
  const fixed = $('#fixed-rows');
  overview = makeRow(fixed, 'overview', 26, [el('span', null, 'Genome')]);
  ruler = makeRow(fixed, 'ruler', 26, [el('span', null, '')]);
  ruler.label.firstChild.id = 'span-label';
  overview.canvas.classList.add('nav');
  ruler.canvas.classList.add('plot');

  const main = $('#tracks');
  const io = new IntersectionObserver(entries => {
    for (const en of entries) {
      const t = objByCanvas.get(en.target.querySelector('canvas'));
      if (!t) continue;
      t.onScreen = en.isIntersecting;
      if (t.onScreen && t.dirty) requestRedraw();
    }
  });
  for (const m of state.manifest.tracks) {
    const kind = m.type === 'gwas' ? 'gwas' : 'line';
    const name = el('div', 'name', displayName(m));
    const sub = el('div', 'sub', kind === 'gwas' ? gwasSub(m)
      : `${m.value_name}${m.preset ? ' · ' + m.preset : ''}`);
    const hide = el('button', 'hide', '×');
    hide.title = 'Hide track';
    const t = makeRow(main, kind, kind === 'gwas' ? state.height : Math.round(state.height * 0.75), [name, sub, hide]);
    Object.assign(t, { meta: m, kind, onScreen: true, dirty: true });
    t.label.title = kind === 'gwas' ? gwasTitle(m) : (m.source || '');
    t.canvas.classList.add('plot');
    hide.addEventListener('click', () => setHidden(m.id, true));
    tracks.push(t);
    io.observe(t.row);
  }
  geneRow = makeRow($('#gene-row'), 'genes', GENE_ROW_H, [el('div', 'name', 'Genes'), el('div', 'sub', '')]);
  geneRow.canvas.classList.add('plot');
  geneRow.row.style.display = state.manifest.genes ? '' : 'none';

  const list = $('#track-list');
  for (const t of tracks) {
    const lab = el('label');
    const cb = el('input');
    cb.type = 'checkbox';
    cb.dataset.id = t.meta.id;
    lab.append(cb, document.createTextNode(t.meta.label || t.meta.id));
    cb.addEventListener('change', () => setHidden(t.meta.id, !cb.checked));
    list.append(lab);
  }
  applyHidden();
}

function applyHidden() {
  for (const t of tracks) t.row.style.display = state.hidden.has(t.meta.id) ? 'none' : '';
  for (const cb of document.querySelectorAll('#track-list input')) cb.checked = !state.hidden.has(cb.dataset.id);
}

function setHidden(id, hidden) {
  if (hidden) state.hidden.add(id); else state.hidden.delete(id);
  applyHidden();
  requestRedraw();
}

function applyHeight() {
  for (const t of tracks) t.cssH = t.kind === 'gwas' ? state.height : Math.round(state.height * 0.75);
  requestRedraw();
}

// ---------- drawing ----------

let rafId = 0;
function requestRedraw() {
  if (!rafId) rafId = requestAnimationFrame(drawAll);
}

function drawAll() {
  rafId = 0;
  if (!state.layout) return;
  sizeCanvas(ruler);
  const w = ruler.w;
  if (w < AXIS + PAD_R + 20) return;
  const PW = w - AXIS - PAD_R;
  const span = state.ge - state.gs;
  geom = { w, PW, gs: state.gs, ge: state.ge, span, bpPerPx: span / PW };
  drawOverview();
  drawRuler();
  for (const t of tracks) {
    if (state.hidden.has(t.meta.id)) continue;
    sizeCanvas(t);
    if (!t.onScreen || !t.loaded) { t.dirty = true; continue; }
    t.dirty = false;
    if (t.kind === 'gwas') drawGwas(t); else drawLine(t);
  }
  if (state.manifest.genes) drawGenes();
  updateLocusBox();
  scheduleHash();
}

const xOf = g => AXIS + (g - geom.gs) / geom.bpPerPx;
const gAtX = x => geom.gs + (Math.min(Math.max(x, AXIS), AXIS + geom.PW) - AXIS) * geom.bpPerPx;

function clipPlot(ctx, h) {
  ctx.save();
  ctx.beginPath();
  ctx.rect(AXIS, 0, geom.PW, h);
  ctx.clip();
}

function drawBackdrop(o, top, bot) {
  const { ctx } = o;
  const vcs = core.visibleContigs(state.layout, geom.gs, geom.ge);
  if (vcs.length < 2) return vcs;
  for (const vc of vcs) {
    const x0 = xOf(vc.entry.offset), x1 = xOf(vc.entry.offset + vc.entry.length);
    if (vc.entry.index % 2) {
      ctx.fillStyle = palette.band;
      ctx.fillRect(x0, top, x1 - x0, bot - top);
    }
  }
  return vcs;
}

function drawSelectionMarks(o, top, bot) {
  const { ctx } = o;
  if (state.band) {
    const x0 = xOf(Math.min(...state.band)), x1 = xOf(Math.max(...state.band));
    ctx.fillStyle = palette.accent;
    ctx.globalAlpha = 0.12;
    ctx.fillRect(x0, top, x1 - x0, bot - top);
    ctx.globalAlpha = 1;
  }
  if (state.sel && state.layout.byName.has(state.sel.contig)) {
    const x = xOf(core.toGenome(state.layout, state.sel.contig, state.sel.pos) + 0.5);
    ctx.strokeStyle = palette.hit;
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(Math.round(x) + 0.5, top);
    ctx.lineTo(Math.round(x) + 0.5, bot);
    ctx.stroke();
    ctx.setLineDash([]);
  }
}

function yMaxFor(t) {
  if (state.yMax) return state.yMax;
  let base;
  if (state.yMode === 'shared') {
    base = 0;
    for (const o of tracks) if (o.kind === 'gwas' && !state.hidden.has(o.meta.id)) base = Math.max(base, o.meta.max_logp);
  } else if (state.yMode === 'view') {
    let mx = 0;
    const bs = state.binSize;
    for (const vc of core.visibleContigs(state.layout, geom.gs, geom.ge)) {
      const c = t.meta.contigs[vc.entry.name];
      if (!c) continue;
      const k1 = Math.min(c.bin_n - 1, Math.floor((vc.posHi - 1) / bs));
      for (let k = Math.floor((vc.posLo - 1) / bs); k <= k1; k++) mx = Math.max(mx, t.bins[c.bin_off + k]);
    }
    base = mx / state.scale;
  } else {
    base = t.meta.max_logp;
  }
  return Math.max(base, t.meta.sig || 0, 2) * 1.08;
}

function getImage(t, W, H) {
  if (!t.off || t.off.width !== W || t.off.height !== H) {
    t.off = document.createElement('canvas');
    t.off.width = W; t.off.height = H;
    t.offCtx = t.off.getContext('2d');
    t.img = t.offCtx.createImageData(W, H);
    t.buf32 = new Uint32Array(t.img.data.buffer);
  }
  t.buf32.fill(0);
  return t.buf32;
}

// Square dots of side `d` device px for sorted points [lo, hi). Hot loop: whole-chromosome views
// push millions of points per frame, so it avoids divisions and Math calls and skips points whose
// centre pixel is already painted (most of them in dense columns).
function rasterPoints(buf, W, H, pos, val, lo, hi, off, gs, bpd, yBase, k, vmax, color, d) {
  const r = d >> 1, inv = 1 / bpd, base = off - 1 - gs, yb = yBase + 0.5;
  for (let i = lo; i < hi; i++) {
    const xf = (base + pos[i]) * inv;
    if (xf < 0 || xf >= W) continue;
    const x = xf | 0;
    const v = val[i] > vmax ? vmax : val[i];
    const y = (yb - v * k) | 0;
    if (y < 0 || y >= H || buf[y * W + x] === color) continue;
    const x0 = Math.max(0, x - r), x1 = Math.min(W - 1, x - r + d - 1);
    const y0 = Math.max(0, y - r), y1 = Math.min(H - 1, y - r + d - 1);
    for (let yy = y0; yy <= y1; yy++) {
      const row = yy * W;
      for (let xx = x0; xx <= x1; xx++) buf[row + xx] = color;
    }
  }
}

function drawGwas(t) {
  const { ctx, w, h, dpr } = t;
  const m = t.meta, S = state.scale;
  ctx.clearRect(0, 0, w, h);
  const top = 9, bot = h - 4, ph = bot - top;
  const ymax = yMaxFor(t);
  t.ymax = ymax; t.top = top; t.bot = bot;
  const yOf = v => bot - Math.min(v, ymax) / ymax * ph;

  clipPlot(ctx, h);
  const vcs = drawBackdrop(t, top, bot);
  drawPeaks(t, top, bot);

  // raster layer (device pixels, plot area only)
  const W = Math.max(1, Math.round(geom.PW * dpr)), H = Math.max(1, Math.round(h * dpr));
  const buf = getImage(t, W, H);
  const bpd = geom.bpPerPx / dpr;
  const vmax = Math.round(ymax * S);
  const k = ph * dpr / vmax;
  const yBase = bot * dpr;
  const dot = Math.max(2, Math.round(1.6 * dpr));
  const full = geom.span <= FULL_SPAN;
  const cutoff = Math.round(m.overview_cutoff * S);
  const cols = new Uint16Array(W);
  for (const vc of vcs) {
    const name = vc.entry.name, c = m.contigs[name];
    if (!c) continue;
    const odd = vc.entry.index % 2 === 1;
    const color = odd ? palette.dotB : palette.dotA;
    const off = vc.entry.offset;
    const data = full ? getFull(t, name) : null;
    if (data) {
      const lo = core.lowerBound(data.pos, vc.posLo), hi = core.upperBound(data.pos, vc.posHi);
      rasterPoints(buf, W, H, data.pos, data.val, lo, hi, off, geom.gs, bpd, yBase, k, vmax, color, dot);
    } else {
      cols.fill(0);
      core.binColumnMax(t.bins, c.bin_off, c.bin_n, state.binSize, off, geom.gs, bpd, W, cols);
      const mass = odd ? palette.massB : palette.massA;
      for (let x = 0; x < W; x++) {
        if (!cols[x]) continue;
        const v = Math.min(cols[x], cutoff, vmax);
        const y0 = Math.max(0, Math.round(yBase - v * k)), y1 = Math.min(H, Math.round(yBase));
        for (let y = y0; y < y1; y++) buf[y * W + x] = mass;
      }
      const lo = core.lowerBound(t.ovPos, vc.posLo, c.ov_off, c.ov_off + c.ov_n);
      const hi = core.upperBound(t.ovPos, vc.posHi, c.ov_off, c.ov_off + c.ov_n);
      rasterPoints(buf, W, H, t.ovPos, t.ovVal, lo, hi, off, geom.gs, bpd, yBase, k, vmax, color, dot);
    }
  }
  t.offCtx.putImageData(t.img, 0, 0);
  ctx.drawImage(t.off, AXIS, 0, geom.PW, h);

  // thresholds
  ctx.lineWidth = 1;
  ctx.setLineDash([5, 4]);
  for (const [v, col] of [[m.sug, palette.sug], [m.sig, palette.sig]]) {
    if (v == null || v > ymax) continue;
    const y = Math.round(yOf(v)) + 0.5;
    ctx.strokeStyle = col;
    ctx.beginPath(); ctx.moveTo(AXIS, y); ctx.lineTo(AXIS + geom.PW, y); ctx.stroke();
  }
  ctx.setLineDash([]);

  // SNPs above the suggestive line as circles
  if (m.sug != null) {
    const sugS = Math.round(m.sug * S), sigS = m.sig != null ? Math.round(m.sig * S) : Infinity;
    for (const vc of vcs) {
      const c = m.contigs[vc.entry.name];
      if (!c) continue;
      const lo = core.lowerBound(t.ovPos, vc.posLo, c.ov_off, c.ov_off + c.ov_n);
      const hi = core.upperBound(t.ovPos, vc.posHi, c.ov_off, c.ov_off + c.ov_n);
      const base = vc.entry.index % 2 ? palette.dotBCss : palette.dotACss;
      for (let i = lo; i < hi; i++) {
        const v = t.ovVal[i];
        if (v < sugS) continue;
        ctx.beginPath();
        ctx.arc(xOf(vc.entry.offset + t.ovPos[i] - 0.5), yOf(v / S), 2.6, 0, 2 * Math.PI);
        ctx.fillStyle = v >= sigS ? palette.sig : base;
        ctx.fill();
      }
    }
  }

  drawSelectionMarks(t, top, bot);
  if (state.sel && state.layout.byName.has(state.sel.contig)) {
    const s = snpAt(fullCache.get(`${m.id}/${state.sel.contig}`), state.sel.pos);
    if (s) {
      ctx.strokeStyle = palette.hit;
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(xOf(core.toGenome(state.layout, state.sel.contig, state.sel.pos) + 0.5), yOf(s.logp), 5, 0, 2 * Math.PI);
      ctx.stroke();
    }
  }
  ctx.restore();
  drawYAxis(t, 0, ymax, yOf, top, bot);
}

function drawPeaks(t, top, bot) {
  const { ctx } = t;
  const L = state.layout;
  const draw = (list, fill, label) => {
    for (const p of list) {
      if (!L.byName.has(p.contig)) continue;
      let x0 = xOf(core.toGenome(L, p.contig, p.start)), x1 = xOf(core.toGenome(L, p.contig, p.end) + 1);
      if (x1 < AXIS || x0 > AXIS + geom.PW) continue;
      if (x1 - x0 < 2) { const c = (x0 + x1) / 2; x0 = c - 1; x1 = c + 1; }
      ctx.fillStyle = fill;
      ctx.fillRect(x0, top - 9, x1 - x0, bot - top + 9);
      if (!label) continue;
      if (p.lead) {
        const xl = xOf(core.toGenome(L, p.contig, p.lead) + 0.5);
        ctx.fillStyle = palette.peakEdge;
        ctx.beginPath(); ctx.moveTo(xl - 3, 0); ctx.lineTo(xl + 3, 0); ctx.lineTo(xl, 4); ctx.fill();
      }
      if (x1 - x0 > 14 || geom.span < 5e7) {
        ctx.fillStyle = palette.peakEdge;
        ctx.font = '600 9.5px system-ui, sans-serif';
        ctx.textBaseline = 'top';
        ctx.fillText('P' + p.id, Math.max(AXIS + 1, x0) + 2, 0);
      }
    }
  };
  if (state.showUnfiltered) draw(t.meta.peaks_unfiltered || [], palette.peakFaint, false);
  draw(t.meta.peaks || [], palette.peak, true);
}

function drawYAxis(o, lo, hi, yOf, top, bot, fmt) {
  const { ctx, h } = o;
  ctx.fillStyle = palette.panel;
  ctx.fillRect(0, 0, AXIS, h);
  ctx.strokeStyle = palette.line;
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(AXIS - 0.5, top); ctx.lineTo(AXIS - 0.5, bot); ctx.stroke();
  const ticks = core.niceTicks(lo, hi, Math.max(1, Math.floor((bot - top) / 26)));
  const step = ticks.length > 1 ? ticks[1] - ticks[0] : 1;
  ctx.fillStyle = palette.muted;
  ctx.font = '10px system-ui, sans-serif';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  for (const v of ticks) {
    const y = Math.round(yOf(v)) + 0.5;
    if (y < top - 1 || y > bot + 1) continue;
    ctx.beginPath(); ctx.moveTo(AXIS - 4, y); ctx.lineTo(AXIS, y); ctx.stroke();
    ctx.fillText((fmt || fmtAxis)(v, step), AXIS - 6, Math.min(Math.max(y, 6), o.h - 6));
  }
  ctx.textAlign = 'left';
}

function drawLine(t) {
  const { ctx, w, h } = t;
  const m = t.meta;
  ctx.clearRect(0, 0, w, h);
  const vcs = core.visibleContigs(state.layout, geom.gs, geom.ge);
  const Wc = Math.max(1, Math.floor(geom.PW));
  const point = m.kind === 'point';
  // aggregate first: tracks without a fixed range scale to the data in view
  const aggs = [];
  let vlo = Infinity, vhi = -Infinity;
  for (const vc of vcs) {
    const c = m.contigs[vc.entry.name];
    if (!c) continue;
    const i0 = core.lowerBound(t.starts, vc.posLo - t.maxLen, c.off, c.off + c.n);
    const i1 = core.upperBound(t.starts, vc.posHi, c.off, c.off + c.n);
    const a = core.windowColumns(t.starts, t.ends, t.vals, i0, i1, vc.entry.offset, geom.gs, geom.bpPerPx, Wc);
    for (let x = 0; x < Wc; x++) if (a.n[x]) { vlo = Math.min(vlo, a.min[x]); vhi = Math.max(vhi, a.max[x]); }
    aggs.push({ vc, a, i0, i1 });
  }
  const fixed = m.ymin != null || m.ymax != null;
  let lo, hi;
  if (fixed) { lo = m.ymin ?? m.vmin; hi = m.ymax ?? m.vmax; } else if (vhi >= vlo) { lo = vlo; hi = vhi; } else { lo = m.vmin; hi = m.vmax; }
  if (!(hi > lo)) { lo = (lo || 0) - 0.5; hi = lo + 1; }
  if (!fixed) { const pad = (hi - lo) * 0.06; if (lo !== 0) lo -= pad; hi += pad; }
  const top = 6, bot = h - 5, ph = bot - top;
  const yOf = v => bot - (Math.min(Math.max(v, lo), hi) - lo) / (hi - lo) * ph;
  t.yOf = yOf;

  clipPlot(ctx, h);
  drawBackdrop(t, top, bot);
  if (m.zero_line && lo < 0 && hi > 0) {
    ctx.strokeStyle = palette.line;
    ctx.lineWidth = 1;
    const y = Math.round(yOf(0)) + 0.5;
    ctx.beginPath(); ctx.moveTo(AXIS, y); ctx.lineTo(AXIS + geom.PW, y); ctx.stroke();
  }
  for (const { vc, a, i0, i1 } of aggs) {
    ctx.fillStyle = palette.lineEnv;
    for (let x = 0; x < Wc; x++) {
      if (a.n[x] < 2) continue;
      const y0 = yOf(a.max[x]), y1 = yOf(a.min[x]);
      if (y1 - y0 >= 1) ctx.fillRect(AXIS + x, y0, 1, y1 - y0);
    }
    ctx.strokeStyle = palette.lineMean;
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    let open = false;
    for (let x = 0; x < Wc; x++) {
      if (!a.n[x]) { if (!point) open = false; continue; }
      const y = yOf(a.sum[x] / a.n[x]);
      if (open) ctx.lineTo(AXIS + x + 0.5, y); else ctx.moveTo(AXIS + x + 0.5, y);
      open = true;
    }
    ctx.stroke();
    if (point && i1 - i0 < geom.PW / 4) {
      ctx.fillStyle = palette.lineMean;
      for (let i = i0; i < i1; i++) {
        if (t.vals[i] !== t.vals[i]) continue;
        ctx.beginPath();
        ctx.arc(xOf(vc.entry.offset + t.starts[i] - 0.5), yOf(t.vals[i]), 2, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
  }
  drawSelectionMarks(t, top, bot);
  ctx.restore();
  drawYAxis(t, lo, hi, yOf, top, bot);
}

function drawOverview() {
  sizeCanvas(overview);
  const { ctx, w, h } = overview;
  ctx.clearRect(0, 0, w, h);
  const L = state.layout, PW = geom.PW;
  const xg = g => AXIS + g / L.total * PW;
  ctx.font = '10px system-ui, sans-serif';
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'center';
  for (const c of L.list) {
    const x0 = xg(c.offset), x1 = xg(c.offset + c.length);
    ctx.fillStyle = c.index % 2 ? palette.ovB : palette.ovA;
    ctx.fillRect(x0, 4, Math.max(0.5, x1 - x0 - 0.5), h - 8);
    const name = c.name.replace(/^chr/i, '');
    if (x1 - x0 > ctx.measureText(name).width + 4) {
      ctx.fillStyle = palette.muted;
      ctx.fillText(name, (x0 + x1) / 2, h / 2);
    }
  }
  // filtered peaks of visible tracks
  ctx.fillStyle = palette.sig;
  for (const t of tracks) {
    if (t.kind !== 'gwas' || state.hidden.has(t.meta.id)) continue;
    for (const p of t.meta.peaks || []) {
      const g = core.toGenome(L, p.contig, p.lead || p.start);
      if (g !== null) ctx.fillRect(Math.round(xg(g)) - 0.5, h - 6, 1.5, 3);
    }
  }
  ctx.textAlign = 'left';
  let x0 = xg(geom.gs), x1 = xg(geom.ge);
  if (x1 - x0 < 3) { const c = (x0 + x1) / 2; x0 = c - 1.5; x1 = c + 1.5; }
  ctx.fillStyle = palette.accent;
  ctx.globalAlpha = 0.15;
  ctx.fillRect(x0, 1, x1 - x0, h - 2);
  ctx.globalAlpha = 1;
  ctx.strokeStyle = palette.accent;
  ctx.lineWidth = 1.5;
  ctx.strokeRect(x0, 1.5, x1 - x0, h - 3);
  if (overview.band) {
    const [a, b] = overview.band;
    ctx.fillStyle = palette.accent;
    ctx.globalAlpha = 0.25;
    ctx.fillRect(Math.min(a, b), 1, Math.abs(b - a), h - 2);
    ctx.globalAlpha = 1;
  }
}

function drawRuler() {
  const { ctx, w, h } = ruler;
  ctx.clearRect(0, 0, w, h);
  clipPlot(ctx, h);
  const vcs = core.visibleContigs(state.layout, geom.gs, geom.ge);
  const multi = vcs.length > 1;
  ctx.font = '11px system-ui, sans-serif';
  ctx.textBaseline = 'top';
  ctx.strokeStyle = palette.line;
  ctx.lineWidth = 1;
  for (const vc of vcs) {
    const off = vc.entry.offset;
    const x0 = Math.max(AXIS, xOf(off + vc.posLo - 1)), x1 = Math.min(AXIS + geom.PW, xOf(off + vc.posHi));
    const segW = x1 - x0;
    const cs = xOf(off);
    if (cs >= AXIS) {
      ctx.beginPath(); ctx.moveTo(Math.round(cs) + 0.5, 0); ctx.lineTo(Math.round(cs) + 0.5, h); ctx.stroke();
    }
    const name = vc.entry.name;
    if (multi && (geom.span > 5e7 || segW < 160)) {
      const label = ctx.measureText(name).width + 6 < segW ? name : name.replace(/^chr/i, '');
      if (ctx.measureText(label).width + 4 < segW) {
        ctx.fillStyle = palette.fg;
        ctx.textAlign = 'center';
        ctx.fillText(label, (x0 + x1) / 2, 7);
        ctx.textAlign = 'left';
      }
      continue;
    }
    const ticks = core.niceTicks(vc.posLo, vc.posHi, Math.max(2, Math.floor(segW / 110)));
    const step = ticks.length > 1 ? ticks[1] - ticks[0] : 1;
    ctx.fillStyle = palette.muted;
    for (const v of ticks) {
      const x = Math.round(xOf(off + v - 1)) + 0.5;
      if (x < x0 || x > x1) continue;
      ctx.beginPath(); ctx.moveTo(x, h - 6); ctx.lineTo(x, h); ctx.stroke();
      ctx.fillText(fmtTick(v, step), x + 3, 3);
    }
    if (multi) {
      ctx.fillStyle = palette.fg;
      ctx.fillText(name, x0 + 3, 14);
    }
  }
  drawSelectionMarks(ruler, 0, h);
  ctx.restore();
  ctx.strokeStyle = palette.line;
  ctx.beginPath(); ctx.moveTo(0, h - 0.5); ctx.lineTo(w, h - 0.5); ctx.stroke();
  $('#span-label').textContent = core.formatBp(geom.span);
}

const labelWidths = new Map();
function labelWidth(ctx, s) {
  let v = labelWidths.get(s);
  if (v === undefined) { v = ctx.measureText(s).width; labelWidths.set(s, v); }
  return v;
}

function drawGenes() {
  const o = geneRow;
  sizeCanvas(o);
  const { ctx, w, h } = o;
  ctx.clearRect(0, 0, w, h);
  o.hits = [];
  clipPlot(ctx, h);
  const vcs = drawBackdrop(o, 0, h);
  const models = geom.span <= GENE_MODEL_SPAN;
  o.label.querySelector('.sub').textContent = models ? (geom.span <= GENE_LABEL_SPAN ? 'gene models' : 'gene models (zoom in for labels)') : 'gene density (zoom < 10 Mb for models)';
  if (!models) {
    const Wc = Math.max(1, Math.floor(geom.PW));
    const counts = new Uint32Array(Wc);
    for (const vc of vcs) {
      const idx = geneIndexByContig.get(vc.entry.name);
      if (!idx) continue;
      const i0 = core.lowerBound(idx.starts, vc.posLo), i1 = core.upperBound(idx.starts, vc.posHi);
      for (let i = i0; i < i1; i++) {
        const x = Math.floor((vc.entry.offset + idx.starts[i] - 1 - geom.gs) / geom.bpPerPx);
        if (x >= 0 && x < Wc) counts[x]++;
      }
    }
    // smooth to a few pixels so density reads at genome scale
    const r = 2, sm = new Float32Array(Wc);
    let mx = 0;
    for (let x = 0; x < Wc; x++) {
      let s = 0;
      for (let d = -r; d <= r; d++) s += counts[Math.min(Wc - 1, Math.max(0, x + d))];
      sm[x] = s; mx = Math.max(mx, s);
    }
    ctx.fillStyle = palette.geneNc;
    const ph = h - 14;
    for (let x = 0; x < Wc; x++) {
      if (!sm[x]) continue;
      const bh = sm[x] / mx * ph;
      ctx.fillRect(AXIS + x, h - 4 - bh, 1, bh);
    }
  } else {
    const labels = geom.span <= GENE_LABEL_SPAN;
    const rowH = labels ? 25 : 13;
    const maxRows = Math.max(1, Math.floor((h - 6) / rowH));
    ctx.font = '10.5px system-ui, sans-serif';
    const items = [];
    for (const vc of vcs) {
      const genes = getGenes(vc.entry.name);
      if (!genes) continue;
      const off = vc.entry.offset;
      for (const g of genes) {
        if (g.e < vc.posLo || g.s > vc.posHi) continue;
        const x0 = xOf(off + g.s - 1), x1 = Math.max(xOf(off + g.e), x0 + 1);
        const lw = labels ? labelWidth(ctx, g.n || g.id) : 0;
        items.push({ g, off, contig: vc.entry.name, x0, x1, ext: [Math.max(x0, AXIS), Math.max(x1, Math.max(x0, AXIS) + lw + 2)] });
      }
    }
    items.sort((a, b) => a.ext[0] - b.ext[0]);
    const rows = core.packRows(items.map(i => i.ext), maxRows, 6);
    items.forEach((it, k) => drawGeneModel(ctx, it, 4 + rows[k] * rowH, labels, o));
  }
  drawSelectionMarks(o, 0, h);
  ctx.restore();
  ctx.fillStyle = palette.panel;
  ctx.fillRect(0, 0, AXIS, h);
}

function drawGeneModel(ctx, it, y, labels, o) {
  const { g, off, x0, x1 } = it;
  const coding = g.b === 'protein_coding';
  const color = coding ? palette.gene : palette.geneNc;
  const yc = y + 5;
  ctx.fillStyle = color;
  ctx.strokeStyle = color;
  ctx.fillRect(x0, yc - 0.5, x1 - x0, 1);
  const X = p => xOf(off + p - 1);
  if (g.x && g.x.length) {
    for (let i = 0; i < g.x.length; i += 2) {
      const a = g.x[i], b = g.x[i + 1];
      const ex0 = X(a), ex1 = Math.max(X(b + 1), ex0 + 1);
      ctx.fillRect(ex0, yc - 3, ex1 - ex0, 6);
      if (g.cs != null) {
        const ca = Math.max(a, g.cs), cb = Math.min(b, g.ce);
        if (cb >= ca) {
          const c0 = X(ca), c1 = Math.max(X(cb + 1), c0 + 1);
          ctx.fillRect(c0, yc - 5, c1 - c0, 10);
        }
      }
    }
  } else {
    ctx.fillRect(x0, yc - 4, x1 - x0, 8);
  }
  // strand chevrons along introns
  if (x1 - x0 > 30 && (g.st === '+' || g.st === '-')) {
    const dir = g.st === '+' ? 1 : -1;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let x = Math.max(x0, AXIS) + 12; x < Math.min(x1, AXIS + geom.PW) - 6; x += 24) {
      ctx.moveTo(x - 2 * dir, yc - 3); ctx.lineTo(x + 1 * dir, yc); ctx.lineTo(x - 2 * dir, yc + 3);
    }
    ctx.stroke();
  }
  if (labels) {
    ctx.fillStyle = palette.fg;
    ctx.textBaseline = 'top';
    ctx.fillText(g.n || g.id, Math.max(x0, AXIS + 1), yc + 7);
  }
  o.hits.push({ x0: Math.min(x0, x1 - 3), x1: Math.max(x1, x0 + 3), y0: y - 1, y1: y + (labels ? 23 : 11), g, contig: it.contig });
}

function updateLocusBox() {
  const box = $('#locus');
  if (document.activeElement === box) return;
  const a = core.fromGenome(state.layout, geom.gs), b = core.fromGenome(state.layout, Math.max(geom.gs, geom.ge - 1));
  if (!a) return;
  if (geom.gs <= 0 && geom.ge >= state.layout.total) box.value = 'All';
  else if (a.contig === b.contig) {
    box.value = (a.pos <= 1 && b.pos >= a.entry.length) ? a.contig : `${a.contig}:${core.formatInt(a.pos)}-${core.formatInt(b.pos)}`;
  } else box.value = `${a.contig}:${core.formatInt(a.pos)} – ${b.contig}:${core.formatInt(b.pos)}`;
}

// ---------- hit testing ----------

function hitGwas(t, x, y) {
  if (!geom || t.ymax == null) return null;
  const S = state.scale;
  const gC = gAtX(x), dg = HIT_R * geom.bpPerPx;
  const yOf = v => t.bot - Math.min(v, t.ymax) / t.ymax * (t.bot - t.top);
  let best = null, bd = HIT_R * HIT_R;
  for (const vc of core.visibleContigs(state.layout, gC - dg, gC + dg + 1)) {
    const name = vc.entry.name, c = t.meta.contigs[name];
    if (!c) continue;
    const full = geom.span <= FULL_SPAN ? fullCache.get(`${t.meta.id}/${name}`) : null;
    const pos = full ? full.pos : t.ovPos, val = full ? full.val : t.ovVal;
    const a = full ? 0 : c.ov_off, b = full ? full.pos.length : c.ov_off + c.ov_n;
    const lo = core.lowerBound(pos, vc.posLo, a, b), hi = core.upperBound(pos, vc.posHi, a, b);
    for (let i = lo; i < hi; i++) {
      const dx = xOf(vc.entry.offset + pos[i] - 0.5) - x, dy = yOf(val[i] / S) - y;
      const d = dx * dx + dy * dy;
      if (d <= bd) { bd = d; best = { contig: name, pos: pos[i], logp: val[i] / S }; }
    }
  }
  return best;
}

function hitGene(x, y) {
  for (const hbox of geneRow.hits || []) {
    if (x >= hbox.x0 - 2 && x <= hbox.x1 + 2 && y >= hbox.y0 && y <= hbox.y1) return hbox;
  }
  return null;
}

function tooltipFor(o, x, y) {
  if (x < AXIS) return null;
  if (o.kind === 'gwas') {
    const h = hitGwas(o, x, y);
    if (!h) return null;
    let s = `${o.meta.label || o.meta.id}\n${h.contig}:${core.formatInt(h.pos)} · −log10 p ${h.logp.toFixed(2)}`;
    const det = o.detail?.get(h.contig);
    if (det && det.resolved) {
      const d = detailAt(det.resolved, h.pos);
      if (d) s += `\n${d.rs ?? ''}  p = ${d.p}`;
    }
    return s;
  }
  if (o.kind === 'line') {
    const p = core.fromGenome(state.layout, gAtX(x));
    const v = lineValueAt(o, p.contig, p.pos);
    if (!v) return null;
    return `${o.meta.label || o.meta.id}: ${fmtNum(v.value)}\n${p.contig}:${core.formatInt(v.start)}${v.end !== v.start ? '-' + core.formatInt(v.end) : ''}`;
  }
  if (o === geneRow) {
    const hbox = hitGene(x, y);
    if (!hbox) return null;
    const g = hbox.g;
    return `${g.n || g.id}${g.d ? ' — ' + g.d : ''}\n${hbox.contig}:${core.formatInt(g.s)}-${core.formatInt(g.e)} (${g.st})`;
  }
  return null;
}

// ---------- selection drawer ----------

function openDrawer(nodes) {
  const d = $('#drawer');
  d.replaceChildren();
  const close = el('button', 'close', '×');
  close.title = 'Close (Esc)';
  close.addEventListener('click', clearSelection);
  d.append(close, ...nodes);
  d.hidden = false;
}

function clearSelection() {
  state.sel = null;
  $('#drawer').hidden = true;
  requestRedraw();
}

function geneBox(rel, contig) {
  const box = el('div', 'gene');
  if (!rel) { box.textContent = 'No annotated gene on this contig'; return box; }
  const g = rel.gene;
  const head = el('div');
  head.append(document.createTextNode(core.describeRelation(rel) + ' '));
  const go = el('button', 'link', '→ view gene');
  go.addEventListener('click', () => goGene({ contig, start: g.s, end: g.e }));
  head.append(go);
  box.append(head);
  if (g.d) box.append(el('div', 'desc', g.d));
  return box;
}

async function selectSnp(contig, pos) {
  const sel = { contig, pos };
  state.sel = sel;
  requestRedraw();
  scheduleHash();
  openDrawer([el('h2', null, `${contig}:${core.formatInt(pos)}`), el('div', 'sub', 'Loading…')]);
  const gw = tracks.filter(t => t.kind === 'gwas' && !state.hidden.has(t.meta.id));
  const [genes, fulls, dets] = await Promise.all([
    loadGenes(contig),
    Promise.all(gw.map(t => loadFull(t, contig))),
    Promise.all(gw.map(t => loadDetail(t, contig))),
  ]);
  if (state.sel !== sel) return;
  gw.forEach((t, i) => { const p = t.detail.get(contig); p.resolved = dets[i]; });
  requestRedraw();

  let rs = null, alleles = '';
  const rows = gw.map((t, i) => {
    const s = snpAt(fulls[i], pos);
    const d = detailAt(dets[i], pos);
    if (d && !rs) { rs = d.rs; if (d.allele1 != null) alleles = `${d.allele1}/${d.allele0}`; }
    return { t, s, d };
  });
  const hasField = f => rows.some(r => r.d && r.d[f] != null);
  const cols = ['af', 'beta', 'se'].filter(hasField);
  const table = el('table');
  const hr = el('tr');
  for (const c of ['Track', '−log10 p', 'p', ...cols]) hr.append(el('th', null, c === 'beta' ? 'β' : c));
  table.append(hr);
  for (const { t, s, d } of rows) {
    const tr = el('tr');
    tr.append(el('td', 'track', displayName(t.meta)));
    const lp = el('td', s ? 'strong' : null, s ? s.logp.toFixed(2) : '—');
    if (s && t.meta.sig != null && s.logp >= t.meta.sig) lp.classList.add('sig');
    tr.append(lp);
    tr.append(el('td', null, d ? d.p : s ? '≈ ' + fmtP(s.logp) : '—'));
    for (const c of cols) tr.append(el('td', null, d && d[c] != null ? fmtNum(Number(d[c])) : '—'));
    table.append(tr);
  }
  const nodes = [el('h2', null, `${contig}:${core.formatInt(pos)}`),
    el('div', 'sub', [rs, alleles].filter(Boolean).join(' · ') || 'not in any detail file (−log10 p < 2 everywhere)'),
    geneBox(core.nearestGene(genes, pos), contig), table];
  const lt = tracks.filter(t => t.kind === 'line' && !state.hidden.has(t.meta.id));
  if (lt.length) {
    nodes.push(el('div', 'section', 'Line tracks'));
    const t2 = el('table');
    for (const t of lt) {
      const v = lineValueAt(t, contig, pos);
      const tr = el('tr');
      tr.append(el('td', 'track', displayName(t.meta)));
      let txt = '—';
      if (v) txt = fmtNum(v.value) + (t.meta.kind === 'point' ? ` (${core.formatBp(v.dist)} away)` : '');
      tr.append(el('td', null, txt));
      t2.append(tr);
    }
    nodes.push(t2);
  }
  openDrawer(nodes);
}

async function showGene(g, contig) {
  state.sel = null;
  requestRedraw();
  const title = el('h2', null, g.n || g.id);
  const sub = el('div', 'sub', `${contig}:${core.formatInt(g.s)}-${core.formatInt(g.e)} (${g.st}) · ${core.formatBp(g.e - g.s + 1)} · ${g.b || ''}`);
  const box = el('div', 'gene');
  if (g.d) box.append(el('div', null, g.d));
  box.append(el('div', 'desc', `${g.id || ''}${g.t ? ' · ' + g.t : ''}${g.x ? ` · ${g.x.length / 2} exons` : ''}`));
  const go = el('button', 'link', '→ zoom to gene');
  go.addEventListener('click', () => goGene({ contig, start: g.s, end: g.e }));
  box.append(go);
  openDrawer([title, sub, box, el('div', 'sub', 'Loading SNPs…')]);
  const gw = tracks.filter(t => t.kind === 'gwas' && !state.hidden.has(t.meta.id));
  const fulls = await Promise.all(gw.map(t => loadFull(t, contig)));
  const flank = 10000;
  const table = el('table');
  const hr = el('tr');
  for (const c of ['Track', `max −log10 p (gene ± ${core.formatBp(flank)})`, 'at']) hr.append(el('th', null, c));
  table.append(hr);
  gw.forEach((t, i) => {
    const f = fulls[i];
    let best = -1;
    if (f) {
      const lo = core.lowerBound(f.pos, g.s - flank), hi = core.upperBound(f.pos, g.e + flank);
      for (let k = lo; k < hi; k++) if (best < 0 || f.val[k] > f.val[best]) best = k;
    }
    const tr = el('tr');
    tr.append(el('td', 'track', displayName(t.meta)));
    const v = best >= 0 ? f.val[best] / state.scale : null;
    const td = el('td', v != null ? 'strong' : null, v != null ? v.toFixed(2) : '—');
    if (v != null && t.meta.sig != null && v >= t.meta.sig) td.classList.add('sig');
    tr.append(td);
    const at = el('td');
    if (best >= 0) {
      const b = el('button', 'link', core.formatInt(f.pos[best]));
      const p = f.pos[best];
      b.addEventListener('click', () => selectSnp(contig, p));
      at.append(b);
    } else at.textContent = '—';
    tr.append(at);
    table.append(tr);
  });
  openDrawer([title, sub, box, table]);
}

// ---------- interaction ----------

const tooltip = $('#tooltip');
function showTooltip(text, cx, cy) {
  if (!text) { tooltip.hidden = true; return; }
  tooltip.textContent = text;
  tooltip.hidden = false;
  const r = tooltip.getBoundingClientRect();
  let x = cx + 14, y = cy + 14;
  if (x + r.width > innerWidth - 8) x = cx - r.width - 14;
  if (y + r.height > innerHeight - 8) y = cy - r.height - 14;
  tooltip.style.left = x + 'px';
  tooltip.style.top = y + 'px';
}

function localXY(e) {
  const r = e.currentTarget.getBoundingClientRect();
  return [e.clientX - r.left, e.clientY - r.top];
}

function handleClick(o, x, y) {
  if (x < AXIS) return;
  if (o.kind === 'gwas') {
    const h = hitGwas(o, x, y);
    if (h) { selectSnp(h.contig, h.pos); return; }
  } else if (o === geneRow) {
    const hbox = hitGene(x, y);
    if (hbox) { showGene(hbox.g, hbox.contig); return; }
  }
}

let hoverRaf = 0, hoverEvt = null;
function onPlotMove(e) {
  const o = objByCanvas.get(e.currentTarget);
  const d = state.drag;
  if (d && d.canvas === e.currentTarget) {
    const dx = e.clientX - d.cx;
    if (Math.abs(dx) > 3) d.moved = true;
    if (!d.moved) return;
    showTooltip(null);
    if (d.shift) {
      const [x] = localXY(e);
      state.band = [gAtX(d.x), gAtX(x)];
      requestRedraw();
    } else {
      const shift = -dx * d.bpPerPx;
      setView(d.gs + shift, d.ge + shift);
    }
    return;
  }
  hoverEvt = { o, xy: localXY(e), cx: e.clientX, cy: e.clientY };
  if (!hoverRaf) {
    hoverRaf = requestAnimationFrame(() => {
      hoverRaf = 0;
      if (!hoverEvt) return;
      const { o, xy, cx, cy } = hoverEvt;
      showTooltip(tooltipFor(o, xy[0], xy[1]), cx, cy);
    });
  }
}

function onPlotDown(e) {
  if (e.button !== 0) return;
  const [x] = localXY(e);
  e.currentTarget.setPointerCapture(e.pointerId);
  e.currentTarget.classList.add('dragging');
  state.drag = { canvas: e.currentTarget, cx: e.clientX, x, gs: state.gs, ge: state.ge, bpPerPx: geom.bpPerPx, shift: e.shiftKey, moved: false };
}

function onPlotUp(e) {
  const d = state.drag;
  if (!d || d.canvas !== e.currentTarget) return;
  e.currentTarget.classList.remove('dragging');
  state.drag = null;
  if (!d.moved) {
    const [x, y] = localXY(e);
    handleClick(objByCanvas.get(e.currentTarget), x, y);
  } else if (d.shift && state.band) {
    const [a, b] = state.band;
    state.band = null;
    setView(Math.min(a, b), Math.max(a, b));
  }
}

function onWheel(e) {
  e.preventDefault();
  if (!geom) return;
  let dx = e.deltaX, dy = e.deltaY;
  if (e.deltaMode === 1) { dx *= 16; dy *= 16; }
  if (Math.abs(dx) > Math.abs(dy)) {
    const d = dx * geom.bpPerPx;
    setView(state.gs + d, state.ge + d);
    return;
  }
  const [x] = localXY(e);
  zoomBy(Math.exp(dy * (e.ctrlKey ? 0.01 : 0.0015)), gAtX(x));
}

function wirePlot(canvas) {
  canvas.addEventListener('pointerdown', onPlotDown);
  canvas.addEventListener('pointermove', onPlotMove);
  canvas.addEventListener('pointerup', onPlotUp);
  canvas.addEventListener('pointercancel', onPlotUp);
  canvas.addEventListener('pointerleave', () => { hoverEvt = null; showTooltip(null); });
  canvas.addEventListener('wheel', onWheel, { passive: false });
}

function wireOverview() {
  const c = overview.canvas;
  const gOf = x => Math.min(Math.max(x - AXIS, 0), geom.PW) / geom.PW * state.layout.total;
  let start = null;
  c.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    c.setPointerCapture(e.pointerId);
    start = localXY(e)[0];
  });
  c.addEventListener('pointermove', e => {
    const [x] = localXY(e);
    if (start == null) {
      const p = core.fromGenome(state.layout, gOf(x));
      showTooltip(p ? `${p.contig} · click to view, drag to select` : null, e.clientX, e.clientY);
      return;
    }
    if (Math.abs(x - start) > 3) { overview.band = [start, x]; requestRedraw(); }
  });
  c.addEventListener('pointerleave', () => showTooltip(null));
  c.addEventListener('pointerup', e => {
    if (start == null) return;
    const [x] = localXY(e);
    if (overview.band) {
      const a = gOf(Math.min(start, x)), b = gOf(Math.max(start, x));
      overview.band = null;
      setView(a, b);
    } else {
      const p = core.fromGenome(state.layout, gOf(x));
      if (p) setView(p.entry.offset, p.entry.offset + p.entry.length);
    }
    start = null;
  });
}

function wireControls() {
  const box = $('#locus'), sug = $('#suggest');
  let active = -1, results = [];
  const render = () => {
    sug.replaceChildren();
    results.forEach((e, i) => {
      const li = el('li', i === active ? 'active' : null);
      li.append(el('b', null, e.name), el('span', null, `${e.contig}:${core.formatInt(e.start)}${e.description ? ' · ' + e.description : ''}`));
      li.addEventListener('mousedown', ev => { ev.preventDefault(); pick(i); });
      sug.append(li);
    });
    sug.hidden = !results.length;
  };
  const pick = i => {
    const e = results[i];
    results = []; render();
    box.blur();
    if (e) goGene(e);
  };
  box.addEventListener('focus', () => box.select());
  box.addEventListener('input', () => {
    const q = box.value;
    results = core.parseLocus(q, state.contigMap) ? [] : core.searchGenes(geneEntries, q, 15);
    active = -1;
    render();
  });
  box.addEventListener('keydown', e => {
    if (e.key === 'ArrowDown' && results.length) { active = (active + 1) % results.length; render(); e.preventDefault(); }
    else if (e.key === 'ArrowUp' && results.length) { active = (active - 1 + results.length) % results.length; render(); e.preventDefault(); }
    else if (e.key === 'Enter') {
      e.preventDefault();
      if (active >= 0) { pick(active); return; }
      const q = box.value.trim();
      results = []; render();
      if (/^(all|genome)$/i.test(q)) setView(0, state.layout.total);
      else if (q) goText(q);
      box.blur();
    } else if (e.key === 'Escape') { results = []; render(); box.blur(); }
  });
  box.addEventListener('blur', () => { results = []; render(); requestRedraw(); });

  $('#btn-genome').addEventListener('click', () => setView(0, state.layout.total));
  $('#btn-in').addEventListener('click', () => zoomBy(0.5));
  $('#btn-out').addEventListener('click', () => zoomBy(2));
  const ymode = $('#ymode'), ymax = $('#ymax'), height = $('#height'), ups = $('#ups'), unf = $('#unfiltered');
  ymode.value = state.yMode;
  ymax.value = state.yMax || '';
  height.value = String(state.height);
  if (height.value !== String(state.height)) { state.height = 90; height.value = '90'; }
  ups.checked = state.showUps;
  unf.checked = state.showUnfiltered;
  ymode.addEventListener('change', () => { state.yMode = ymode.value; requestRedraw(); });
  ymax.addEventListener('change', () => { state.yMax = Number(ymax.value) > 0 ? Number(ymax.value) : null; requestRedraw(); });
  height.addEventListener('change', () => { state.height = Number(height.value); applyHeight(); });
  ups.addEventListener('change', () => setUps(ups.checked));
  unf.addEventListener('change', () => { state.showUnfiltered = unf.checked; requestRedraw(); });
  $('#show-all').addEventListener('click', () => { state.hidden.clear(); applyHidden(); requestRedraw(); });
  $('#hide-all').addEventListener('click', () => { for (const t of tracks) state.hidden.add(t.meta.id); applyHidden(); requestRedraw(); });
  document.addEventListener('click', e => {
    const menu = $('#track-menu');
    if (menu.open && !menu.contains(e.target)) menu.open = false;
  });

  window.addEventListener('keydown', e => {
    if (e.target.closest && e.target.closest('input, select, textarea')) return;
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === '+' || e.key === '=') zoomBy(0.5);
    else if (e.key === '-' || e.key === '_') zoomBy(2);
    else if (e.key === 'ArrowLeft') panBy(-0.2);
    else if (e.key === 'ArrowRight') panBy(0.2);
    else if (e.key === '0' || e.key === 'Home') setView(0, state.layout.total);
    else if (e.key === 'Escape') clearSelection();
    else if (e.key === '/') { e.preventDefault(); $('#locus').focus(); }
    else return;
    e.preventDefault();
  });
  window.addEventListener('resize', requestRedraw);
  new ResizeObserver(requestRedraw).observe($('#tracks'));
  matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => { readPalette(); requestRedraw(); });
  window.addEventListener('hashchange', () => {
    const h = new URLSearchParams(location.hash.slice(1));
    const loc = h.get('loc');
    if (loc && !/^chr/i.test(loc)) setUps(true);
    const v = core.decodeView(state.layout, loc);
    if (v) setView(v[0], v[1]);
    applySelHash(h.get('sel'));
  });
}

// ---------- start ----------

async function init() {
  readPalette();
  setStatus('Loading manifest…');
  let m;
  try {
    m = await fetchJSON('manifest.json');
  } catch (err) {
    fatal(`Could not load <code>${DATA}manifest.json</code> (${err.message}).<br><br>Prepare the data first:<br>
      <code>python3 code/00_common/gwas_browser/prepare_tracks.py code/00_common/gwas_browser/tracks.json</code><br><br>
      and serve the repository root with <code>bash code/00_common/gwas_browser/serve.sh</code>.`);
    return;
  }
  state.manifest = m;
  state.scale = m.logp_scale;
  state.binSize = m.bin_size;
  for (const c of m.contigs) state.contigMap.set(c.name, c);
  // contigs on the axis come from GWAS tracks (pixy has windows on every scaffold); line-only fallback
  const axisTracks = m.tracks.some(t => t.type === 'gwas') ? m.tracks.filter(t => t.type === 'gwas') : m.tracks;
  for (const t of axisTracks) for (const c of Object.keys(t.contigs || {})) state.dataContigs.add(c);
  readHash();
  buildLayout();
  state.gs = 0; state.ge = state.layout.total;
  buildRows();
  wireControls();
  for (const o of [ruler, geneRow, ...tracks]) wirePlot(o.canvas);
  wireOverview();
  setStatus(`Loading ${m.tracks.length} tracks…`);
  const t0 = performance.now();
  const results = await Promise.allSettled([...tracks.map(t => loadTrack(t).catch(err => {
    t.label.querySelector('.sub').textContent = 'failed to load, see console';
    throw err;
  })), loadGeneIndex()]);
  const failed = results.filter(r => r.status === 'rejected');
  failed.forEach(r => console.error(r.reason));
  if (state.pendingLoc) {
    const v = core.decodeView(state.layout, state.pendingLoc);
    if (v) setView(v[0], v[1]);
  }
  updateStatus();
  if (failed.length) setStatus(`${failed.length} track file(s) failed to load, see console`);
  console.info(`GWAS browser: loaded in ${Math.round(performance.now() - t0)} ms`);
  requestRedraw();
  applySelHash(state.pendingSel);
}

function applySelHash(text) {
  const mm = /^([^:]+):(\d+)$/.exec(text || '');
  if (mm && ensureContig(mm[1])) selectSnp(mm[1], Number(mm[2]));
}

// Console handle for debugging and timing, e.g. gwasBrowser.setView(...), gwasBrowser.drawAll().
window.gwasBrowser = {
  state, tracks, setView, drawAll: () => drawAll(), goText, selectSnp,
  setTheme: theme => { document.documentElement.dataset.theme = theme; readPalette(); requestRedraw(); },
};

init();
