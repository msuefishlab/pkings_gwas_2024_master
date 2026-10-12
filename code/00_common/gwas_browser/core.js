// Pure helpers for the GWAS browser (no DOM). Tested by test_core.mjs.
// Genome coordinates: g = contig offset + (pos - 1), so a contig of length L covers [offset, offset + L).

// Lay out contigs end to end. `contigs` is [{name, length}] in display order.
export function buildLayout(contigs) {
  const list = [];
  const byName = new Map();
  let offset = 0;
  for (const c of contigs) {
    const entry = { name: c.name, length: c.length, offset, index: list.length };
    list.push(entry);
    byName.set(c.name, entry);
    offset += c.length;
  }
  return { list, byName, total: offset };
}

export function toGenome(layout, contig, pos) {
  const c = layout.byName.get(contig);
  return c ? c.offset + pos - 1 : null;
}

// Contig and 1-based position at genome coordinate g (clamped to the layout).
export function fromGenome(layout, g) {
  const list = layout.list;
  if (!list.length) return null;
  let lo = 0, hi = list.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (list[mid].offset <= g) lo = mid; else hi = mid - 1;
  }
  const c = list[lo];
  const pos = Math.min(c.length, Math.max(1, Math.floor(g - c.offset) + 1));
  return { contig: c.name, pos, entry: c };
}

// Contigs overlapping genome range [gs, ge), each with the visible 1-based position range.
export function visibleContigs(layout, gs, ge) {
  const out = [];
  for (const c of layout.list) {
    const cs = c.offset, ce = c.offset + c.length;
    if (ce <= gs) continue;
    if (cs >= ge) break;
    out.push({
      entry: c,
      posLo: Math.max(1, Math.floor(gs - cs) + 1),
      posHi: Math.min(c.length, Math.ceil(ge - cs) + 1),
    });
  }
  return out;
}

// First index i in sorted arr[lo, hi) with arr[i] >= x.
export function lowerBound(arr, x, lo = 0, hi = arr.length) {
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] < x) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// First index i in sorted arr[lo, hi) with arr[i] > x.
export function upperBound(arr, x, lo = 0, hi = arr.length) {
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (arr[mid] <= x) lo = mid + 1; else hi = mid;
  }
  return lo;
}

// Contig names as typed by a person: "16" -> "chr16" when only the prefixed name exists.
export function resolveContig(name, contigNames) {
  if (contigNames.has(name)) return name;
  const lower = name.toLowerCase();
  for (const n of contigNames) if (n.toLowerCase() === lower) return n;
  if (/^\d+$/.test(name) && contigNames.has('chr' + name)) return 'chr' + name;
  return null;
}

// Parse "chr16", "chr16:10,800,000-11,100,000" or "chr16:10916883".
// Returns {contig, start, end} (1-based inclusive) or null. A single position is centred in `pad` bp.
export function parseLocus(text, contigs, pad = 25000) {
  const m = /^\s*([A-Za-z0-9_.|-]+?)\s*(?::\s*([\d,]+)\s*(?:(?:-|\.\.)\s*([\d,]+))?)?\s*$/.exec(text);
  if (!m) return null;
  const names = contigs instanceof Map ? new Set(contigs.keys()) : new Set(contigs.map(c => c.name));
  const contig = resolveContig(m[1], names);
  if (!contig) return null;
  const length = contigs instanceof Map ? contigs.get(contig).length
    : contigs.find(c => c.name === contig).length;
  if (m[2] === undefined) return { contig, start: 1, end: length };
  const a = Number(m[2].replace(/,/g, ''));
  if (m[3] === undefined) {
    return { contig, start: Math.max(1, a - pad), end: Math.min(length, a + pad) };
  }
  const b = Number(m[3].replace(/,/g, ''));
  const start = Math.max(1, Math.min(a, b)), end = Math.min(length, Math.max(a, b));
  return end >= start ? { contig, start, end } : null;
}

export function formatInt(n) {
  return Math.round(n).toLocaleString('en-US');
}

export function formatBp(bp) {
  const a = Math.abs(bp);
  if (a >= 1e6) return (bp / 1e6).toFixed(a >= 1e7 ? 1 : 2).replace(/\.?0+$/, '') + ' Mb';
  if (a >= 1e3) return (bp / 1e3).toFixed(a >= 1e4 ? 0 : 1).replace(/\.0$/, '') + ' kb';
  return Math.round(bp) + ' bp';
}

// Roughly n "nice" tick values covering [min, max].
export function niceTicks(min, max, n = 4) {
  if (!(max > min)) return [min];
  const raw = (max - min) / Math.max(1, n);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const ticks = [];
  for (let v = Math.ceil(min / step) * step; v <= max + step * 1e-9; v += step) {
    ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  }
  return ticks;
}

// Max value per pixel column for sorted points pos[lo, hi) (values val, scaled ints, 0 = none).
// Column of a point: floor((offset + pos - 1 - gs) / bpPerPx).
export function columnMax(pos, val, lo, hi, offset, gs, bpPerPx, width, out) {
  const cols = out || new Uint16Array(width);
  for (let i = lo; i < hi; i++) {
    const x = Math.floor((offset + pos[i] - 1 - gs) / bpPerPx);
    if (x < 0 || x >= width) continue;
    if (val[i] > cols[x]) cols[x] = val[i];
  }
  return cols;
}

// Per-column max of fixed-size bins (bin k covers positions k*binSize+1 .. (k+1)*binSize).
export function binColumnMax(bins, binOff, binN, binSize, offset, gs, bpPerPx, width, out) {
  const cols = out || new Uint16Array(width);
  const kLo = Math.max(0, Math.floor((gs - offset) / binSize));
  const kHi = Math.min(binN, Math.ceil((gs + width * bpPerPx - offset) / binSize));
  for (let k = kLo; k < kHi; k++) {
    const v = bins[binOff + k];
    if (!v) continue;
    const g0 = offset + k * binSize, g1 = g0 + binSize;
    let x0 = Math.floor((g0 - gs) / bpPerPx), x1 = Math.floor((g1 - 1 - gs) / bpPerPx);
    if (x1 < 0 || x0 >= width) continue;
    x0 = Math.max(0, x0); x1 = Math.min(width - 1, x1);
    for (let x = x0; x <= x1; x++) if (v > cols[x]) cols[x] = v;
  }
  return cols;
}

// Aggregate windows [start, end] with values into per-column min/max/sum/count (NaN skipped).
export function windowColumns(starts, ends, vals, lo, hi, offset, gs, bpPerPx, width, agg) {
  const a = agg || {
    min: new Float32Array(width).fill(Infinity),
    max: new Float32Array(width).fill(-Infinity),
    sum: new Float64Array(width),
    n: new Uint32Array(width),
  };
  for (let i = lo; i < hi; i++) {
    const v = vals[i];
    if (v !== v) continue;
    let x0 = Math.floor((offset + starts[i] - 1 - gs) / bpPerPx);
    let x1 = Math.floor((offset + ends[i] - 1 - gs) / bpPerPx);
    if (x1 < 0 || x0 >= width) continue;
    x0 = Math.max(0, x0); x1 = Math.min(width - 1, x1);
    for (let x = x0; x <= x1; x++) {
      if (v < a.min[x]) a.min[x] = v;
      if (v > a.max[x]) a.max[x] = v;
      a.sum[x] += v;
      a.n[x]++;
    }
  }
  return a;
}

// Nearest gene to a position, strand-aware. genes: [{s, e, st, n, x}] sorted by s.
// Returns {gene, where: 'exon'|'intron'|'gene'|'upstream'|'downstream', dist} or null.
export function nearestGene(genes, pos) {
  let inside = null, insideRank = 9, best = null, bestDist = Infinity;
  for (const g of genes) {
    if (pos >= g.s && pos <= g.e) {
      let where = 'gene';
      if (g.x && g.x.length) {
        where = 'intron';
        for (let i = 0; i < g.x.length; i += 2) {
          if (pos >= g.x[i] && pos <= g.x[i + 1]) { where = 'exon'; break; }
        }
      }
      const rank = where === 'exon' ? 0 : where === 'intron' ? 1 : 2;
      if (rank < insideRank || (rank === insideRank && g.e - g.s < inside.gene.e - inside.gene.s)) {
        inside = { gene: g, where, dist: 0 };
        insideRank = rank;
      }
    } else if (!inside) {
      const before = pos < g.s;
      const dist = before ? g.s - pos : pos - g.e;
      if (dist < bestDist) {
        bestDist = dist;
        const upstream = (before && g.st !== '-') || (!before && g.st === '-');
        best = { gene: g, where: upstream ? 'upstream' : 'downstream', dist };
      }
    }
  }
  return inside || best;
}

export function describeRelation(rel) {
  if (!rel) return 'no gene on this contig';
  const name = rel.gene.n || rel.gene.id;
  if (rel.where === 'exon') return `in an exon of ${name}`;
  if (rel.where === 'intron') return `in an intron of ${name}`;
  if (rel.where === 'gene') return `inside ${name}`;
  const side = rel.where === 'upstream' ? "5′ (upstream)" : "3′ (downstream)";
  return `${formatBp(rel.dist)} ${side} of ${name}`;
}

// Rank gene index entries [{name, description, ...}] for a query; best first.
export function searchGenes(entries, query, limit = 12) {
  const q = query.trim().toLowerCase();
  if (!q) return [];
  const scored = [];
  for (const e of entries) {
    const name = e.nameLower;
    let score = -1;
    if (name === q) score = 0;
    else if (name.startsWith(q)) score = 1;
    else if (name.includes(q)) score = 2;
    else if (q.length >= 3 && e.descLower) {
      const i = e.descLower.indexOf(q);
      if (i === 0 || (i > 0 && /[^a-z0-9]/.test(e.descLower[i - 1]))) score = 3;
      else if (i > 0) score = 4;
    }
    if (score >= 0) scored.push([score, name.length, e]);
  }
  scored.sort((a, b) => a[0] - b[0] || a[1] - b[1] || (a[2].name < b[2].name ? -1 : 1));
  return scored.slice(0, limit).map(s => s[2]);
}

// Greedy row packing of [x0, x1] extents (sorted by x0). Returns the row of each item.
export function packRows(extents, maxRows, gap = 4) {
  const rowEnd = [];
  return extents.map(([x0, x1]) => {
    let r = rowEnd.findIndex(end => end + gap <= x0);
    if (r < 0) {
      if (rowEnd.length < maxRows) { r = rowEnd.length; rowEnd.push(-Infinity); }
      else r = rowEnd.indexOf(Math.min(...rowEnd));
    }
    rowEnd[r] = Math.max(rowEnd[r], x1);
    return r;
  });
}

// Parse genes/index.tsv text into entries grouped by contig.
export function parseGeneIndex(text) {
  const entries = [];
  const lines = text.split('\n');
  const header = lines[0].split('\t');
  const col = name => header.indexOf(name);
  const iName = col('name'), iContig = col('contig'), iStart = col('start'), iEnd = col('end'),
    iStrand = col('strand'), iBiotype = col('biotype'), iDesc = col('description');
  for (let k = 1; k < lines.length; k++) {
    const line = lines[k];
    if (!line) continue;
    const f = line.split('\t');
    const description = iDesc >= 0 ? (f[iDesc] || '') : '';
    entries.push({
      name: f[iName], nameLower: f[iName].toLowerCase(), contig: f[iContig],
      start: Number(f[iStart]), end: Number(f[iEnd]),
      strand: iStrand >= 0 ? f[iStrand] : '.', biotype: iBiotype >= 0 ? f[iBiotype] : '',
      description, descLower: description.toLowerCase(),
    });
  }
  return entries;
}

// Locus state for the URL hash: "chr16:10800000-11100000" spanning possibly several contigs is
// stored as genome coordinates relative to contigs, so it survives layout changes.
export function encodeView(layout, gs, ge) {
  const a = fromGenome(layout, gs), b = fromGenome(layout, Math.max(gs, ge - 1));
  if (!a || !b) return '';
  if (a.contig === b.contig) return `${a.contig}:${a.pos}-${b.pos}`;
  return `${a.contig}:${a.pos}-${b.contig}:${b.pos}`;
}

export function decodeView(layout, text) {
  const m = /^([^:]+):(\d+)-(?:([^:]+):)?(\d+)$/.exec(text || '');
  if (!m) return null;
  const gs = toGenome(layout, m[1], Number(m[2]));
  const ge = toGenome(layout, m[3] || m[1], Number(m[4]));
  if (gs === null || ge === null || ge < gs) return null;
  return [gs, ge + 1];
}
