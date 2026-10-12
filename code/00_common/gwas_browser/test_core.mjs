// Run: node --test code/00_common/gwas_browser/test_core.mjs
import test from 'node:test';
import assert from 'node:assert/strict';
import * as core from './core.js';

const contigs = [{ name: 'chr1', length: 1000 }, { name: 'chr2', length: 500 }, { name: 'ups3', length: 100 }];
const layout = core.buildLayout(contigs);

test('layout offsets and round trip', () => {
  assert.equal(layout.total, 1600);
  assert.equal(core.toGenome(layout, 'chr2', 1), 1000);
  assert.deepEqual(core.fromGenome(layout, 1000).pos, 1);
  assert.equal(core.fromGenome(layout, 1000).contig, 'chr2');
  assert.equal(core.fromGenome(layout, 999).contig, 'chr1');
  assert.equal(core.fromGenome(layout, 999).pos, 1000);
  assert.equal(core.fromGenome(layout, 99999).contig, 'ups3');
});

test('visible contigs', () => {
  const v = core.visibleContigs(layout, 900, 1100);
  assert.deepEqual(v.map(c => c.entry.name), ['chr1', 'chr2']);
  assert.equal(v[0].posLo, 901);
  assert.equal(v[1].posHi, 101);
});

test('bounds', () => {
  const a = Uint32Array.from([1, 3, 3, 7]);
  assert.equal(core.lowerBound(a, 3), 1);
  assert.equal(core.upperBound(a, 3), 3);
  assert.equal(core.lowerBound(a, 8), 4);
  assert.equal(core.lowerBound(a, 0), 0);
});

test('parse locus', () => {
  assert.deepEqual(core.parseLocus('chr1', contigs), { contig: 'chr1', start: 1, end: 1000 });
  assert.deepEqual(core.parseLocus('chr1:100-200', contigs), { contig: 'chr1', start: 100, end: 200 });
  assert.deepEqual(core.parseLocus(' 2:1,00-4,00 ', contigs), { contig: 'chr2', start: 100, end: 400 });
  assert.deepEqual(core.parseLocus('chr1:500', contigs, 50), { contig: 'chr1', start: 450, end: 550 });
  assert.deepEqual(core.parseLocus('chr2:300-100', contigs), { contig: 'chr2', start: 100, end: 300 });
  assert.equal(core.parseLocus('chr9:1-2', contigs), null);
  assert.equal(core.parseLocus('ikbip', contigs), null);
});

test('format', () => {
  assert.equal(core.formatBp(10916883), '10.9 Mb');
  assert.equal(core.formatBp(1500000), '1.5 Mb');
  assert.equal(core.formatBp(12345), '12 kb');
  assert.equal(core.formatBp(1500), '1.5 kb');
  assert.equal(core.formatBp(12), '12 bp');
  assert.deepEqual(core.niceTicks(0, 9.7, 3), [0, 5]);
  assert.deepEqual(core.niceTicks(-1, 1, 4), [-1, -0.5, 0, 0.5, 1]);
});

test('column max of points and bins', () => {
  const pos = Uint32Array.from([1, 5, 6, 10]);
  const val = Uint16Array.from([10, 30, 20, 40]);
  // chr2 points, view starts at chr2:1, 5 bp per pixel, 2 columns
  const cols = core.columnMax(pos, val, 0, 4, 1000, 1000, 5, 2);
  assert.deepEqual([...cols], [30, 40]);
  const bins = Uint16Array.from([0, 7, 9, 3]);
  // bins of 10 bp on chr1 (offset 0); view [0, 40) at 20 bp per pixel -> columns {bins 0,1} {bins 2,3}
  assert.deepEqual([...core.binColumnMax(bins, 0, 4, 10, 0, 0, 20, 2)], [7, 9]);
  // zoomed in: one bin spans several columns
  assert.deepEqual([...core.binColumnMax(bins, 0, 4, 10, 0, 10, 5, 4)], [7, 7, 9, 9]);
});

test('window columns skip NaN', () => {
  const s = Uint32Array.from([1, 11, 21]), e = Uint32Array.from([10, 20, 30]);
  const v = Float32Array.from([1, NaN, 3]);
  const a = core.windowColumns(s, e, v, 0, 3, 0, 0, 10, 3);
  assert.deepEqual([...a.n], [1, 0, 1]);
  assert.equal(a.max[2], 3);
});

test('nearest gene', () => {
  const genes = [
    { s: 100, e: 200, st: '+', n: 'a', x: [100, 120, 180, 200] },
    { s: 500, e: 600, st: '-', n: 'b', x: [] },
  ];
  assert.equal(core.nearestGene(genes, 110).where, 'exon');
  assert.equal(core.nearestGene(genes, 150).where, 'intron');
  assert.equal(core.nearestGene(genes, 550).where, 'gene');
  const up = core.nearestGene(genes, 90);
  assert.equal(up.gene.n, 'a'); assert.equal(up.where, 'upstream'); assert.equal(up.dist, 10);
  const minus = core.nearestGene(genes, 650);
  assert.equal(minus.gene.n, 'b'); assert.equal(minus.where, 'upstream');
  assert.equal(core.nearestGene(genes, 260).where, 'downstream');
  assert.match(core.describeRelation(up), /10 bp 5′/);
  assert.equal(core.nearestGene([], 5), null);
});

test('gene search ranking', () => {
  const entries = core.parseGeneIndex(
    'name\tcontig\tstart\tend\tstrand\tbiotype\tdescription\n' +
    'kcnj12\tchr1\t1\t2\t+\tprotein_coding\tpotassium channel\n' +
    'kcnj2a\tchr1\t3\t4\t+\tprotein_coding\tpotassium inwardly rectifying\n' +
    'scn4aa\tchr2\t5\t6\t-\tprotein_coding\tsodium channel, voltage-gated\n' +
    'kcn\tchr2\t7\t8\t+\tprotein_coding\tx\n');
  assert.deepEqual(core.searchGenes(entries, 'kcn').map(e => e.name), ['kcn', 'kcnj12', 'kcnj2a']);
  assert.deepEqual(core.searchGenes(entries, 'sodium').map(e => e.name), ['scn4aa']);
  assert.deepEqual(core.searchGenes(entries, 'channel').map(e => e.name), ['kcnj12', 'scn4aa']);
});

test('row packing', () => {
  assert.deepEqual(core.packRows([[0, 10], [5, 20], [16, 30], [40, 50]], 3), [0, 1, 0, 0]);
  assert.deepEqual(core.packRows([[0, 10], [1, 10], [2, 10]], 2), [0, 1, 0]);
});

test('view hash round trip', () => {
  const t = core.encodeView(layout, 950, 1100);
  assert.equal(t, 'chr1:951-chr2:100');
  assert.deepEqual(core.decodeView(layout, t), [950, 1100]);
  assert.deepEqual(core.decodeView(layout, 'chr2:1-500'), [1000, 1500]);
  assert.equal(core.decodeView(layout, 'nope'), null);
});
