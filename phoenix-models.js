/* phoenix-models.js — substitution models, transition probabilities and
   pairwise evolutionary distances (DNADIST / PROTDIST / RESTDIST core). */
(function (root, factory) {
  const api = factory(typeof require === 'function' ? require('./phoenix-core.js') : root.PCore);
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PModels = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (PCore) {
  'use strict';

  const DNA_INDEX = { A: 0, C: 1, G: 2, T: 3 };
  const AA_STATES = PCore.AA;

  const DNA_MODELS = {
    jc: { label: 'Jukes–Cantor (equal rates)', needs: [] },
    f81: { label: 'F81 (empirical frequencies)', needs: [] },
    k2p: { label: 'Kimura 2-parameter (+ κ)', needs: ['kappa'] },
    f84: { label: 'F84 (+ κ)', needs: ['kappa'] },
    hky: { label: 'HKY85 (+ κ)', needs: ['kappa'] },
    tn93: { label: 'Tamura–Nei (+ κ1, κ2)', needs: ['kappa', 'kappa2'] },
    gtr: { label: 'General time reversible (+ 6 rates)', needs: ['rates'] },
    cat: { label: 'Discrete-gamma rates (+ α, k)', needs: ['alpha', 'categories'] }
  };

  const PROTEIN_MODELS = {
    poisson: { label: 'Poisson (equal rates)', needs: [] },
    'equal-input': { label: 'Equal-input (empirical frequencies)', needs: [] }
  };

  function uniformFreqs(n) {
    return new Array(n).fill(1 / n);
  }

  function normalizeFreqs(f) {
    const s = f.reduce((a, b) => a + b, 0);
    return s > 0 ? f.map(x => x / s) : uniformFreqs(f.length);
  }

  /** Empirical base/state frequencies from encoded per-site vectors. */
  function empiricalFrequencies(chars, nStates) {
    const counts = new Array(nStates).fill(0);
    for (const col of chars) {
      for (const v of col) {
        let total = 0;
        for (let k = 0; k < nStates; k++) total += v[k];
        if (!total) continue;
        for (let k = 0; k < nStates; k++) counts[k] += v[k] / total;
      }
    }
    const sum = counts.reduce((a, b) => a + b, 0);
    if (!sum) return uniformFreqs(nStates);
    return normalizeFreqs(counts);
  }

  /** Relative exchangeability matrix (off-diagonal) for a DNA model. */
  function dnaExchange(model, p) {
    const r = [1, 1, 1, 1].map(() => [1, 1, 1, 1]);
    const k = p.kappa || 2;
    const k2 = p.kappa2 || 1;
    const set = (i, j, v) => { r[i][j] = v; r[j][i] = v; };
    if (model === 'k2p' || model === 'f84' || model === 'hky') {
      set(0, 2, k); set(1, 3, k);
    } else if (model === 'tn93') {
      set(0, 2, k); set(1, 3, k2);
    } else if (model === 'gtr') {
      const g = p.rates || [2, 2, 1, 1, 1, 1];
      set(0, 2, g[0]); set(1, 3, g[1]);
      set(0, 1, g[2]); set(0, 3, g[3]); set(1, 2, g[4]); set(2, 3, g[5]);
    }
    for (let i = 0; i < 4; i++) r[i][i] = 0;
    return r;
  }

  /** Build a normalized rate generator Q (expected substitutions = 1). */
  function buildQ(opts) {
    const kind = opts.kind || 'dna';
    const n = kind === 'dna' ? 4 : kind === 'protein' ? 20 : 2;
    let pi = opts.freqs ? normalizeFreqs(opts.freqs.slice()) : uniformFreqs(n);
    if (n !== pi.length) pi = uniformFreqs(n);
    let r;
    if (kind === 'protein') {
      r = pi.map((_, i) => pi.map((_, j) => (i === j ? 0 : 1)));
    } else if (kind === 'sites') {
      r = [[0, 1], [1, 0]];
    } else {
      r = dnaExchange(opts.model || 'jc', opts);
      if ((opts.model === 'jc' || opts.model === 'k2p') && opts.freqs === undefined) pi = uniformFreqs(4);
    }
    const Q = pi.map((_, i) => pi.map((_, j) => (i === j ? 0 : r[i][j] * pi[j])));
    for (let i = 0; i < n; i++) {
      const rowSum = Q[i].reduce((a, b) => a + b, 0);
      Q[i][i] = -rowSum;
    }
    let rate = 0;
    for (let i = 0; i < n; i++) rate += -pi[i] * Q[i][i];
    if (rate > 0) for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Q[i][j] /= rate;
    return { Q, pi, n, kind };
  }

  function identity(n) {
    const I = [];
    for (let i = 0; i < n; i++) { const row = new Array(n).fill(0); row[i] = 1; I.push(row); }
    return I;
  }

  function matMul(A, B, n) {
    const out = new Array(n);
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      const ai = A[i];
      for (let k = 0; k < n; k++) {
        const a = ai[k];
        if (a === 0) continue;
        const bk = B[k];
        for (let j = 0; j < n; j++) row[j] += a * bk[j];
      }
      out[i] = row;
    }
    return out;
  }

  /** P(t) = exp(Qt) by scaling, Taylor series and repeated squaring. */
  function transitionMatrix(Q, t, n) {
    if (!(t > 0)) return identity(n);
    let norm = 0;
    for (let i = 0; i < n; i++) norm = Math.max(norm, -Q[i][i]);
    let scaled = t, squarings = 0;
    while (norm * scaled > 0.25 && squarings < 50) { scaled /= 2; squarings++; }
    const A = Q.map(row => row.map(v => v * scaled));
    let P = identity(n);
    let term = identity(n);
    for (let order = 1; order <= 14; order++) {
      term = matMul(term, A, n);
      const inv = 1 / order;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) { term[i][j] *= inv; P[i][j] += term[i][j]; }
      let big = 0;
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) big = Math.max(big, Math.abs(term[i][j]));
      if (big < 1e-18) break;
    }
    for (let k = 0; k < squarings; k++) P = matMul(P, P, n);
    return P;
  }

  // ---------------------------------------------------------- base counting

  function baseCounts(a, b, alphaIndex) {
    // Pairs where both taxa have an unambiguous state.
    const n = a.length;
    const counts = new Array(n * n).fill(0);
    let comparable = 0;
    for (let s = 0; s < a.length; s++) {
      const x = alphaIndex[a[s]], y = alphaIndex[b[s]];
      if (x === undefined || y === undefined) continue;
      counts[x * n + y]++;
      comparable++;
    }
    return { counts, comparable, n };
  }

  const DNA_LOOKUP = { A: 0, C: 1, G: 2, T: 3, U: 3 };

  function pDistance(a, b) {
    let compared = 0, diff = 0;
    for (let s = 0; s < a.length; s++) {
      const x = DNA_LOOKUP[a[s]], y = DNA_LOOKUP[b[s]];
      if (x === undefined || y === undefined) continue;
      compared++;
      if (x !== y) diff++;
    }
    if (!compared) return { comparable: 0, p: NaN };
    return { comparable: compared, p: diff / compared };
  }

  function jukesCantor(p, states) {
    const b = (states - 1) / states;
    const arg = 1 - p / b;
    if (arg <= 0) return { distance: Infinity, saturated: true };
    return { distance: -b * Math.log(arg) };
  }

  function kimura2p(pTrans, pTransv) {
    const P = pTrans, Q = pTransv;
    const a = 1 - 2 * P - Q, b = 1 - 2 * Q;
    if (a <= 0 || b <= 0) return { distance: Infinity, saturated: true };
    return { distance: 0.5 * Math.log(1 / a) + 0.25 * Math.log(1 / b) };
  }

  /** Tamura–Nei (1993) closed form. p1/p2 are purine/pyrimidine
      transition proportions, pQ the transversion proportion. */
  function tn93(pi, p1, p2, pQ) {
    const gR = pi[0] + pi[2], gY = pi[1] + pi[3];
    const h = 2 * gR * gY;
    if (!(h > 0)) return { distance: Infinity, saturated: true };
    const c1 = 1 - p1 / gR - pQ / h;
    const c2 = 1 - p2 / gY - pQ / h;
    const c3 = 1 - pQ / h;
    if (c1 <= 0 || c2 <= 0 || c3 <= 0) return { distance: Infinity, saturated: true };
    const d = -h * Math.log(c3) - (gR * h) * Math.log(c1) - (gY * h) * Math.log(c2) + 0;
    return { distance: Math.max(0, d) };
  }

  // ------------------------------------------------------- ML pair distance

  function pairLogLik(a, b, Q, pi, n) {
    // Pruning on a two-taxon tree: sum over states of pi * P(ab) * ...
    return t => {
      const P = transitionMatrix(Q, t, n);
      let ll = 0, sites = 0;
      for (let s = 0; s < a.length; s++) {
        const x = DNA_LOOKUP[a[s]], y = DNA_LOOKUP[b[s]];
        if (x === undefined || y === undefined) continue;
        sites++;
        let p = 0;
        for (let i = 0; i < n; i++) if (i === x) p += pi[i] * P[i][y];
        if (p > 0) ll += Math.log(p);
      }
      if (!sites) return -Infinity;
      return { ll, sites, perSite: ll / sites };
    };
  }

  function goldenMax(fn, lo, hi, iters = 60) {
    const gr = (Math.sqrt(5) - 1) / 2;
    let a = lo, b = hi;
    let c = b - gr * (b - a), d = a + gr * (b - a);
    let fc = fn(c), fd = fn(d);
    for (let k = 0; k < iters; k++) {
      if (fc < fd) { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = fn(d); }
      else { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = fn(c); }
    }
    return { x: (a + b) / 2, fx: fn((a + b) / 2) };
  }

  function mlDistance(a, b, opts) {
    const model = buildQ(Object.assign({ kind: 'dna' }, opts));
    const objective = pairLogLik(a, b, model.Q, model.pi, model.n);
    if (objective(1e-9) === -Infinity) return { comparable: 0, distance: NaN };
    const maximize = t => {
      const r = objective(t);
      return r === -Infinity ? -Infinity : r.perSite;
    };
    const fit = goldenMax(maximize, 1e-7, 10);
    const atZero = maximize(1e-7), atBig = maximize(10);
    if (!isFinite(fit.fx)) return { comparable: 0, distance: NaN, saturated: true };
    void atZero; void atBig;
    return { distance: fit.x, perSite: fit.fx };
  }

  // ------------------------------------------------------------ DNA distance

  /** Base composition pooled over one pair of sequences. */
  function empiricalPairFreqs(a, b) {
    const counts = [0, 0, 0, 0];
    for (const seq of [a, b]) {
      for (const ch of seq) {
        const i = DNA_LOOKUP[ch];
        if (i !== undefined) counts[i]++;
      }
    }
    const total = counts.reduce((x, y) => x + y, 0);
    return total ? counts.map(c => c / total) : [0.25, 0.25, 0.25, 0.25];
  }

  function dnaDistance(a, b, opts = {}) {
    const model = opts.model || 'jc';
    const { p, comparable } = pDistance(a, b);
    if (!comparable) return { comparable: 0, distance: NaN };
    if (model === 'p') return { comparable, distance: p, p };
    if (model === 'ml') return Object.assign({ comparable, p }, mlDistance(a, b, opts));
    if (model === 'jc') return Object.assign({ comparable, p }, jukesCantor(p, 4));

    // Base composition counts for the analytic corrections.
    const pi = opts.freqs ? normalizeFreqs(opts.freqs) : [0.25, 0.25, 0.25, 0.25];
    const { counts } = baseCounts(a, b, DNA_LOOKUP);
    const at = (i, j) => counts[i * 4 + j] / comparable;
    const pAG = at(0, 2), pGA = at(2, 0), pCT = at(1, 3), pTC = at(3, 1);
    const pR = pAG + pGA, pY = pCT + pTC;
    const pQ = (p - pR - pY);
    if (model === 'k2p') return Object.assign({ comparable, p }, kimura2p(pR, pQ));
    if (model === 'f84' || model === 'hky' || model === 'tn93' || model === 'gtr' || model === 'cat') {
      const freqs = opts.freqs || empiricalPairFreqs(a, b);
      return Object.assign({ comparable, p }, mlDistance(a, b, Object.assign({}, opts, { kind: 'dna', freqs })));
    }
    return Object.assign({ comparable, p }, jukesCantor(p, 4));
  }

  // -------------------------------------------------------- protein distance

  function proteinDistance(a, b, opts = {}) {
    const model = opts.model || 'poisson';
    let compared = 0, diff = 0;
    for (let s = 0; s < a.length; s++) {
      const x = AA_STATES.indexOf(a[s]), y = AA_STATES.indexOf(b[s]);
      if (x < 0 || y < 0) continue;
      compared++;
      if (x !== y) diff++;
    }
    if (!compared) return { comparable: 0, distance: NaN };
    const p = diff / compared;
    if (model === 'p') return { comparable: compared, p, distance: p };
    if (model === 'ml') {
      return Object.assign({ comparable: compared, p }, mlDistanceProtein(a, b, opts));
    }
    if (model === 'poisson') {
      const arg = 1 - p;
      return { comparable: compared, p, distance: arg > 0 ? -Math.log(arg) : Infinity, saturated: arg <= 0 };
    }
    const pi = opts.freqs ? normalizeFreqs(opts.freqs) : uniformFreqs(20);
    const b0 = 1 - pi.reduce((s, x) => s + x * x, 0);
    const arg = 1 - p / b0;
    return { comparable: compared, p, distance: arg > 0 ? -b0 * Math.log(arg) : Infinity, saturated: arg <= 0 };
  }

  function mlDistanceProtein(a, b, opts) {
    const model = buildQ({ kind: 'protein', freqs: opts.freqs, model: 'poisson' });
    const objective = t => {
      const P = transitionMatrix(model.Q, t, model.n);
      let ll = 0, sites = 0;
      for (let s = 0; s < a.length; s++) {
        const x = AA_STATES.indexOf(a[s]), y = AA_STATES.indexOf(b[s]);
        if (x < 0 || y < 0) continue;
        sites++;
        let pr = 0;
        for (let i = 0; i < model.n; i++) if (i === x) pr += model.pi[i] * P[i][y];
        if (pr > 0) ll += Math.log(pr);
      }
      return sites ? ll / sites : -Infinity;
    };
    const fit = goldenMax(objective, 1e-7, 20);
    return { distance: fit.x, perSite: fit.fx };
  }

  // -------------------------------------------------- restriction-site model

  function siteDistance(a, b, opts = {}) {
    const model = opts.model || 'fraction';
    let compared = 0, diff = 0;
    for (let s = 0; s < a.length; s++) {
      const x = a[s], y = b[s];
      if (!'01'.includes(x) || !'01'.includes(y)) continue;
      compared++;
      if (x !== y) diff++;
    }
    if (!compared) return { comparable: 0, distance: NaN };
    const p = diff / compared;
    if (model === 'fraction') return { comparable: compared, p, distance: p };
    if (model === 'ml') {
      const m = buildQ({ kind: 'sites' });
      const objective = t => {
        const P = transitionMatrix(m.Q, t, 2);
        let ll = 0, sites = 0;
        for (let s = 0; s < a.length; s++) {
          const x = '01'.indexOf(a[s]), y = '01'.indexOf(b[s]);
          if (x < 0 || y < 0) continue;
          sites++;
          let pr = 0;
          for (let i = 0; i < 2; i++) if (i === x) pr += m.pi[i] * P[i][y];
          if (pr > 0) ll += Math.log(pr);
        }
        return sites ? ll / sites : -Infinity;
      };
      const fit = goldenMax(objective, 1e-7, 5);
      return { comparable: compared, p, distance: fit.x };
    }
    const arg = 1 - 2 * p;
    return { comparable: compared, p, distance: arg > 0 ? -0.5 * Math.log(arg) : Infinity, saturated: arg <= 0 };
  }

  /** Discrete gamma rate multipliers (mean 1). */
  function gammaRates(alpha, k) {
    if (!(alpha > 0) || !(k >= 1)) return [1];
    const rates = [];
    for (let i = 0; i < k; i++) {
      const q = (i + 0.5) / k;
      rates.push(approxGammaQuantile(q, alpha) * alpha);
    }
    const mean = rates.reduce((a, b) => a + b, 0) / rates.length;
    return rates.map(r => r / mean);
  }

  function approxGammaQuantile(q, alpha) {
    // Wilson–Hilferty approximation, adequate for rate categories.
    const z = inverseNormal(q);
    const t = 1 - 1 / (9 * alpha) + z / Math.sqrt(9 * alpha);
    return alpha * Math.max(t, 1e-6) ** 3;
  }

  function inverseNormal(p) {
    const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
      1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
    const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
      6.680131188771972e+01, -1.328068155288572e+01];
    const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
      -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
    const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    const pl = 0.02425;
    if (p < pl) {
      const q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    if (p > 1 - pl) {
      const q = Math.sqrt(-2 * Math.log(1 - p));
      return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    }
    const q = p - 0.5, r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }

  return {
    DNA_MODELS, PROTEIN_MODELS, buildQ, transitionMatrix, empiricalFrequencies,
    normalizeFreqs, uniformFreqs, pDistance, jukesCantor, kimura2p, tn93,
    dnaDistance, proteinDistance, siteDistance, mlDistance, gammaRates, goldenMax
  };
});
