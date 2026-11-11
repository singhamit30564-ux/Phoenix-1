/* phoenix-continuous.js — Brownian-motion restricted maximum likelihood for
   continuous characters and gene frequencies (the CONTML method).
   Reference: Felsenstein (1973, 1981c); Cavalli-Sforza & Edwards (1967). */
(function (root, factory) {
  const api = factory(
    typeof require === 'function' ? require('./phoenix-core.js') : root.PCore,
    typeof require === 'function' ? require('./phoenix-infer.js') : root.PInfer,
    typeof require === 'function' ? require('./phoenix-linalg.js') : root.PLinAlg,
    typeof require === 'function' ? require('./phoenix-models.js') : root.PModels
  );
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PContinuous = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function (C, I, L, M) {
  'use strict';

  /**
   * Covariance matrix implied by a tree under Brownian motion:
   * C[i][j] = distance from the root to the most recent common ancestor of i and j,
   * C[i][i] = distance from the root to tip i.
   */
  function treeCovariance(tree) {
    const names = C.leaves(tree).map(l => l.name);
    const ctx = I.treeIndices(tree);
    const depth = new Map();
    C.preorder(tree, n => {
      if (n === tree) depth.set(n, 0);
      else depth.set(n, depth.get(ctx.parent.get(n)) + (n.length || 0));
    });
    const n = names.length;
    const K = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      const li = ctx.leafByName.get(names[i]);
      for (let j = i; j < n; j++) {
        const lj = ctx.leafByName.get(names[j]);
        const a = I.mrca(li, lj, ctx);
        const v = Math.max(depth.get(a), 1e-9);
        K[i][j] = v;
        K[j][i] = v;
      }
    }
    return { names, C: K };
  }

  /**
   * Restricted (REML) profile log-likelihood with the root state and the rate
   * sigma^2 both profiled out, so the value is invariant to where the tree is
   * rooted and to the overall scale of the branch lengths. Only relative branch
   * lengths are identifiable; the scale is fixed afterwards.
   */
  function remlLogLik(Cov, X) {
    const n = Cov.length;
    const d = X[0].length;
    const Lc = L.cholesky(Cov);
    if (!Lc) return null;
    const Ainv = L.inverse(Lc);
    const ones = new Float64Array(n).fill(1);
    const denom = L.quadratic(Ainv, ones);
    if (!(denom > 0)) return null;
    let rss = 0;
    for (let k = 0; k < d; k++) {
      const x = new Float64Array(n);
      for (let i = 0; i < n; i++) x[i] = X[i][k];
      let cross = 0;
      for (let i = 0; i < n; i++) {
        let row = 0;
        for (let j = 0; j < n; j++) row += Ainv[i][j] * x[j];
        cross += ones[i] * row;
      }
      const mu = cross / denom;
      const shifted = new Float64Array(n);
      for (let i = 0; i < n; i++) shifted[i] = x[i] - mu;
      rss += L.quadratic(Ainv, shifted);
    }
    if (!(rss > 0)) return null;
    const dof = (n - 1) * d;
    const sigma2 = rss / dof;
    const logDet = L.logDet(Lc);
    // Restricted likelihood: the ln|1' C^-1 1| term is what makes the value
    // invariant to the position of the root once the mean is profiled out.
    const ll = -0.5 * (
      n * d * Math.log(2 * Math.PI)
      + d * logDet
      + d * Math.log(denom)
      + dof * Math.log(sigma2)
      + dof
    );
    return { logLikelihood: ll, logDet, rss, sigma2, dof };
  }

  function covarianceFor(tree) {
    const { names, C: K } = treeCovariance(tree);
    return { names, K };
  }

  function allEdges(tree) {
    const out = [];
    const parent = C.parentMap(tree);
    for (const [child] of parent) out.push(child);
    return out;
  }

  /** Coordinate-wise golden-section maximisation of the restricted likelihood. */
  function optimiseBranchLengths(tree, X, sweeps) {
    const edges = allEdges(tree);
    let best = null;
    const objective = () => {
      const { K } = covarianceFor(tree);
      return remlLogLik(K, X);
    };
    for (let round = 0; round < (sweeps || 3); round++) {
      for (const edge of edges) {
        const saved = edge.length || 0.1;
        let bestT = Math.max(saved, 1e-4);
        edge.length = bestT;
        let bestFit = objective();
        if (!bestFit) { edge.length = saved; continue; }
        const f = t => {
          edge.length = t;
          const r = objective();
          return r ? r.logLikelihood : -Infinity;
        };
        const fit = M.goldenMax(f, 1e-5, Math.max(4, bestT * 4), 22);
        if (fit.fx > bestFit.logLikelihood) {
          edge.length = fit.x;
          best = objective();
        } else {
          edge.length = bestT;
          best = bestFit;
        }
      }
    }
    return best;
  }

  function euclideanDistances(points) {
    const n = points.length;
    const D = points.map(() => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let sum = 0;
        for (let k = 0; k < points[i].length; k++) sum += (points[i][k] - points[j][k]) ** 2;
        D[i][j] = D[j][i] = Math.sqrt(sum);
      }
    }
    return D;
  }

  /**
   * Full CONTML-style analysis: start from a neighbour-joining tree on
   * Euclidean distances, fit branch lengths by restricted maximum likelihood,
   * then refine the topology by nearest-neighbour interchange.
   */
  function fit(points, names, opts = {}) {
    const n = names.length;
    if (n < 2) throw new Error('CONTML needs at least two taxa.');
    const d = points[0].length;
    if (points.some(r => r.length !== d)) throw new Error('All rows must have the same number of characters.');
    const D = euclideanDistances(points);
    let tree = n === 2
      ? { name: null, children: [{ name: names[0], length: 0.2 }, { name: names[1], length: 0.2 }] }
      : I.neighborJoining(D, names);
    I.assignIds(tree);
    for (const [child] of C.parentMap(tree)) if (!(child.length > 0)) child.length = 0.1;

    let fitResult = optimiseBranchLengths(tree, points, opts.sweeps || 3);
    if (!fitResult) throw new Error('The covariance matrix is not positive definite; check for identical taxa or zero-length branches.');

    let improved = true;
    let rounds = 0;
    while (improved && rounds < (opts.maxNni || 6)) {
      improved = false;
      const current = fitResult.logLikelihood;
      const candidates = nniNeighbours(tree);
      for (const candidate of candidates) {
        I.assignIds(candidate);
        for (const [child] of C.parentMap(candidate)) if (!(child.length > 0)) child.length = 0.1;
        const trial = optimiseBranchLengths(candidate, points, 2);
        if (trial && trial.logLikelihood > current + 1e-6) {
          tree = candidate;
          fitResult = optimiseBranchLengths(tree, points, 2) || trial;
          improved = true;
          break;
        }
      }
      rounds++;
    }
    fitResult = optimiseBranchLengths(tree, points, 3) || fitResult;
    // Fix the unidentifiable scale by rescaling so that sigma^2 = 1; the
    // profile likelihood is unchanged by this rescaling.
    const scale = fitResult.sigma2;
    C.preorder(tree, node => { if (node !== tree) node.length = (node.length || 0) * scale; });
    return {
      tree,
      logLikelihood: fitResult.logLikelihood,
      // After rescaling the rate is exactly one; only relative branch lengths
      // are identifiable from the data.
      sigma2: 1,
      rawSigma2: scale,
      dof: fitResult.dof,
      dimensions: d,
      taxa: n
    };
  }

  function nniNeighbours(tree) {
    const out = [];
    const internals = [];
    C.postorder(tree, n => { if (n.children && n.children.length === 2) internals.push(n.id); });
    for (const id of internals) {
      for (let which = 0; which < 2; which++) {
        for (let sub = 0; sub < 2; sub++) {
          const trial = C.cloneTree(tree);
          const node = I.findById(trial, id);
          if (!node || !node.children) continue;
          const A = node.children[which];
          if (!A || !A.children || A.children.length !== 2) continue;
          const sibling = node.children[1 - which];
          const moved = A.children[sub];
          A.children[sub] = sibling;
          node.children[1 - which] = moved;
          out.push(trial);
        }
      }
    }
    return out;
  }

  /** Cavalli-Sforza & Edwards (1967) variance-stabilising transform. */
  function arcsineSqrt(p) {
    const clamped = Math.min(Math.max(p, 0), 1);
    return Math.asin(Math.sqrt(clamped));
  }

  return { treeCovariance, remlLogLik, fit, optimiseBranchLengths, arcsineSqrt, euclideanDistances };
});
