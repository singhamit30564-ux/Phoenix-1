/* phoenix-linalg.js — small dense linear algebra used by the multivariate
   likelihood code. Cholesky based, no external dependencies. */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.PLinAlg = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';

  /** Cholesky factorisation of a symmetric positive definite matrix. */
  function cholesky(A) {
    const n = A.length;
    const L = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = A[i][j];
        for (let k = 0; k < j; k++) sum -= L[i][k] * L[j][k];
        if (i === j) {
          if (!(sum > 0)) return null;
          L[i][i] = Math.sqrt(sum);
        } else {
          L[i][j] = sum / L[j][j];
        }
      }
    }
    return L;
  }

  function logDet(L) {
    let total = 0;
    for (let i = 0; i < L.length; i++) total += Math.log(L[i][i]);
    return 2 * total;
  }

  /** Solve L Lᵀ x = b for a single right-hand side. */
  function solve(L, b) {
    const n = L.length;
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let sum = b[i];
      for (let k = 0; k < i; k++) sum -= L[i][k] * y[k];
      y[i] = sum / L[i][i];
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
      let sum = y[i];
      for (let k = i + 1; k < n; k++) sum -= L[k][i] * x[k];
      x[i] = sum / L[i][i];
    }
    return x;
  }

  /** Inverse of a symmetric positive definite matrix from its Cholesky factor. */
  function inverse(L) {
    const n = L.length;
    const out = Array.from({ length: n }, () => new Float64Array(n));
    for (let j = 0; j < n; j++) {
      const e = new Float64Array(n);
      e[j] = 1;
      const col = solve(L, e);
      for (let i = 0; i < n; i++) out[i][j] = col[i];
    }
    return out;
  }

  /** Quadratic form xᵀ A⁻¹ x given A⁻¹. */
  function quadratic(Ainv, x) {
    let total = 0;
    for (let i = 0; i < x.length; i++) {
      let row = 0;
      for (let j = 0; j < x.length; j++) row += Ainv[i][j] * x[j];
      total += x[i] * row;
    }
    return total;
  }

  /** Determinant of a small matrix by LU with partial pivoting. */
  function determinant(M) {
    const n = M.length;
    const a = M.map(r => Array.from(r));
    let det = 1;
    for (let i = 0; i < n; i++) {
      let pivot = i;
      for (let r = i + 1; r < n; r++) if (Math.abs(a[r][i]) > Math.abs(a[pivot][i])) pivot = r;
      if (Math.abs(a[pivot][i]) < 1e-300) return 0;
      if (pivot !== i) { const t = a[pivot]; a[pivot] = a[i]; a[i] = t; det = -det; }
      det *= a[i][i];
      for (let r = i + 1; r < n; r++) {
        const f = a[r][i] / a[i][i];
        for (let c = i; c < n; c++) a[r][c] -= f * a[i][c];
      }
    }
    return det;
  }

  return { cholesky, logDet, solve, inverse, quadratic, determinant };
});
