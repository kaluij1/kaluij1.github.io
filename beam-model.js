/*
 * Cantilever beam, Euler-Bernoulli finite elements.
 *
 * Each node has two degrees of freedom: vertical displacement v (up positive) and rotation θ.
 * Each element has the 4 × 4 stiffness matrix
 *
 *            EI   |  12    6l   -12    6l  |
 *     k  =  ----  |  6l   4l²   -6l   2l²  |        l = element length
 *            l³   | -12   -6l    12   -6l  |
 *                 |  6l   2l²   -6l   4l²  |
 *
 * The global matrix K is assembled from the elements, the two degrees of freedom at the fixed
 * end are removed, and K u = F is solved by Gaussian elimination with partial pivoting.
 * Bending moment comes from the solution: M = EI v'' with v'' from the element shape functions.
 * Bending stress: σ = -M y / I.
 *
 * Linear elastic, small displacements, no shear deformation, no self-weight.
 * Works in the browser (window.BeamModel) and in Node (require).
 */
(function (root) {
  "use strict";
  var MATERIALS = {
    steel:    { name: "Structural steel",   E: 200e9,   nu: 0.30,  yield: 250e6 },
    aluminum: { name: "Aluminum 6061-T6",   E: 68.9e9,  nu: 0.33,  yield: 276e6 },
    titanium: { name: "Titanium Ti-6Al-4V", E: 113.8e9, nu: 0.342, yield: 880e6 }
  };

  function elementK(EI, l) {
    var c = EI / (l * l * l);
    return [[12 * c, 6 * l * c, -12 * c, 6 * l * c],
            [6 * l * c, 4 * l * l * c, -6 * l * c, 2 * l * l * c],
            [-12 * c, -6 * l * c, 12 * c, -6 * l * c],
            [6 * l * c, 2 * l * l * c, -6 * l * c, 4 * l * l * c]];
  }

  // p: { L, b, h (m), E (Pa), F (N, positive = downward), n (elements) }
  function solve(p) {
    var n = Math.round(p.n);
    if (!(p.L > 0 && p.b > 0 && p.h > 0 && p.E > 0 && n >= 1) || !isFinite(p.F)) return { error: "Length, width, height, modulus and element count must all be greater than zero." };
    var I = p.b * Math.pow(p.h, 3) / 12, EI = p.E * I, l = p.L / n, N = 2 * (n + 1), i, j, e;
    var K = new Float64Array(N * N), f = new Float64Array(N), ke = elementK(EI, l);
    for (e = 0; e < n; e++) for (i = 0; i < 4; i++) for (j = 0; j < 4; j++) K[(2 * e + i) * N + 2 * e + j] += ke[i][j];
    f[N - 2] = -p.F;                                   // point load on the last node's displacement

    // fixed end: v = 0 and θ = 0 at node 0, so rows and columns 0 and 1 are removed
    var m = N - 2, A = new Float64Array(m * m), bvec = new Float64Array(m);
    for (i = 0; i < m; i++) { bvec[i] = f[i + 2]; for (j = 0; j < m; j++) A[i * m + j] = K[(i + 2) * N + j + 2]; }
    var scale = 0; for (i = 0; i < m; i++) scale = Math.max(scale, Math.abs(A[i * m + i]));
    var minPivot = Infinity;
    for (var col = 0; col < m; col++) {
      var piv = col, big = Math.abs(A[col * m + col]);
      for (i = col + 1; i < m; i++) { var a = Math.abs(A[i * m + col]); if (a > big) { big = a; piv = i; } }
      if (!(big > 1e-13 * scale)) return { error: "The stiffness matrix is singular: the beam is not held in place." };
      if (big / scale < minPivot) minPivot = big / scale;
      if (piv !== col) { for (j = 0; j < m; j++) { var t = A[col * m + j]; A[col * m + j] = A[piv * m + j]; A[piv * m + j] = t; } var tb = bvec[col]; bvec[col] = bvec[piv]; bvec[piv] = tb; }
      for (i = col + 1; i < m; i++) {
        var fac = A[i * m + col] / A[col * m + col]; if (fac === 0) continue;
        for (j = col; j < m; j++) A[i * m + j] -= fac * A[col * m + j];
        bvec[i] -= fac * bvec[col];
      }
    }
    var u = new Float64Array(N);
    for (i = m - 1; i >= 0; i--) { var sum = bvec[i]; for (j = i + 1; j < m; j++) sum -= A[i * m + j] * u[j + 2]; u[i + 2] = sum / A[i * m + i]; }

    // residual of the free equations and reactions at the wall, both from the original K
    var r2 = 0, f2 = 0, R = [0, 0];
    for (i = 0; i < N; i++) { var s = 0; for (j = 0; j < N; j++) s += K[i * N + j] * u[j]; if (i < 2) R[i] = s - f[i]; else { r2 += (s - f[i]) * (s - f[i]); f2 += f[i] * f[i]; } }

    var x = [], v = [], th = [], Mn = [];
    for (i = 0; i <= n; i++) { x.push(i * l); v.push(u[2 * i]); th.push(u[2 * i + 1]); }
    // moment at any point, from the element's shape functions
    function at(xq) {
      var el = Math.min(n - 1, Math.max(0, Math.floor(xq / l))), z = Math.min(1, Math.max(0, xq / l - el));
      var d1 = u[2 * el], d2 = u[2 * el + 1], d3 = u[2 * el + 2], d4 = u[2 * el + 3], z2 = z * z, z3 = z2 * z;
      var vv = (1 - 3 * z2 + 2 * z3) * d1 + l * (z - 2 * z2 + z3) * d2 + (3 * z2 - 2 * z3) * d3 + l * (z3 - z2) * d4;
      var tt = ((-6 * z + 6 * z2) * d1 + l * (1 - 4 * z + 3 * z2) * d2 + (6 * z - 6 * z2) * d3 + l * (3 * z2 - 2 * z) * d4) / l;
      var curv = ((-6 + 12 * z) * d1 + l * (-4 + 6 * z) * d2 + (6 - 12 * z) * d3 + l * (-2 + 6 * z) * d4) / (l * l);
      return { v: vv, theta: tt, M: EI * curv, sigmaTop: -p.E * curv * p.h / 2 };
    }
    for (i = 0; i <= n; i++) Mn.push(at(i * l).M);
    var vMax = 0, iMax = 0, Mmax = 0, xM = 0;
    for (i = 0; i <= n; i++) { if (Math.abs(v[i]) > Math.abs(vMax)) { vMax = v[i]; iMax = i; } if (Math.abs(Mn[i]) > Math.abs(Mmax)) { Mmax = Mn[i]; xM = x[i]; } }
    var sigMax = Math.abs(Mmax) * (p.h / 2) / I;
    var delta = p.F * Math.pow(p.L, 3) / (3 * EI);      // analytical tip deflection, positive in the load direction
    return {
      I: I, EI: EI, l: l, n: n, ndof: N, nfree: m, ke: ke, u: u, x: x, v: v, theta: th, M: Mn, at: at,
      vTip: v[n], thetaTip: th[n], vMax: vMax, xVmax: x[iMax], Mmax: Mmax, xMmax: xM, sigmaMax: sigMax,
      deltaExact: delta, errPct: delta !== 0 ? Math.abs(-v[n] - delta) / Math.abs(delta) * 100 : 0,
      exact: function (xq) { return -p.F * xq * xq * (3 * p.L - xq) / (6 * EI); },
      reactionForce: R[0], reactionMoment: R[1], residual: f2 > 0 ? Math.sqrt(r2 / f2) : 0, minPivot: minPivot
    };
  }

  var api = { solve: solve, elementK: elementK, MATERIALS: MATERIALS };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.BeamModel = api;
})(this);
