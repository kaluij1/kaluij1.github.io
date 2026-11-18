/*
 * Viscous flow past a circular cylinder: the numerical model behind the "real air" mode.
 *
 * 2D incompressible Navier-Stokes in dimensionless form (cylinder diameter D = 1, free-stream
 * speed U = 1, density = 1), so the only parameter is the Reynolds number Re = U D / ν.
 *
 * Grid:      uniform staggered (MAC) grid, finite volume.
 * Advection: van Leer flux limiter. Diffusion: central differences. Time: explicit Euler,
 *            step set by the Courant and diffusion limits.
 * Pressure:  projection method; the Poisson equation is solved exactly every step
 *            (sine or cosine transform in y, tridiagonal solve in x).
 * Cylinder:  immersed boundary by direct forcing. Each velocity point carries the fraction of
 *            its cell that lies inside the cylinder, and the velocity there is driven to zero.
 *            The force needed to do that, summed over the cylinder, is the fluid's force on it.
 * Boundaries: left = uniform inflow. Right = open. Top and bottom = open (default) or free-slip walls.
 *
 * Works in the browser (window.WakeModel) and in Node (require).
 */
(function (root) {
  "use strict";

  function create(o) {
    o = o || {};
    var s = { N: o.N || 10, nx: o.nx || 180, ny: (o.ny || 80) & ~1, Re: o.Re || 100, open: o.open !== false,
              xc: o.xc != null ? o.xc : 4.5, CFL: 0.4, DIFF: 0.4, t: 0, steps: 0, dt: 0, cfl: 0, dnum: 0, Fx: 0, Fy: 0, cd: 0, cl: 0 };
    var nx = s.nx, ny = s.ny, h = s.h = 1 / s.N, W = nx + 1, i, j, k;
    s.Lx = nx * h; s.Ly = ny * h; s.yc = o.yc != null ? o.yc : s.Ly / 2 + 0.31 * h;   // a fraction of a cell off centre
    s.u = new Float64Array(W * ny).fill(1); s.us = new Float64Array(W * ny);
    s.v = new Float64Array(nx * (ny + 1)); s.vs = new Float64Array(nx * (ny + 1));
    s.p = new Float64Array(nx * ny); s.rhs = new Float64Array(nx * ny); s.hat = new Float64Array(nx * ny);
    s.tb = new Float64Array(nx); s.td = new Float64Array(nx);
    // which velocity points lie inside the cylinder (more than half their cell, by 5 x 5 sampling)
    s.au = new Float32Array(W * ny); s.av = new Float32Array(nx * (ny + 1));
    function frac(cx, cy) { var n = 0; for (var a = 0; a < 5; a++) for (var b = 0; b < 5; b++) { var x = cx + ((a + 0.5) / 5 - 0.5) * h - s.xc, y = cy + ((b + 0.5) / 5 - 0.5) * h - s.yc; if (x * x + y * y < 0.25) n++; } return o.mask === 'smooth' ? n / 25 : (n >= 13 ? 1 : 0); }
    s.solid = [];
    for (j = 0; j < ny; j++) for (i = 0; i <= nx; i++) { s.au[i + j * W] = frac(i * h, (j + 0.5) * h); }
    for (j = 0; j <= ny; j++) for (i = 0; i < nx; i++) { s.av[i + j * nx] = frac((i + 0.5) * h, j * h); }
    s.i0 = Math.max(1, Math.floor((s.xc - 0.5) / h) - 1); s.i1 = Math.min(nx - 1, Math.ceil((s.xc + 0.5) / h) + 1);
    s.j0 = Math.max(1, Math.floor((s.yc - 0.5) / h) - 1); s.j1 = Math.min(ny - 1, Math.ceil((s.yc + 0.5) / h) + 1);
    // transform basis for the pressure solve
    s.bas = new Float64Array(ny * ny); s.lam = new Float64Array(ny); s.nrm = new Float64Array(ny);
    for (k = 0; k < ny; k++) {
      var kk = s.open ? k + 1 : k;
      s.lam[k] = -(2 - 2 * Math.cos(Math.PI * kk / ny)) / (h * h);
      s.nrm[k] = (s.open ? kk === ny : kk === 0) ? 1 / ny : 2 / ny;
      for (j = 0; j < ny; j++) s.bas[k * ny + j] = s.open ? Math.sin(Math.PI * kk * (j + 0.5) / ny) : Math.cos(Math.PI * kk * (j + 0.5) / ny);
    }
    s.sgn = new Int8Array(ny); s.fs = new Float64Array(nx * (ny >> 1)); s.fd = new Float64Array(nx * (ny >> 1));
    for (k = 0; k < ny; k++) s.sgn[k] = s.bas[k * ny] * s.bas[k * ny + ny - 1] > 0 ? 1 : -1;
    // a small sideways kick behind the cylinder so shedding starts without a long wait
    for (j = 0; j <= ny; j++) for (i = 0; i < nx; i++) { var dx = (i + 0.5) * h - s.xc - 1.2, dy = j * h - s.yc; s.v[i + j * nx] = 0.25 * Math.exp(-(dx * dx + dy * dy) / 1.5); }
    // inflow: the ideal-flow velocity around the cylinder, so the inlet does not squeeze the flow
    s.uin = new Float64Array(ny); s.vin = new Float64Array(ny + 1);
    for (j = 0; j < ny; j++) { var yy = (j + 0.5) * h - s.yc, xx = -s.xc, r2 = xx * xx + yy * yy; s.uin[j] = o.inlet === 'uniform' ? 1 : 1 - 0.25 * (xx * xx - yy * yy) / (r2 * r2); }
    for (j = 0; j <= ny; j++) { var y2 = j * h - s.yc, x2 = -s.xc, q2 = x2 * x2 + y2 * y2; s.vin[j] = o.inlet === 'uniform' ? 0 : -0.5 * x2 * y2 / (q2 * q2); }
    for (j = 0; j < ny; j++) s.u[j * W] = s.uin[j];
    // top and bottom edges in the 'far field' setting: the ideal-flow vertical velocity
    s.vbot = new Float64Array(nx); s.vtop = new Float64Array(nx);
    if (!s.open && o.far !== 'wall') for (i = 0; i < nx; i++) {
      var x3 = (i + 0.5) * h - s.xc, yb = -s.yc, yt = s.Ly - s.yc, rb = x3 * x3 + yb * yb, rt = x3 * x3 + yt * yt;
      s.vbot[i] = -0.5 * x3 * yb / (rb * rb); s.vtop[i] = -0.5 * x3 * yt / (rt * rt);
    }
    s.hist = { t: [], cd: [], cl: [] };
    project(s, 1, true);
    return s;
  }

  function fv(a, b, c) { var d = c - b; if (d === 0) return b; var r = (b - a) / d; return r <= 0 ? b : b + d * r / (1 + r); }

  // solve the pressure equation for the velocity in (us, vs) and write the divergence-free result to (u, v)
  function project(s, dt, fromCurrent) {
    var nx = s.nx, ny = s.ny, W = nx + 1, h = s.h, i, j, k, us = fromCurrent ? s.u : s.us, vs = fromCurrent ? s.v : s.vs, u = s.u, v = s.v;
    var rhs = s.rhs, hat = s.hat, p = s.p, bas = s.bas, tb = s.tb, td = s.td, r2 = 1 / (h * h), sc = 1 / (dt * h), acc, row;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) rhs[i + j * nx] = sc * (us[i + 1 + j * W] - us[i + j * W] + vs[i + (j + 1) * nx] - vs[i + j * nx]);
    // the basis functions are symmetric or antisymmetric about mid-height, so fold each column first: half the work
    var hf = ny >> 1, S = s.fs, Dm = s.fd, jm, src;
    for (j = 0; j < hf; j++) { row = j * nx; jm = (ny - 1 - j) * nx; for (i = 0; i < nx; i++) { S[row + i] = rhs[row + i] + rhs[jm + i]; Dm[row + i] = rhs[row + i] - rhs[jm + i]; } }
    for (k = 0; k < ny; k++) {
      for (i = 0; i < nx; i++) td[i] = 0;
      src = s.sgn[k] > 0 ? S : Dm;
      for (j = 0; j < hf; j++) { acc = bas[k * ny + j] * s.nrm[k]; row = j * nx; for (i = 0; i < nx; i++) td[i] += src[row + i] * acc; }
      var lam = s.lam[k];
      for (i = 0; i < nx; i++) tb[i] = -2 * r2 + lam;
      tb[0] = -r2 + lam; tb[nx - 1] = -3 * r2 + lam;
      for (i = 1; i < nx; i++) { var w = r2 / tb[i - 1]; tb[i] -= w * r2; td[i] -= w * td[i - 1]; }
      hat[nx - 1 + k * nx] = td[nx - 1] / tb[nx - 1];
      for (i = nx - 2; i >= 0; i--) hat[i + k * nx] = (td[i] - r2 * hat[i + 1 + k * nx]) / tb[i];
    }
    for (j = 0; j < hf; j++) {
      row = j * nx; for (i = 0; i < nx; i++) { S[row + i] = 0; Dm[row + i] = 0; }
      for (k = 0; k < ny; k++) { acc = bas[k * ny + j]; var kr = k * nx; src = s.sgn[k] > 0 ? S : Dm; for (i = 0; i < nx; i++) src[row + i] += hat[kr + i] * acc; }
      jm = (ny - 1 - j) * nx; for (i = 0; i < nx; i++) { p[row + i] = S[row + i] + Dm[row + i]; p[jm + i] = S[row + i] - Dm[row + i]; }
    }
    var f = dt / h;
    for (j = 0; j < ny; j++) {
      u[j * W] = s.uin[j];
      for (i = 1; i < nx; i++) u[i + j * W] = us[i + j * W] - f * (p[i + j * nx] - p[i - 1 + j * nx]);
      u[nx + j * W] = us[nx + j * W] + f * 2 * p[nx - 1 + j * nx];
    }
    for (i = 0; i < nx; i++) {
      v[i] = s.open ? vs[i] - f * 2 * p[i] : s.vbot[i];
      for (j = 1; j < ny; j++) v[i + j * nx] = vs[i + j * nx] - f * (p[i + j * nx] - p[i + (j - 1) * nx]);
      v[i + ny * nx] = s.open ? vs[i + ny * nx] + f * 2 * p[i + (ny - 1) * nx] : s.vtop[i];
    }
  }

  function step(s) {
    var nx = s.nx, ny = s.ny, W = nx + 1, h = s.h, i, j, idx, u = s.u, v = s.v, us = s.us, vs = s.vs, nu = 1 / s.Re;
    var um = 1e-9, vm = 1e-9, q;
    for (i = 0; i < u.length; i++) { q = u[i] < 0 ? -u[i] : u[i]; if (q > um) um = q; }
    for (i = 0; i < v.length; i++) { q = v[i] < 0 ? -v[i] : v[i]; if (q > vm) vm = q; }
    var dtc = s.CFL * h / (um + vm), dtd = s.DIFF * h * h / (4 * nu), dt = Math.min(dtc, dtd);
    s.dt = dt; s.cfl = dt * (um + vm) / h; s.dnum = 4 * nu * dt / (h * h); s.limit = dtc < dtd ? "CFL" : "diffusion";
    var a = dt / h, d = nu * dt / (h * h), uc, uE, uW, uN, uS, fe, fw, fn, fs, vel, vc, vN, vS, vE, vW;

    for (j = 0; j < ny; j++) {
      us[j * W] = s.uin[j];
      for (i = 1; i < nx; i++) {
        idx = i + j * W; uc = u[idx]; uE = u[idx + 1]; uW = u[idx - 1];
        vel = 0.5 * (uc + uE); fe = vel >= 0 ? vel * fv(uW, uc, uE) : vel * fv(i + 2 <= nx ? u[idx + 2] : uE, uE, uc);
        vel = 0.5 * (uW + uc); fw = vel >= 0 ? vel * fv(i >= 2 ? u[idx - 2] : uW, uW, uc) : vel * fv(uE, uc, uW);
        vel = 0.5 * (v[i - 1 + (j + 1) * nx] + v[i + (j + 1) * nx]);
        if (j === ny - 1) { fn = vel * uc; uN = uc; }
        else { uN = u[idx + W]; fn = vel >= 0 ? vel * fv(j > 0 ? u[idx - W] : uc, uc, uN) : vel * fv(j + 2 < ny ? u[idx + 2 * W] : uN, uN, uc); }
        vel = 0.5 * (v[i - 1 + j * nx] + v[i + j * nx]);
        if (j === 0) { fs = vel * uc; uS = uc; }
        else { uS = u[idx - W]; fs = vel >= 0 ? vel * fv(j > 1 ? u[idx - 2 * W] : uS, uS, uc) : vel * fv(j + 1 < ny ? u[idx + W] : uc, uc, uS); }
        us[idx] = uc - a * (fe - fw + fn - fs) + d * (uE + uW + uN + uS - 4 * uc);
      }
      idx = nx + j * W; us[idx] = u[idx] - a * Math.max(u[idx], 0) * (u[idx] - u[idx - 1]);
    }
    for (j = 1; j < ny; j++) for (i = 0; i < nx; i++) {
      idx = i + j * nx; vc = v[idx]; vN = v[idx + nx]; vS = v[idx - nx];
      vel = 0.5 * (vc + vN); fn = vel >= 0 ? vel * fv(vS, vc, vN) : vel * fv(j + 2 <= ny ? v[idx + 2 * nx] : vN, vN, vc);
      vel = 0.5 * (vS + vc); fs = vel >= 0 ? vel * fv(j >= 2 ? v[idx - 2 * nx] : vS, vS, vc) : vel * fv(vN, vc, vS);
      vel = 0.5 * (u[i + 1 + (j - 1) * W] + u[i + 1 + j * W]);
      if (i === nx - 1) { fe = vel * vc; vE = vc; }
      else { vE = v[idx + 1]; fe = vel >= 0 ? vel * fv(i > 0 ? v[idx - 1] : vc, vc, vE) : vel * fv(i + 2 < nx ? v[idx + 2] : vE, vE, vc); }
      if (i === 0) { fw = 0.5 * (u[(j - 1) * W] + u[j * W]) * s.vin[j]; vW = 2 * s.vin[j] - vc; }
      else { vW = v[idx - 1]; vel = 0.5 * (u[i + (j - 1) * W] + u[i + j * W]); fw = vel >= 0 ? vel * fv(i > 1 ? v[idx - 2] : vW, vW, vc) : vel * fv(i + 1 < nx ? v[idx + 1] : vc, vc, vW); }
      vs[idx] = vc - a * (fe - fw + fn - fs) + d * (vE + vW + vN + vS - 4 * vc);
    }
    for (i = 0; i < nx; i++) { vs[i] = s.open ? vs[i + nx] : s.vbot[i]; vs[i + ny * nx] = s.open ? vs[i + (ny - 1) * nx] : s.vtop[i]; }

    // immersed boundary: stop the fluid inside the cylinder and add up the force it took
    var Fx = 0, Fy = 0, al;
    for (j = s.j0 - 1; j <= s.j1; j++) for (i = s.i0; i <= s.i1 + 1; i++) { al = s.au[i + j * W]; if (al > 0) { Fx += al * us[i + j * W]; us[i + j * W] *= 1 - al; } }
    for (j = s.j0; j <= s.j1 + 1; j++) for (i = s.i0 - 1; i <= s.i1; i++) { al = s.av[i + j * nx]; if (al > 0) { Fy += al * vs[i + j * nx]; vs[i + j * nx] *= 1 - al; } }
    s.Fx = Fx * h * h / dt; s.Fy = Fy * h * h / dt;          // force per unit span, in units of ρ U² D
    project(s, dt, false);
    s.t += dt; s.steps++;
    s.cd = 2 * s.Fx; s.cl = 2 * s.Fy;
    if (s.steps % 2 === 0) { var H = s.hist; H.t.push(s.t); H.cd.push(s.cd); H.cl.push(s.cl); if (H.t.length > 6000) { H.t.splice(0, 1000); H.cd.splice(0, 1000); H.cl.splice(0, 1000); } }
  }

  // Strouhal number, mean drag and lift amplitude from the last few shedding cycles
  function stats(s, window) {
    var H = s.hist, n = H.t.length, t0 = s.t - (window || 30), i, k0 = 0;
    while (k0 < n && H.t[k0] < t0) k0++;
    if (n - k0 < 20) return { St: 0, cd: s.cd, clAmp: 0, cycles: 0 };
    var mean = 0, cdm = 0, lo = Infinity, hi = -Infinity, cdlo = Infinity, cdhi = -Infinity;
    for (i = k0; i < n; i++) { mean += H.cl[i]; cdm += H.cd[i]; }
    mean /= n - k0; cdm /= n - k0;
    var cross = [];
    for (i = k0 + 1; i < n; i++) {
      var a = H.cl[i - 1] - mean, b = H.cl[i] - mean;
      if (a < 0 && b >= 0) cross.push(H.t[i - 1] + (H.t[i] - H.t[i - 1]) * (-a) / (b - a));
    }
    var St = 0, cyc = cross.length - 1;
    if (cyc >= 2) {
      St = cyc / (cross[cyc] - cross[0]);
      // average over whole cycles only
      var ta = cross[0], tb = cross[cyc], m2 = 0, c2 = 0, cnt = 0;
      for (i = k0; i < n; i++) if (H.t[i] >= ta && H.t[i] <= tb) { m2 += H.cl[i]; c2 += H.cd[i]; cnt++; if (H.cl[i] < lo) lo = H.cl[i]; if (H.cl[i] > hi) hi = H.cl[i]; if (H.cd[i] < cdlo) cdlo = H.cd[i]; if (H.cd[i] > cdhi) cdhi = H.cd[i]; }
      cdm = c2 / cnt; mean = m2 / cnt;
    }
    var amp = isFinite(hi) ? (hi - lo) / 2 : 0;
    return { St: amp > 0.01 ? St : 0, cd: cdm, clAmp: amp, clMean: mean, cdAmp: isFinite(cdhi) ? (cdhi - cdlo) / 2 : 0, cycles: Math.max(0, cyc) };
  }
  // Williamson's fit to experiments for the shedding frequency, valid for about 49 < Re < 180
  function strouhalRef(Re) { return Re < 49 ? 0 : -3.3265 / Re + 0.1816 + 1.6e-4 * Re; }
  function maxDivergence(s) {
    var nx = s.nx, ny = s.ny, W = nx + 1, m = 0;
    for (var j = 0; j < ny; j++) for (var i = 0; i < nx; i++) { var dv = Math.abs(s.u[i + 1 + j * W] - s.u[i + j * W] + s.v[i + (j + 1) * nx] - s.v[i + j * nx]) / s.h; if (dv > m) m = dv; }
    return m;
  }
  var api = { create: create, step: step, stats: stats, strouhalRef: strouhalRef, maxDivergence: maxDivergence };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.WakeModel = api;
})(this);
