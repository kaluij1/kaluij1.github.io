/*
 * Forced convection over a heated aluminum plate: the numerical model.
 *
 * Air:    2D incompressible Navier-Stokes + advection-diffusion of temperature, constant properties.
 *         Staggered (MAC) grid, finite volume. Advection: van Leer flux limiter (second order, bounded).
 *         Diffusion: central differences. Time: explicit Euler, step limited by CFL and diffusion numbers.
 *         Pressure: projection method. The Poisson equation is solved exactly each step
 *         (cosine transform in y, tridiagonal solve in x), so the velocity is divergence free to round-off.
 * Plate:  thin aluminum plate with a uniform heat flux on its underside (Biot number << 1, so one
 *         temperature through the thickness). 1D conduction along the plate, implicit in time,
 *         coupled to the air through the wall heat flux.
 *
 * Boundaries: left = uniform inflow at T_in. Right and top = open (p = 0, zero gradient).
 *             Bottom = no-slip on the plate; free-slip and insulated ahead of and behind it,
 *             so the plate acts like a thin plate in a free stream (upper half shown).
 * Not modeled: buoyancy, radiation, turbulence, property variation with temperature, the third dimension.
 *
 * Works in the browser (window.ConvectionModel) and in Node (require).
 */
(function (root) {
  "use strict";
  var P_ATM = 101325, R_AIR = 287.05, CP_AIR = 1006;

  // Dry air at 1 atm: ideal gas density, Sutherland viscosity and conductivity.
  function airProps(Tc) {
    var T = Tc + 273.15;
    var rho = P_ATM / (R_AIR * T);
    var mu = 1.716e-5 * Math.pow(T / 273.15, 1.5) * (273.15 + 110.4) / (T + 110.4);
    var k = 0.0241 * Math.pow(T / 273.15, 1.5) * (273.15 + 194) / (T + 194);
    return { rho: rho, mu: mu, k: k, cp: CP_AIR, nu: mu / rho, alpha: k / (rho * CP_AIR), Pr: mu * CP_AIR / k };
  }
  var PLATE = { name: "Aluminum 6061-T6", k: 167, rho: 2700, cp: 896, t: 0.003 };

  function create(o) {
    o = o || {};
    var s = {
      nx: o.nx || 100, ny: o.ny || 48, L: o.L || 0.20, H: o.H || 0.03,
      xLead: o.xLead != null ? o.xLead : 0.03, plateLen: o.plateLen || 0.12, depth: 0.10,
      U: o.U != null ? o.U : 0.5, Tin: o.Tin != null ? o.Tin : 20, q: o.q != null ? o.q : 500,
      accel: o.accel || 1,                 // divides the plate's heat capacity (1 = real aluminum)
      plate: { k: o.plateK != null ? o.plateK : PLATE.k, rho: PLATE.rho, cp: PLATE.cp, t: PLATE.t },
      CFL: 0.4, DIFF: 0.4, t: 0, steps: 0, dt: 0, cfl: 0, dnum: 0, div: 0
    };
    s.dx = s.L / s.nx; s.dy = s.H / s.ny;
    s.i0 = Math.round(s.xLead / s.dx); s.i1 = s.i0 + Math.round(s.plateLen / s.dx); s.np = s.i1 - s.i0;
    var nx = s.nx, ny = s.ny;
    s.u = new Float64Array((nx + 1) * ny).fill(s.U); s.us = new Float64Array((nx + 1) * ny);
    s.v = new Float64Array(nx * (ny + 1)); s.vs = new Float64Array(nx * (ny + 1));
    s.T = new Float64Array(nx * ny).fill(s.Tin); s.Tn = new Float64Array(nx * ny);
    s.p = new Float64Array(nx * ny); s.rhs = new Float64Array(nx * ny); s.hat = new Float64Array(nx * ny);
    s.Fx = new Float64Array((nx + 1) * ny); s.Fy = new Float64Array(nx * (ny + 1));
    s.Tw = new Float64Array(s.np).fill(s.Tin); s.qw = new Float64Array(s.np);
    s.ta = new Float64Array(Math.max(nx, s.np)); s.tb = new Float64Array(Math.max(nx, s.np)); s.tc = new Float64Array(Math.max(nx, s.np)); s.td = new Float64Array(Math.max(nx, s.np));
    // cosine basis for the pressure solve: Neumann at the bottom wall, p = 0 at the open top
    s.cos = new Float64Array(ny * ny); s.lam = new Float64Array(ny);
    for (var k = 0; k < ny; k++) {
      s.lam[k] = -(2 - 2 * Math.cos(Math.PI * (k + 0.5) / ny)) / (s.dy * s.dy);
      for (var j = 0; j < ny; j++) s.cos[k * ny + j] = Math.cos(Math.PI * (k + 0.5) * (j + 0.5) / ny);
    }
    s.air = airProps(s.Tin);
    return s;
  }

  // van Leer limited face value: a = far upwind, b = upwind, c = downwind
  function fv(a, b, c) {
    var d = c - b;
    if (d === 0) return b;
    var r = (b - a) / d;
    return r <= 0 ? b : b + d * r / (1 + r);
  }

  function step(s) {
    var nx = s.nx, ny = s.ny, W = nx + 1, dx = s.dx, dy = s.dy, i, j, k, c;
    var u = s.u, v = s.v, T = s.T, a = s.air = airProps(s.Tin);
    var nu = a.nu, al = a.alpha, rho = a.rho;

    /* ---- time step from the stability limits ---- */
    var um = 1e-9, vm = 1e-9;
    for (i = 0; i < u.length; i++) { var q1 = u[i] < 0 ? -u[i] : u[i]; if (q1 > um) um = q1; }
    for (i = 0; i < v.length; i++) { var q2 = v[i] < 0 ? -v[i] : v[i]; if (q2 > vm) vm = q2; }
    var amax = Math.max(nu, al), g = 1 / (dx * dx) + 1 / (dy * dy);
    var dtc = s.CFL / (um / dx + vm / dy), dtd = s.DIFF / (2 * amax * g);
    var dt = Math.min(dtc, dtd);
    s.dt = dt; s.cfl = dt * (um / dx + vm / dy); s.dnum = 2 * amax * g * dt; s.limit = dtc < dtd ? "CFL" : "diffusion";

    /* ---- plate: implicit 1D conduction, heated from below, cooled by the air above ---- */
    var np = s.np, i0 = s.i0, hc = a.k / (dy / 2);
    var C = s.plate.rho * s.plate.cp * s.plate.t / s.accel, kap = s.plate.k * s.plate.t / (dx * dx);
    var A = s.ta, B = s.tb, Cc = s.tc, D = s.td;
    for (i = 0; i < np; i++) {
      var nb = (i > 0 ? 1 : 0) + (i < np - 1 ? 1 : 0);
      A[i] = i > 0 ? -kap : 0; Cc[i] = i < np - 1 ? -kap : 0;
      B[i] = C / dt + hc + nb * kap;
      D[i] = C / dt * s.Tw[i] + s.q + hc * T[i0 + i];
    }
    for (i = 1; i < np; i++) { var m = A[i] / B[i - 1]; B[i] -= m * Cc[i - 1]; D[i] -= m * D[i - 1]; }
    s.Tw[np - 1] = D[np - 1] / B[np - 1];
    for (i = np - 2; i >= 0; i--) s.Tw[i] = (D[i] - Cc[i] * s.Tw[i + 1]) / B[i];
    for (i = 0; i < np; i++) s.qw[i] = hc * (s.Tw[i] - T[i0 + i]);      // W/m2 from plate to air

    /* ---- air temperature: limited advection + diffusion ---- */
    var Fx = s.Fx, Fy = s.Fy, Tin = s.Tin, uf, vf;
    for (j = 0; j < ny; j++) {
      Fx[j * W] = u[j * W] * Tin;
      for (i = 1; i < nx; i++) {
        uf = u[i + j * W]; c = i + j * nx;
        Fx[i + j * W] = uf >= 0 ? uf * fv(i > 1 ? T[c - 2] : T[c - 1], T[c - 1], T[c]) : uf * fv(i < nx - 1 ? T[c + 1] : T[c], T[c], T[c - 1]);
      }
      uf = u[nx + j * W]; Fx[nx + j * W] = uf >= 0 ? uf * T[nx - 1 + j * nx] : uf * Tin;
    }
    for (i = 0; i < nx; i++) {
      Fy[i] = 0;
      for (j = 1; j < ny; j++) {
        vf = v[i + j * nx]; c = i + j * nx;
        Fy[c] = vf >= 0 ? vf * fv(j > 1 ? T[c - 2 * nx] : T[c - nx], T[c - nx], T[c]) : vf * fv(j < ny - 1 ? T[c + nx] : T[c], T[c], T[c - nx]);
      }
      vf = v[i + ny * nx]; Fy[i + ny * nx] = vf >= 0 ? vf * T[i + (ny - 1) * nx] : vf * Tin;
    }
    var Tn = s.Tn, rdx2 = 1 / (dx * dx), rdy2 = 1 / (dy * dy), rcp = rho * a.cp;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++) {
      c = i + j * nx;
      var Tc = T[c];
      var dif = al * (((i < nx - 1 ? T[c + 1] : Tc) - 2 * Tc + (i > 0 ? T[c - 1] : Tin * 2 - Tc)) * rdx2
                    + ((j < ny - 1 ? T[c + nx] : Tc) - Tc) * rdy2 + (j > 0 ? (T[c - nx] - Tc) * rdy2 : 0));
      if (j === 0 && i >= i0 && i < s.i1) dif += s.qw[i - i0] / (rcp * dy);
      Tn[c] = Tc + dt * (-(Fx[i + 1 + j * W] - Fx[i + j * W]) / dx - (Fy[c + nx] - Fy[c]) / dy + dif);
    }

    /* ---- momentum predictor ---- */
    var us = s.us, vs = s.vs, idx, uc, uE, uW, fe, fw, fn, fs, vel;
    for (j = 0; j < ny; j++) {
      us[j * W] = s.U;
      for (i = 1; i < nx; i++) {
        idx = i + j * W; uc = u[idx]; uE = u[idx + 1]; uW = u[idx - 1];
        vel = 0.5 * (uc + uE); fe = vel >= 0 ? vel * fv(uW, uc, uE) : vel * fv(i + 2 <= nx ? u[idx + 2] : uE, uE, uc);
        vel = 0.5 * (uW + uc); fw = vel >= 0 ? vel * fv(i >= 2 ? u[idx - 2] : uW, uW, uc) : vel * fv(uE, uc, uW);
        vel = 0.5 * (v[i - 1 + (j + 1) * nx] + v[i + (j + 1) * nx]);
        if (j === ny - 1) fn = vel * uc;
        else { var uN = u[idx + W]; fn = vel >= 0 ? vel * fv(j > 0 ? u[idx - W] : uc, uc, uN) : vel * fv(j + 2 < ny ? u[idx + 2 * W] : uN, uN, uc); }
        if (j === 0) fs = 0;
        else { var uS = u[idx - W]; vel = 0.5 * (v[i - 1 + j * nx] + v[i + j * nx]); fs = vel >= 0 ? vel * fv(j > 1 ? u[idx - 2 * W] : uS, uS, uc) : vel * fv(j + 1 < ny ? u[idx + W] : uc, uc, uS); }
        var uNn = j < ny - 1 ? u[idx + W] : uc, uSs = j > 0 ? u[idx - W] : (i >= i0 && i <= s.i1 ? -uc : uc);
        us[idx] = uc + dt * (-(fe - fw) / dx - (fn - fs) / dy + nu * ((uE - 2 * uc + uW) * rdx2 + (uNn - 2 * uc + uSs) * rdy2));
      }
      idx = nx + j * W; us[idx] = u[idx] - dt * Math.max(u[idx], 0) * (u[idx] - u[idx - 1]) / dx;
    }
    var vc, vN, vS;
    for (i = 0; i < nx; i++) vs[i] = 0;
    for (j = 1; j < ny; j++) for (i = 0; i < nx; i++) {
      idx = i + j * nx; vc = v[idx]; vN = v[idx + nx]; vS = v[idx - nx];
      vel = 0.5 * (vc + vN); fn = vel >= 0 ? vel * fv(vS, vc, vN) : vel * fv(j + 2 <= ny ? v[idx + 2 * nx] : vN, vN, vc);
      vel = 0.5 * (vS + vc); fs = vel >= 0 ? vel * fv(j >= 2 ? v[idx - 2 * nx] : vS, vS, vc) : vel * fv(vN, vc, vS);
      vel = 0.5 * (u[i + 1 + (j - 1) * W] + u[i + 1 + j * W]);
      if (i === nx - 1) fe = vel * vc;
      else { var vE = v[idx + 1]; fe = vel >= 0 ? vel * fv(i > 0 ? v[idx - 1] : vc, vc, vE) : vel * fv(i + 2 < nx ? v[idx + 2] : vE, vE, vc); }
      if (i === 0) fw = 0;
      else { var vW = v[idx - 1]; vel = 0.5 * (u[i + (j - 1) * W] + u[i + j * W]); fw = vel >= 0 ? vel * fv(i > 1 ? v[idx - 2] : vW, vW, vc) : vel * fv(i + 1 < nx ? v[idx + 1] : vc, vc, vW); }
      var vEe = i < nx - 1 ? v[idx + 1] : vc, vWw = i > 0 ? v[idx - 1] : -vc;
      vs[idx] = vc + dt * (-(fe - fw) / dx - (fn - fs) / dy + nu * ((vEe - 2 * vc + vWw) * rdx2 + (vN - 2 * vc + vS) * rdy2));
    }
    for (i = 0; i < nx; i++) vs[i + ny * nx] = vs[i + (ny - 1) * nx];

    /* ---- pressure: solve the Poisson equation exactly, then project ---- */
    var rhs = s.rhs, hat = s.hat, p = s.p, cs = s.cos, sc = rho / dt;
    for (j = 0; j < ny; j++) for (i = 0; i < nx; i++)
      rhs[i + j * nx] = sc * ((us[i + 1 + j * W] - us[i + j * W]) / dx + (vs[i + (j + 1) * nx] - vs[i + j * nx]) / dy);
    var acc, row, tb = s.tb, td = s.td;
    for (k = 0; k < ny; k++) {
      for (i = 0; i < nx; i++) td[i] = 0;
      for (j = 0; j < ny; j++) { acc = cs[k * ny + j] * 2 / ny; row = j * nx; for (i = 0; i < nx; i++) td[i] += rhs[row + i] * acc; }
      var lam = s.lam[k];
      for (i = 0; i < nx; i++) tb[i] = -2 * rdx2 + lam;
      tb[0] = -rdx2 + lam; tb[nx - 1] = -3 * rdx2 + lam;
      for (i = 1; i < nx; i++) { var w = rdx2 / tb[i - 1]; tb[i] -= w * rdx2; td[i] -= w * td[i - 1]; }
      hat[nx - 1 + k * nx] = td[nx - 1] / tb[nx - 1];
      for (i = nx - 2; i >= 0; i--) hat[i + k * nx] = (td[i] - rdx2 * hat[i + 1 + k * nx]) / tb[i];
    }
    for (j = 0; j < ny; j++) {
      row = j * nx; for (i = 0; i < nx; i++) p[row + i] = 0;
      for (k = 0; k < ny; k++) { acc = cs[k * ny + j]; var kr = k * nx; for (i = 0; i < nx; i++) p[row + i] += hat[kr + i] * acc; }
    }
    var f = dt / rho;
    for (j = 0; j < ny; j++) {
      u[j * W] = s.U;
      for (i = 1; i < nx; i++) u[i + j * W] = us[i + j * W] - f * (p[i + j * nx] - p[i - 1 + j * nx]) / dx;
      u[nx + j * W] = us[nx + j * W] + f * 2 * p[nx - 1 + j * nx] / dx;
    }
    for (i = 0; i < nx; i++) {
      v[i] = 0;
      for (j = 1; j < ny; j++) v[i + j * nx] = vs[i + j * nx] - f * (p[i + j * nx] - p[i + (j - 1) * nx]) / dy;
      v[i + ny * nx] = vs[i + ny * nx] + f * 2 * p[i + (ny - 1) * nx] / dy;
    }
    s.T = Tn; s.Tn = T;
    s.t += dt; s.steps++;
  }

  function maxDivergence(s) {
    var nx = s.nx, ny = s.ny, W = nx + 1, m = 0;
    for (var j = 0; j < ny; j++) for (var i = 0; i < nx; i++) {
      var d = Math.abs((s.u[i + 1 + j * W] - s.u[i + j * W]) / s.dx + (s.v[i + (j + 1) * nx] - s.v[i + j * nx]) / s.dy);
      if (d > m) m = d;
    }
    return m;
  }

  // Everything the page reports, computed from the fields.
  function diagnostics(s) {
    var nx = s.nx, ny = s.ny, W = nx + 1, a = s.air, i, j, np = s.np, Lp = np * s.dx;
    var tmax = -1e9, tmin = 1e9, tsum = 0, qsum = 0;
    for (i = 0; i < np; i++) { if (s.Tw[i] > tmax) tmax = s.Tw[i]; if (s.Tw[i] < tmin) tmin = s.Tw[i]; tsum += s.Tw[i]; qsum += s.qw[i]; }
    var tmean = tsum / np;
    var Qin = s.q * Lp * s.depth, Qconv = qsum * s.dx * s.depth;
    // outlet: flow weighted (mixing cup) temperature, and the hottest air leaving
    var mt = 0, mf = 0, tpk = -1e9, Hout = 0;
    for (j = 0; j < ny; j++) { var uo = s.u[nx + j * W], To = s.T[nx - 1 + j * nx]; mt += uo * To; mf += uo; if (To > tpk) tpk = To; Hout += uo * (To - s.Tin) * s.dy; }
    for (i = 0; i < nx; i++) { var vt = s.v[i + ny * nx]; Hout += vt * ((vt >= 0 ? s.T[i + (ny - 1) * nx] : s.Tin) - s.Tin) * s.dx; }
    var Qair = a.rho * a.cp * Hout * s.depth;
    var hbar = qsum / np / Math.max(tmean - s.Tin, 1e-9);
    var ReL = s.U * Lp / a.nu, pr3 = Math.cbrt(a.Pr);
    return {
      TwMax: tmax, TwMin: tmin, TwMean: tmean, ToutMix: mt / mf, ToutPeak: tpk,
      Qin: Qin, Qconv: Qconv, Qair: Qair, hbar: hbar, NuL: hbar * Lp / a.k,
      ReL: ReL, Pr: a.Pr, Lp: Lp,
      hUHF: 0.6795 * Math.sqrt(ReL) * pr3 * a.k / Lp,     // laminar, uniform flux, from mean wall temperature
      hUWT: 0.664 * Math.sqrt(ReL) * pr3 * a.k / Lp,      // laminar, uniform wall temperature
      Ri: 9.81 / (s.Tin + 273.15) * Math.max(tmean - s.Tin, 0) * Lp / (s.U * s.U),
      Bi: hbar * s.plate.t / s.plate.k,
      tauPlate: s.plate.rho * s.plate.cp * s.plate.t / s.accel / Math.max(hbar, 1e-9)
    };
  }
  // local laminar correlations at distance x from the leading edge (W/m2K)
  function hLocal(s, x, coef) { var a = s.air; return coef * Math.sqrt(s.U * x / a.nu) * Math.cbrt(a.Pr) * a.k / x; }

  var api = { create: create, step: step, diagnostics: diagnostics, airProps: airProps, hLocal: hLocal, maxDivergence: maxDivergence, PLATE: PLATE };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  else root.ConvectionModel = api;
})(this);
