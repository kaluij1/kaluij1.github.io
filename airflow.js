/*
 * Air flow past a cylinder, for the home page hero.
 *
 * The model is ideal (potential) flow: a steady stream moving left to right,
 * plus a doublet that makes the air part around a cylinder. The cylinder
 * follows the cursor. When the cylinder moves, the air it sees is the stream
 * minus its own velocity, so a fast sweep pushes air ahead and pulls it in behind.
 *
 *   velocity = U + (a^2 / r^2) * [ W - 2 (W . r^) r^ ]      W = U - V_cylinder
 *
 * a is the cylinder radius, r the vector from its centre. Ideal flow has no
 * friction, so there is no wake, and the fastest air (2x the stream) is at the
 * top and bottom of the cylinder.
 *
 * Tweak the numbers in CONFIG to change the feel.
 */
(function () {
  "use strict";

  var CONFIG = {
    stream: 52,        // px/s, speed of the undisturbed air (left to right)
    radius: 46,        // px, cylinder radius
    density: 1 / 1500, // tracer particles per square px
    maxParticles: 950, // hard cap for performance
    trailPoints: 16,   // points kept per tracer
    trailStep: 0.035,  // s between trail points (trail lasts trailPoints * trailStep)
    life: [6, 14],     // s, each tracer respawns somewhere new after this long
    lineWidth: 1.1,    // px
    alpha: 0.5,        // opacity of the strongest streaks
    follow: 9,         // 1/s, how tightly the cylinder follows the cursor
    maxCylSpeed: 260,  // px/s, cap on the cylinder speed fed to the model
    home: [0.72, 0.5]  // where the cylinder rests, as a fraction of the hero (wide screens)
  };

  var canvas = document.getElementById("flow");
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext("2d");
  var hero = canvas.parentElement;
  var motionQuery = window.matchMedia("(prefers-reduced-motion: reduce)");

  var W = 0, H = 0, dpr = 1, rgb = "70,96,127";
  var P = [];                                  // tracers
  var cyl = { x: 0, y: 0, vx: 0, vy: 0 };      // cylinder centre and velocity
  var target = { x: 0, y: 0, held: false };    // where the cylinder is heading
  var raf = 0, running = false, inView = true, last = 0, trailClock = 0;
  var NB = 4;                                  // opacity buckets

  function readColor() {
    var v = getComputedStyle(document.documentElement).getPropertyValue("--field-rgb").trim();
    if (v) rgb = v.replace(/\s+/g, ",");
  }
  function homePos() {
    return W > 900 ? { x: W * CONFIG.home[0], y: H * CONFIG.home[1] } : { x: W * 0.5, y: H - CONFIG.radius - 104 };
  }

  /* The flow model: air velocity at (x, y), in px/s. */
  function velocity(x, y, out) {
    var U = CONFIG.stream, a = CONFIG.radius;
    var rx = x - cyl.x, ry = y - cyl.y, r2 = rx * rx + ry * ry;
    if (r2 < a * a) r2 = a * a;
    var r = Math.sqrt(r2), nx = rx / r, ny = ry / r;
    var wx = U - cyl.vx, wy = -cyl.vy;          // stream as the cylinder sees it
    var k = (a * a) / r2, dot = wx * nx + wy * ny;
    out.x = U + k * (wx - 2 * dot * nx);
    out.y = k * (wy - 2 * dot * ny);
    return out;
  }

  function spawn(p, anywhere) {
    p.x = anywhere ? Math.random() * W : -8;
    p.y = Math.random() * H;
    var dx = p.x - cyl.x, dy = p.y - cyl.y, a = CONFIG.radius + 2;
    if (dx * dx + dy * dy < a * a) p.y = cyl.y + (dy < 0 ? -a : a);
    p.age = 0;
    p.life = CONFIG.life[0] + Math.random() * (CONFIG.life[1] - CONFIG.life[0]);
    p.n = 0; p.speed = CONFIG.stream;
  }

  function build() {
    var rect = hero.getBoundingClientRect();
    W = rect.width; H = rect.height;
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(W * dpr);
    canvas.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    var h = homePos();
    cyl.x = h.x; cyl.y = h.y; cyl.vx = cyl.vy = 0;
    if (!target.held) { target.x = h.x; target.y = h.y; }
    var n = Math.min(CONFIG.maxParticles, Math.round(W * H * CONFIG.density));
    P = [];
    for (var i = 0; i < n; i++) {
      var p = { tx: new Float32Array(CONFIG.trailPoints), ty: new Float32Array(CONFIG.trailPoints) };
      spawn(p, true);
      p.age = Math.random() * p.life;
      P.push(p);
    }
  }

  var v = { x: 0, y: 0 };
  function step(dt) {
    // the cylinder eases toward the cursor; its velocity feeds the model
    var k = 1 - Math.exp(-CONFIG.follow * dt);
    var nx = cyl.x + (target.x - cyl.x) * k, ny = cyl.y + (target.y - cyl.y) * k;
    var vx = (nx - cyl.x) / dt, vy = (ny - cyl.y) / dt;
    var sp = Math.sqrt(vx * vx + vy * vy), cap = CONFIG.maxCylSpeed;
    if (sp > cap) { vx *= cap / sp; vy *= cap / sp; }
    cyl.vx += (vx - cyl.vx) * 0.18; cyl.vy += (vy - cyl.vy) * 0.18;
    cyl.x = nx; cyl.y = ny;

    trailClock += dt;
    var record = trailClock >= CONFIG.trailStep;
    if (record) trailClock = 0;
    var a = CONFIG.radius, N = CONFIG.trailPoints;

    for (var i = 0; i < P.length; i++) {
      var p = P[i];
      p.age += dt;
      if (p.age > p.life) { spawn(p, true); continue; }
      // midpoint (RK2) step along the flow
      velocity(p.x, p.y, v);
      var mx = p.x + v.x * dt * 0.5, my = p.y + v.y * dt * 0.5;
      velocity(mx, my, v);
      p.x += v.x * dt; p.y += v.y * dt;
      p.speed = Math.sqrt(v.x * v.x + v.y * v.y);
      // air cannot pass through the cylinder
      var rx = p.x - cyl.x, ry = p.y - cyl.y, r2 = rx * rx + ry * ry;
      if (r2 < a * a) { var r = Math.sqrt(r2) || 1; p.x = cyl.x + rx / r * a; p.y = cyl.y + ry / r * a; }
      if (p.x > W + 12 || p.x < -40 || p.y < -40 || p.y > H + 40) { spawn(p, p.x <= W + 12); continue; }
      if (record) {
        if (p.n < N) { p.tx[p.n] = p.x; p.ty[p.n] = p.y; p.n++; }
        else { p.tx.copyWithin(0, 1); p.ty.copyWithin(0, 1); p.tx[N - 1] = p.x; p.ty[N - 1] = p.y; }
      }
    }
  }

  function drawCylinder() {
    ctx.beginPath();
    ctx.arc(cyl.x, cyl.y, CONFIG.radius, 0, 6.2832);
    ctx.strokeStyle = "rgba(" + rgb + ",.6)";
    ctx.lineWidth = 1.2;
    ctx.stroke();
  }

  function draw() {
    ctx.clearRect(0, 0, W, H);
    ctx.lineWidth = CONFIG.lineWidth;
    ctx.lineCap = "round"; ctx.lineJoin = "round";
    var U = CONFIG.stream, b, i;
    // sort tracers into a few opacity levels so each level is one stroke
    for (i = 0; i < P.length; i++) {
      var p = P[i];
      var fade = Math.min(1, p.age / 0.7, (p.life - p.age) / 0.7);
      var e = fade * (0.45 + 0.55 * Math.min(1, p.speed / (2 * U)));
      p.b = Math.max(0, Math.min(NB - 1, Math.round(e * NB - 0.5)));
    }
    for (b = 0; b < NB; b++) {
      ctx.beginPath();
      for (i = 0; i < P.length; i++) {
        var q = P[i];
        if (q.b !== b || q.n < 2) continue;
        ctx.moveTo(q.tx[0], q.ty[0]);
        for (var j = 1; j < q.n; j++) ctx.lineTo(q.tx[j], q.ty[j]);
        ctx.lineTo(q.x, q.y);
      }
      ctx.strokeStyle = "rgba(" + rgb + "," + (CONFIG.alpha * (b + 1) / NB).toFixed(3) + ")";
      ctx.stroke();
    }
    drawCylinder();
  }

  /* Still picture for people who ask for reduced motion: streamlines past a resting cylinder. */
  function drawStatic() {
    var h = homePos();
    cyl.x = h.x; cyl.y = h.y; cyl.vx = cyl.vy = 0;
    ctx.clearRect(0, 0, W, H);
    ctx.lineWidth = CONFIG.lineWidth;
    ctx.strokeStyle = "rgba(" + rgb + "," + (CONFIG.alpha * 0.7).toFixed(3) + ")";
    ctx.beginPath();
    for (var y0 = 13; y0 < H; y0 += 26) {
      var x = 0, y = y0 + 0.01, guard = 0;
      ctx.moveTo(x, y);
      while (x < W && guard++ < 4000) {
        velocity(x, y, v);
        var s = Math.sqrt(v.x * v.x + v.y * v.y) || 1;
        x += v.x / s * 4; y += v.y / s * 4;
        ctx.lineTo(x, y);
      }
    }
    ctx.stroke();
    drawCylinder();
  }

  function warmUp() { for (var i = 0; i < 90; i++) step(1 / 30); }

  function frame(now) {
    var dt = Math.min(0.05, (now - last) / 1000 || 0.016);
    last = now;
    step(dt);
    draw();
    raf = window.requestAnimationFrame(frame);
  }
  function start() {
    if (running || motionQuery.matches || !inView || document.hidden) return;
    running = true; last = performance.now();
    raf = window.requestAnimationFrame(frame);
  }
  function stop() { running = false; window.cancelAnimationFrame(raf); }

  /* Pointer: the cylinder goes to the cursor while it is over the hero, and home otherwise. */
  function onMove(e) {
    var rect = hero.getBoundingClientRect();
    var x = e.clientX - rect.left, y = e.clientY - rect.top;
    if (x >= 0 && x <= rect.width && y >= 0 && y <= rect.height) { target.x = x; target.y = y; target.held = true; }
    else goHome();
  }
  function goHome() { var h = homePos(); target.x = h.x; target.y = h.y; target.held = false; }

  var resizeTimer = 0;
  function onResize() {
    window.clearTimeout(resizeTimer);
    resizeTimer = window.setTimeout(function () {
      build();
      if (motionQuery.matches) drawStatic(); else { warmUp(); draw(); }
    }, 120);
  }
  function onVisibility() { if (document.hidden) stop(); else start(); }
  function onMotionChange() { if (motionQuery.matches) { stop(); drawStatic(); } else { warmUp(); start(); } }

  readColor();
  build();
  if (motionQuery.matches) drawStatic(); else { warmUp(); draw(); }

  window.addEventListener("pointermove", onMove, { passive: true });
  document.documentElement.addEventListener("pointerleave", goHome);
  window.addEventListener("pointerup", function (e) { if (e.pointerType === "touch") goHome(); });
  document.addEventListener("visibilitychange", onVisibility);
  if (motionQuery.addEventListener) motionQuery.addEventListener("change", onMotionChange);

  var ro = null, io = null;
  if ("ResizeObserver" in window) { ro = new ResizeObserver(onResize); ro.observe(hero); }
  else window.addEventListener("resize", onResize);
  if ("IntersectionObserver" in window) {
    io = new IntersectionObserver(function (entries) {
      inView = entries[0].isIntersecting;
      if (inView) start(); else stop();
    });
    io.observe(hero);
  }
  start();

  // window.airflow.config holds the live settings; change a value in the console to try it.
  window.airflow = {
    config: CONFIG,
    rebuild: function () { build(); warmUp(); },
    destroy: function () {
      stop();
      window.removeEventListener("pointermove", onMove);
      document.documentElement.removeEventListener("pointerleave", goHome);
      document.removeEventListener("visibilitychange", onVisibility);
      if (ro) ro.disconnect();
      if (io) io.disconnect();
      ctx.clearRect(0, 0, W, H);
    }
  };
})();
