// Òrbites en un potencial central: separació radial/angular amb el moment
// angular, potencial efectiu i petites oscil·lacions entorn de l'òrbita
// circular (Goldstein, cap. 3). Unitats: m = k = 1.
(function () {
  'use strict';

  const m = 1, k = 1;
  const COL = {
    verd: '#009639', verdFosc: '#00722b', zona: 'rgba(0,150,57,0.11)',
    blau: '#33469B', taronja: '#e67e22', gris: '#8a929e',
    eix: '#3C4250', reixa: '#e6ece8', text: '#4a5260'
  };

  // ---- Paràmetres (els controla l'usuari) ----
  const p = { model: 'kepler', n: -1.5, h: 0.02, l: 1, fracE: 0.1, vel: 1, enMarxa: true, lineal: true, comp: true };

  // ---- Magnituds derivades (es recalculen a calcula()) ----
  const d = {};
  // ---- Estat dinàmic ----
  let s = null;          // { r, pr, th, t }
  let traca = [];        // punts { r, th, t }
  let finalitzat = null; // { motiu, quan }
  let E0num = 0;

  // =====================================================================
  //  Física
  // =====================================================================
  const esLog = () => Math.abs(p.n + 1) < 1e-9;

  function V(r) {
    switch (p.model) {
      case 'kepler': return -k / r;
      case 'potencia': return esLog() ? k * Math.log(r) : k * Math.pow(r, p.n + 1) / (p.n + 1);
      default: return -k / r - d.h / (r * r * r);
    }
  }
  // Força radial f(r) = -dV/dr
  function f(r) {
    switch (p.model) {
      case 'kepler': return -k / (r * r);
      case 'potencia': return -k * Math.pow(r, p.n);
      default: return -k / (r * r) - 3 * d.h / Math.pow(r, 4);
    }
  }
  function dfdr(r) {
    switch (p.model) {
      case 'kepler': return 2 * k / (r * r * r);
      case 'potencia': return -k * p.n * Math.pow(r, p.n - 1);
      default: return 2 * k / (r * r * r) + 12 * d.h / Math.pow(r, 5);
    }
  }
  const Vcf = r => p.l * p.l / (2 * m * r * r);
  const Vef = r => V(r) + Vcf(r);
  const d2Vef = r => -dfdr(r) + 3 * p.l * p.l / (m * Math.pow(r, 4));

  function biseca(g, a, b) {
    let ga = g(a);
    for (let i = 0; i < 200; i++) {
      const c = 0.5 * (a + b), gc = g(c);
      if ((gc > 0) === (ga > 0)) { a = c; ga = gc; } else b = c;
    }
    return 0.5 * (a + b);
  }

  function calcula() {
    const L2m = p.l * p.l / m;

    // Correcció -h/r^3: cal l^4 > 12 m^2 k h perquè hi hagi òrbita circular
    d.h = 0; d.avisH = '';
    if (p.model === 'pert') {
      const hMax = L2m * L2m / (12 * k);
      d.h = Math.min(p.h, 0.95 * hMax);
      if (p.h > d.h) d.avisH = `Amb aquest l no hi ha òrbita circular si h ≥ ${hMax.toFixed(4)}. S'usa h = ${d.h.toFixed(4)}.`;
    }

    // Òrbita circular r0 (mínim de V_ef) i, si n'hi ha, màxim rb
    d.rb = 0;
    if (p.model === 'kepler') d.r0 = L2m / k;
    else if (p.model === 'potencia') d.r0 = Math.pow(L2m / k, 1 / (p.n + 3));
    else {
      const disc = Math.sqrt(L2m * L2m - 12 * k * d.h);
      d.r0 = (L2m + disc) / (2 * k);
      d.rb = (L2m - disc) / (2 * k);
    }
    const r0 = d.r0;

    d.Vmin = Vef(r0);
    d.S = p.l * p.l / (2 * m * r0 * r0);        // escala d'energia
    d.dE = p.fracE * d.S;
    d.E = d.Vmin + d.dE;
    d.Vb = d.rb > 0 ? Vef(d.rb) : Infinity;

    // Petites oscil·lacions
    d.k2 = d2Vef(r0);
    d.wt = p.l / (m * r0 * r0);
    d.wr = Math.sqrt(d.k2 / m);
    d.beta = d.wr / d.wt;
    d.A = Math.sqrt(2 * d.dE / d.k2);

    // Punts de retorn: V_ef(r) = E
    const g = r => Vef(r) - d.E;
    d.circular = d.dE < 1e-9 * d.S;
    if (d.circular) { d.rMin = r0; d.rMax = r0; }
    else {
      const rIn = d.rb > 0 ? d.rb : r0 * 1e-4;
      d.rMin = g(rIn) > 0 ? biseca(g, rIn, r0) : null;
      let R = 2 * r0;
      while (g(R) < 0 && R < 1e5 * r0) R *= 2;
      d.rMax = g(R) >= 0 ? biseca(g, Math.max(r0, R / 2), R) : null;
    }

    d.Rview = d.rMax !== null ? Math.min(Math.max(1.15 * d.rMax, 2 * r0), 8 * r0) : 4 * r0;

    // Angle apsidal i període radial exactes (quadratura amb r = c + a sin σ)
    d.psi = null; d.Tr = null;
    if (d.circular) { d.psi = Math.PI / d.beta; d.Tr = 2 * Math.PI / d.wr; }
    else if (d.rMin !== null && d.rMax !== null) {
      const c = 0.5 * (d.rMax + d.rMin), a = 0.5 * (d.rMax - d.rMin), N = 4000;
      let sPsi = 0, sT = 0;
      for (let i = 0; i < N; i++) {
        const sg = -Math.PI / 2 + (i + 0.5) * Math.PI / N;
        const r = c + a * Math.sin(sg), dr = a * Math.cos(sg);
        const K = d.E - Vef(r);
        if (K <= 0) continue;
        const arrel = Math.sqrt(2 * m * K);
        sPsi += p.l / (r * r) * dr / arrel;
        sT += m * dr / arrel;
      }
      d.psi = sPsi * Math.PI / N;
      d.Tr = 2 * sT * Math.PI / N;
    }
  }

  // ---- Integració: r, p_r, θ amb RK4 i pas adaptatiu ----
  function deriv(r, pr) {
    return [pr / m, p.l * p.l / (m * r * r * r) + f(r), p.l / (m * r * r)];
  }
  function rk4(st, h) {
    const k1 = deriv(st.r, st.pr);
    const k2 = deriv(st.r + 0.5 * h * k1[0], st.pr + 0.5 * h * k1[1]);
    const k3 = deriv(st.r + 0.5 * h * k2[0], st.pr + 0.5 * h * k2[1]);
    const k4 = deriv(st.r + h * k3[0], st.pr + h * k3[1]);
    st.r += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    st.pr += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    st.th += h / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    st.t += h;
  }
  const pasAdaptiu = st => 0.004 * st.r / Math.hypot(st.pr / m, p.l / (m * st.r));
  const energia = st => st.pr * st.pr / (2 * m) + Vef(st.r);

  function reinicia() {
    traca = []; finalitzat = null;
    const sortida = 1.45 * d.Rview;
    if (d.rMin !== null) {
      // Comença al pericentre (t = 0, θ = 0)
      s = { r: d.rMin, pr: 0, th: 0, t: 0 };
      if (d.rMax === null) {
        // No lligada: integra enrere fins que la partícula és lluny
        for (let i = 0; i < 400000 && s.r < sortida; i++) rk4(s, -pasAdaptiu(s));
      }
    } else if (d.rMax !== null) {
      s = { r: d.rMax, pr: 0, th: 0, t: 0 };            // cau al centre des de l'apocentre
    } else {
      s = { r: sortida, pr: -Math.sqrt(2 * m * Math.max(d.E - Vef(sortida), 0)), th: 0, t: 0 };
    }
    s.t0 = s.t; s.th0 = s.th;
    E0num = energia(s);
    afegeixPunt(true);
  }

  function afegeixPunt(forca) {
    const u = traca[traca.length - 1];
    const Tth = 2 * Math.PI / d.wt;
    if (forca || !u || Math.abs(s.th - u.th) > 0.004 || s.t - u.t > Tth / 600) {
      traca.push({ r: s.r, th: s.th, t: s.t });
      if (traca.length > 40000) traca.splice(0, 10000);
    }
  }

  function avanca(ts) {
    if (!p.enMarxa || finalitzat) return;
    const objectiu = s.t + p.vel * (2 * Math.PI / d.wt) / 300;
    for (let i = 0; i < 20000 && s.t < objectiu; i++) {
      rk4(s, Math.min(pasAdaptiu(s), objectiu - s.t + 1e-12));
      afegeixPunt(false);
      if (s.r < 0.02 * d.r0) { finalitzat = { motiu: 'centre', quan: ts }; break; }
      if (d.rMax === null && s.pr > 0 && s.r > 1.5 * d.Rview) { finalitzat = { motiu: 'escapa', quan: ts }; break; }
    }
  }

  // =====================================================================
  //  Dibuix
  // =====================================================================
  const cvPot = document.getElementById('cv_pot');
  const cvOrb = document.getElementById('cv_orb');
  const cvRt = document.getElementById('cv_rt');

  function prepara(cv) {
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    const W = Math.round(w * dpr), H = Math.round(h * dpr);
    if (cv.width !== W || cv.height !== H) { cv.width = W; cv.height = H; }
    const ctx = cv.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.font = '11px system-ui, sans-serif';
    return { ctx, w, h };
  }

  function marques(a, b, n) {
    const pas0 = (b - a) / n, mag = Math.pow(10, Math.floor(Math.log10(pas0))), q = pas0 / mag;
    const pas = (q < 1.5 ? 1 : q < 3 ? 2 : q < 7 ? 5 : 10) * mag;
    const dec = Math.max(0, -Math.floor(Math.log10(pas) + 1e-9));
    const llista = [];
    for (let v = Math.ceil(a / pas) * pas; v <= b + 1e-9 * pas; v += pas) llista.push(Math.abs(v) < 1e-12 ? 0 : v);
    return { llista, dec };
  }

  const retalla = y => Math.max(-1e4, Math.min(1e4, y));

  function eixos(ctx, x0, x1, y0, y1, X, Y, xa, xb, ya, yb, etX, etY) {
    ctx.strokeStyle = COL.reixa; ctx.lineWidth = 1; ctx.fillStyle = COL.text;
    const mx = marques(xa, xb, 6), my = marques(ya, yb, 5);
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    mx.llista.forEach(v => {
      const x = X(v);
      ctx.beginPath(); ctx.moveTo(x, y0); ctx.lineTo(x, y1); ctx.stroke();
      ctx.fillText(v.toFixed(mx.dec), x, y1 + 4);
    });
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    my.llista.forEach(v => {
      const y = Y(v);
      ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x1, y); ctx.stroke();
      ctx.fillText(v.toFixed(my.dec), x0 - 5, y);
    });
    ctx.strokeStyle = COL.eix; ctx.lineWidth = 1.2;
    ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
    ctx.font = 'italic 12px serif'; ctx.fillStyle = COL.eix;
    ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
    ctx.fillText(etX, x1, y1 + 30);
    ctx.save(); ctx.translate(12, (y0 + y1) / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(etY, 0, 0); ctx.restore();
    ctx.font = '11px system-ui, sans-serif';
  }

  function corba(ctx, fn, X, Y, a, b, N) {
    ctx.beginPath();
    for (let i = 0; i <= N; i++) {
      const r = a + (b - a) * i / N, y = retalla(Y(fn(r)));
      if (i === 0) ctx.moveTo(X(r), y); else ctx.lineTo(X(r), y);
    }
    ctx.stroke();
  }

  let transfPot = null;

  function dibuixaPot() {
    const { ctx, w, h } = prepara(cvPot);
    const x0 = 52, x1 = w - 12, y0 = 12, y1 = h - 34;
    const rMaxP = d.Rview, rMinP = rMaxP / 800;
    const yMin = d.Vmin - 0.4 * d.S, yMax = Math.max(d.Vmin + 2.4 * d.S, d.E + 0.6 * d.S);
    const X = r => x0 + r / rMaxP * (x1 - x0);
    const Y = v => y1 - (v - yMin) / (yMax - yMin) * (y1 - y0);
    transfPot = { y0, y1, yMin, yMax };

    eixos(ctx, x0, x1, y0, y1, X, Y, 0, rMaxP, yMin, yMax, 'r', 'V_ef(r)');

    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();

    // V = 0
    if (yMin < 0 && yMax > 0) {
      ctx.strokeStyle = COL.gris; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0, Y(0)); ctx.lineTo(x1, Y(0)); ctx.stroke();
    }

    // Zona permesa: entre V_ef i E
    if (!d.circular) {
      const a = d.rMin !== null ? d.rMin : rMinP, b = d.rMax !== null ? Math.min(d.rMax, rMaxP) : rMaxP;
      ctx.fillStyle = COL.zona;
      ctx.beginPath(); ctx.moveTo(X(a), Y(d.E));
      for (let i = 0; i <= 300; i++) {
        const r = a + (b - a) * i / 300;
        ctx.lineTo(X(r), retalla(Y(Math.min(Vef(r), d.E))));
      }
      ctx.lineTo(X(b), Y(d.E)); ctx.closePath(); ctx.fill();
    }

    // Components V(r) i l²/2mr²
    if (p.comp) {
      ctx.setLineDash([5, 4]); ctx.lineWidth = 1.3; ctx.strokeStyle = COL.gris;
      corba(ctx, V, X, Y, rMinP, rMaxP, 600);
      corba(ctx, Vcf, X, Y, rMinP, rMaxP, 600);
      ctx.setLineDash([]);
    }

    // Aproximació parabòlica
    if (p.lineal) {
      ctx.setLineDash([6, 4]); ctx.lineWidth = 1.8; ctx.strokeStyle = COL.taronja;
      corba(ctx, r => d.Vmin + 0.5 * d.k2 * (r - d.r0) * (r - d.r0), X, Y, rMinP, rMaxP, 400);
      ctx.setLineDash([]);
    }

    // V_ef
    ctx.lineWidth = 2.6; ctx.strokeStyle = COL.verd;
    corba(ctx, Vef, X, Y, rMinP, rMaxP, 800);

    // Energia
    ctx.strokeStyle = COL.blau; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x0, Y(d.E)); ctx.lineTo(x1, Y(d.E)); ctx.stroke();
    ctx.fillStyle = COL.blau; ctx.font = 'italic 13px serif'; ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
    ctx.fillText('E', x1 - 4, Y(d.E) - 3);
    ctx.font = '11px system-ui, sans-serif';

    // Punts de retorn i r0
    ctx.setLineDash([3, 3]); ctx.lineWidth = 1;
    const vertical = (r, col, et) => {
      if (r === null || r > rMaxP) return;
      ctx.strokeStyle = col;
      ctx.beginPath(); ctx.moveTo(X(r), retalla(Y(Vef(r)))); ctx.lineTo(X(r), y1); ctx.stroke();
      ctx.fillStyle = col; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
      ctx.fillText(et, X(r), y1 - 2);
    };
    vertical(d.r0, COL.gris, 'r₀');
    if (!d.circular) { vertical(d.rMin, COL.blau, 'r_min'); vertical(d.rMax, COL.blau, 'r_max'); }
    ctx.setLineDash([]);

    // Partícula: segment = ½ m ṙ²
    if (s && s.r < rMaxP) {
      const xp = X(s.r), yp = Y(Vef(s.r));
      ctx.strokeStyle = COL.blau; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.moveTo(xp, yp); ctx.lineTo(xp, Y(d.E)); ctx.stroke();
      ctx.fillStyle = COL.verdFosc; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(xp, yp, 6.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    }
    ctx.restore();
  }

  function dibuixaOrb() {
    const { ctx, w, h } = prepara(cvOrb);
    const cx = w / 2, cy = h / 2, sc = (Math.min(w, h) / 2 - 8) / d.Rview;
    const P = (r, th) => [cx + sc * r * Math.cos(th), cy - sc * r * Math.sin(th)];

    // Zona permesa
    if (!d.circular) {
      ctx.fillStyle = COL.zona;
      ctx.beginPath();
      ctx.arc(cx, cy, sc * (d.rMax !== null ? d.rMax : 3 * d.Rview), 0, 2 * Math.PI);
      if (d.rMin !== null) { ctx.moveTo(cx + sc * d.rMin, cy); ctx.arc(cx, cy, sc * d.rMin, 0, 2 * Math.PI); }
      ctx.fill('evenodd');
    }

    // Eixos i òrbita circular
    ctx.strokeStyle = COL.reixa; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(0, cy); ctx.lineTo(w, cy); ctx.moveTo(cx, 0); ctx.lineTo(cx, h); ctx.stroke();
    ctx.strokeStyle = COL.gris; ctx.setLineDash([2, 4]);
    ctx.beginPath(); ctx.arc(cx, cy, sc * d.r0, 0, 2 * Math.PI); ctx.stroke();
    ctx.setLineDash([]);

    // Traça exacta
    if (traca.length > 1) {
      ctx.strokeStyle = COL.verd; ctx.lineWidth = 2;
      ctx.beginPath();
      traca.forEach((q, i) => { const [x, y] = P(q.r, q.th); if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y); });
      const [xs, ys] = P(s.r, s.th); ctx.lineTo(xs, ys);
      ctx.stroke();
    }

    // Aproximació lineal en u = 1/r: u = u0 ± a cos βθ
    if (p.lineal && s && traca.length) {
      const u0 = 1 / d.r0, a = d.A / (d.r0 * d.r0), signe = d.rMin !== null ? 1 : -1;
      const thA = traca[0].th, thB = s.th, N = Math.min(20000, Math.ceil(Math.abs(thB - thA) / 0.01) + 1);
      ctx.strokeStyle = COL.taronja; ctx.lineWidth = 1.8; ctx.setLineDash([6, 4]);
      ctx.beginPath();
      let dins = false;
      for (let i = 0; i <= N; i++) {
        const th = thA + (thB - thA) * i / N, u = u0 + signe * a * Math.cos(d.beta * th);
        if (u <= 0 || 1 / u > 2 * d.Rview) { dins = false; continue; }
        const [x, y] = P(1 / u, th);
        if (dins) ctx.lineTo(x, y); else ctx.moveTo(x, y);
        dins = true;
      }
      ctx.stroke(); ctx.setLineDash([]);
    }

    // Centre de forces
    ctx.fillStyle = COL.eix;
    ctx.beginPath(); ctx.arc(cx, cy, 4, 0, 2 * Math.PI); ctx.fill();

    // Partícula i radi vector
    if (s) {
      const [x, y] = P(s.r, s.th);
      ctx.strokeStyle = COL.gris; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(x, y); ctx.stroke();
      ctx.fillStyle = COL.verdFosc; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(x, y, 6.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
    }

    if (finalitzat) {
      ctx.fillStyle = '#b03a2e'; ctx.font = 'bold 13px system-ui, sans-serif';
      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
      ctx.fillText(finalitzat.motiu === 'centre' ? 'Ha caigut al centre' : 'Ha escapat', cx, 8);
    }
  }

  function dibuixaRt() {
    const { ctx, w, h } = prepara(cvRt);
    const x0 = 52, x1 = w - 12, y0 = 10, y1 = h - 34;
    const lligada = d.rMin !== null && d.rMax !== null;
    const W = lligada ? 2.5 * d.Tr : Math.max(2 * Math.PI / d.wt, s.t - traca[0].t);
    const tA = Math.max(traca[0].t, s.t - W), tB = tA + W;
    const rTop = 1.05 * d.Rview;
    const X = t => x0 + (t - tA) / W * (x1 - x0);
    const Y = r => y1 - r / rTop * (y1 - y0);

    eixos(ctx, x0, x1, y0, y1, X, Y, tA, tB, 0, rTop, 't', 'r');

    ctx.save();
    ctx.beginPath(); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.clip();

    ctx.setLineDash([2, 4]); ctx.strokeStyle = COL.gris; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(x0, Y(d.r0)); ctx.lineTo(x1, Y(d.r0)); ctx.stroke();
    ctx.setLineDash([]);

    if (p.lineal && lligada) {
      ctx.strokeStyle = COL.taronja; ctx.lineWidth = 1.8; ctx.setLineDash([6, 4]);
      corba(ctx, t => d.r0 - d.A * Math.cos(d.wr * t), X, Y, tA, tB, 600);
      ctx.setLineDash([]);
    }

    ctx.strokeStyle = COL.verd; ctx.lineWidth = 2;
    ctx.beginPath();
    let primer = true;
    for (const q of traca) {
      if (q.t < tA) continue;
      const y = Y(q.r);
      if (primer) { ctx.moveTo(X(q.t), y); primer = false; } else ctx.lineTo(X(q.t), y);
    }
    ctx.lineTo(X(s.t), Y(s.r));
    ctx.stroke();
    ctx.restore();
  }

  // =====================================================================
  //  Resultats
  // =====================================================================
  const tex = t => (window.katex ? katex.renderToString(t, { throwOnError: false }) : t);
  const fmt = (x, n = 3) => (x === null || x === undefined) ? '—' : !isFinite(x) ? '∞' : x.toFixed(n);
  const graus = x => x === null ? '—' : `${(x * 180 / Math.PI).toFixed(2)}° (${(x / Math.PI).toFixed(4)}π)`;

  function taula(id, files) {
    document.getElementById(id).innerHTML = files
      .map(([clau, et]) => `<tr><td>${et}</td><td id="o_${clau}"></td></tr>`).join('');
  }
  taula('t_orbita', [
    ['r0', tex('r_0') + ' (circular)'],
    ['Vmin', tex('V_\\text{ef,min}')],
    ['E', tex('E')],
    ['e', tex('e') + ' (Kepler)'],
    ['rmin', tex('r_\\text{min}')],
    ['rmax', tex('r_\\text{max}')],
    ['r', tex('r(t)')],
    ['rdot', tex('\\dot r')],
    ['thdot', tex('\\dot\\theta = l/mr^2')],
    ['deriva', 'Error numèric en ' + tex('E')]
  ]);
  taula('t_lineal', [
    ['wt', tex('\\omega_\\theta = l/mr_0^2')],
    ['wr', tex("\\omega_r = \\sqrt{V_\\text{ef}''(r_0)/m}")],
    ['beta', tex('\\beta = \\omega_r/\\omega_\\theta')],
    ['A', tex('A/r_0')],
    ['psil', tex('\\psi_\\text{lin} = \\pi/\\beta')],
    ['psie', tex('\\psi') + ' exacte'],
    ['dosPsi', tex('2\\psi') + ' (entre pericentres)'],
    ['prec', tex('\\Delta\\varpi = 2\\psi - 2\\pi')]
  ]);
  const posa = (clau, v) => { document.getElementById('o_' + clau).textContent = v; };

  function tipusOrbita() {
    if (d.circular) return 'Circular';
    if (p.model === 'kepler') {
      if (Math.abs(d.E) < 1e-9 * d.S) return 'Parabòlica';
      return d.E < 0 ? 'El·líptica' : 'Hiperbòlica';
    }
    if (d.rMin === null) return 'Cau al centre (no hi ha prou barrera centrífuga)';
    if (d.rMax === null) return 'No lligada: escapa a l\'infinit';
    const q = d.psi / Math.PI;
    for (let den = 1; den <= 8; den++) {
      const num = Math.round(q * den);
      if (Math.abs(q - num / den) < 1e-3) return `Lligada i tancada (ψ = ${num === 1 ? '' : num}π${den === 1 ? '' : '/' + den})`;
    }
    return 'Lligada, no tancada (roseta)';
  }

  function actualitzaEstatics() {
    document.getElementById('o_tipus').textContent = tipusOrbita();
    posa('r0', fmt(d.r0)); posa('Vmin', fmt(d.Vmin)); posa('E', fmt(d.E));
    posa('e', p.model === 'kepler' ? fmt(Math.sqrt(Math.max(0, 1 + 2 * d.E * p.l * p.l / (m * k * k)))) : '—');
    posa('rmin', d.rMin === null ? '— (cau)' : fmt(d.rMin));
    posa('rmax', d.rMax === null ? '∞' : fmt(d.rMax));
    posa('wt', fmt(d.wt)); posa('wr', fmt(d.wr)); posa('beta', fmt(d.beta, 4));
    posa('A', fmt(d.A / d.r0)); posa('psil', graus(Math.PI / d.beta));
    posa('psie', graus(d.psi));
    posa('dosPsi', d.psi === null ? '—' : `${(2 * d.psi * 180 / Math.PI).toFixed(2)}°`);
    // La precessió només té sentit per a òrbites quasi keplerianes (β ≈ 1)
    posa('prec', d.psi === null || Math.abs(d.beta - 1) > 0.5 ? '—' : `${((2 * d.psi - 2 * Math.PI) * 180 / Math.PI).toFixed(2)}°`);

    document.getElementById('v_E').textContent = d.E.toFixed(3);
    document.getElementById('v_l').textContent = p.l.toFixed(2);
    document.getElementById('v_n').textContent = p.n.toFixed(2);
    document.getElementById('v_h').textContent = p.h.toFixed(3);
    document.getElementById('v_vel').textContent = p.vel.toFixed(1) + '×';
    const nota = document.getElementById('nota_h');
    nota.textContent = d.avisH; nota.classList.toggle('avis', !!d.avisH);
  }

  function actualitzaDinamics() {
    posa('r', fmt(s.r)); posa('rdot', fmt(s.pr / m)); posa('thdot', fmt(p.l / (m * s.r * s.r)));
    posa('deriva', (Math.abs(energia(s) - E0num) / d.S).toExponential(1));
  }

  // =====================================================================
  //  Controls
  // =====================================================================
  const $ = id => document.getElementById(id);
  const sliders = { s_n: 'n', s_h: 'h', s_l: 'l', s_E: 'fracE' };

  function sincronitza() {
    document.querySelectorAll('input[name=model]').forEach(r => { r.checked = r.value === p.model; });
    $('ctrl_n').hidden = p.model !== 'potencia';
    $('ctrl_h').hidden = p.model !== 'pert';
    Object.entries(sliders).forEach(([id, clau]) => { $(id).value = p[clau]; });
    $('s_vel').value = p.vel;
    $('c_lineal').checked = p.lineal; $('c_comp').checked = p.comp;
    $('b_marxa').textContent = p.enMarxa ? 'Pausa' : 'Continua';
  }

  function canvi(desmarcaExemple = true) {
    if (desmarcaExemple) document.querySelectorAll('.preset').forEach(b => b.classList.remove('actiu'));
    calcula(); reinicia(); sincronitza(); actualitzaEstatics();
  }

  document.querySelectorAll('input[name=model]').forEach(r =>
    r.addEventListener('change', () => { p.model = r.value; canvi(); }));
  Object.entries(sliders).forEach(([id, clau]) =>
    $(id).addEventListener('input', () => { p[clau] = parseFloat($(id).value); canvi(); }));
  $('s_vel').addEventListener('input', () => { p.vel = parseFloat($('s_vel').value); actualitzaEstatics(); });
  $('c_lineal').addEventListener('change', () => { p.lineal = $('c_lineal').checked; });
  $('c_comp').addEventListener('change', () => { p.comp = $('c_comp').checked; });
  $('b_marxa').addEventListener('click', () => { p.enMarxa = !p.enMarxa; sincronitza(); });
  $('b_reinicia').addEventListener('click', () => reinicia());

  const EXEMPLES = {
    circular:    { model: 'kepler', l: 1, fracE: 0 },
    pertorbada:  { model: 'kepler', l: 1, fracE: 0.08 },
    eliptica:    { model: 'kepler', l: 1, fracE: 0.36 },   // e = 0.6 → fracE = e²
    parabolica:  { model: 'kepler', l: 1, fracE: 1 },      // E = 0
    hiperbolica: { model: 'kepler', l: 1, fracE: 2.25 },   // e = 1.5
    harmonic:    { model: 'potencia', n: 1, l: 1, fracE: 0.6 },
    roseta:      { model: 'potencia', n: -1.5, l: 1, fracE: 0.5 },
    ondulada:    { model: 'potencia', n: 6, l: 1, fracE: 0.15 },     // β = 3
    relativista: { model: 'pert', h: 0.03, l: 1, fracE: 0.25 }
  };
  document.querySelectorAll('.preset').forEach(b => b.addEventListener('click', () => {
    Object.assign(p, EXEMPLES[b.dataset.preset]);
    document.querySelectorAll('.preset').forEach(x => x.classList.toggle('actiu', x === b));
    canvi(false);
  }));

  // Arrossegar la línia d'energia a la gràfica A
  let arrossegant = false;
  function energiaDesDe(ev) {
    if (!transfPot) return;
    const rect = cvPot.getBoundingClientRect(), y = ev.clientY - rect.top;
    const { y0, y1, yMin, yMax } = transfPot;
    const E = yMin + (y1 - y) / (y1 - y0) * (yMax - yMin);
    p.fracE = Math.round(Math.max(0, Math.min(3, (E - d.Vmin) / d.S)) * 1000) / 1000;
    canvi();
  }
  cvPot.addEventListener('pointerdown', ev => { arrossegant = true; cvPot.setPointerCapture(ev.pointerId); energiaDesDe(ev); });
  cvPot.addEventListener('pointermove', ev => { if (arrossegant) energiaDesDe(ev); });
  cvPot.addEventListener('pointerup', () => { arrossegant = false; });
  cvPot.addEventListener('pointercancel', () => { arrossegant = false; });

  // =====================================================================
  //  Inici
  // =====================================================================
  if (window.renderMathInElement) {
    renderMathInElement(document.body, {
      delimiters: [{ left: '$$', right: '$$', display: true }, { left: '$', right: '$', display: false }],
      throwOnError: false
    });
  }

  // ?exemple=nom obre directament un dels exemples
  const exInicial = document.querySelector(`.preset[data-preset="${new URLSearchParams(location.search).get('exemple')}"]`);
  if (exInicial) exInicial.click(); else canvi();
  let comptador = 0;
  function bucle(ts) {
    avanca(ts);
    if (finalitzat && ts - finalitzat.quan > 1500) reinicia();
    dibuixaPot(); dibuixaOrb(); dibuixaRt();
    if (comptador++ % 6 === 0) actualitzaDinamics();
    requestAnimationFrame(bucle);
  }
  requestAnimationFrame(bucle);
})();
