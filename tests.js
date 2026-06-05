/*
 * Node tests for the Oman Evidence OS engine. Run: node tests.js
 * Independently hand-computed expectations — NOT a re-derivation of the engine.
 */
const { NMA, runMarkovTrace, mulberry32 } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// ---------------------------------------------------------------------------
// normalCDF — A&S 7.1.26 erf form, Phi(x)=0.5(1+sign*y).
// Reference points: Phi(0)=0.5, Phi(1.96)=0.975, Phi(-1.96)=0.025, symmetry.
// ---------------------------------------------------------------------------
ok('normalCDF(0) ~ 0.5', close(NMA.normalCDF(0), 0.5, 1e-6), 'got ' + NMA.normalCDF(0));
ok('normalCDF(1.96) ~ 0.975', close(NMA.normalCDF(1.96), 0.975, 1e-3), 'got ' + NMA.normalCDF(1.96));
ok('normalCDF(-1.96) ~ 0.025', close(NMA.normalCDF(-1.96), 0.025, 1e-3), 'got ' + NMA.normalCDF(-1.96));
ok('normalCDF symmetric: Phi(z)+Phi(-z)=1', close(NMA.normalCDF(1.3) + NMA.normalCDF(-1.3), 1, 1e-9));
ok('normalCDF(1.645) ~ 0.95', close(NMA.normalCDF(1.645), 0.95, 1e-3), 'got ' + NMA.normalCDF(1.645));
ok('normalCDF monotone increasing', NMA.normalCDF(0.5) < NMA.normalCDF(0.6));
ok('normalCDF in [0,1]', NMA.normalCDF(5) <= 1 && NMA.normalCDF(-5) >= 0);

// ---------------------------------------------------------------------------
// calcLogOR — logOR = ln( (t2Ev*(t1Tot-t1Ev)) / ((t2Tot-t2Ev)*t1Ev) ),
// variance = 1/a + 1/b + 1/c + 1/d. Hand case: (10/100) vs (20/100):
//   logOR = ln((20*90)/(80*10)) = ln(2.25) = 0.8109302162
//   var   = 1/10 + 1/90 + 1/20 + 1/80 = 0.1736111111
// ---------------------------------------------------------------------------
const lor = NMA.calcLogOR(10, 100, 20, 100);
ok('calcLogOR logOR = ln(2.25)', close(lor.logOR, Math.log(2.25), 1e-12), 'got ' + lor.logOR);
ok('calcLogOR variance = sum of reciprocals', close(lor.variance, 1 / 10 + 1 / 90 + 1 / 20 + 1 / 80, 1e-12), 'got ' + lor.variance);
// Zero-cell -> Haldane-Anscombe 0.5 continuity correction (no NaN/Inf).
const lorZero = NMA.calcLogOR(0, 100, 5, 100);
ok('calcLogOR zero-cell finite (0.5 CC)', Number.isFinite(lorZero.logOR) && Number.isFinite(lorZero.variance), 'got ' + lorZero.logOR);
// CC formula check: a=0.5,b=100.5,c=5.5,d=95.5 -> ln((5.5*100.5)/(95.5*0.5))
ok('calcLogOR zero-cell uses 0.5 adj', close(lorZero.logOR, Math.log((5.5 * 100.5) / (95.5 * 0.5)), 1e-12), 'got ' + lorZero.logOR);
// No correction when all cells non-zero.
ok('calcLogOR no CC when full cells', NMA.calcLogOR(10, 100, 20, 100).variance === (1 / 10 + 1 / 90 + 1 / 20 + 1 / 80));

// ---------------------------------------------------------------------------
// invertMatrix — Gauss-Jordan with partial pivoting.
//   A=[[4,3],[6,3]], det = 12-18 = -6
//   A^-1 = (1/det)[[3,-3],[-6,4]] = [[-0.5,0.5],[1,-2/3]]
// ---------------------------------------------------------------------------
const inv = NMA.invertMatrix([[4, 3], [6, 3]]);
ok('invertMatrix [0][0] = -0.5', close(inv[0][0], -0.5, 1e-12), 'got ' + inv[0][0]);
ok('invertMatrix [0][1] = 0.5', close(inv[0][1], 0.5, 1e-12), 'got ' + inv[0][1]);
ok('invertMatrix [1][0] = 1', close(inv[1][0], 1, 1e-12), 'got ' + inv[1][0]);
ok('invertMatrix [1][1] = -2/3', close(inv[1][1], -2 / 3, 1e-12), 'got ' + inv[1][1]);
// A * A^-1 = I (property test).
const I = [[0, 0], [0, 0]];
const A = [[4, 3], [6, 3]];
for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) for (let k = 0; k < 2; k++) I[i][j] += A[i][k] * inv[k][j];
ok('invertMatrix A*Ainv = I', close(I[0][0], 1) && close(I[1][1], 1) && close(I[0][1], 0, 1e-12) && close(I[1][0], 0, 1e-12));
// Identity inverts to identity.
const invI = NMA.invertMatrix([[1, 0], [0, 1]]);
ok('invertMatrix identity -> identity', invI[0][0] === 1 && invI[1][1] === 1 && invI[0][1] === 0);
// Singular matrix throws.
let threw = false;
try { NMA.invertMatrix([[1, 2], [2, 4]]); } catch (e) { threw = true; }
ok('invertMatrix singular throws', threw);

// ---------------------------------------------------------------------------
// isConnected — DFS over adjacency.
// ---------------------------------------------------------------------------
ok('isConnected: linear chain connected', NMA.isConnected([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) === true);
ok('isConnected: isolated node disconnected', NMA.isConnected([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) === false);
ok('isConnected: empty -> false', NMA.isConnected([]) === false);

// ---------------------------------------------------------------------------
// buildAdjacency — symmetric edge counts from contrasts.
// ---------------------------------------------------------------------------
const adj = NMA.buildAdjacency(['a', 'b', 'c'], [
    { treatment1: 'a', treatment2: 'b' },
    { treatment1: 'a', treatment2: 'b' },
    { treatment1: 'b', treatment2: 'c' }
]);
ok('buildAdjacency symmetric', adj[0][1] === adj[1][0]);
ok('buildAdjacency counts parallel edges', adj[0][1] === 2 && adj[1][2] === 1);
ok('buildAdjacency no self-edge for absent pair', adj[0][2] === 0);

// ---------------------------------------------------------------------------
// pickReference — treatment appearing in the most multi-arm comparisons wins;
// preferred is honoured if present.
// ---------------------------------------------------------------------------
const arms = [
    { study: 'S1', treatment: 'x' }, { study: 'S1', treatment: 'y' },
    { study: 'S2', treatment: 'x' }, { study: 'S2', treatment: 'z' }
];
ok('pickReference honours preferred', NMA.pickReference(arms, ['x', 'y', 'z'], 'y') === 'y');
ok('pickReference picks hub (x) by score', NMA.pickReference(arms, ['x', 'y', 'z'], null) === 'x');

// ---------------------------------------------------------------------------
// computeNMA — single 2-treatment study: WLS must recover the direct logOR
// exactly (one contrast, one parameter). Control t4=10/100, drug t1=20/100.
//   logOR(t1 vs t4) = ln(2.25) = 0.8109302162
//   se = sqrt(1/10+1/90+1/20+1/80) = sqrt(0.1736111) = 0.4166666667
//   OR = 2.25; z = 0.81093/0.41667 = 1.94623; Phi(z) = 0.97419
//   pScore(t1) = Phi(z)*100 = 97.419; pScore(t4) = (1-Phi(z))*100 = 2.581
// ---------------------------------------------------------------------------
const nma = NMA.computeNMA({
    effects: { t1: { name: 'Drug' }, t4: { name: 'Control' } },
    trials: [{ id: 'S1', drug: 't1', ctrl: 't4', ev_tx: 20, n_tx: 100, ev_ctrl: 10, n_ctrl: 100 }]
});
ok('computeNMA no error', !nma.error, nma.error);
ok('computeNMA ref = t4', nma.ref === 't4');
ok('computeNMA t1 logOR = ln(2.25)', close(nma.effects.t1.logOR, Math.log(2.25), 1e-9), 'got ' + nma.effects.t1.logOR);
ok('computeNMA t1 OR = 2.25', close(nma.effects.t1.or, 2.25, 1e-9), 'got ' + nma.effects.t1.or);
ok('computeNMA t1 se = sqrt(var)', close(nma.effects.t1.se, Math.sqrt(1 / 10 + 1 / 90 + 1 / 20 + 1 / 80), 1e-9), 'got ' + nma.effects.t1.se);
ok('computeNMA ref OR = 1', nma.effects.t4.or === 1);
ok('computeNMA p<=1 in CIs (OR_low<OR<OR_high)', nma.effects.t1.ci_low < nma.effects.t1.or && nma.effects.t1.or < nma.effects.t1.ci_high);
// pScore: t1 favoured (higher OR = more events; harmful but ranks first by raw effect).
const pT1 = nma.rankings.find(r => r.key === 't1').pScore;
const pT4 = nma.rankings.find(r => r.key === 't4').pScore;
ok('computeNMA pScore(t1) ~ 97.419', close(pT1, NMA.normalCDF(Math.log(2.25) / 0.4166666667) * 100, 1e-2), 'got ' + pT1);
ok('computeNMA pScores sum to ~100 (k=2)', close(pT1 + pT4, 100, 1e-6), 'got ' + (pT1 + pT4));
ok('computeNMA rankings sorted desc', nma.rankings[0].pScore >= nma.rankings[1].pScore);

// computeNMA guards: too little data returns an error object, not a throw.
ok('computeNMA empty -> error', !!NMA.computeNMA({ effects: {}, trials: [] }).error);
ok('computeNMA single-treatment -> error', !!NMA.computeNMA({
    effects: { t1: {} },
    trials: [{ id: 'S', drug: 't1', ctrl: 't1', ev_tx: 5, n_tx: 50, ev_ctrl: 4, n_ctrl: 50 }]
}).error);

// ---------------------------------------------------------------------------
// runMarkovTrace — deterministic 10-year cohort.
//   No death, no events, annual cost 100, util 1, discount 0:
//   cost = 10*100 = 1000; QALY = 10*1 = 10.
// ---------------------------------------------------------------------------
const mk = runMarkovTrace({ costAnnual: 100, costEvent: 0, utilState: 1, probEvent: 0, probDeath: 0, discount: 0 });
ok('runMarkovTrace undiscounted cost = 1000', close(mk.c, 1000, 1e-9), 'got ' + mk.c);
ok('runMarkovTrace undiscounted QALY = 10', close(mk.e, 10, 1e-9), 'got ' + mk.e);
// Discounting reduces totals.
const mkD = runMarkovTrace({ costAnnual: 100, costEvent: 0, utilState: 1, probEvent: 0, probDeath: 0, discount: 0.035 });
ok('runMarkovTrace discounting lowers cost', mkD.c < mk.c && mkD.c > 0, 'got ' + mkD.c);
// Year-0 factor is 1, so discounted >= 9 years of cost added at <1 factor: 100..~870ish.
ok('runMarkovTrace discounted cost in plausible range', mkD.c > 800 && mkD.c < 1000, 'got ' + mkD.c);
// With death probability, cohort depletes -> totals strictly less than the
// no-death case for cost (fewer alive person-years).
const mkDeath = runMarkovTrace({ costAnnual: 100, costEvent: 0, utilState: 1, probEvent: 0, probDeath: 0.1, discount: 0 });
ok('runMarkovTrace mortality reduces person-years', mkDeath.e < 10 && mkDeath.e > 0, 'got ' + mkDeath.e);

// ---------------------------------------------------------------------------
// mulberry32 — deterministic, seed-stable, in [0,1).
// ---------------------------------------------------------------------------
const r1 = mulberry32(12345);
const r2 = mulberry32(12345);
const a = r1(), b = r1();
ok('mulberry32 deterministic per seed', r2() === a && r2() === b);
ok('mulberry32 in [0,1)', a >= 0 && a < 1 && b >= 0 && b < 1);
ok('mulberry32 advances state', a !== b);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
