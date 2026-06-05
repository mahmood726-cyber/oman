/*
 * Oman Evidence OS — pure statistical / deterministic core.
 *
 * Extracted VERBATIM from the inline script of oman.html so the deterministic
 * math (frequentist network meta-analysis, logOR contrasts, Gaussian-elimination
 * matrix inversion, normal CDF, P-score ranking, and the deterministic Markov
 * cost-effectiveness trace) is a single source of truth, importable under Node
 * for testing.
 *
 * Browser: this file is loaded with <script src="engine.js"> BEFORE the inline
 * script, exposing NMA / runMarkovTrace / mulberry32 as globals. Node: the same
 * symbols are exported via module.exports at the bottom.
 *
 * Method: frequentist inverse-variance NMA on the logOR scale (fixed-effect
 * contrast-synthesis via weighted least squares, (XᵀWX)⁻¹XᵀWy), Haldane–Anscombe
 * 0.5 continuity correction only when a cell is zero, P-score ranking from the
 * normal CDF of pairwise z-statistics. Faithful to the shipped app — no
 * methodology changes.
 */

// Deterministic PRNG factory (mulberry32). Pure given a seed; returns a closure.
function mulberry32(a) { return function() { var t = a += 0x6D2B79F5; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; } }

// Deterministic 10-year Markov cohort trace (no RNG). Pure function of params.
function runMarkovTrace(params) {
    const { costAnnual, costEvent, utilState, probEvent, probDeath, discount } = params;
    let totalCost = 0; let totalQALY = 0; let stateHealth = 1.0; let stateDead = 0.0;
    for (let year = 0; year < 10; year++) {
        const discFactor = 1 / Math.pow(1 + discount, year);
        const died = stateHealth * probDeath;
        const hadEvent = stateHealth * probEvent;
        const cycleCost = (stateHealth * costAnnual) + (hadEvent * costEvent);
        const cycleUtil = (stateHealth * utilState) - (hadEvent * 0.05);
        totalCost += cycleCost * discFactor; totalQALY += cycleUtil * discFactor;
        stateDead += died; stateHealth = 1.0 - stateDead;
    }
    return { c: totalCost, e: totalQALY };
}

// Frequentist NMA engine. Extracted verbatim from App.engine in oman.html;
// methods call each other through `this`, so they are kept on a single object.
const NMA = {
    calcLogOR(t1Events, t1Total, t2Events, t2Total) {
        const a0 = t1Events, b0 = t1Total - t1Events;
        const c0 = t2Events, d0 = t2Total - t2Events;
        const useCC = (a0 === 0 || b0 === 0 || c0 === 0 || d0 === 0);
        const adj = useCC ? 0.5 : 0;
        const a = a0 + adj, b = b0 + adj, c = c0 + adj, d = d0 + adj;
        const logOR = Math.log((c * b) / (d * a));
        const variance = 1 / a + 1 / b + 1 / c + 1 / d;
        return { logOR, variance };
    },
    invertMatrix(A) {
        const n = A.length;
        const augmented = A.map((row, i) => {
            const out = [...row];
            for (let j = 0; j < n; j++) out.push(i === j ? 1 : 0);
            return out;
        });
        for (let i = 0; i < n; i++) {
            let maxRow = i;
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) maxRow = k;
            }
            [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];
            const pivot = augmented[i][i];
            if (Math.abs(pivot) < 1e-12) throw new Error('Singular matrix');
            for (let j = 0; j < 2 * n; j++) augmented[i][j] /= pivot;
            for (let k = 0; k < n; k++) {
                if (k === i) continue;
                const factor = augmented[k][i];
                for (let j = 0; j < 2 * n; j++) augmented[k][j] -= factor * augmented[i][j];
            }
        }
        return augmented.map(row => row.slice(n));
    },
    normalCDF(z) {
        const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        const sign = z < 0 ? -1 : 1;
        const x = Math.abs(z) / Math.sqrt(2);
        const t = 1 / (1 + p * x);
        const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
        return 0.5 * (1 + sign * y);
    },
    pickReference(arms, treatments, preferred) {
        if (preferred && treatments.includes(preferred)) return preferred;
        const byStudy = new Map();
        arms.forEach(a => {
            if (!byStudy.has(a.study)) byStudy.set(a.study, new Set());
            byStudy.get(a.study).add(a.treatment);
        });
        const scores = {};
        for (const set of byStudy.values()) {
            const list = Array.from(set);
            list.forEach((t) => {
                scores[t] = (scores[t] || 0) + (list.length - 1);
            });
        }
        return Object.entries(scores).sort((a, b) => b[1] - a[1])[0]?.[0] || treatments[0];
    },
    buildContrasts(arms, ref) {
        const byStudy = new Map();
        arms.forEach(a => {
            if (!byStudy.has(a.study)) byStudy.set(a.study, []);
            byStudy.get(a.study).push(a);
        });
        const contrasts = [];
        for (const studyArms of byStudy.values()) {
            if (studyArms.length < 2) continue;
            const baseline = studyArms.find(a => a.treatment === ref) || studyArms[0];
            studyArms.forEach((arm) => {
                if (arm === baseline || arm.treatment === baseline.treatment) return;
                const res = this.calcLogOR(baseline.events, baseline.n, arm.events, arm.n);
                if (!Number.isFinite(res.variance) || res.variance <= 0) return;
                contrasts.push({
                    treatment1: baseline.treatment,
                    treatment2: arm.treatment,
                    effect: res.logOR,
                    variance: res.variance
                });
            });
        }
        return contrasts;
    },
    buildAdjacency(treatments, contrasts) {
        const n = treatments.length;
        const matrix = Array(n).fill(null).map(() => Array(n).fill(0));
        contrasts.forEach(c => {
            const i = treatments.indexOf(c.treatment1);
            const j = treatments.indexOf(c.treatment2);
            if (i < 0 || j < 0) return;
            matrix[i][j] += 1;
            matrix[j][i] += 1;
        });
        return matrix;
    },
    isConnected(adjMatrix) {
        if (!adjMatrix.length) return false;
        const n = adjMatrix.length;
        const visited = new Array(n).fill(false);
        const stack = [0];
        while (stack.length) {
            const node = stack.pop();
            if (visited[node]) continue;
            visited[node] = true;
            for (let i = 0; i < n; i++) {
                if (adjMatrix[node][i] > 0 && !visited[i]) stack.push(i);
            }
        }
        return visited.every(Boolean);
    },
    computeNMA(data) {
        const effects = data.effects || {};
        const treatKeys = Object.keys(effects);
        const trials = data.trials || [];
        const arms = [];
        const warnings = [];

        trials.forEach((t, idx) => {
            const nTx = Number(t.n_tx || 0);
            const nCtrl = Number(t.n_ctrl || 0);
            if (nTx <= 0 || nCtrl <= 0) return;
            const baseId = String(t.id || t.author || `Study_${idx + 1}`);
            const studyId = `${baseId}-${idx + 1}`;
            const txKey = String(t.drug || 't1');
            const ctrlKey = String(t.ctrl || 't4');
            arms.push({ study: studyId, treatment: ctrlKey, events: t.ev_ctrl, n: nCtrl });
            arms.push({ study: studyId, treatment: txKey, events: t.ev_tx, n: nTx });
        });

        if (treatKeys.length < 2 || arms.length < 2) {
            return { error: 'Not enough trial data for NMA.', warnings, treatments: treatKeys };
        }

        const usedTreatments = [...new Set(arms.map(a => a.treatment))];
        if (usedTreatments.length < 2) {
            return { error: 'Network requires at least two treatments.', warnings, treatments: treatKeys };
        }

        const preferredRef = treatKeys.includes('t4') ? 't4' : null;
        const ref = this.pickReference(arms, usedTreatments, preferredRef);
        const treatments = [ref, ...treatKeys.filter(k => usedTreatments.includes(k) && k !== ref)];
        const missingTreatments = treatKeys.filter(k => !usedTreatments.includes(k));

        const contrasts = this.buildContrasts(arms, ref);
        if (contrasts.length < Math.max(1, treatments.length - 1)) {
            return { error: 'Insufficient contrasts to fit NMA.', warnings, treatments, missingTreatments, ref };
        }

        const adjacency = this.buildAdjacency(treatments, contrasts);
        if (!this.isConnected(adjacency)) {
            return { error: 'Network is not fully connected.', warnings, treatments, missingTreatments, ref };
        }

        const nParams = treatments.length - 1;
        const X = contrasts.map(c => {
            const row = new Array(nParams).fill(0);
            const idx1 = treatments.indexOf(c.treatment1);
            const idx2 = treatments.indexOf(c.treatment2);
            if (idx1 > 0) row[idx1 - 1] = -1;
            if (idx2 > 0) row[idx2 - 1] = 1;
            return row;
        });
        const y = contrasts.map(c => c.effect);
        const V = contrasts.map(c => c.variance);

        const XtWX = Array(nParams).fill(null).map(() => Array(nParams).fill(0));
        const XtWy = Array(nParams).fill(0);
        for (let i = 0; i < nParams; i++) {
            for (let j = 0; j < nParams; j++) {
                for (let k = 0; k < contrasts.length; k++) {
                    XtWX[i][j] += X[k][i] * X[k][j] / V[k];
                }
            }
            for (let k = 0; k < contrasts.length; k++) {
                XtWy[i] += X[k][i] * y[k] / V[k];
            }
        }

        let varCov;
        try {
            varCov = this.invertMatrix(XtWX);
        } catch (e) {
            return { error: 'NMA matrix inversion failed.', warnings, treatments, missingTreatments, ref };
        }

        const d = new Array(nParams).fill(0);
        for (let i = 0; i < nParams; i++) {
            for (let j = 0; j < nParams; j++) d[i] += varCov[i][j] * XtWy[j];
        }

        const effectsByKey = {};
        effectsByKey[ref] = { logOR: 0, se: 0, or: 1, ci_low: 1, ci_high: 1 };
        treatments.slice(1).forEach((key, idx) => {
            const se = Math.sqrt(Math.max(0, varCov[idx][idx]));
            const logOR = d[idx];
            const ciLow = logOR - 1.96 * se;
            const ciHigh = logOR + 1.96 * se;
            effectsByKey[key] = {
                logOR,
                se,
                or: Math.exp(logOR),
                ci_low: Math.exp(ciLow),
                ci_high: Math.exp(ciHigh)
            };
        });
        missingTreatments.forEach((key) => { effectsByKey[key] = { missing: true }; });

        const nTreat = treatments.length;
        const fullVar = Array(nTreat).fill(null).map(() => Array(nTreat).fill(0));
        for (let i = 1; i < nTreat; i++) {
            for (let j = 1; j < nTreat; j++) {
                fullVar[i][j] = varCov[i - 1][j - 1];
            }
        }

        const effectsVec = treatments.map((key, idx) => (idx === 0 ? 0 : d[idx - 1]));
        const league = {};
        const pairKey = (a, b) => `${a}|${b}`;
        for (let i = 0; i < nTreat; i++) {
            for (let j = 0; j < nTreat; j++) {
                if (i === j) {
                    league[pairKey(treatments[i], treatments[j])] = { self: true };
                    continue;
                }
                const diff = effectsVec[i] - effectsVec[j];
                const varDiff = fullVar[i][i] + fullVar[j][j] - 2 * fullVar[i][j];
                const seDiff = Math.sqrt(Math.max(0.0001, varDiff));
                const ciLow = diff - 1.96 * seDiff;
                const ciHigh = diff + 1.96 * seDiff;
                league[pairKey(treatments[i], treatments[j])] = {
                    or: Math.exp(diff),
                    ci_low: Math.exp(ciLow),
                    ci_high: Math.exp(ciHigh)
                };
            }
        }

        const rankings = [];
        for (let i = 0; i < nTreat; i++) {
            let sumP = 0;
            for (let j = 0; j < nTreat; j++) {
                if (i === j) continue;
                const diff = effectsVec[i] - effectsVec[j];
                const varDiff = fullVar[i][i] + fullVar[j][j] - 2 * fullVar[i][j];
                const seDiff = Math.sqrt(Math.max(0.0001, varDiff));
                const z = diff / seDiff;
                sumP += this.normalCDF(z);
            }
            rankings.push({ key: treatments[i], pScore: (sumP / (nTreat - 1)) * 100 });
        }
        rankings.sort((a, b) => b.pScore - a.pScore);

        return {
            ref,
            treatments,
            effects: effectsByKey,
            league,
            rankings,
            missingTreatments,
            warnings
        };
    }
};

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { NMA, runMarkovTrace, mulberry32 };
}
