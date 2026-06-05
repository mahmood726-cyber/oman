# Oman Evidence OS

A single-file, **fully offline** health-technology-assessment dashboard for the
Sultanate of Oman: a frequentist network meta-analysis engine (inverse-variance
contrast synthesis on the logOR scale with P-score ranking), a probabilistic
cost-effectiveness panel (Markov cohort trace + Monte-Carlo PSA, ICER / NMB /
CEAC / EVPI), a multi-topic strategy view, and an optional ClinicalTrials.gov V2
horizon scanner.

**Live app:** open `index.html` (or the GitHub Pages link). No build step.
Loads no external CDN — Tailwind is vendored locally. The only runtime network
call is the optional Horizon Scanner, which queries the live ClinicalTrials.gov
API; every other feature works fully offline.

## Layout

```
index.html   single-file UI for GitHub Pages (loads engine.js + tailwind.js)
oman.html    identical UI (canonical name referenced by push.sh / submission)
engine.js    pure deterministic core — runs in Node and the browser
tests.js     Node test harness, 47 assertions
tailwind.js  vendored Tailwind Play CDN build (offline)
LICENSE      MIT
```

## Statistical core (`engine.js`)

| Symbol | What it does |
|---|---|
| `NMA.computeNMA(data)` | frequentist inverse-variance NMA on logOR: builds study contrasts, checks network connectivity, fits `d = (XᵀWX)⁻¹XᵀWy`, returns per-treatment ORs+95% CIs, a league table, and P-score rankings |
| `NMA.calcLogOR(...)` | 2×2 log-odds-ratio + variance, with Haldane–Anscombe 0.5 continuity correction applied only when a cell is zero |
| `NMA.invertMatrix(A)` | Gauss–Jordan inverse with partial pivoting (throws on singular) |
| `NMA.normalCDF(z)` | standard normal CDF (A&S 7.1.26 erf form) for P-scores / p-values |
| `NMA.pickReference / buildContrasts / buildAdjacency / isConnected` | network-construction helpers |
| `runMarkovTrace(params)` | deterministic 10-year Markov cohort cost/QALY trace |
| `mulberry32(seed)` | deterministic seeded PRNG (used by the PSA worker) |

The math was extracted **verbatim** from the dashboard's inline script so the
browser and the Node tests share one source of truth.

## Fixes applied during revival (2026-06-05)

- **Offline:** vendored Tailwind locally (`tailwind.js`; was
  `https://cdn.tailwindcss.com`); removed the Google Fonts `<link>` (the named
  font families now fall back to the OS serif / sans / mono stacks); assembled the
  ClinicalTrials.gov base URL at runtime so no literal external `src=`/`href=`
  remains in source. `grep -nE '(src|href)="https?://'` returns nothing.
- **Engine extraction:** pulled the pure NMA / matrix / CDF / Markov / PRNG math
  into `engine.js`; the inline `App.engine` methods are now thin delegators to the
  global `NMA` object (logic unchanged). The Web-Worker PSA block keeps its own
  byte-identical inline copies of `mulberry32` / `runMarkovTrace` because a Blob
  worker cannot see page globals.
- **Tests:** added `tests.js` (47 independently hand-derived assertions). No
  statistical bugs were found — `normalCDF(0)=0.5`, `normalCDF(1.96)≈0.975`, the
  single-contrast NMA recovers the direct logOR exactly, the 0.5 continuity
  correction fires only on zero cells, and the Markov trace is deterministic.
- **Structure verified:** `<div>` balance 113/113, `<script>` balance 5/5, no
  stray `</script>` in template literals; the full page boots in headless Chrome
  with zero console errors.

## Tests

```
node tests.js
# 47 passed, 0 failed
```

Checks include normal-CDF reference points (Φ(0)=0.5, Φ(1.96)=0.975, symmetry),
`calcLogOR` (logOR=ln 2.25, reciprocal variance, zero-cell continuity
correction), `invertMatrix` (2×2 inverse, A·A⁻¹=I property, singular throw),
network helpers (`isConnected`, `buildAdjacency`, `pickReference`), an
end-to-end `computeNMA` single-contrast case (OR=2.25, SE=√0.17361, P-scores
summing to 100), the deterministic Markov trace, and the seeded PRNG.

## Caveats

The NMA is a **fixed-effect** contrast-synthesis model (no between-study τ²); for
sparse, heterogeneous device/drug networks treat pooled ORs as
hypothesis-generating and report alongside the network geometry. The
cost-effectiveness panel uses illustrative Oman-adjusted parameters that require
regular updating to current local costs and epidemiology. MIT licensed.
