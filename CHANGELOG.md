# Changelog

## 1.0.1 — 2026-10-16

Restore public-entry numerical option checks lost during cleanup: refinement/backtracking/halving limits, finite steps and tolerances, singular-value cutoffs, deduplication bounds, finite delay grouping, and the supported weak-loss ratio. Valid-input numerical formulas remain unchanged.

Restore independent point-ray spreading/caustic and Gaussian-beam analytical benchmarks, integrate them into CTest/pytest, and allow beam numerical checks without visualization packages.

Validation: 11 CTest targets and 35 numerical Python tests pass; 2,716 point analytical/diagnostic checks and 23 beam analytical/convergence checks pass. All 1,340 pre-cleanup/current spreading, caustic, and convergence CSV rows match exactly. Seventeen regression cases reproduced the missing guards before repair.
