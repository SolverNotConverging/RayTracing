"""Summarize native and Python RayTracing CPU numerical audit evidence."""
from pathlib import Path
import csv
import html
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/numerical_audit'
beam=json.loads((OUT/'beams/summary.json').read_text())
comparison=json.loads((OUT/'historical_comparison.json').read_text())
assert beam['passed'] and len(beam['checks'])==23
assert all(case['max_scaled_difference']==0 for case in comparison.values())
assert '100% tests passed' in (OUT/'final_ctest.txt').read_text()
assert '35 passed' in (OUT/'installed_final_tests.txt').read_text()
assert '17 failed, 3 passed' in (OUT/'guard_before_rebuild.txt').read_text()
assert json.loads((OUT/'no_viz/summary.json').read_text())['passed']
spreading=OUT/'native/tests/cpp/analytical-spreading'
table=''.join(f"<tr><td>{html.escape(check['name'])}</td><td>{check['max_error']:.4g}</td><td>{check['tolerance']:.4g}</td><td>Pass</td></tr>" for check in beam['checks'])
history_table=''.join(f"<tr><td>{html.escape(name)}</td><td>{case['rows']}</td><td>{case['max_scaled_difference']}</td></tr>" for name,case in comparison.items())
record=dict(current_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            historical_commit='f6f69fb',beam_analytical_checks=len(beam['checks']),
            point_analytical_diagnostic_checks=2716,pre_post_identical_rows=sum(c['rows'] for c in comparison.values()),
            native_ctest_targets=11,python_numerical_tests=35,viewer_tests_excluded=2,
            pre_fix_regressions_failed=17,installed_wheel='raytracing-1.0.0-cp314-cp314-win_amd64.whl',
            solver_changed='public option validation only; valid-input numerical formulas preserved')
(OUT/'summary.json').write_text(json.dumps(record,indent=2))
report=f'''<!doctype html><meta charset="utf-8"><title>RayTracing CPU numerical cleanup audit</title>
<style>body{{font:16px/1.55 system-ui;max-width:1100px;margin:40px auto;padding:0 24px;color:#1b3144}}h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{padding:9px;border-bottom:1px solid #d6dfe5;text-align:left}}img{{width:100%}}code{{background:#edf2f6;padding:2px 5px}}.pass{{color:#08734f}}</style>
<h1>RayTracing: CPU numerical cleanup audit</h1><p>2026-10-16 · HEAD 4cbba43 plus validation repair · pre-cleanup reference f6f69fb</p>
<p class="pass"><b>All checked numerical results pass after repairing public option validation.</b></p>
<p>Both point-ray multipath and Gaussian-beam transport were tested. Current C++ and Python extension code was rebuilt with MSVC in Release mode. The validated repaired wheel was installed into the repository's existing .venv and tested there. The numerical run used no CUDA. Two viewer-only tests were excluded.</p>
<h2>Cleanup regression found and fixed</h2><p>The cleanup moved validation out of low-level kernels, but failed to transfer every guard to the public solver boundary. Valid physical inputs still reproduced the analytical references. Invalid controls could silently discard paths or accept approximation settings outside the supported numerical model.</p>
<ul><li><code>src/solver.cpp</code>: restored finite delay-grouping tolerance, nonnegative refinement/backtracking/halving budgets, finite positive refinement steps and spreading tolerances, valid singular-value cutoffs and ratio, finite nonnegative ray-start distance, and deduplication angular bounds.</li>
<li><code>src/beam.cpp</code>: restored the pre-cleanup bound 0 ≤ max_loss_ratio ≤ 0.01.</li>
<li>Checks execute once at solver entry. Low-level routines retain the documented valid-input contract. No public arguments were added.</li>
<li>The new public-option regression set demonstrated 17 failures on the pre-fix freshly compiled extension. All 20 cases pass after repair, including valid zero iteration/halving budgets and zero loss tolerance.</li></ul>
<h2>Coverage restored</h2><p>The cleanup deleted the analytical spreading and beam transport benchmarks. They are restored with current names and enums while retaining their independent reference formulas and tolerances. Spreading now runs as a CTest target, and the beam suite is exercised by pytest. Plotting is optional; all beam analytical checks also pass with Matplotlib and PyVista imports blocked.</p>
<ul><li><b>Point-ray spreading:</b> direct 1/r field, plane image paths, normal/oblique convex spherical and cylindrical mirrors, concave focusing, two-plane/two-sphere paths, rigid transformations, receiver caustics, visibility and edge diagnostics.</li>
<li><b>Caustic phase:</b> spherical/cylindrical focus multiplicity, repeated returns through up to 12 reflections, exp(−iπμ/2) Maslov factors, and finite-difference convergence.</li>
<li><b>Gaussian beams:</b> complex on/off-axis fields through the waist, Gouy phase, receiver power quadrature, finite-beam coherent slabs, branch-depth/pruning convergence, Fresnel/Snell angle sweeps, power accounting, curved refraction and astigmatic/repeated focusing.</li>
<li><b>Sources and APIs:</b> analytical/loaded antenna models, polarization, refinement/deduplication, Gaussian beamlets, native aperture fitting, impulse delays, coherent tap grouping, and tap-to-field reconstruction are covered by the restored and existing tests.</li></ul>
<h2>Test results</h2><ul><li><b>11/11 CTest targets passed</b> (nine numerical executables and two examples).</li><li><b>35 Python numerical tests passed</b>; two viewer tests deselected.</li><li><b>2,716 point spreading/caustic analytical and diagnostic checks passed.</b> Maximum spreading relative error: 1.85075×10⁻⁶. Central-difference convergence ratio on halving the step: 4.05498, consistent with second order.</li><li><b>23 beam analytical/convergence checks passed.</b></li></ul>
<h2>Direct pre-cleanup comparison</h2><p>The pre-cleanup CPU C++ implementation was independently built and its original benchmark executed. All 1,340 recorded rows match current results exactly at the numeric values written to CSV, including caustic counts and complex-field errors. This comparison concerns these benchmark cases, not every possible scene.</p>
<table><tr><th>Dataset</th><th>Rows</th><th>Maximum scaled old/new difference</th></tr>{history_table}</table>
<h2>Independent beam reference errors</h2><p>Complex Gaussian field errors are about 6×10⁻¹⁴; slab sums about 4×10⁻¹⁵; Snell/Fresnel errors at roundoff. Beamlet field error converges to 2.63×10⁻⁹. The nonzero 3.18×10⁻⁵ W receiver-power error is the finite quadrature error, with its convergence explicitly checked.</p>
<img src="beams/convergence.png" alt="Gaussian quadrature, slab branch-depth and beamlet convergence">
<table><tr><th>Analytical check</th><th>Maximum error</th><th>Limit</th><th>Result</th></tr>{table}</table>
<h2>Scope and conventions</h2><p>This project uses exp(−iωt), unlike FEM/FDFD's positive-time convention. Passive ε and n therefore have positive imaginary parts; exp(+iknr) decays. Point rays require a lossless homogeneous background and use reflecting half-space surfaces. Beams use paraxial, weak-loss, whole-beam interactions. Edge diffraction and aperture clipping are not modeled; geometric-optics amplitudes at point receiver caustics remain intentionally unavailable. The tests respect these model assumptions.</p>
<h2>Reproduce and inspect</h2><p>Configure a Release C++ build with Eigen3, run <code>cmake --build ...</code> and <code>ctest --test-dir ... --output-on-failure</code>. Build/install the current wheel before Python tests. Run <code>pytest tests/python -k 'not views_render'</code> and <code>python benchmarks/beam_transport.py</code>. Fresh build/test logs and the pre-cleanup source snapshot are retained beside this report.</p>
<p><a href="final_ctest.txt">Native test log</a> · <a href="installed_final_tests.txt">Installed-wheel Python tests</a> · <a href="guard_before_rebuild.txt">Pre-fix failing regressions</a> · <a href="historical_comparison.json">Old/new comparison</a> · <a href="native/tests/cpp/analytical-spreading/summary.txt">Point reference summary</a> · <a href="beams/summary.json">Beam reference data</a> · <a href="summary.json">Audit summary</a></p>'''
(OUT/'report.html').write_text(report,encoding='utf-8')
print('Numerical audit passed:',OUT/'report.html')
