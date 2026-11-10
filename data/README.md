# Reference equilibria

The five required physical inputs, copied byte for byte from `landreman/alpha_analysis`.
[manifest.json](manifest.json) records their hashes, sizes, field periods and the legacy
estimates of the radially global extrema of B. Keep the original filenames. Together they
occupy about 5 MB; use ordinary Git.

The diagnostic matrix is five files times six levels, λₙ = 0.05, 0.1, 0.5, 0.8, 0.9, 0.95,
with b = B_min + λₙ (B_max − B_min) and the extrema over all radii. The levels probe
typical trapping classes; they do not replace the pitch integral and are not six
independent equilibria. Production runs record their own extrema and the actual b.

[legacy_probes.json](legacy_probes.json) retains two contour seeds and one generic
transition from the old numerical reports (actual b, root lifts, source hashes), copied
from `alpha_analysis_2`. They are comparison inputs for milestone M2, not acceptance
results.

The 24 MB W7-X reference boozmn, the VMEC wout file, old meshes, plots and run outputs
were not copied.
