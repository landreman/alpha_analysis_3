# M0 bootstrap checks

Performed on 2026-11-10 on the researcher's laptop in a clean `.venv` created from the
conda `20220806-03` interpreter (Python 3.10.5).

- The five boozmn files were copied from `landreman/alpha_analysis` at revision
  `95519166d4eeefc25f1eca331c8d706961174f39`; their SHA-256 hashes in
  `data/manifest.json` match the hashes recorded independently in `alpha_analysis_2`.
- The original note was copied from the researcher's Box folder; its hash is in
  `docs/references/provenance.json` and is verified by `test_reference_documents_unchanged`.
- `src/alpha_analysis_3/boozer_field.py` is byte-identical to the original repository's
  reader at the same revision (origin hash in `provenance.json`).
- `test_five_real_inputs_and_fourier_conventions` checks each file's hash, field-period
  count, the Fourier evaluation against a direct mode sum, periodicity in ζ, and the
  along-line derivative against a finite difference.
- `test_error_target_and_pitch_conventions` pins the accepted tolerance, the six λₙ levels,
  the source profile evaluated at ρ, and the direction sign.
- Results of `make check` and `make smoke`, and the first GitHub `Tests` run, are recorded
  in the bootstrap commit message and in STATUS once CI has run.

No numerical machinery beyond the field reader exists at M0. Nothing here is evidence
about f.
