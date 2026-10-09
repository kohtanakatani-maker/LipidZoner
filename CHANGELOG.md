# Changelog

## 1.1.0 (2026-10-08) – manuscript release

- Peak detection + alignment from raw mzML (pyOpenMS worker) integrated into
  the GUI (Detection tab), with automatic parameter estimation, gap filling
  and RT drift correction.
- Annotation pipeline: Match → Conflict search → IS filter → RT outlier
  filter → Adduct ion filter → Coherence filter → Quant ion selection.
- Positive and negative modes handled side by side; cross-mode RT shift
  estimated robustly (median + MAD).
- Coherence filter with per-class / linkage-aware RT regression, leverage-
  based confidence and loose attribution for single-scored conflicts.
- Export of quant-ion tables (height and area), non-target feature tables
  with confidence columns, reports and figures; LipidQuant relay export.
- Session save/load (JSON v2 format + data sidecar).
- Absolute quantification tab deferred to a future release (placeholder).

## 1.0.0 (2026-05)

- First public-ready version: annotation and filtering of MS-DIAL alignment
  tables with LipidQuant relay export.
