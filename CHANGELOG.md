# Changelog

## 1.1.0 (2026-10-21) – manuscript release

Supersedes the 2026-10-08 build of 1.1.0, which was never published. The
annotation results change; see *Filtering and isotope grouping* below.

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

### Filtering and isotope grouping

- **IS adduct fingerprints are now taken from the row the IS filter kept.**
  Previously the fingerprint was built from the first matching row
  (smallest |Δm/z|). When a class had several IS candidates, that row could
  be one the IS filter had already rejected, which shifted the fingerprint's
  reference RT by up to 1.6 min and dropped the genuine minor adducts. A
  class whose only expected adduct is shared with another class then lost
  every conflict in the adduct ion filter.
- **Fingerprint entries with no usable intensity are no longer marked as
  expected.** A weight of 0 made the ratio method estimate a contribution of
  exactly 0 while the count method still voted, which is inconsistent.
- **The adduct ion filter no longer rejects candidates it could not score.**
  Only classes that actually received a score compete; a class with no
  fingerprint is left undecided and handed to the coherence filter instead
  of being dropped.
- **New: minimum IS intensity for fingerprints** (Advanced Setting ③, two
  absolute-count thresholds, both 0 = disabled). *IS fp adduct min* drops
  individual adducts below the threshold; *IS fp class min* drops the whole
  class if its strongest remaining adduct is below it. Classes without a
  fingerprint are left undecided rather than rejected. Stored in the session
  as `is_fp_min_adduct_int` / `is_fp_min_class_int`.
- **Isotope grouping now compares against the root of the isotope ladder**
  (the true M+0) instead of the immediate neighbour. A neighbour that was
  itself flagged as co-eluting still carries the intensity of the isotopes
  below it, so predicting from it overestimates and the "intensity above the
  predicted isotope" rule never fired. On the verification data it fired for
  0 of 38 known co-eluting species before and 27 of 38 after. This is what
  separates a saturated species (X n:0, +2.01565 Da) from the M+2 isotope of
  X n:1 (+2.00671 Da) — only 8.94 mDa apart, inside the 0.006 Da m/z
  tolerance at these masses.

### Compatibility

- `Alignment_raw_worker.py` now works with both pyOpenMS API shapes.
  pyOpenMS ≥ 3.6 (pybind11) returns results, while earlier builds
  (Cython/autowrap) write into an output argument, and the two raise
  different exceptions when called the wrong way (`TypeError` vs a bare
  `Exception`). The worker tries the new shape, falls back to the old one,
  and re-raises the original error if both fail, so a genuine failure is not
  hidden by the version probe. Which shape was used is logged once per
  function.
