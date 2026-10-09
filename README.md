# LipidZoner

**LipidZoner** is a desktop GUI for processing SFC/MS (and LC/MS) untargeted
lipidomics data: peak detection and alignment from raw mzML files, library-based
lipid annotation, a chain of filters that removes false annotations, selection
of the quantitative ion per lipid class, and export of per-class peak tables
(including the input format of [LipidQuant](https://holcapek.upce.cz/lipidquant)).

Version 1.1.0 is the version used in the accompanying manuscript.

> Note on the *Quantification* tab: absolute quantification inside LipidZoner is
> planned for a future release and the tab is currently a placeholder. Export the
> annotated peak tables (Utility → Reports / Export, or the LipidQuant relay
> export) and quantify externally.

## Features

- **Detection tab** – raw mzML → centroiding → mass-trace / elution-peak
  detection (pyOpenMS, run in a worker process) → cross-sample alignment,
  isotope grouping, gap filling and optional RT drift correction.
  Parameters (mass accuracy, noise, S/N, FWHM) can be estimated automatically
  from the data.
- **Annotation tab** – positive and negative modes side by side:
  1. Match (library m/z search, ppm tolerance)
  2. Conflict search (one feature, several candidates)
  3. Internal-standard (IS) filter – class RT windows anchored on the IS
  4. RT outlier filter
  5. Adduct ion filter – IS-derived adduct fingerprints
  6. Coherence filter – intra-class RT regression (carbon number / double
     bonds / linkage) to resolve inter-class conflicts, with confidence scores
  7. Quant ion selection – choose pos or neg per class
  - interactive RT × m/z scatter plots, EIC viewer, candidate picker, manual
    accept/reject with full provenance (`final_status` per candidate)
- **Utility tab** – session save/load (JSON + data sidecar), import of
  externally produced tables (MS-DIAL per-sample peak lists or alignment
  results), reports (filter summaries, adduct heatmaps, scatter PNGs),
  quant-ion table export (height/area per sample), export of *all* features
  including non-annotated ones, and the LipidQuant relay export.

## Requirements

- Python 3.10 or newer (tested with 3.10 and 3.13 on Windows 11 and Linux)
- Packages listed in `requirements.txt`:
  PySide6, numpy, pandas, matplotlib, openpyxl, chardet, and **pyopenms**
  (needed only for raw mzML processing; without it the external-import route
  still works).

```bash
python -m pip install -r requirements.txt
```

## Running

```bash
python LipidZoner.py
```

`Alignment_raw_worker.py` must stay in the same directory as `LipidZoner.py`.
pyOpenMS ships its own Qt5 libraries and cannot live in the same process as
PySide6 (Qt6), so peak detection runs in a separate Python process started
from the GUI.

## Typical workflow

1. **Detection** → Step 1: add the mzML files of one polarity (positive and
   negative are handled as two columns). Step 2: peak-detection parameters
   (or *Estimate*). Step 3: alignment / isotope-grouping parameters. Run.
   Send the aligned table to Annotation.
2. **Annotation** → load the lipid library (`.xlsm`/`.xlsx`, see below), set
   the ppm tolerance, then *Run All* (①–⑥). Review conflicts in the scatter
   plot, pick the quant ion per class (⑦).
3. **Utility → Reports / Export** → export the quant-ion table (xlsx), the
   non-target feature table, PNG figures, or the LipidQuant relay files.
4. **Utility → Session** → save the session (JSON) to resume later; the
   session stores parameters, IS choices, manual decisions and the aligned
   data.

## Lipid library format

One worksheet per lipid class (sheet name = class). Row 2 is the header, data
start at row 3. Columns (0-based):

| col | content |
|----:|---------|
| 0 | compound name (e.g. `PC 34:1`) |
| 1–8 | element counts C, H, O, N, P, S, D, ¹³C |
| 11 | positive-mode adduct for the class (header cell holds the adduct, e.g. `[M+H]+`; cell = 1 to enable) |
| 12 | negative-mode adduct (same convention) |
| 13 | `1` if the row is an internal standard |
| 14 | name of the IS assigned to this species (optional) |

An example library is provided in `examples/`.

## Repository layout

```
LipidZoner.py              main application (single file)
Alignment_raw_worker.py    pyOpenMS worker (peak detection, EIC, gap filling)
LipidZoner_icon.ico        application icon
requirements.txt
examples/                  example lipid library
CITATION.cff               how to cite
LICENSE                    GPL-3.0-or-later
```

## Citing

If you use LipidZoner, please cite the Zenodo record (DOI in `CITATION.cff`)
and the manuscript once published.

## License

LipidZoner is released under the **GNU General Public License v3.0 or later**
(see `LICENSE`). You are free to use, study, modify and redistribute it,
provided that derivative works distributed to others are also released under
the GPL with source code.

**Commercial licensing.** If you wish to incorporate LipidZoner into a
product or otherwise use it under terms other than the GPL, a separate
commercial license is available from the copyright holder. Please contact
the author.
