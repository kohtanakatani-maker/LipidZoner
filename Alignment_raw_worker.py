#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# LipidZoner — SFC/MS lipidomics data processing GUI
# Copyright (C) 2026 Kohta Nakatani
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Commercial licensing (use outside the terms of the GPL) is available
# from the copyright holder.
"""Alignment raw mzML peak detection worker (subprocess).

Invoked by LipidZoner.py via subprocess to process .mzML files
with pyOpenMS. Runs in its own Python process so that pyOpenMS's bundled
Qt5 DLLs do NOT conflict with PySide6's Qt6 DLLs in the main GUI process.

v0.2.4 fix: pyOpenMS's MzMLFile().load() cannot read paths containing
non-ASCII characters (e.g., Japanese) on Windows due to a C++-layer
limitation. As a workaround, this worker copies the input mzML to an
ASCII-only tempfile path before loading.

v0.5 fix (2026-05-31): centroiding stage added before MassTraceDetection.
  Root cause of the "peak splitting" problem: the input mzML is in
  PROFILE mode (cvParam MS:1000128). When profile spectra are fed directly
  to MassTraceDetection, the TOF profile sampling points (~5 mDa apart at
  m/z 753) are picked up as separate mass traces, so a single chemical
  entity is split into several m/z traces. Confirmed against real data:
  PC 33:1 D7 (IS, theor. 753.6134) at RT 2.36 min produced 4 traces
  (753.6093 / .6143 / .6193 / .6244) from profile-direct detection, but
  PeakPickerHiRes centroiding collapses them to a single peak.
  Fix: run PeakPickerHiRes (profile -> centroid) before MassTraceDetection.
  This matches the standard pyOpenMS metabolomics workflow and what MS-DIAL
  does internally. Controlled by --centroid {auto,on,off} (default auto).
  Spectrum-type detection uses spec.getType() (reliable, cvParam-based);
  PeakTypeEstimator.estimateType() is NOT used (it misjudged profile data
  as centroid in testing).

Usage (called from main process):
  python Alignment_raw_worker.py \
      --mzml <input.mzML> --out <output.csv> \
      --mass-error-ppm 20 --noise 1000 --snr 3 --fwhm 5 --centroid auto

Or to just check pyOpenMS availability:
  python Alignment_raw_worker.py --check

Gap filling (v0.6, 2026-10-03):
  python Alignment_raw_worker.py --fill \
      --mzml <input.mzML> --targets <targets.csv> --out <filled.csv> \
      --fill-ppm 10 --fill-rt-half 0.035

  targets.csv : feature_id, mz, rt[, rt_half]
  filled.csv  : feature_id, height, area, n_points, rt_apex, source

  MS-DIAL の "Gap filling by compulsion"(既定 ON)に相当する。極大が無くても
  指定窓で積分する。既存の --eic と違い mzML のロードは 1 回だけなので、
  数千ターゲットでも 1 ファイル 2 秒程度で終わる。
"""
import sys
import argparse
import traceback


def check_pyopenms_available():
    """Try importing pyOpenMS. Print version on success, error on failure."""
    try:
        import pyopenms as oms
        print(oms.__version__)
        return 0
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1


# Spectrum type enum (pyOpenMS SpectrumSettings::SpectrumType):
#   0 = UNKNOWN, 1 = CENTROID, 2 = PROFILE
_SPEC_TYPE_NAME = {0: 'UNKNOWN', 1: 'CENTROID', 2: 'PROFILE'}


def _first_ms1_spectrum_type(exp):
    """Return the integer spectrum type of the first MS1 spectrum, or None."""
    for i in range(exp.size()):
        s = exp[i]
        if s.getMSLevel() == 1:
            return int(s.getType())
    return None


def maybe_centroid(exp, mode):
    """Centroid profile-mode spectra with PeakPickerHiRes.

    Run BEFORE MassTraceDetection. Profile-mode mzML (e.g. SFC-QTOF raw
    converted by ProteoWizard) must be centroided first, otherwise the
    profile sampling points are detected as separate mass traces and one
    chemical entity is split into several m/z traces.

    mode:
      'off'  -> do nothing, return exp unchanged.
      'on'   -> force-pick every spectrum with PeakPickerHiRes.
      'auto' -> (default) skip if spectra are already centroided; otherwise
                let PeakPickerHiRes pick only PROFILE spectra
                (check_spectrum_type=True, which reads spec.getType()).

    Returns (exp_out, did_centroid: bool).
    """
    import pyopenms as oms

    if mode == 'off':
        print("Centroiding: mode=off (skipped)", file=sys.stderr)
        return exp, False

    first_type = _first_ms1_spectrum_type(exp)
    type_name = _SPEC_TYPE_NAME.get(first_type, str(first_type))
    print(f"Centroiding: mode={mode}, first MS1 spectrum type={type_name}",
          file=sys.stderr)

    if mode == 'auto' and first_type == 1:
        # Already centroided -> nothing to do.
        print("  Spectra already centroided; skipping PeakPickerHiRes.",
              file=sys.stderr)
        return exp, False

    pp = oms.PeakPickerHiRes()
    pp_params = pp.getDefaults()
    # Disable S/N filtering at the picking stage; downstream MassTraceDetection
    # handles noise via noise_threshold_int. This keeps low-intensity real
    # peaks instead of dropping them during picking.
    try:
        pp_params.setValue('signal_to_noise', 0.0)
    except Exception:
        pass
    pp.setParameters(pp_params)

    centroided = oms.MSExperiment()
    # check_spectrum_type=True for 'auto': PeakPickerHiRes inspects each
    # spectrum's type and only picks PROFILE spectra (passes centroid ones
    # through unchanged). For 'on': force-pick all (check=False).
    check_type = (mode == 'auto')
    pp.pickExperiment(exp, centroided, check_type)
    print(f"  PeakPickerHiRes done: {exp.size()} -> {centroided.size()} "
          f"spectra (profile -> centroid).", file=sys.stderr)
    return centroided, True


def run_peak_detection(args):
    """Run pyOpenMS centroiding + Stage 1+2 peak detection on a single mzML."""
    try:
        import pyopenms as oms
        import pandas as pd
    except Exception as e:
        print(f"Import error: {e}", file=sys.stderr)
        return 1

    # Load mzML.
    # pyOpenMS on Windows cannot handle non-ASCII paths (e.g., Japanese).
    # Workaround: if the source path contains non-ASCII characters,
    # copy to an ASCII-only tempfile first.
    import os as _os
    import shutil as _shutil
    import tempfile as _tempfile
    src_path = args.mzml
    use_temp = False
    try:
        src_path.encode('ascii')
    except UnicodeEncodeError:
        use_temp = True

    if use_temp:
        # Copy to ASCII tempfile (in system temp dir, ASCII-only path)
        tmp_dir = _tempfile.gettempdir()
        tmp_basename = "alignment_mzml_temp_%d.mzML" % _os.getpid()
        load_path = _os.path.join(tmp_dir, tmp_basename)
        print(f"  Source path has non-ASCII characters; "
              f"copying to ASCII tempfile: {load_path}", file=sys.stderr)
        _shutil.copy2(args.mzml, load_path)
    else:
        load_path = args.mzml

    try:
        print(f"Loading mzML: {load_path}", file=sys.stderr)
        exp = oms.MSExperiment()
        oms.MzMLFile().load(load_path, exp)
        print(f"  {exp.size()} spectra loaded", file=sys.stderr)
    finally:
        if use_temp:
            try:
                _os.unlink(load_path)
                print(f"  Cleaned up tempfile: {load_path}", file=sys.stderr)
            except Exception as _e:
                print(f"  Warning: failed to clean tempfile: {_e}",
                      file=sys.stderr)

    # Stage 0: Centroiding (profile -> centroid) if needed.
    exp, did_centroid = maybe_centroid(exp, args.centroid)

    # Stage 1: MassTraceDetection
    print(f"Stage 1: MassTraceDetection (mass_error={args.mass_error_ppm} ppm, "
          f"noise={args.noise})", file=sys.stderr)
    mtd = oms.MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue('mass_error_ppm', args.mass_error_ppm)
    mtd_params.setValue('noise_threshold_int', args.noise)
    mtd_params.setValue('chrom_peak_snr', args.snr)
    # v0.5.16 fix: pyOpenMS default min_trace_length=5.0s discards SFC's narrow
    #   (~4s) peaks -> whole classes (ceramides, ether lipids) were lost.
    mtd_params.setValue('min_trace_length', args.min_trace_length)
    mtd.setParameters(mtd_params)
    mass_traces = []
    mtd.run(exp, mass_traces, 1000000)
    print(f"  {len(mass_traces)} mass traces", file=sys.stderr)

    # Stage 2: ElutionPeakDetection
    print(f"Stage 2: ElutionPeakDetection (chrom_fwhm={args.fwhm} sec)",
          file=sys.stderr)
    epd = oms.ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue('chrom_fwhm', args.fwhm)
    epd_params.setValue('chrom_peak_snr', args.snr)
    epd.setParameters(epd_params)
    out_traces = []
    epd.detectPeaks(mass_traces, out_traces)
    print(f"  {len(out_traces)} elution peaks", file=sys.stderr)

    # Compute FWHM
    for mt in out_traces:
        mt.estimateFWHM(False)

    # Build DataFrame
    records = []
    for i, mt in enumerate(out_traces):
        cmz = float(mt.getCentroidMZ())
        # v0.5.9: per-trace m/z scatter (instrument mass precision), in ppm.
        #   getCentroidSD() = stdev of the trace's m/z across its scans.
        try:
            mz_sd = float(mt.getCentroidSD())
        except Exception:
            mz_sd = 0.0
        mz_sd_ppm = (mz_sd / cmz * 1.0e6) if cmz > 0 else 0.0
        records.append({
            'Peak ID':       i,
            'Precursor m/z': cmz,
            'RT (min)':      float(mt.getCentroidRT()) / 60.0,
            'Height':        float(mt.getMaxIntensity(False)),
            'Area':          float(mt.computePeakArea()),
            'fwhm_sec':      float(mt.getFWHM()),
            'mz_sd_ppm':     mz_sd_ppm,
            'Isotope':       0,
            'Adduct':        '',
        })
    df = pd.DataFrame(records)
    df.to_csv(args.out, index=False)
    print(f"Wrote {len(df)} peaks to {args.out} "
          f"(centroided={did_centroid})", file=sys.stderr)
    return 0


def extract_eic(args):
    """Extract XIC (Extracted Ion Chromatogram) for a given m/z window.

    Note: EIC is computed directly from the loaded spectra (profile or
    centroid) by taking the max intensity in the m/z window per MS1 scan,
    so it does NOT require centroiding.
    """
    try:
        import pyopenms as oms
        import pandas as pd
    except Exception as e:
        print(f"Import error: {e}", file=sys.stderr)
        return 1
    import os as _os, shutil as _shutil, tempfile as _tempfile

    # Same ASCII-path workaround as run_peak_detection
    src_path = args.mzml
    use_temp = False
    try:
        src_path.encode('ascii')
    except UnicodeEncodeError:
        use_temp = True
    if use_temp:
        tmp_dir = _tempfile.gettempdir()
        tmp_basename = "alignment_eic_temp_%d.mzML" % _os.getpid()
        load_path = _os.path.join(tmp_dir, tmp_basename)
        _shutil.copy2(args.mzml, load_path)
    else:
        load_path = args.mzml
    try:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(load_path, exp)
        print(f"  Loaded {exp.size()} spectra for EIC at m/z={args.eic_mz} "
              f"+/-{args.eic_mz_tol} Da", file=sys.stderr)

        # Extract intensities at the target m/z (+/- tol) for each MS1 scan
        records = []
        mz_lo = args.eic_mz - args.eic_mz_tol
        mz_hi = args.eic_mz + args.eic_mz_tol
        for i in range(exp.size()):
            spec = exp[i]
            if spec.getMSLevel() != 1:
                continue
            rt_min = spec.getRT() / 60.0
            mzs, intensities = spec.get_peaks()
            mask = (mzs >= mz_lo) & (mzs <= mz_hi)
            if mask.any():
                inten = float(intensities[mask].max())
            else:
                inten = 0.0
            records.append({'rt_min': rt_min, 'intensity': inten})

        df = pd.DataFrame(records)
        df.to_csv(args.out, index=False)
        print(f"  Wrote {len(df)} EIC points to {args.out}", file=sys.stderr)
    finally:
        if use_temp:
            try: _os.unlink(load_path)
            except Exception: pass
    return 0


def fill_gaps(args):
    """ターゲット一覧の m/z x RT 窓を 1 パスで積分する。

    MS-DIAL の "Gap filling by compulsion" に相当。クロマトグラムに極大が
    無くても積分する。

    height の採り方は「窓内で RT が最も近い極大」。極大が 1 つも無いときだけ
    窓内の最大値に落とす（source 列で 'apex' / 'forced' と区別できる）。

    単純に窓内の最大を採ると、10〜40 倍強い隣のピークが 0.05 分ほど離れて
    並んでいる場合にその裾を拾ってしまう（実測 3 件で真値の 2.1〜3.1 倍）。
    最近接の極大にすると 0.7〜1.0 倍に収まり、埋めた値が真値の 3 倍を超える
    件数は 2,251 件中 14 → 2 件になった。

    targets.csv の列:
      feature_id : 呼び出し側の ID(そのまま返す)
      mz         : 中心 m/z
      rt         : 中心 RT(分)
      rt_half    : 片側の窓幅(分)。省略時は --fill-rt-half

    m/z の許容幅は --fill-ppm(ppm)。アノテーション側と同じ単位にしてある。
    絶対 Da 固定にすると高 m/z で窓が狭くなりすぎる(本体のアノテーション側と同じ理由)。

    窓幅について: MS-DIAL は「その feature を持つ検体の平均ピーク幅」を使う。
    呼び出し側がそれを rt_half として渡す想定。検証データ(chrom_fwhm 2.2 秒
    = 0.037 分)では ±0.035 分が最適だった。±0.05 分まで広げると、RT が
    0.05 分差で並ぶ DG の sn 位置異性体で隣のピークを拾い始める。
    """
    try:
        import pyopenms as oms
        import numpy as np
        import pandas as pd
    except Exception as e:
        print(f"Import error: {e}", file=sys.stderr)
        return 1
    import os as _os, shutil as _shutil, tempfile as _tempfile, time as _time

    try:
        tg = pd.read_csv(args.targets)
    except Exception as e:
        print(f"Cannot read targets: {e}", file=sys.stderr)
        return 1
    for c in ('feature_id', 'mz', 'rt'):
        if c not in tg.columns:
            print(f"targets.csv needs column '{c}'", file=sys.stderr)
            return 1
    if 'rt_half' not in tg.columns:
        tg['rt_half'] = float(args.fill_rt_half)
    tg['rt_half'] = (tg['rt_half'].astype(float)
                     .fillna(float(args.fill_rt_half)))
    tg.loc[tg['rt_half'] <= 0, 'rt_half'] = float(args.fill_rt_half)

    # ASCII-path workaround (run_peak_detection と同じ)
    src_path = args.mzml
    use_temp = False
    try:
        src_path.encode('ascii')
    except UnicodeEncodeError:
        use_temp = True
    if use_temp:
        tmp_dir = _tempfile.gettempdir()
        tmp_basename = "alignment_fill_temp_%d.mzML" % _os.getpid()
        load_path = _os.path.join(tmp_dir, tmp_basename)
        _shutil.copy2(args.mzml, load_path)
    else:
        load_path = args.mzml

    try:
        t0 = _time.time()
        exp = oms.MSExperiment()
        oms.MzMLFile().load(load_path, exp)
        print(f"  Loaded {exp.size()} spectra in {_time.time()-t0:.1f}s",
              file=sys.stderr)

        # MS1 を 1 回だけデコードして保持(ターゲットごとの再デコードを避ける)
        t0 = _time.time()
        rts = []
        mz_arrays = []
        int_arrays = []
        for i in range(exp.size()):
            sp = exp[i]
            if sp.getMSLevel() != 1:
                continue
            m, it = sp.get_peaks()
            rts.append(sp.getRT() / 60.0)
            mz_arrays.append(m)
            int_arrays.append(it)
        if not rts:
            print("No MS1 spectra found", file=sys.stderr)
            pd.DataFrame(columns=['feature_id', 'height', 'area',
                                  'n_points', 'rt_apex',
                                  'source']).to_csv(
                args.out, index=False)
            return 0
        rts = np.asarray(rts, dtype=float)
        order = np.argsort(rts)
        rts = rts[order]
        mz_arrays = [mz_arrays[k] for k in order]
        int_arrays = [int_arrays[k] for k in order]
        print(f"  Decoded {len(rts)} MS1 scans in {_time.time()-t0:.1f}s",
              file=sys.stderr)

        t0 = _time.time()
        fid = tg['feature_id'].to_numpy()
        tmz = tg['mz'].to_numpy(dtype=float)
        trt = tg['rt'].to_numpy(dtype=float)
        thalf = tg['rt_half'].to_numpy(dtype=float)
        ppm = float(args.fill_ppm)
        n = len(tg)
        out_h = np.zeros(n)
        out_a = np.zeros(n)
        out_n = np.zeros(n, dtype=int)
        out_rt = np.full(n, np.nan)
        out_src = np.empty(n, dtype=object)
        out_src[:] = 'none'
        PAD = 2      # 窓の端の点も極大か判定できるよう左右に 2 スキャン足す
        n_scan = len(rts)
        for t in range(n):
            tol = tmz[t] * ppm * 1e-6
            w0 = int(np.searchsorted(rts, trt[t] - thalf[t]))
            w1 = int(np.searchsorted(rts, trt[t] + thalf[t]))
            if w1 <= w0:
                continue
            k0 = max(0, w0 - PAD)
            k1 = min(n_scan, w1 + PAD)
            ys = np.empty(k1 - k0, dtype=float)
            lo = tmz[t] - tol
            hi = tmz[t] + tol
            for q, k in enumerate(range(k0, k1)):
                m = mz_arrays[k]
                j0 = np.searchsorted(m, lo)
                j1 = np.searchsorted(m, hi)
                ys[q] = int_arrays[k][j0:j1].max() if j1 > j0 else 0.0
            # 窓そのものの範囲(pad を除いた部分)
            a = w0 - k0
            b = w1 - k0
            out_n[t] = b - a
            if b <= a:
                continue
            # 窓内の極大(両隣以上、かつ 0 より大きい)
            seg = ys[a:b]
            left = ys[a - 1:b - 1] if a >= 1 else np.r_[-1.0, ys[a:b - 1]]
            right = (ys[a + 1:b + 1] if b + 1 <= len(ys)
                     else np.r_[ys[a + 1:b], -1.0])
            is_max = (seg > 0) & (seg >= left) & (seg >= right)
            if is_max.any():
                cand = np.flatnonzero(is_max)
                # RT が目的値に最も近いものを選ぶ
                pick = cand[int(np.argmin(np.abs(rts[k0 + a + cand] - trt[t])))]
                out_h[t] = float(seg[pick])
                out_rt[t] = float(rts[k0 + a + pick])
                out_src[t] = 'apex'
            else:
                # MS-DIAL の compulsion: 極大が無くても窓内の最大で埋める
                pick = int(np.argmax(seg))
                out_h[t] = float(seg[pick])
                out_rt[t] = float(rts[k0 + a + pick])
                out_src[t] = 'forced'
            # 面積は秒換算の台形積分(検出側の computePeakArea と桁を揃える)
            if b - a > 1:
                try:
                    out_a[t] = float(
                        np.trapezoid(seg, rts[k0 + a:k0 + b]) * 60.0)
                except AttributeError:
                    out_a[t] = float(np.trapz(seg, rts[k0 + a:k0 + b]) * 60.0)
        print(f"  Filled {n} targets in {_time.time()-t0:.1f}s "
              f"({int((out_h > 0).sum())} non-zero, "
              f"{int((out_src == 'forced').sum())} without a local maximum)",
              file=sys.stderr)

        pd.DataFrame({'feature_id': fid, 'height': out_h, 'area': out_a,
                      'n_points': out_n, 'rt_apex': out_rt,
                      'source': out_src}).to_csv(
            args.out, index=False)
        print(f"  Wrote {n} rows to {args.out}", file=sys.stderr)
    finally:
        if use_temp:
            try:
                _os.unlink(load_path)
            except Exception:
                pass
    return 0


def extract_tic(args):
    """Extract Total Ion Chromatogram (sum intensity per scan)."""
    try:
        import pyopenms as oms
        import pandas as pd
    except Exception as e:
        print(f"Import error: {e}", file=sys.stderr)
        return 1
    import os as _os, shutil as _shutil, tempfile as _tempfile
    src_path = args.mzml
    use_temp = False
    try: src_path.encode('ascii')
    except UnicodeEncodeError: use_temp = True
    if use_temp:
        tmp_dir = _tempfile.gettempdir()
        tmp_basename = "alignment_tic_temp_%d.mzML" % _os.getpid()
        load_path = _os.path.join(tmp_dir, tmp_basename)
        _shutil.copy2(args.mzml, load_path)
    else:
        load_path = args.mzml
    try:
        exp = oms.MSExperiment()
        oms.MzMLFile().load(load_path, exp)
        records = []
        for i in range(exp.size()):
            spec = exp[i]
            if spec.getMSLevel() != 1:
                continue
            rt_min = spec.getRT() / 60.0
            mzs, intensities = spec.get_peaks()
            tic = float(intensities.sum()) if len(intensities) > 0 else 0.0
            records.append({'rt_min': rt_min, 'tic': tic})
        df = pd.DataFrame(records)
        df.to_csv(args.out, index=False)
        print(f"  Wrote {len(df)} TIC points to {args.out}", file=sys.stderr)
    finally:
        if use_temp:
            try: _os.unlink(load_path)
            except Exception: pass
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true',
                        help="Check pyOpenMS availability and print version")
    parser.add_argument('--eic', action='store_true',
                        help="Extract EIC (XIC) at a given m/z instead of "
                             "peak detection")
    parser.add_argument('--tic', action='store_true',
                        help="Extract Total Ion Chromatogram instead of "
                             "peak detection")
    parser.add_argument('--fill', action='store_true',
                        help="Gap filling: integrate a list of m/z x RT "
                             "windows in one pass (MS-DIAL's 'Gap filling "
                             "by compulsion')")
    parser.add_argument('--targets',
                        help="--fill: CSV with feature_id, mz, rt[, rt_half]")
    parser.add_argument('--fill-ppm', type=float, default=10.0,
                        help="--fill: m/z window in ppm. Default 10")
    parser.add_argument('--fill-rt-half', type=float, default=0.035,
                        help="--fill: half RT window in minutes, used when "
                             "targets.csv has no rt_half column. Default "
                             "0.035 (about one chromatographic peak width "
                             "for the SFC method)")
    parser.add_argument('--eic-mz', type=float,
                        help="Target m/z for EIC extraction")
    parser.add_argument('--eic-mz-tol', type=float, default=0.005,
                        help="m/z tolerance for EIC (Da). Default 0.005")
    parser.add_argument('--mzml', help="Input .mzML path")
    parser.add_argument('--out', help="Output CSV path")
    parser.add_argument('--mass-error-ppm', type=float, default=20.0)
    parser.add_argument('--noise', type=float, default=1000.0)
    parser.add_argument('--snr', type=float, default=3.0)
    parser.add_argument('--fwhm', type=float, default=5.0)
    parser.add_argument('--min-trace-length', type=float, default=3.0)
    parser.add_argument('--centroid', choices=['auto', 'on', 'off'],
                        default='auto',
                        help="Centroid PROFILE-mode spectra with "
                             "PeakPickerHiRes before mass-trace detection. "
                             "auto (default): pick only PROFILE spectra "
                             "(skip if already centroided); on: force-pick "
                             "all; off: disable centroiding.")
    args = parser.parse_args()

    if args.check:
        return check_pyopenms_available()
    if args.eic:
        if not args.mzml or not args.out or args.eic_mz is None:
            parser.error("--eic requires --mzml, --out, and --eic-mz")
        return extract_eic(args)
    if args.tic:
        if not args.mzml or not args.out:
            parser.error("--tic requires --mzml and --out")
        return extract_tic(args)
    if args.fill:
        if not args.mzml or not args.out or not args.targets:
            parser.error("--fill requires --mzml, --targets and --out")
        return fill_gaps(args)
    if not args.mzml or not args.out:
        parser.error("--mzml and --out are required "
                     "(or use --check / --eic / --fill)")
    return run_peak_detection(args)


if __name__ == '__main__':
    sys.exit(main())
