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
"""LipidZoner — SFC/MS lipidomics data processing GUI.

Peak detection and alignment (raw mzML, via pyOpenMS in a worker process),
library-based lipid annotation, multi-stage filtering (internal standard,
adduct ion, RT outlier, coherence), quant-ion selection
and export of per-class peak tables.

Pipeline (Annotation tab):
    1 Match  ->  2 Conflict search  ->  3 IS filter  ->  4 RT outlier filter
    ->  5 Adduct ion filter  ->  6 Coherence filter  ->  7 Select quant ion
    ->  Export

Run:  python LipidZoner.py
Requires Alignment_raw_worker.py in the same directory for raw mzML
processing (pyOpenMS bundles Qt5 and cannot share a process with PySide6).

License: GPL-3.0-or-later (see LICENSE).
"""


import json
import logging
import os
import sys
import re
import glob
import subprocess as _subprocess
import tempfile as _tempfile
import time as _time
import datetime
from pathlib import Path

import chardet
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.colors as mcolors
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import (QAction, QBrush, QColor, QCursor, QFont,
                            QIcon, QKeySequence, QPixmap)
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog,
    QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QFrame,
    QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget,
    QListWidgetItem, QMainWindow,
    QMenu, QMenuBar, QMessageBox, QProgressDialog, QPushButton,
    QScrollArea, QSizePolicy,
    QSpinBox, QRadioButton, QSplitter, QTabWidget, QTableWidget,
    QTableWidgetItem,
    QTextEdit, QToolButton, QVBoxLayout, QGridLayout, QWidget,
)

# ヒートマップ共通ウィジェット(LipidIonAudit と共通)
# DEFAULT_ADDUCT_COLUMNS は前方参照される(_compute_is_fingerprints 等)ため
# ここで先に定義。AdductHeatmapWidget クラスは後段(AdductPatternsDialog
# の直前)で定義する。
DEFAULT_ADDUCT_COLUMNS = (
    '[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+', '[M-H2O+H]+',
    '[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+CH3OCOO]-',
)

log = logging.getLogger("lipidzoner")

matplotlib.rcParams['font.family'] = ['Arial', 'DejaVu Sans']
matplotlib.rcParams['font.size']   = 9

__version__ = "1.1.0"

# ─── MS-DIAL 列分類キーワード ─────────────────────────────────────────
META_KEYWORDS = [
    "alignment id", "average rt", "average rt(min)", "average mz",
    "metabolite name", "adduct type", "post curation result",
    "fill %", "ms/ms assigned", "reference rt", "reference m/z",
    "formula", "ontology", "inchikey", "smiles",
    "annotation tag", "rt matched", "m/z matched", "ms/ms matched",
    "comment", "manually modified", "isotope tracking parent id",
    "isotope tracking weight number", "total score",
    "rt similarity", "m/z similarity",
    "simple dot product", "weighted dot product",
    "reverse dot product", "matched peaks count",
    "matched peaks percentage", "fragment presence %",
    "s/n average", "spectrum reference file name",
    "ms1 isotopic spectrum", "ms/ms spectrum",
    "snr", "s/n", "signal/noise",
]
STAT_KEYWORDS = ["average", "stdev", "cv"]
HEADER_MARKER_KEYWORDS = [
    "alignment id", "average rt", "average rt(min)", "average mz",
    "metabolite name", "adduct type", "sample", "area", "height",
]

# ─── 脂質ライブラリ定数 ───────────────────────────────────────────────
MONO_MASSES = {
    'C': 12.000000000, 'H': 1.00782503207, 'D': 2.01410177785,
    'N': 14.0030740048, 'O': 15.99491461956, 'P': 30.97376163,
    'S': 31.97207100,  'F': 18.99840322,   'Cl': 34.96885271,
    'Br': 78.91833710, 'Na': 22.98976928,  'K':  38.96370668,
    '13C': 13.003354835, '18O': 17.99915961286,
}
ADDUCT_OFFSETS = {
    # pos
    '[M+H]+':         1.007276,
    '[M+Na]+':       22.989218,
    '[M+NH4]+':      18.034164,
    '[M+K]+':        38.963158,
    '[M-H2O+H]+':   -17.003289,
    # neg
    '[M-H]-':        -1.007276,
    '[M+HCOO]-':     44.998201,
    '[M+CH3COO]-':   59.013304,
    '[M+CH3OCOO]-':  75.008768,   # methyl carbonate adduct
}

# Adduct Ion Filter で扱う 9 アダクトの集合
ADDUCT_SET_POS_V54 = ('[M+H]+', '[M+Na]+', '[M+NH4]+', '[M+K]+', '[M-H2O+H]+')
ADDUCT_SET_NEG_V54 = ('[M-H]-', '[M+HCOO]-', '[M+CH3COO]-', '[M+CH3OCOO]-')
ADDUCT_SET_V54 = ADDUCT_SET_POS_V54 + ADDUCT_SET_NEG_V54
# ⑤ の強度比判定のパラメータ
#   RATIO_WIN_FRAC   … 寄与率がこれ以上なら勝者（未満は混合として未決）
#   RATIO_CONSIST_*  … 推定寄与の合計が実測ピーク強度から外れすぎたら
#                      強度比モデルが成立していないとみなし、本数ルールへ
# pos↔neg RT シフトの σ の下限。IS の実測ばらつきが小さく出ても、
# アライメントのゆらぎを考えるとクロスモード窓をこれ以下にはしない。
POS_NEG_SHIFT_SIGMA_FLOOR = 0.005
ADDUCT_RATIO_WIN_FRAC   = 0.8
ADDUCT_RATIO_CONSIST_LO = 0.1
ADDUCT_RATIO_CONSIST_HI = 10.0
# ライブラリ内の表記ゆれ → 正規化
ADDUCT_ALIAS = {
    '[M+OAc]-':    '[M+CH3COO]-',
    '[M+H-H2O]+':  '[M-H2O+H]+',
}


# ─── タブラベル ────────────────────────────────────────
# AP のトップレベルタブ名を 1 箇所に集約する。日本語表記にしたい場合は
# ここだけ書き換えればよい(例: 'annot': "アノテーション")。
#   align : ピーク検出 + アライメント(表示名は "Detection")
#   annot : アノテーション (旧 Analysis タブ)
#   quant : 定量(将来追加予定のプレースホルダ)
#   util  : Session / Reports / relay export を束ねるサブタブ container
TAB_LABELS = {
    'align': "Detection",
    'annot': "Annotation",
    'quant': "Quantification (TBA)",
    'util':  "Utility",
}
# Utility 内サブタブ
SUBTAB_LABELS = {
    'session': "Session",
    'import':  "External data import",
    'reports': "Reports / Export",
    'relay':   "LipidQuant relay export",
}


# ─── 複数ラベル方式の status 値 ────────────────────────────
# match_df の各行は 1 つの (ライブラリエントリ, 候補ピーク) ペアを表す。
# 4 つの status 列(is_filter / rt_outlier / coherence / manual)が独立して
# 立つ rejection をトラッキングし、final_status はそれらの AND で決まる。
# match_library_to_peaks が multi_label=True で多候補を返す。各 status 列は
# デフォルト値で初期化され、IS Filter / RT Outlier / Coherence の mark で
# 更新される。

STATUS_KEPT          = 'kept'
STATUS_NA            = 'n/a'
STATUS_REJ_WINDOW    = 'rejected_window'      # IS Filter
STATUS_REJ_ADDUCT    = 'rejected_adduct'      # Adduct Ion Filter
STATUS_REJ_OUTLIER   = 'rejected_outlier'     # RT Outlier Filter
STATUS_REJ_RESIDUAL  = 'rejected_residual'    # Coherence Filter intra-compound
STATUS_REJ_LOSER     = 'rejected_loser'       # Coherence Filter inter-class
STATUS_REJ_MANUAL    = 'rejected_manual'      # 手動
# ノンターゲット表で「そもそもライブラリに当たらなかった」印。
# STATUS_NA(= 当たったが判定前)とは意味が違うので別の語にする。
NONTARGET_UNANNOTATED = 'not annotated'
# final_status は kept / rejected / n/a の 3 値しか持たない。
# 理由はステップごとの status 列に分かれているので、表に出すときは
# ここで名前に直す(ノンターゲット表の rejected_by 列)。
REJECTION_REASON_COLUMNS = (
    ('is_filter_status',     STATUS_REJ_WINDOW,   'IS Filter'),
    ('rt_outlier_status',    STATUS_REJ_OUTLIER,  'RT Outlier Filter'),
    ('adduct_filter_status', STATUS_REJ_ADDUCT,   'Adduct Ion Filter'),
    ('coherence_status',     STATUS_REJ_RESIDUAL,
     'Coherence Filter (residual)'),
    ('coherence_status',     STATUS_REJ_LOSER,
     'Coherence Filter (inter-class)'),
    ('manual_status',        STATUS_REJ_MANUAL,   'manual'),
)


def rejection_reasons(row) -> str:
    """1 行がどのフィルタで落ちたかを ' + ' 区切りで返す。
    落ちていなければ空文字。"""
    out = []
    for col, bad, label in REJECTION_REASON_COLUMNS:
        try:
            v = row.get(col)
        except AttributeError:
            v = None
        if v is not None and str(v) == bad and label not in out:
            out.append(label)
    return ' + '.join(out)

ANNOTATION_CONFIDENCE_COLUMNS = (
    'confidence', 'sigma', 'decided_by', 'competitors')


def annotation_confidence(entry: dict | None, row=None,
                          sigma_threshold: float = 3.0) -> tuple:
    """帰属 1 件の信頼度を (ラベル, σ, 判定経路, 対立候補) にまとめる。

    エクスポートに出すための要約。ノンターゲット解析では
    「間違いが混ざっても帰属しておいて、信頼度で絞れる」方が使いやすい、
    という方針に合わせた。

    Parameters
    ----------
    entry : dict | None
        `_build_merged_attribution()` の 1 エントリ。None は
        「その座標で衝突していない」= 競合なし。
    row : match_df の 1 行 | None
        手動帰属かどうかを見るために使う。manual_status が kept なら
        人が決めたものとして最優先で扱う。
    sigma_threshold : float
        high / medium の境界（⑥ の σ しきい値と同じ値を渡す）。

    戻り値
    ------
    (label, sigma, decided_by, competitors)
      label       'high' / 'medium' / 'low' / 'none' / 'no conflict'
      sigma       ⑥ の σ（float）。⑤ 決着・手動・衝突なしでは NaN。
                  Excel で並べ替え・絞り込みができるよう数値で返す。
      decided_by  'adduct:ratio' / 'adduct:count' / 'coherence' /
                  'coherence:link-borrowed' / 'coherence:single' /
                  'manual' / 'undecided' / ''
      competitors '化合物[クラス]' を '; ' で連ねたもの
    """
    def _losers(e) -> str:
        if not e:
            return ''
        out = []
        for l in (e.get('losers') or []):
            c = str(l.get('compound', '') or '')
            k = str(l.get('class', '') or '')
            if c or k:
                out.append(f'{c}[{k}]' if k else c)
        if not out:
            c = str(e.get('loser_compound', '') or '')
            k = str(e.get('loser_class', '') or '')
            if c or k:
                out.append(f'{c}[{k}]' if k else c)
        return '; '.join(out[:8]) + (' ...' if len(out) > 8 else '')

    # 手動帰属は人が決めたものなので最優先
    try:
        if row is not None and str(row.get('manual_status', '')) == STATUS_KEPT:
            return ('high', float('nan'), 'manual', _losers(entry))
    except AttributeError:
        pass

    if not entry:
        return ('no conflict', float('nan'), '', '')

    method = str(entry.get('method', '') or '')
    comp = _losers(entry)

    # 判定 C で救済された行(勝者が落ちたので次点が残った)。
    # 比較で勝ったわけではないので低信頼度として出す。
    try:
        if (row is not None
                and str(row.get('final_status', '')) == STATUS_KEPT
                and entry.get('winner_compound')
                and str(row.get('compound', ''))
                != str(entry.get('winner_compound'))):
            return ('low', float('nan'), 'coherence:fallback', comp)
    except AttributeError:
        pass
    try:
        sig = float(entry.get('confidence'))
    except (TypeError, ValueError):
        sig = float('nan')

    if method == 'undecided':
        return ('none', float('nan'), 'undecided', comp)

    if method == 'adduct':
        am = str(entry.get('adduct_method', '') or '')
        if am == 'ratio':
            return ('high', float('nan'), 'adduct:ratio', comp)
        if am == 'count':
            return ('medium', float('nan'), 'adduct:count', comp)
        return ('medium', float('nan'), 'adduct', comp)

    # ⑥ Coherence
    _s = float(sig) if np.isfinite(sig) else float('nan')
    if entry.get('single_scored'):
        # 比較していない（相手が採点できなかった）
        return ('low', _s, 'coherence:single', comp)
    _fm = (str(entry.get('winner_fit_mode', '') or '')
           + '|' + str(entry.get('loser_fit_mode', '') or ''))
    if 'link_borrowed' in _fm:
        return ('low', _s, 'coherence:link-borrowed', comp)
    if not np.isfinite(sig):
        return ('low', _s, 'coherence', comp)
    if sig >= float(sigma_threshold):
        return ('high', _s, 'coherence', comp)
    if sig >= 1.0:
        return ('medium', _s, 'coherence', comp)
    return ('low', _s, 'coherence', comp)


# multi_label のデフォルト挙動切替
MULTI_LABEL_DEFAULT = True


# ─── キュレーション済み主モード+主アダクト ────────────────
#
# 各脂質クラスの「主モード+主アダクト」の内蔵知識。
# 各クラスについて1つの理論 m/z で統一的にマッチングする際に参照される。
#
# 編集方針:
#   - LipidZoner.py を直接編集(コードの一部、git で管理)
#   - LipidIonAudit が出力する snippet をコピー&ペーストで反映
#   - LipidQuant ライブラリの adduct フラグとは独立(マクロ専用と割り切る)
#
# 形式:
#   {class_name: {'mode': 'pos'|'neg', 'adduct': adduct文字列}}
#
# 未登録クラスは内蔵知識ベースのマッチングからスキップされる(起動時に警告ログ出力)。

CURATED_ADDUCTS: dict[str, dict[str, str]] = {
    # Glycerophospholipids
    'PC':       {'mode': 'pos', 'adduct': '[M+H]+'},
    'PE':       {'mode': 'pos', 'adduct': '[M+H]+'},
    'PS':       {'mode': 'neg', 'adduct': '[M-H]-'},
    'PG':       {'mode': 'neg', 'adduct': '[M-H]-'},
    'PI':       {'mode': 'neg', 'adduct': '[M-H]-'},
    'PA':       {'mode': 'neg', 'adduct': '[M-H]-'},
    'BMP':      {'mode': 'neg', 'adduct': '[M-H]-'},
    # Lyso
    'LPC':      {'mode': 'pos', 'adduct': '[M+H]+'},
    'LPE':      {'mode': 'pos', 'adduct': '[M+H]+'},
    'LPS':      {'mode': 'neg', 'adduct': '[M-H]-'},
    'LPG':      {'mode': 'neg', 'adduct': '[M-H]-'},
    'LPI':      {'mode': 'neg', 'adduct': '[M-H]-'},
    'LPA':      {'mode': 'neg', 'adduct': '[M-H]-'},
    # Sphingolipids
    'SM':       {'mode': 'pos', 'adduct': '[M+H]+'},
    'Cer':      {'mode': 'pos', 'adduct': '[M+H]+'},
    'CerPE':    {'mode': 'pos', 'adduct': '[M+H]+'},
    'HexCer':   {'mode': 'pos', 'adduct': '[M+H]+'},
    'Hex2Cer':  {'mode': 'pos', 'adduct': '[M+H]+'},
    'Hex3Cer':  {'mode': 'pos', 'adduct': '[M+H]+'},
    'Hex4Cer':  {'mode': 'pos', 'adduct': '[M+H]+'},
    'SHexCer':  {'mode': 'neg', 'adduct': '[M-H]-'},
    'GM3':      {'mode': 'neg', 'adduct': '[M-H]-'},
    'SPB':      {'mode': 'pos', 'adduct': '[M+H]+'},
    'SPBP':     {'mode': 'pos', 'adduct': '[M+H]+'},
    # Neutral lipids
    'TG':       {'mode': 'pos', 'adduct': '[M+NH4]+'},
    'DG':       {'mode': 'pos', 'adduct': '[M+NH4]+'},
    'MG':       {'mode': 'pos', 'adduct': '[M+NH4]+'},
    'SE':       {'mode': 'pos', 'adduct': '[M+NH4]+'},
    'CE':       {'mode': 'pos', 'adduct': '[M+NH4]+'},
    'ST':       {'mode': 'pos', 'adduct': '[M-H2O+H]+'},
    'FA':       {'mode': 'neg', 'adduct': '[M-H]-'},
    'CAR':      {'mode': 'pos', 'adduct': '[M+H]+'},
}


def _validate_curated_adducts() -> list[str]:
    """CURATED_ADDUCTS の整合性をチェックし、警告メッセージのリストを返す。

    起動時に呼ばれ、不正なエントリがあれば警告を表示する(該当クラスは
    マッチング対象からスキップされる)。
    """
    warnings: list[str] = []
    for cls, entry in CURATED_ADDUCTS.items():
        if not isinstance(entry, dict):
            warnings.append(f"CURATED_ADDUCTS['{cls}']: not a dict")
            continue
        mode = entry.get('mode')
        adu = entry.get('adduct')
        if mode not in ('pos', 'neg'):
            warnings.append(
                f"CURATED_ADDUCTS['{cls}']: invalid mode {mode!r} "
                f"(must be 'pos' or 'neg')")
        resolved = ADDUCT_ALIAS.get(adu, adu)
        if resolved not in ADDUCT_OFFSETS:
            warnings.append(
                f"CURATED_ADDUCTS['{cls}']: unknown adduct {adu!r}")
    return warnings


# ─── キュレーション済みアダクト fingerprint ────────────────
#
# 各脂質クラスの「9 アダクトの期待出方」を記述する内蔵知識。Adduct Ion
# Filter (Step ⑤) が衝突解決時に参照する。
#
#   オプション 2「全空、ランタイム IS 派生のみ」を採用。
#   起動時は空辞書、② IS Filter 後に ② IS の周辺で 9 アダクトを観測して
#   _runtime_is_fingerprints を生成し、これだけで運用する。
#
# 将来の拡張パス:
#   その snippet をここに流し込む運用に移行可能。データ構造はそのまま。
#
# 形式:
#   {
#       class_name: {
#           adduct_name: {
#               'expected':     bool,       # 期待されるアダクトか
#               'weight':       float,      # 判別スコアでの重み
#               'ratio_target': float|None, # 主アダクトに対する相対強度
#                                           #
#           },
#           ...
#       },
#       ...
#   }

CURATED_ADDUCT_FINGERPRINTS: dict[str, dict[str, dict[str, object]]] = {}


def _validate_curated_adduct_fingerprints() -> list[str]:
    """CURATED_ADDUCT_FINGERPRINTS の整合性をチェックし、警告メッセージの
    リストを返す。

    起動時に呼ばれる。空辞書(オプション 2 のデフォルト)は何も警告しない。
    エントリがある場合は、各クラスごとに以下を検証する:
      - 値が dict であること
      - 各アダクト名が ADDUCT_SET_V54 に含まれること
      - expected が bool、weight が float、ratio_target が float または None
    """
    warnings: list[str] = []
    for cls, entry in CURATED_ADDUCT_FINGERPRINTS.items():
        if not isinstance(entry, dict):
            warnings.append(
                f"CURATED_ADDUCT_FINGERPRINTS['{cls}']: not a dict")
            continue
        for adu, spec in entry.items():
            if adu not in ADDUCT_SET_V54:
                warnings.append(
                    f"CURATED_ADDUCT_FINGERPRINTS['{cls}']['{adu}']: "
                    f"unknown adduct (expected one of {ADDUCT_SET_V54})")
                continue
            if not isinstance(spec, dict):
                warnings.append(
                    f"CURATED_ADDUCT_FINGERPRINTS['{cls}']['{adu}']: "
                    f"not a dict")
                continue
            if not isinstance(spec.get('expected'), bool):
                warnings.append(
                    f"CURATED_ADDUCT_FINGERPRINTS['{cls}']['{adu}']: "
                    f"'expected' must be bool")
            w = spec.get('weight')
            if not isinstance(w, (int, float)):
                warnings.append(
                    f"CURATED_ADDUCT_FINGERPRINTS['{cls}']['{adu}']: "
                    f"'weight' must be float")
            rt = spec.get('ratio_target')
            if rt is not None and not isinstance(rt, (int, float)):
                warnings.append(
                    f"CURATED_ADDUCT_FINGERPRINTS['{cls}']['{adu}']: "
                    f"'ratio_target' must be float or None")
    return warnings


# ════════════════════════════════════════════════════════════════════
#  ライブラリ関連ユーティリティ
# ════════════════════════════════════════════════════════════════════

def _class_colors(n: int) -> list[str]:
    # matplotlib >= 3.9 removed matplotlib.cm.get_cmap; use the registry.
    cmap = matplotlib.colormaps['tab20'].resampled(max(n, 20))
    return [mcolors.to_hex(cmap(i % cmap.N)) for i in range(n)]


def fe_intensity_metric(fe) -> str:
    """fe.df のサンプル列が height か area かを返す。

    height_s* には **Detection Step 3 で選んだ指標** が入る。Area を
    選んでいれば height_s* の中身は面積である。列名からは分からないので
    align_params['height_col'](session にも保存されている)を見る。

    古いデータで align_params が無い場合は area_s* の有無から推定する
    (area_s* は height_col='Height' のときにだけ書かれるため)。
    それでも分からなければ 'intensity' を返す。
    """
    hc = str(((getattr(fe, 'align_params', None) or {})
              .get('height_col') or '')).strip().lower()
    if hc in ('height', 'area'):
        return hc
    try:
        cols = list(getattr(fe, '_aligned_df').columns)
    except Exception:
        cols = []
    if any(str(c).startswith('area_s') for c in cols):
        return 'height'
    return 'intensity'


def fe_intensity_tables(fe) -> list:
    """[(指標名, サンプル値の DataFrame), ...] を返す。

    1 つ目は fe.df のサンプル列(= Detection で選んだ指標)。
    2 つ目は area_s*。**Height で走らせたときだけ** 存在する
    (Area で走らせると area_s* は書かれず、height 側も残らない)。

    どちらも fe.df と同じ行順・同じ列名(サンプル名)に揃えてある。
    """
    out: list = []
    base = getattr(fe, 'df', None)
    if base is None or getattr(base, 'empty', True):
        return out
    try:
        smp = list(fe.sample_columns())
    except Exception:
        smp = []
    if not smp:
        return out
    primary = fe_intensity_metric(fe)
    out.append((primary, base[smp].copy()))
    if primary == 'area':
        # Area で走らせた場合、height は保持されていない
        return out
    ad = getattr(fe, '_aligned_df', None)
    if ad is None or '_feature_id_internal' not in base.columns:
        return out
    acols = [f'area_s{i}' for i in range(len(smp))]
    if not all(c in ad.columns for c in acols):
        return out
    try:
        idx = ad.drop_duplicates(subset='feature_id').set_index('feature_id')
        sub = idx.loc[base['_feature_id_internal'].values, acols]
        if len(sub) != len(base):
            return out
        area = pd.DataFrame(np.asarray(sub, dtype=float),
                            columns=smp, index=base.index)
    except Exception as e:
        log.warning(f"area_s lookup failed: {e}")
        return out
    out.append(('area', area))
    return out


def fe_filled_table(fe) -> "pd.DataFrame | None":
    """fe.df の行順・サンプル名に揃えた filled フラグ表を返す。

    ギャップフィリングを通していない(または旧データの)エントリでは None。
    """
    base = getattr(fe, 'df', None)
    if base is None or getattr(base, 'empty', True):
        return None
    try:
        smp = list(fe.sample_columns())
    except Exception:
        smp = []
    if not smp:
        return None
    ad = getattr(fe, '_aligned_df', None)
    if ad is None or '_feature_id_internal' not in base.columns:
        return None
    fcols = [f'filled_s{i}' for i in range(len(smp))]
    if not all(c in ad.columns for c in fcols):
        return None
    try:
        idx = ad.drop_duplicates(subset='feature_id').set_index('feature_id')
        sub = idx.loc[base['_feature_id_internal'].values, fcols]
        if len(sub) != len(base):
            return None
        return pd.DataFrame(np.asarray(sub).astype(bool),
                            columns=smp, index=base.index)
    except Exception as e:
        log.warning(f"filled_s lookup failed: {e}")
        return None


def ap_window_title(label: str) -> str:
    """Annotation Pipeline ウィンドウのタイトル。

    トップレベルの .py を固定名で運用するようにしたので、いま動いて
    いるのがどの版かをファイル名から読めない。唯一いつも見えている
    このウィンドウのタイトルに版番号を出す。
    """
    return f"LipidZoner v{__version__} — Annotation Pipeline — {label}"


def feature_uid(mode, rt, mz, rt_digits: int = 3,
                mz_digits: int = 4) -> str:
    """RT と m/z から feature のユニーク ID を作る。

      feature_uid('pos', 2.9703, 786.60112) -> 'pos_2.970_786.6011'

    極性を先頭に付けるのは、pos と neg を 1 つの行列に積んで多変量解析
    するときに偶然の衝突を避けるため。RT 3 桁 / m/z 4 桁は、アライメント
    の再現性(RT 中央値 0.005 min)と m/z 許容値(0.005 Da)より十分細かく、
    かつ表示がぶれない桁として選んだ。

    値が欠けているところは 'na' にする。空文字にすると結合キーとして
    使えなくなるため。
    """
    def _f(v, d):
        try:
            fv = float(v)
        except (TypeError, ValueError):
            return 'na'
        if fv != fv:          # NaN
            return 'na'
        return f"{fv:.{d}f}"
    _m = str(mode or 'na').strip() or 'na'
    return f"{_m}_{_f(rt, rt_digits)}_{_f(mz, mz_digits)}"


def timestamped_filename(name: str) -> str:
    """既定のファイル名の末尾に、いま(= ボタンを押した時刻)の
    日時を挟む。拡張子はそのまま残す。

      quantification_results.xlsx
        → quantification_results_260828_143512.xlsx

    書式は検証データのファイル名の慣習(260227_… / 260818_…)に合わせて
    YYMMDD_HHMMSS。名前順に並べるとそのまま時系列になる。

    export のたびに手で名前を変えないと前のファイルを潰してしまう、
    という事故を防ぐのが目的。ダイアログ上で編集できるので、
    同じ名前で上書きしたいときはその場で消せる。
    """
    _ts = datetime.datetime.now().strftime('%y%m%d_%H%M%S')
    _name = str(name or '')
    stem, dot, ext = _name.rpartition('.')
    if not dot:
        return f"{_name}_{_ts}" if _name else _ts
    return f"{stem}_{_ts}.{ext}"


def _exact_mass(c, h, o, n, p, s, d, c13) -> float:
    m = MONO_MASSES
    return (m['C']*c + m['H']*h + m['O']*o + m['N']*n +
            m['P']*p + m['S']*s + m['D']*d + m['13C']*c13)


def load_lipid_library(path: str | Path, ion_mode: str) -> pd.DataFrame:
    """
    xlsmライブラリを読み込み、指定イオンモードのエントリを返す。
    ion_mode: 'pos' | 'neg' | 'both'

    各エントリに is_IS フラグ（内部標準かどうか）を付与する。
    IS列は13列目（インデックス13、「IS」列ヘッダ）を参照する。
    """
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True)
    records = []
    for sheet in wb.sheetnames:
        ws = wb[sheet]
        rows = list(ws.iter_rows(max_row=2, values_only=True))
        if len(rows) < 2 or len(rows[1]) <= 12:
            continue
        hdr = rows[1]
        adduct_pos = str(hdr[11]) if hdr[11] else None
        adduct_neg = str(hdr[12]) if hdr[12] else None

        if ion_mode == 'pos':
            use = [(adduct_pos, 11)] if adduct_pos else []
        elif ion_mode == 'neg':
            use = [(adduct_neg, 12)] if adduct_neg else []
        else:
            use = []
            if adduct_pos: use.append((adduct_pos, 11))
            if adduct_neg: use.append((adduct_neg, 12))

        for row in ws.iter_rows(min_row=3, values_only=True):
            if not row or not row[0]:
                continue
            try:
                elems = tuple(int(row[i] or 0) for i in range(1, 9))
            except (TypeError, ValueError):
                continue
            em = _exact_mass(*elems)
            # IS列（列13、0-indexed）: 1/True/"1" ならIS
            is_is_flag = row[13] if len(row) > 13 else None
            is_IS = is_is_flag in (1, '1', True)
            for adduct_str, col_idx in use:
                flag = row[col_idx] if len(row) > col_idx else None
                if flag not in (1, '1', True):
                    continue
                resolved = ADDUCT_ALIAS.get(adduct_str, adduct_str)
                if resolved not in ADDUCT_OFFSETS:
                    continue
                records.append({
                    'compound':       row[0],
                    'lipid_class':    sheet,
                    'adduct':         adduct_str,
                    'exact_mass':     em,
                    'theoretical_mz': em + ADDUCT_OFFSETS[resolved],
                    'is_IS':          is_IS,
                    # col14 = 各分子種に割り当てられた IS 名
                    'is_ref_name':    (str(row[14]).strip()
                                       if len(row) > 14 and row[14]
                                       not in (None, '') else None),
                })
    return pd.DataFrame(records)


def _collect_match_candidates(
    th: float,
    obs_mz: np.ndarray,
    obs_rt: np.ndarray,
    ppm_tol: float,
    sample_intensities: np.ndarray | None = None,
    void_rt: float = 0.0,
) -> list[dict]:
    """ppm tol 内の全候補ピークを m/z 誤差絶対値の小→大 順で返す。

    用途: ISPreviewDialog の Candidate Picker 機能。混雑領域で同 m/z の
    ピークが複数ある場合、ユーザーが正しい IS を手動選択できるようにする。

    Parameters
    ----------
    th : float
        理論 m/z
    obs_mz, obs_rt : ndarray
        観測ピークの m/z, RT(マスク済みの 1D 配列)
    ppm_tol : float
        m/z 許容誤差(ppm)
    sample_intensities : ndarray | None
        shape (n_peaks, n_samples) のサンプル強度行列。指定があれば
        各候補に sample_intensities と mean_intensity を付与する。

    戻り値:
        [{peak_idx, obs_rt, obs_mz, delta_ppm,
          [mean_intensity, sample_intensities]}, ...]
        m/z 誤差絶対値小→大でソート済み。候補なしなら空リスト。
    """
    diffs = obs_mz - th
    abs_diffs = np.abs(diffs)
    cand_idx = np.where(abs_diffs <= th * ppm_tol / 1e6)[0]
    # ボイド領域の候補を落とす(Candidate Picker にも出さない)
    if void_rt and void_rt > 0 and len(cand_idx):
        cand_idx = cand_idx[obs_rt[cand_idx] >= void_rt]
    if len(cand_idx) == 0:
        return []
    # m/z 誤差絶対値で昇順ソート(現行の m/z 最小誤差マッチが先頭になる)
    sort_order = np.argsort(abs_diffs[cand_idx])
    cand_idx = cand_idx[sort_order]

    out: list[dict] = []
    for i in cand_idx:
        c = {
            'peak_idx':  int(i),
            'obs_mz':    float(obs_mz[i]),
            'obs_rt':    float(obs_rt[i]),
            'delta_ppm': float((obs_mz[i] - th) / th * 1e6),
        }
        if sample_intensities is not None:
            sv = sample_intensities[i].astype(float)
            c['sample_intensities'] = sv.tolist()
            c['mean_intensity']     = float(np.mean(sv))
        out.append(c)
    return out


def _match_single_entry(
    th: float,
    cls: str,
    obs_mz: np.ndarray,
    obs_rt: np.ndarray,
    ppm_tol: float,
    ref_rt: float | None = None,
    rt_tol: float | None = None,
    excluded_idx: set | None = None,
    void_rt: float = 0.0,
) -> dict:
    """
    単一ライブラリエントリをマッチングする内部ヘルパー。

    ref_rt/rt_tol が指定されていれば (ref_rt ± rt_tol) のRT窓に限定し
    RT誤差最小の候補を選ぶ。それ以外はm/z誤差最小。
    excluded_idx に含まれる観測ピークインデックスは候補から除外する。

    void_rt > 0 なら、それより早く出る観測ピークは候補にしない。
    ボイドボリューム(非保持成分・キャリーオーバー)を弾くため。
    **観測配列そのものは触らない**ので matched_idx の意味は変わらない。
    """
    diffs = np.abs(obs_mz - th)
    cands = np.where(diffs <= th * ppm_tol / 1e6)[0]

    # ボイド領域の候補を落とす(ppm 窓を作った直後に 1 回だけ)
    if void_rt and void_rt > 0 and len(cands):
        cands = cands[obs_rt[cands] >= void_rt]

    if excluded_idx:
        cands = np.array([i for i in cands if i not in excluded_idx])

    if len(cands) == 0:
        return {'obs_mz': np.nan, 'obs_rt': np.nan,
                'delta_ppm': np.nan, 'matched': False,
                'matched_idx': -1}

    if ref_rt is not None and rt_tol is not None:
        rt_diffs = np.abs(obs_rt[cands] - ref_rt)
        in_window = rt_diffs <= rt_tol
        if not np.any(in_window):
            return {'obs_mz': np.nan, 'obs_rt': np.nan,
                    'delta_ppm': np.nan, 'matched': False,
                    'matched_idx': -1}
        sub_cands = cands[in_window]
        sub_rtd   = rt_diffs[in_window]
        best = sub_cands[np.argmin(sub_rtd)]
    else:
        best = cands[np.argmin(diffs[cands])]

    return {'obs_mz':      float(obs_mz[best]),
            'obs_rt':      float(obs_rt[best]),
            'delta_ppm':   float(diffs[best] / th * 1e6),
            'matched':     True,
            'matched_idx': int(best)}


# ─── IS 自動選択 ───────────────────────────────────────
#  複数候補があるときに何を選ぶか。**手動選択(manual_rt)が最優先**で、
#  以下はピン留めが無いときのフォールバック。
#
#  IS_AUTO_MIN_REL_INT: 最強候補に対する相対強度がこれ未満の候補を捨てる。
#    残った候補の中から |Δm/z| 最小を選ぶので、強度が同程度なら現行と同じ挙動。
#
#  既定 OFF の理由(検証データでの実測、S6):
#    複数候補を持つ IS 4 件のうち、強度規則が改善するのは 1 件
#    (pos/HexCer、正解は相対強度 1.00)だが、別の 1 件は逆に悪化する
#    (neg/HexCer、正解は相対強度 0.31 の弱い方)。さらに pos/SE の正解は
#    相対強度 0.43・+7.32 ppm で、最強でも最近接でもない。
#    単純な自動規則では正解に届かないため、既定は従来どおり |Δm/z| 最小とし、
#    強度加味は明示的に有効化するオプションとする。
# ボイドボリューム(非保持成分・キャリーオーバー)の RT 上限。
# ここより早く出る観測ピークはライブラリマッチの候補にしない。
# カラム・流量・メソッドで変わる値なので UI(Advanced Setting ① )で
# 変更でき、session にも保存される。0 にすると無効。
#
# 既定 0.6 min の根拠(検証データ実測):
#   pos は RT 0.35-0.50 に 691 feature(全強度の 11.5%)の密集帯があり、
#   0.55-0.75 が谷、0.80 から本体。本物の最早は SE の IS 0.818。
VOID_RT_CUTOFF_DEFAULT = 0.60

# 画面フォントの大きさ。1700 px のウィンドウに対して 9pt は
# 小さく、各列の右側が余っていた。11pt で余白が埋まり読みやすくなる。
# 環境に合わせて調整したいときはここだけ変えればよい(全タブに効く)。
UI_FONT_POINT_SIZE: int = 11

IS_AUTO_MIN_REL_INT: float = 0.5


def _select_is_peak_idx(
    th: float,
    obs_mz: np.ndarray,
    ppm_tol: float,
    intensities: np.ndarray | None = None,
    use_intensity: bool = False,
    min_rel_int: float = IS_AUTO_MIN_REL_INT,
) -> int:
    """IS の観測ピークを 1 本選ぶ。該当なしなら -1。

    use_intensity=False(既定)なら |Δm/z| 最小(従来の挙動)。
    True なら「最強候補の min_rel_int 倍未満の候補を捨ててから
    |Δm/z| 最小」を選ぶ。
    """
    diffs = np.abs(obs_mz - th)
    cand = np.where(diffs <= th * ppm_tol / 1e6)[0]
    if len(cand) == 0:
        return -1
    if use_intensity and intensities is not None and len(cand) > 1:
        try:
            ii = np.asarray(intensities, dtype=float)[cand]
            ii = np.where(np.isfinite(ii), ii, 0.0)
            if ii.max() > 0:
                keep = cand[ii >= ii.max() * float(min_rel_int)]
                if len(keep):
                    cand = keep
        except Exception:
            pass
    return int(cand[int(np.argmin(diffs[cand]))])


def collect_multi_candidate_is(
    is_choices: dict | None,
    lib_df: 'pd.DataFrame',
    obs_mz: np.ndarray,
    ppm_tol: float,
) -> list[str]:
    """複数候補がありながら手動でピン留めされていない IS を列挙する。

    対話実行時は ISPreviewDialog の has_unconfirmed_multi_candidates() が
    同等の警告を出すが、session を読んで Run All する auto-replay 経路では
    ダイアログごとスキップされるため警告が出ていなかった。
    """
    out: list[str] = []
    if not is_choices:
        return out
    try:
        is_rows = lib_df[lib_df.get('is_IS', False) == True]
    except Exception:
        return out
    for _, r in is_rows.iterrows():
        cls = r.get('lipid_class')
        rows = is_choices.get(cls) or []
        ch = next((x for x in rows
                   if x.get('compound') == r.get('compound')
                   and x.get('adduct') == r.get('adduct')), None)
        if ch is None or not ch.get('use'):
            continue
        if ch.get('manual_rt') is not None:
            continue          # ピン留め済み
        th = float(r['theoretical_mz'])
        n = int((np.abs(obs_mz - th) <= th * ppm_tol / 1e6).sum())
        if n > 1:
            out.append(f"{cls} / {r.get('compound')} {r.get('adduct')} "
                       f"({n} candidates)")
    return out


def _class_rt_tol_override(is_choices, cls: str, default_tol: float) -> float:
    """③ の表で指定されたクラス別 RT tolerance を返す。

    is_choices[cls] の各行に 'rt_tol' があり、0 より大きければそれを使う。
    未指定(None / 0 / キー無し)なら default_tol(Advanced の IS RT
    tolerance)にフォールバックする。RT tol はクラスに 1 つなので、
    最初に見つかった正の値を採る(ダイアログ側で同期済み)。
    """
    if not is_choices:
        return float(default_tol)
    for row in (is_choices.get(cls) or []):
        try:
            v = row.get('rt_tol')
        except AttributeError:
            continue
        if v:
            try:
                fv = float(v)
            except (TypeError, ValueError):
                continue
            if fv > 0:
                return fv
    return float(default_tol)


def match_library_to_peaks(
    lib_df: pd.DataFrame,
    obs_mz: np.ndarray,
    obs_rt: np.ndarray,
    ppm_tol: float,
    class_ref_rt: dict[str, tuple[float, float, str]] | None = None,
    is_tol: float = 0.10,
    apply_is_filter: bool = True,
    is_choices: dict[str, list[dict]] | None = None,
    is_peak_intensities: 'np.ndarray | None' = None,
    is_use_intensity: bool = False,
    multi_label: bool = MULTI_LABEL_DEFAULT,
    void_rt: float = 0.0,
) -> tuple[pd.DataFrame, dict[str, tuple[float, float, str]], list[str]]:
    """
    ライブラリを観測ピークに対してマッチングする。

    apply_is_filter=True の場合、2パスマッチングを実行:

    第1パス: IS保持クラスを優先処理
      - ISエントリをm/z最小誤差でマッチング
      - obs_rtの中央値を ref_rt として、同クラスの全エントリを
        自動Fix RT (ref=ref_rt, tol=is_tol) で再マッチ
      - 手動Fix RT設定 (source="manual") はそちらを優先
      - ISが見つからない場合は通常m/z最小誤差マッチングにフォールバック
      - マッチした観測ピークインデックスを reserved_idx に登録

    第2パス: IS非保持クラス
      - reserved_idx に含まれる候補を除外してm/z最小誤差マッチング
      - 手動Fix RT設定があればそちらを優先（除外は適用しない）

    apply_is_filter=False の場合、通常マッチング:
      - 全エントリをm/z最小誤差でマッチング
      - 手動Fix RT設定があるクラスのみ RT 制約を適用
      - class_ref_rt の auto_is エントリは無視
      - reserved_idx は空

    Parameters
    ----------
    lib_df : DataFrame
        ライブラリエントリ（is_IS列を含む）
    obs_mz, obs_rt : ndarray
        観測ピークのm/z, RT
    ppm_tol : float
        m/z許容誤差 (ppm)
    class_ref_rt : dict[str, (ref_rt, tol_min, source)] | None
        クラス単位のreference RT指定。source は "manual" または "auto_is"。
    is_tol : float
        IS自動Fix RTのtolerance（class_ref_rtに未登録クラス用のデフォルト）
    apply_is_filter : bool
        True の場合 2パスマッチング、False の場合通常マッチング（IS処理なし）

    Returns
    -------
    match_df : DataFrame
        マッチ結果（lib_df + obs_mz/obs_rt/delta_ppm/matched列）
    updated_class_ref_rt : dict
        IS自動算出された ref_rt を反映した class_ref_rt
        apply_is_filter=False の場合は入力の manual エントリのみ保持
    is_missing_classes : list[str]
        ISが検出されずフォールバックしたクラス名のリスト
        apply_is_filter=False の場合は常に空
    """
    if lib_df.empty:
        return lib_df.copy(), dict(class_ref_rt or {}), []

    class_ref_rt = dict(class_ref_rt or {})

    # apply_is_filter=False の場合は auto_is エントリを無視
    if not apply_is_filter:
        class_ref_rt = {
            cls: v for cls, v in class_ref_rt.items()
            if v[2] == "manual"
        }

    classes_arr  = lib_df['lipid_class'].values
    th_arr       = lib_df['theoretical_mz'].values
    is_IS_arr    = lib_df['is_IS'].values if 'is_IS' in lib_df.columns \
                   else np.zeros(len(lib_df), dtype=bool)

    # 通常マッチングモード: 全エントリを単純マッチング
    if not apply_is_filter:
        results: dict[int, dict] = {}
        for idx in range(len(lib_df)):
            cls = classes_arr[idx]
            if cls in class_ref_rt:
                ref_rt = class_ref_rt[cls][0]
                rt_tol = class_ref_rt[cls][1]
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    ref_rt=ref_rt, rt_tol=rt_tol, void_rt=void_rt)
            else:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    void_rt=void_rt)
            results[idx] = r
        result_rows = [results.get(i, {
            'obs_mz': np.nan, 'obs_rt': np.nan,
            'delta_ppm': np.nan, 'matched': False, 'matched_idx': -1
        }) for i in range(len(lib_df))]
        match_df = pd.concat(
            [lib_df.reset_index(drop=True), pd.DataFrame(result_rows)], axis=1)
        # Match Overlay でも multi_label 展開を行う。
        # これにより DG sn 異性体のように同一 lib エントリに対して ppm 内に
        # 複数候補がある場合、すべての候補が matched 行として記録される。
        # IS Filter で新たなスポットが現れる現象を解消する。
        if multi_label:
            match_df = _expand_to_multi_label(
                match_df, obs_mz, obs_rt, ppm_tol,
                class_ref_rt=class_ref_rt,
                is_tol=is_tol,
                apply_is_filter=False,
                reserved_peak_idx=None,
                void_rt=void_rt,
            )
        return match_df, class_ref_rt, []

    # IS保持クラスを特定
    is_bearing_classes = set()
    for cls in lib_df['lipid_class'].unique():
        mask = (classes_arr == cls) & is_IS_arr
        if np.any(mask):
            is_bearing_classes.add(cls)

    # is_choices が指定されている場合、ユーザーが全 IS を拒否したクラスは
    # IS 保持クラスから除外する
    if is_choices is not None:
        filtered = set()
        for cls in is_bearing_classes:
            choices = is_choices.get(cls, [])
            if any(c.get("use") for c in choices):
                filtered.add(cls)
        is_bearing_classes = filtered

    # is_choices を (cls, compound, adduct) → {use, manual_rt} にルックアップ可能に
    choice_lookup: dict[tuple, dict] = {}
    if is_choices is not None:
        for cls, rows in is_choices.items():
            for row in rows:
                key = (cls, row.get("compound", ""), row.get("adduct", ""))
                choice_lookup[key] = row

    compound_arr = lib_df['compound'].values
    adduct_arr   = lib_df['adduct'].values

    # 結果格納（indexで引けるように辞書で持つ）
    results: dict[int, dict] = {}
    reserved_idx: set = set()  # 第1パスで帰属した観測ピークのインデックス
    is_missing_classes: list[str] = []

    # ─── 第1パス: IS保持クラス ─────────────────────────────────────
    for cls in is_bearing_classes:
        cls_mask    = classes_arr == cls
        cls_indices = np.where(cls_mask)[0]

        # 手動Fix RT設定がある場合はそれを使用
        if cls in class_ref_rt and class_ref_rt[cls][2] == "manual":
            ref_rt = class_ref_rt[cls][0]
            rt_tol = class_ref_rt[cls][1]
            for idx in cls_indices:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    ref_rt=ref_rt, rt_tol=rt_tol, void_rt=void_rt)
                results[idx] = r
                if r['matched']:
                    reserved_idx.add(r['matched_idx'])
            continue

        # ISエントリのみ先にマッチング（m/z誤差最小、RT制約なし）
        # is_choices がある場合は Use=True のIS のみ使用
        is_indices = cls_indices[is_IS_arr[cls_indices]]
        is_rts = []
        for idx in is_indices:
            if is_choices is not None:
                key = (cls, compound_arr[idx], adduct_arr[idx])
                choice = choice_lookup.get(key)
                if choice is None or not choice.get("use"):
                    # Use=False のISは無視（マッチング結果なし）
                    results[idx] = {
                        'obs_mz': np.nan, 'obs_rt': np.nan,
                        'delta_ppm': np.nan, 'matched': False,
                        'matched_idx': -1}
                    continue
                # Manual RT があればそれを採用
                manual_rt = choice.get("manual_rt")
                if manual_rt is not None:
                    # Manual RT を obs_rt として「仮想的に」採用
                    # マッチング自体は通常通り行うが、ref_rt算出用に manual_rt を使う
                    # クラス別 RT tol があればそれを使う
                    r = _match_single_entry(
                        th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                        ref_rt=float(manual_rt),
                        rt_tol=_class_rt_tol_override(is_choices, cls, is_tol),
                        void_rt=void_rt)
                    results[idx] = r
                    if r['matched']:
                        is_rts.append(float(manual_rt))
                        reserved_idx.add(r['matched_idx'])
                    else:
                        # Manual RT 近傍にピークがなくても、ref_rt としては採用
                        is_rts.append(float(manual_rt))
                    continue

            # 通常のm/z最小マッチング
            # is_use_intensity=True のときだけ、明らかに弱い候補を
            #   落としてから |Δm/z| 最小を採る。既定は従来どおり。
            if is_use_intensity and is_peak_intensities is not None:
                _pi = _select_is_peak_idx(
                    th_arr[idx], obs_mz, ppm_tol,
                    intensities=is_peak_intensities, use_intensity=True)
                # 強度ベースの自動選択でもボイド領域は候補にしない
                if (_pi >= 0 and void_rt and void_rt > 0
                        and float(obs_rt[_pi]) < void_rt):
                    _pi = -1
                if _pi >= 0:
                    r = {'obs_mz': float(obs_mz[_pi]),
                         'obs_rt': float(obs_rt[_pi]),
                         'delta_ppm': float(
                             (obs_mz[_pi] - th_arr[idx]) / th_arr[idx] * 1e6),
                         'matched': True, 'matched_idx': int(_pi)}
                else:
                    r = {'obs_mz': np.nan, 'obs_rt': np.nan,
                         'delta_ppm': np.nan, 'matched': False,
                         'matched_idx': -1}
            else:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    void_rt=void_rt)
            results[idx] = r
            if r['matched']:
                is_rts.append(r['obs_rt'])
                reserved_idx.add(r['matched_idx'])

        if is_rts:
            # IS検出成功 → 中央値を ref_rt として自動Fix RT適用
            ref_rt = float(np.median(is_rts))
            # ③ の表でクラス別 RT tol が入っていればそれを使う。
            # 1 クラス内で溶出位置が分かれる脂質(HexCer の O2 系 / O3 系は
            # 0.26 min 離れる)では、既定の ±0.1 min だと片方が丸ごと
            # rejected_window で落ちるため。
            rt_tol = _class_rt_tol_override(is_choices, cls, is_tol)
            class_ref_rt[cls] = (ref_rt, rt_tol, "auto_is", "ref_tol")

            # IS以外のエントリにRT制約付きマッチング
            non_is_indices = cls_indices[~is_IS_arr[cls_indices]]
            for idx in non_is_indices:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    ref_rt=ref_rt, rt_tol=rt_tol, void_rt=void_rt)
                results[idx] = r
                if r['matched']:
                    reserved_idx.add(r['matched_idx'])
        else:
            # IS検出失敗 → 通常マッチングにフォールバック
            is_missing_classes.append(cls)
            # 自動IS設定は残さない（ISが復活したら次回再計算されるため）
            if cls in class_ref_rt and class_ref_rt[cls][2] == "auto_is":
                del class_ref_rt[cls]
            non_is_indices = cls_indices[~is_IS_arr[cls_indices]]
            for idx in non_is_indices:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    void_rt=void_rt)
                results[idx] = r
                if r['matched']:
                    reserved_idx.add(r['matched_idx'])

    # ─── 第2パス: IS非保持クラス ───────────────────────────────────
    non_is_classes = set(lib_df['lipid_class'].unique()) - is_bearing_classes
    for cls in non_is_classes:
        cls_mask    = classes_arr == cls
        cls_indices = np.where(cls_mask)[0]

        # 手動Fix RT設定がある場合はそれを使用（除外は適用しない）
        if cls in class_ref_rt and class_ref_rt[cls][2] == "manual":
            ref_rt = class_ref_rt[cls][0]
            rt_tol = class_ref_rt[cls][1]
            for idx in cls_indices:
                r = _match_single_entry(
                    th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                    ref_rt=ref_rt, rt_tol=rt_tol, void_rt=void_rt)
                results[idx] = r
            continue

        # 通常マッチング（reserved_idxを除外）
        for idx in cls_indices:
            r = _match_single_entry(
                th_arr[idx], cls, obs_mz, obs_rt, ppm_tol,
                excluded_idx=reserved_idx, void_rt=void_rt)
            results[idx] = r

    # ─── 結果をDataFrameに変換 ─────────────────────────────────────
    result_rows = [results.get(i, {
        'obs_mz': np.nan, 'obs_rt': np.nan,
        'delta_ppm': np.nan, 'matched': False, 'matched_idx': -1
    }) for i in range(len(lib_df))]
    result_df = pd.DataFrame(result_rows)
    # matched_idx はエクスポート除外にも使うので残す
    match_df = pd.concat(
        [lib_df.reset_index(drop=True), result_df], axis=1)

    # multi_label 展開
    # 既存ロジックで「ベスト 1 候補」のマッチが入っている match_df を、
    # 「全候補(ppm 内)」の多行データフレームに拡張する。
    # 各候補に isomer_rank を付与し、status 列をデフォルト初期化。
    if multi_label:
        match_df = _expand_to_multi_label(
            match_df, obs_mz, obs_rt, ppm_tol,
            class_ref_rt=class_ref_rt,
            is_tol=is_tol,
            apply_is_filter=apply_is_filter,
            reserved_peak_idx=reserved_idx if apply_is_filter else None,
            void_rt=void_rt)
    else:
        # 単一ラベル: status 列だけ追加して初期化
        match_df = _add_status_columns(match_df)

    return match_df, class_ref_rt, is_missing_classes


def _add_status_columns(match_df: pd.DataFrame) -> pd.DataFrame:
    """match_df に status 列をデフォルト値で追加する。"""
    df = match_df.copy()
    n = len(df)
    df['isomer_rank']         = np.where(df['matched'], 1, 0)
    df['is_filter_status']    = STATUS_NA
    df['adduct_filter_status'] = STATUS_NA
    df['rt_outlier_status']   = STATUS_NA
    df['coherence_status']    = STATUS_NA
    df['manual_status']       = STATUS_NA
    df['final_status']        = np.where(df['matched'], STATUS_KEPT, STATUS_NA)
    return df


def _apply_rt_outlier_status(
    match_df: pd.DataFrame,
    outlier_coords: dict[str, set],
) -> pd.DataFrame:
    """outlier_coords を rt_outlier_status 列にマークする。

    削除はせず status 列を更新するだけ。final_status は再計算する。
    """
    if match_df.empty:
        return match_df
    df = match_df.copy()
    if not outlier_coords:
        # 全クリア
        df['rt_outlier_status'] = STATUS_NA
        return _update_final_status_df(df)

    statuses = []
    for _, row in df.iterrows():
        if not row.get('matched'):
            statuses.append(STATUS_NA)
            continue
        cls = row.get('lipid_class')
        rt_r = round(float(row['obs_rt']), 6)
        mz_r = round(float(row['obs_mz']), 6)
        cls_set = outlier_coords.get(cls, set())
        if (rt_r, mz_r) in cls_set:
            statuses.append(STATUS_REJ_OUTLIER)
        else:
            statuses.append(STATUS_KEPT)
    df['rt_outlier_status'] = statuses
    return _update_final_status_df(df)


def _apply_coherence_status(
    match_df: pd.DataFrame,
    coherence_models: dict[tuple[str, str], dict],
    coherence_assignments: dict[tuple, dict],
) -> pd.DataFrame:
    """Coherence による status 列更新。

    判定 A: coherence_assignments が winner を決めた座標では、
            winner と一致しない matched 行を全て
            coherence_status = STATUS_REJ_LOSER にマークする。
            fix48 までは loser_class(次点)だけを見ていたため、
            3 候補以上の衝突で 3 番目以降が kept のまま残っていた。
    判定 B(新規): 同一 compound に複数の生存候補がある場合、
                  pooled 回帰の予測 RT との残差が最小の行を kept、
                  他を STATUS_REJ_RESIDUAL にマーク。
    判定 C: 判定 A の勝者が 判定 B で落ちて、そのスポットに生き残りが
                  いなくなった場合、次点以降を残差順に見て 1 つだけ kept に
                  戻す。手動で決めてあるスポットと、③④⑤ で落ちている候補、
                  その compound が別ピークで採用済みのものは対象外。

    生存候補とは、is_filter_status != REJ_WINDOW(IS Filter を通過)で
    かつ coherence_status != REJ_LOSER(衝突敗者でない)の matched 行。
    """
    if match_df.empty:
        return match_df
    df = match_df.copy()

    # matched 行はデフォルト KEPT に。matched=False は NA のまま。
    df['coherence_status'] = np.where(
        df.get('matched', False) == True, STATUS_KEPT, STATUS_NA)

    # 判定 A: 衝突敗者をマーク(現状の既存挙動を status 列に反映)
    if coherence_assignments:
        rt_r = df['obs_rt'].round(6)
        mz_r = df['obs_mz'].round(6)
        for i in df.index:
            if not df.at[i, 'matched']:
                continue
            coord = (float(rt_r[i]), float(mz_r[i]))
            a = coherence_assignments.get(coord)
            if a is None:
                continue
            # winner 本人だけを残し、それ以外は全部 loser にする。
            # 1 つの観測ピークに帰属できる分子種は 1 つだけ、という
            # ⑥ の建前をそのまま status 列に反映する。
            _is_winner = (
                str(df.at[i, 'lipid_class']) == str(a.get('winner_class'))
                and str(df.at[i, 'compound']) == str(a.get('winner_compound'))
                and str(df.at[i, 'adduct']) == str(a.get('winner_adduct')))
            if not _is_winner:
                df.at[i, 'coherence_status'] = STATUS_REJ_LOSER

    # 判定 B: 同一 compound 内で残差最小選択
    if coherence_models:
        # C, U, link を計算(後で削除)
        df['_C']    = df['compound'].apply(
            lambda s: _parse_c_and_u(str(s) if s else '')[0])
        df['_U']    = df['compound'].apply(
            lambda s: _parse_c_and_u(str(s) if s else '')[1])
        df['_link'] = df['compound'].apply(
            lambda s: _link_kind(str(s) if s else ''))

        # 生存候補のみを対象
        alive_mask = (
            (df.get('matched', False) == True) &
            (df.get('is_filter_status', STATUS_NA) != STATUS_REJ_WINDOW) &
            (df['coherence_status'] != STATUS_REJ_LOSER) &
            df['_C'].notna() & df['_U'].notna()
        )

        for compound, group in df[alive_mask].groupby('compound'):
            if len(group) <= 1:
                continue  # 単一候補なら選択不要

            # 残差最小の 1 本を「系列ごとに」残す。二峰性クラス
            # (DG)は 1 つの分子種が 2 本のピークとして出るので、
            # 「1 compound = 1 ピーク」を系列単位に緩める。単一系列の
            # クラスは series が None になり、従来どおり compound 全体で
            # 1 本だけが残る。
            residuals = []
            for idx, row in group.iterrows():
                cls  = row['lipid_class']
                link = row['_link']
                M    = coherence_models.get((cls, link))
                if M is None or M.get('rmse', 0.0) == 0.0:
                    residuals.append((idx, np.nan, None))
                    continue
                _p, _r, _s = _coherence_best_pred(
                    M, row['_C'], row['_U'], row['obs_rt'])
                residuals.append((idx, _r, _s))

            valid = [(i, r, s) for i, r, s in residuals
                     if r is not None and not np.isnan(r)]
            if not valid:
                continue  # 該当 link のモデルなし → 全保留

            _best_by_series: dict = {}
            for i, r, s in valid:
                _cur = _best_by_series.get(s)
                if _cur is None or r < _cur[1]:
                    _best_by_series[s] = (i, r)
            _keep_idx = {v[0] for v in _best_by_series.values()}
            for idx, _r, _s in residuals:
                if idx not in _keep_idx:
                    df.at[idx, 'coherence_status'] = STATUS_REJ_RESIDUAL

        df = df.drop(columns=['_C', '_U', '_link'])

    # 判定 C: 勝者が 判定 B で落ちたスポットを救済する。
    #   判定 A で勝者にした行が、判定 B で「同じ compound の別ピークの方が
    #   合う」と判断されて REJ_RESIDUAL になることがある。そうなると勝者も
    #   敗者も落ちて、そのスポットのアノテーションが丸ごと消える。
    #   検証データでは RT 2.8349 / m/z 504.3097 がこれで、⑥ の勝者
    #   LPE 20:3 が別ピークに負けて落ち、敗者 LPA 22:4 も REJ_LOSER のまま
    #   残っていた(fix54 までは △ だったので LPA 22:4 が kept で出ていた)。
    #   ノンターゲット解析では「誰も残らない」より「次点を低い信頼度で
    #   残す」方が使いやすいので、次点以降を残差順に 1 つだけ戻す。
    if coherence_assignments and coherence_models:
        _rt_r = df['obs_rt'].round(6)
        _mz_r = df['obs_mz'].round(6)
        _by_coord: dict = {}
        for i in df.index:
            if not bool(df.at[i, 'matched']):
                continue
            _by_coord.setdefault(
                (float(_rt_r[i]), float(_mz_r[i])), []).append(i)
        # 「その compound は別のピークで採用済み」の判定は、⑥ だけでなく
        # 他のフィルタも通っている行に限る。⑤ で落ちている行の compound を
        # 「採用済み」と数えると、本当はどこにも使われていない compound が
        # 救済対象から外れてしまう(検証データの LPA 22:4 がこれだった)。
        def _series_of(_i):
            """その行がどの系列に属するか。二峰性でなければ None。"""
            _c = str(df.at[_i, 'compound'])
            _C2, _U2 = _parse_c_and_u(_c)
            if _C2 is None:
                return None
            _M2 = coherence_models.get(
                (str(df.at[_i, 'lipid_class']), _link_kind(_c)))
            if _M2 is None:
                return None
            return _coherence_best_pred(
                _M2, _C2, _U2, df.at[_i, 'obs_rt'])[2]

        try:
            _avail = (
                (df['coherence_status'] == STATUS_KEPT)
                & (df.get('matched', False) == True)
                & (df.get('is_filter_status', STATUS_NA) != STATUS_REJ_WINDOW)
                & (df.get('rt_outlier_status', STATUS_NA)
                   != STATUS_REJ_OUTLIER)
                & (df.get('adduct_filter_status', STATUS_NA)
                   != STATUS_REJ_ADDUCT)
                & (df.get('manual_status', STATUS_NA) != STATUS_REJ_MANUAL))
            # 「採用済み」を (compound, 系列) 単位で見る。DG の
            # 早い方が既に採用されていても、遅い方の枠はまだ空いている。
            _kept_compounds = set(
                (str(df.at[_i, 'compound']), _series_of(_i))
                for _i in df.loc[_avail].index)
        except Exception:
            _kept_compounds = set()
        _promoted: list[str] = []
        for coord, a in (coherence_assignments or {}).items():
            idxs = _by_coord.get((float(coord[0]), float(coord[1])))
            if not idxs:
                continue
            # そのスポットに生き残りがあるなら何もしない
            if any(str(df.at[i, 'coherence_status']) == STATUS_KEPT
                   for i in idxs):
                continue
            # 手動で決めてあるスポットは触らない
            if any(str(df.at[i, 'manual_status'])
                   in (STATUS_KEPT, STATUS_REJ_MANUAL) for i in idxs):
                continue
            # 候補の選び方は 2 段階。モデルがある候補は RT 残差で選ぶ。
            # 1 つもモデルが無いときは、何も残さないより m/z 誤差が
            # 最小の候補を残す方が使えるので、それを最後の手段にする。
            best = None       # ((tier, score), idx, compound, how, series)
            for i in idxs:
                if (str(df.at[i, 'is_filter_status']) == STATUS_REJ_WINDOW
                        or str(df.at[i, 'rt_outlier_status'])
                        == STATUS_REJ_OUTLIER
                        or str(df.at[i, 'adduct_filter_status'])
                        == STATUS_REJ_ADDUCT):
                    continue
                _cmp = str(df.at[i, 'compound'])
                _ser = _series_of(i)
                if (_cmp, _ser) in _kept_compounds:
                    continue   # その compound × 系列は別ピークで採用済み
                _C, _U = _parse_c_and_u(_cmp)
                _M = None
                if _C is not None:
                    _M = coherence_models.get(
                        (str(df.at[i, 'lipid_class']), _link_kind(_cmp)))
                if _M is not None and _M.get('rmse'):
                    _key = (0, _coherence_best_pred(
                        _M, _C, _U, df.at[i, 'obs_rt'])[1])
                    _how = 'residual'
                else:
                    try:
                        _d = abs(float(df.at[i, 'obs_mz'])
                                 - float(df.at[i, 'theoretical_mz']))
                    except (KeyError, TypeError, ValueError):
                        _d = float('inf')
                    _key = (1, _d)
                    _how = 'mz'
                if best is None or _key < best[0]:
                    best = (_key, i, _cmp, _how, _ser)
            if best is None:
                continue
            df.at[best[1], 'coherence_status'] = STATUS_KEPT
            _kept_compounds.add((best[2], best[4]))
            try:
                a['promoted_compound'] = best[2]
                a['promoted_class'] = str(df.at[best[1], 'lipid_class'])
                a['promoted_by'] = best[3]
            except Exception:
                pass
            _promoted.append(
                f"{coord[0]:.3f}/{coord[1]:.4f}→{best[2]}({best[3]})")
        if _promoted:
            log.info("[Coherence Filter] promoted after residual rejection "
                  f"({len(_promoted)}): " + ", ".join(_promoted[:8])
                  + (" …" if len(_promoted) > 8 else ""))

    return _update_final_status_df(df)


def _update_final_status_df(df: pd.DataFrame) -> pd.DataFrame:
    """status 列群から final_status をベクトル化再計算する。

    優先順位:
      manual_status == STATUS_KEPT     → 'kept'(他のフィルタを上書き)
      manual_status == STATUS_REJ_MANUAL → 'rejected'
      その他、いずれかの rejection ステータス → 'rejected'
      それ以外 → 'kept'(matched=True 限定)
    """
    if df.empty:
        return df
    if 'final_status' not in df.columns:
        df = _add_status_columns(df)
    rejected_any = (
        (df.get('is_filter_status', STATUS_NA) == STATUS_REJ_WINDOW) |
        (df.get('adduct_filter_status', STATUS_NA) == STATUS_REJ_ADDUCT) |
        (df.get('rt_outlier_status', STATUS_NA) == STATUS_REJ_OUTLIER) |
        (df.get('coherence_status', STATUS_NA).isin(
            [STATUS_REJ_RESIDUAL, STATUS_REJ_LOSER]))
    )
    manual_rej  = df.get('manual_status', STATUS_NA) == STATUS_REJ_MANUAL
    manual_kept = df.get('manual_status', STATUS_NA) == STATUS_KEPT
    matched_true = df.get('matched', False) == True
    df['final_status'] = np.where(
        ~matched_true, STATUS_NA,
        np.where(manual_rej, 'rejected',
                 np.where(manual_kept, STATUS_KEPT,
                          np.where(rejected_any, 'rejected', STATUS_KEPT))))
    return df


def _expand_to_multi_label(
    match_df: pd.DataFrame,
    obs_mz: np.ndarray,
    obs_rt: np.ndarray,
    ppm_tol: float,
    class_ref_rt: dict | None = None,
    is_tol: float = 0.10,
    apply_is_filter: bool = True,
    reserved_peak_idx: set | None = None,
    void_rt: float = 0.0,
) -> pd.DataFrame:
    """既存の 1 行/エントリの match_df を、ppm 内全候補の多行 df に展開。

    各エントリについて _collect_match_candidates で全候補を取得し、
    isomer_rank を 1(m/z 最近接、既存挙動と一致)から順に付与。

    apply_is_filter=True の場合: IS RT 窓で先に絞ってから候補展開する
。

    reserved_peak_idx が指定されたとき、IS-bearing クラスが
    第1パスで確保したピーク (= reserved_idx) に展開された
    「非 IS-bearing かつ手動RT 未指定」クラスの候補は STATUS_REJ_WINDOW に
    マークする。これにより SM (IS あり) 帰属ピーク上の CerPE (IS なし) 候補が
    KEPT のまま残り、calc_conflicts で SM vs CerPE の偽の衝突になる問題を回避する。
    """
    if match_df.empty:
        return _add_status_columns(match_df)

    expanded: list[dict] = []
    class_ref_rt = class_ref_rt or {}
    reserved_peak_idx = reserved_peak_idx or set()

    for _, row in match_df.iterrows():
        cls = row.get('lipid_class')
        th  = row.get('theoretical_mz')
        if pd.isna(th) or th is None or th <= 0:
            # 無効エントリはそのまま保持
            d = dict(row)
            d.update({
                'isomer_rank': 0,
                'is_filter_status': STATUS_NA,
                'rt_outlier_status': STATUS_NA,
                'coherence_status': STATUS_NA,
                'manual_status': STATUS_NA,
                'final_status': STATUS_NA,
            })
            expanded.append(d)
            continue

        # 全候補取得
        cands = _collect_match_candidates(
            float(th), obs_mz, obs_rt, ppm_tol, void_rt=void_rt)

        # IS Filter 用の窓判定(削除はせずマークする)
        is_filter_active = apply_is_filter and cls in class_ref_rt
        if is_filter_active:
            ref_rt = float(class_ref_rt[cls][0])
            rt_tol_cls = float(class_ref_rt[cls][1])

        if not cands:
            # 候補なし → 1 行残す(matched=False)
            d = dict(row)
            d.update({
                'obs_mz':            np.nan,
                'obs_rt':            np.nan,
                'delta_ppm':         np.nan,
                'matched':           False,
                'matched_idx':       -1,
                'isomer_rank':       0,
                'is_filter_status':  STATUS_NA,
                'rt_outlier_status': STATUS_NA,
                'coherence_status':  STATUS_NA,
                'manual_status':     STATUS_NA,
                'final_status':      STATUS_NA,
            })
            expanded.append(d)
            continue

        # 非 IS-bearing かつ手動RT 未指定クラスの判定。
        # cls が class_ref_rt に存在しない (= auto_is でも manual でもない)
        # 場合、IS-bearing クラスに reserved されたピークへの展開は無効化する。
        block_at_reserved = (
            apply_is_filter
            and bool(reserved_peak_idx)
            and (cls not in class_ref_rt)
        )

        # 全候補を展開(rank 1 = m/z 最近接)
        for rank, c in enumerate(cands, start=1):
            d = dict(row)

            # IS Filter mark
            if is_filter_active:
                if abs(c['obs_rt'] - ref_rt) <= rt_tol_cls:
                    is_status = STATUS_KEPT
                    final_st  = STATUS_KEPT
                else:
                    is_status = STATUS_REJ_WINDOW
                    final_st  = 'rejected'
            elif block_at_reserved and c['peak_idx'] in reserved_peak_idx:
                # IS-bearing クラスが確保したピーク上の
                # 非 IS-bearing 候補は IS Filter で reject 扱いとする。
                is_status = STATUS_REJ_WINDOW
                final_st  = 'rejected'
            else:
                is_status = STATUS_NA
                final_st  = STATUS_KEPT

            d.update({
                'obs_mz':            c['obs_mz'],
                'obs_rt':            c['obs_rt'],
                'delta_ppm':         c['delta_ppm'],
                'matched':           True,
                'matched_idx':       c['peak_idx'],
                'isomer_rank':       rank,
                'is_filter_status':  is_status,
                'rt_outlier_status': STATUS_NA,
                'coherence_status':  STATUS_NA,
                'manual_status':     STATUS_NA,
                'final_status':      final_st,
            })
            expanded.append(d)

    return pd.DataFrame(expanded)


def calc_rt_outliers(
    match_df: pd.DataFrame,
    mad_k: float,
    min_n: int = 5,
    exempt_is: bool = True,
    is_exempt_classes: set | None = None,
    class_ref_rt: dict | None = None,
    verbose: bool = True,
) -> dict[str, set]:
    """
    クラスごとにRT方向のMAD外れ値を計算する。

    MAD (Median Absolute Deviation) ベースの外れ値検出:
      |RT - median| > mad_k * MAD * 1.4826 の点を外れ値とする。

    IQR法では、同一方向に複数の外れ値が集まる場合に Q1/Q3 が引き寄せられ
    外れ値を検出できなくなる（masking問題）。MAD法は中央値ベースのため
    この問題に対してロバストである。

    Parameters
    ----------
    exempt_is : bool
        True(既定)なら **IS を持つクラスには ④ を適用しない**。

        ④ の目的は「**IS の無い脂質クラス**の極端な外れ値を除く」こと。
        IS のあるクラスは ③ IS Filter が class_ref_rt で RT 窓を
        既に決めているので、④ が重ねて RT を見る必要がない。

        fix22 では「IS の**行だけ**フラグしない」という行単位の実装に
        していたが、これは指示の読み違いだった。行単位だと、下記の
        事故が起きる:

          検証 session(2026-08-18)実測。pos の HexCer 47 行は
          O2 系(非水酸化、中央値 2.050)と O3 系(水酸化、中央値 2.306)
          の 2 群からなる。IS は O2 系(2.063)なので ③ が
          class_ref_rt=(2.063, ±0.1) を置き、O3 系 27 行を
          rejected_window で落とす。ところが ④ は match_df の
          **matched 行すべて**(= ③ で落とした行も含む)で中央値を取るので
          2.2755 になり、③ が残した O2 系(2.019〜2.038)を
          「中央値から 0.24〜0.26 離れている」という理由で消してしまう。
          結果、HexCer は IS ともう 1 点だけになり画面から消えた。

        クラス単位の免除にすれば、この ③ と ④ の食い違いは
        IS のあるクラスでは起きなくなる。
        (IS の無いクラスでは ④ が唯一の RT 判定なので従来どおり働く)

        免除の条件を「IS を持つ」から「**③ が実際にその
        クラスの RT を拘束している**」に変えた。上の根拠は
        「③ が RT 窓を決めているから ④ は要らない」なので、③ を
        切ったクラスでは前提が成り立たず、④ も効かないと RT 方向の
        判定が **1 つも無くなる**(DG の事例)。

        ③ が効いているかは、呼び出し側から渡された class_ref_rt
        (③ が置いた RT 窓の辞書)で判定する。渡されない場合は
        match_df の is_filter_status に **kept** があるかで代用する
        (kept が付くのは ③ が窓を置いたクラスだけ。reserved ピーク
         由来の rejected_window は窓とは無関係なので数えない)。

          - ③ が窓を置いたクラス   → class_ref_rt に入る / kept が付く
          - ③ 未実行               → 入らない / 全行 n/a
          - IS 未検出              → 入らない
          - IS の Use を全部外した → 入らない
        つまり ③ の Preview で IS のチェックを外した場合も自動的に
        ④ の対象に戻る。

    is_exempt_classes : set | None
        Advanced Setting → Filter Classes → ③ IS Filter で
        チェックを外したクラス(= ③ の rejected を剥がすクラス)。
        ここに入っているクラスは ③ が効いていないものとして扱い、
        ④ を適用する。

    class_ref_rt : dict | None
        ③ が置いた {class: (ref_rt, tol, kind, ...)}。
        渡されればこれを ③ の有効判定に使う(最も直接的な指標)。

    verbose : bool
        True(既定)なら、IS のあるクラスについて ④ を適用したか
        免除したかを 1 行で表示する。

    戻り値: {lipid_class: set of (rt_r6, mz_r6)}
    """
    outliers: dict[str, set] = {}
    matched = match_df[match_df['matched']]
    has_is_col = 'is_IS' in matched.columns
    has_isf_col = 'is_filter_status' in matched.columns
    ex_cls = set(is_exempt_classes or ())
    _exempted: list[str] = []      # IS があり ③ も効いている → ④ 免除
    _applied: list[str] = []       # IS はあるが ③ が効いていない → ④ 適用
    for cls, grp in matched.groupby('lipid_class'):
        if len(grp) < min_n:
            outliers[cls] = set()
            continue
        # IS を持つクラスは丸ごと対象外。③ が class_ref_rt で
        # RT 窓を決めているので、④ が重ねて判定する必要がない。
        # ただし **③ が実際にそのクラスに効いているとき限り**。
        #   ③ を切ったクラス(Filter Classes の除外 / IS の Use を外した /
        #   IS 未検出 / ③ 未実行)では RT 窓が無いので、④ を適用しないと
        #   RT 方向の判定が 1 つも無くなる。
        if exempt_is and has_is_col and bool(
                grp['is_IS'].fillna(False).to_numpy(dtype=bool).any()):
            if class_ref_rt is not None:
                _win = cls in class_ref_rt
            else:
                # 代用指標: is_filter_status に kept があるか。
                # kept が付くのは ③ が窓を置いたクラスだけ。
                # reserved ピーク由来の rejected_window は窓と無関係
                # なので、!= n/a では判定できない。
                _win = has_isf_col and bool(
                    (grp['is_filter_status'].astype(str)
                     == STATUS_KEPT).any())
            _is_filter_on = (cls not in ex_cls) and _win
            if _is_filter_on:
                _exempted.append(str(cls))
                outliers[cls] = set()
                continue
            _applied.append(str(cls))
        rts  = grp['obs_rt'].values
        med  = np.median(rts)
        mad  = np.median(np.abs(rts - med))
        if mad == 0.0:
            # 全エントリが同一RTの場合は外れ値なし
            outliers[cls] = set()
            continue
        threshold = mad_k * mad * 1.4826   # 正規分布換算の一貫性係数
        flag = np.abs(grp['obs_rt'].values - med) > threshold
        # IS 行は外れ値にしない。fix30 で ③ を切ったクラスもここに
        # 来るようになったが、IS を消すとそのクラスの定量が丸ごと
        # 壊れるので、IS 行の保護は維持する。
        # (DG のように sn 異性体で IS が 2 本ある場合、
        #  どちらも本物なので消してはいけない)
        if exempt_is and has_is_col:
            flag = flag & ~(grp['is_IS'].fillna(False).to_numpy(dtype=bool))
        out_rows  = grp[flag]
        outliers[cls] = set(
            zip(out_rows['obs_rt'].round(6), out_rows['obs_mz'].round(6)))
    if verbose and (_exempted or _applied):
        # どのクラスに ④ が効いたのかを必ず出す。
        # 「IS があるのに ④ が適用された」ことに気付けるようにするため。
        _msg = f"[RT Outlier Filter] IS-bearing classes: exempt (③ active) {len(_exempted)}"
        if _applied:
            _msg += (f" / ③ disabled so ④ applied to {len(_applied)} "
                     f"({', '.join(sorted(_applied))})")
        log.info(_msg)
    return outliers


def _parse_carbon_count(compound: str) -> int | None:
    """脂質化合物名から総炭素数 C を抽出する。

    対応形式:
        "PC 34:1"         → 34
        "SM 34:1;O2"      → 34
        "CerPE 30:1 O2"   → 30
        "DG 34:1 "        → 34（末尾スペース対応）

    パース不能な場合は None を返す。
    """
    if not compound:
        return None
    m = re.match(r'^[\w()]+\s+(\d+):\d+', compound.strip())
    if m:
        return int(m.group(1))
    return None


def calc_conflicts(
    match_df: pd.DataFrame,
) -> tuple[set, dict[tuple, list[dict]]]:
    """
    同一観測ピークに2つ以上の異なる脂質クラスが割り当てられた座標を返す。
    同一クラス内での重複（複数エントリが同一ピークにマッチ）は対象外。

    final_status='kept' の行のみで衝突を検出する。
    IS Filter で REJ_WINDOW された(IS RT 窓外で帰属が信頼できない)
    候補は衝突から除外し、HexCer が PC/PE の RT 領域に紛れ込むような
    化学的にあり得ない conflict を排除する。

    戻り値:
      - conflict_coords: set of (rt_r6, mz_r6)
      - conflict_map:    {(rt_r6, mz_r6): [{'lipid_class', 'compound', 'adduct'}, ...]}
                         同一座標に帰属された全エントリ情報
    """
    matched = match_df[match_df['matched']]
    # final_status が存在すれば KEPT 行のみで衝突判定。
    # 存在しない(まだ IS Filter 前)場合は matched=True 全行で従来通り。
    if 'final_status' in matched.columns:
        matched = matched[matched['final_status'] == STATUS_KEPT]
    if matched.empty:
        return set(), {}

    coords = list(zip(
        matched['obs_rt'].round(6),
        matched['obs_mz'].round(6)))
    tmp = matched[['lipid_class', 'compound', 'adduct']].copy()
    tmp['coord'] = pd.Series(coords, index=matched.index)

    # 異なるクラスが2つ以上ある座標のみ抽出
    n_classes = tmp.groupby('coord')['lipid_class'].nunique()
    conflict_coords = set(n_classes[n_classes > 1].index)

    conflict_map: dict[tuple, list[dict]] = {}
    for coord, grp in tmp[tmp['coord'].isin(conflict_coords)].groupby('coord'):
        conflict_map[coord] = [
            {'lipid_class': r['lipid_class'],
             'compound':    r['compound'],
             'adduct':      r['adduct']}
            for _, r in grp.iterrows()
        ]
    return conflict_coords, conflict_map


def calc_adduct_filter(
    conflict_coords: set,
    conflict_map: dict,
    fingerprints: dict,
    pos_rt: np.ndarray | None,
    pos_mz: np.ndarray | None,
    neg_rt: np.ndarray | None,
    neg_mz: np.ndarray | None,
    ppm_tol: float,
    rt_win: float,
    int_floor: float = 0.0,
    vote_threshold: float = 1.0,
    conflict_mode: str = 'pos',
    cross_mode_rt_shift: float = 0.0,
    cross_mode_rt_win: float | None = None,
    pos_int: np.ndarray | None = None,
    neg_int: np.ndarray | None = None,
    ratio_threshold: float = ADDUCT_RATIO_WIN_FRAC,
) -> dict:
    """Adduct Ion Filter のスコアリングロジック本体。

    各衝突 (RT, m/z) について、IS 由来 fingerprint を持つ候補クラスを
    対象に重み付きスコアリングを行い、勝者を決定する。

    引数:
      conflict_coords : set of (rt_r6, mz_r6)
      conflict_map    : {coord: [{'lipid_class', 'compound', 'adduct'}, ...]}
      fingerprints    : {class_name: {adduct: {'expected': bool, 'weight': float}}}
                        内蔵 CURATED + ランタイム IS 派生 のマージ
      pos_rt / pos_mz : pos モードのピーク表(np.ndarray)
      neg_rt / neg_mz : neg モードのピーク表(np.ndarray、無ければ None)
      ppm_tol         : アダクト観測判定の m/z 許容誤差
      rt_win          : 共溶出 RT 窓
      int_floor       : 観測判定の強度しきい値(現バージョン未使用、将来拡張用)
      vote_threshold  : スコア差で勝者を確定するための閾値

    戻り値:
      {
        coord: {
          'winner':       str | None,
          'losers':       set[str],
          'scores':       {class: float},
          'observations': {adduct: bool},
          'method':       str,   # 'ratio' / 'count' / 'undecided' / 'skipped'
          'fracs':        {class: float},  # 強度比から推定した寄与率
          'est':          {class: float},  # ピーク強度スケールの推定寄与
          'est_from':     {class: adduct}, # 推定に使った判別アダクト
          'peak_intensity': float,
          'unjudgeable':  set[str],        # fingerprint が無く採点できない候補
          'skipped':      bool,            # 採点可能な候補が 2 つ未満
        },
        ...
      }

    仕様:
      - fingerprint を持つ候補だけで採点する。持たない候補が
        混じっていても座標ごと放棄しない(採点可能が 2 つ未満のときだけ
        skipped=True で ⑥ に流す)。
      - 採点は「判別アダクトの実測強度を IS 由来の強度比で割り
        戻し、各クラスの寄与を推定する」方式(method='ratio')。強度が
        取れない / 整合しない座標は従来の本数ルール(method='count')。
      - 各 adduct について(本数ルール):
          全員が期待 (共通 adduct) → 判別寄与なし、skip
          誰も期待しない → 判別寄与なし、skip
          一部のみ期待 (= 判別 adduct):
            期待するクラスは 観測あり → +1、観測なし → -1
            期待しないクラスは影響なし
      - これで「判別 adduct の観測 = 強い証拠、不在 = 反対の証拠」のシンプル
        な化学的判定になる。weight は不使用。
      - スコア差 >= vote_threshold で勝者確定。未満なら判定保留(winner=None)
    """
    pos_set = set(ADDUCT_SET_POS_V54)
    results: dict = {}

    for coord in conflict_coords:
        rt, mz = coord
        entries = conflict_map.get(coord, [])
        if not entries:
            continue
        all_classes = sorted(set(e['lipid_class'] for e in entries))
        if len(all_classes) < 2:
            continue  # 衝突候補が 1 つなら判別不要

        # fingerprint を持つ候補だけを採点対象にする。
        #   fix49 までは「1 つでも持たない候補があれば座標ごと放棄」
        #   だったため、PA のようにその極性で IS が出ないクラスが 1 つ
        #   混じるだけで、PC/PE という判定可能な組まで ⑤ が黙っていた
        #   (検証データで 100 座標。うち 64 は採点可能だった)。
        candidate_classes = [c for c in all_classes if c in fingerprints]
        unjudgeable = [c for c in all_classes if c not in fingerprints]
        if len(candidate_classes) < 2:
            results[coord] = {
                'winner':       None,
                'losers':       set(),
                'scores':       {},
                'observations': {},
                'skipped':      True,
                'unjudgeable':  set(unjudgeable),
                'method':       'skipped',
            }
            continue

        # 候補クラスごとに自分の matched adduct から neutral
        # mass を計算し、その neutral mass で各 adduct の観測を判定する。
        # 異なる adduct で衝突するケース(例: PC[M+H]+ vs PI[M+Na]+)で、
        # 各クラスの discriminator adduct を正しい m/z で探せるようにする。
        class_to_adduct: dict[str, str] = {}
        for e in entries:
            class_to_adduct.setdefault(
                e['lipid_class'],
                ADDUCT_ALIAS.get(e['adduct'], e['adduct']))
        # アダクトが未知のクラスは採点対象から外すだけにする
        # (座標ごと放棄はしない)
        _unknown = [c for c in candidate_classes
                    if class_to_adduct.get(c) not in ADDUCT_OFFSETS]
        if _unknown:
            unjudgeable = sorted(set(unjudgeable) | set(_unknown))
            candidate_classes = [c for c in candidate_classes
                                 if c not in _unknown]
        if len(candidate_classes) < 2:
            results[coord] = {
                'winner': None, 'losers': set(),
                'scores': {}, 'observations': {}, 'skipped': True,
                'unjudgeable': set(unjudgeable), 'method': 'skipped',
            }
            continue
        # 各クラスの neutral mass
        class_neutrals = {
            c: mz - ADDUCT_OFFSETS[class_to_adduct[c]]
            for c in candidate_classes
        }

        # 観測判定: (class, adduct) ごとにそのクラスの neutral mass で判定
        cls_observations: dict[tuple, bool] = {}
        # 強度比の計算にクラス別の実測強度が要る
        cls_intensities: dict[tuple, float] = {}
        # ログ用: 全クラス共通 neutral(代表)で 9 adduct の素朴な観測も保持
        observations: dict[str, bool] = {a: False for a in ADDUCT_SET_V54}
        # 観測ピーク強度(adduct → max intensity across classes)
        observation_intensities: dict[str, float] = {
            a: 0.0 for a in ADDUCT_SET_V54}
        for c in candidate_classes:
            neutral = class_neutrals[c]
            for adduct in ADDUCT_SET_V54:
                offset = ADDUCT_OFFSETS.get(adduct)
                if offset is None:
                    cls_observations[(c, adduct)] = False
                    continue
                mz_theo = neutral + offset
                target_is_pos = (adduct in pos_set)
                if target_is_pos:
                    rt_arr, mz_arr = pos_rt, pos_mz
                    int_arr = pos_int
                else:
                    rt_arr, mz_arr = neg_rt, neg_mz
                    int_arr = neg_int
                if rt_arr is None or mz_arr is None or len(rt_arr) == 0:
                    cls_observations[(c, adduct)] = False
                    continue

                conflict_in_pos = (conflict_mode == 'pos')
                cross_mode = (conflict_in_pos != target_is_pos)
                if cross_mode:
                    if conflict_in_pos:
                        target_rt = rt + cross_mode_rt_shift
                    else:
                        target_rt = rt - cross_mode_rt_shift
                    use_tol = (cross_mode_rt_win
                               if cross_mode_rt_win is not None
                               else rt_win)
                else:
                    target_rt = rt
                    use_tol = rt_win

                ppm_err = np.abs(mz_arr - mz_theo) / mz_theo * 1e6
                rt_err = np.abs(rt_arr - target_rt)
                mask = (ppm_err <= ppm_tol) & (rt_err <= use_tol)
                is_obs = bool(np.any(mask))
                cls_observations[(c, adduct)] = is_obs
                if is_obs:
                    observations[adduct] = True
                    # 観測ピーク群の最大 intensity を記録
                    if int_arr is not None:
                        try:
                            cand_idx = np.where(mask)[0]
                            cand_int = np.asarray(int_arr)[cand_idx]
                            if cand_int.size:
                                peak_int = float(np.nanmax(cand_int))
                                cls_intensities[(c, adduct)] = peak_int
                                if peak_int > observation_intensities[adduct]:
                                    observation_intensities[adduct] = peak_int
                        except Exception:
                            pass

        # ── 強度比による寄与推定 ────────────────────────
        # IS の fingerprint は「そのクラスが出すアダクトの強度比」を
        # 持っている(weight、最強 = 1.00)。判別アダクトの実測強度を
        # その比で割り戻すと「この観測ピークのうちクラス c に由来する
        # 分」が推定できる。
        #     est_c = I(a_d) × w_c(a0) / w_c(a_d)
        #   a0  … そのクラスが衝突ピークに当たっているアダクト
        #   a_d … c だけが期待する判別アダクト(複数あれば最大の推定値)
        # 検証データでは Σ est / 実測ピーク強度 = 0.97(中央値)、
        # ⑥ の RT 残差との一致 93%、相関 +0.80。
        disc_map: dict[str, list[str]] = {}
        for adduct in ADDUCT_SET_V54:
            _ex = [c for c in candidate_classes
                   if fingerprints[c].get(adduct, {}).get('expected', False)]
            if 0 < len(_ex) < len(candidate_classes):
                disc_map[adduct] = _ex

        est: dict[str, float] = {c: 0.0 for c in candidate_classes}
        est_from: dict[str, str] = {}
        for c in candidate_classes:
            a0 = class_to_adduct.get(c)
            try:
                w0 = float(fingerprints[c].get(a0, {}).get('weight', 0.0) or 0.0)
            except (TypeError, ValueError):
                w0 = 0.0
            for adduct, _ex in disc_map.items():
                if c not in _ex:
                    continue
                if not cls_observations.get((c, adduct), False):
                    continue
                try:
                    wd = float(fingerprints[c].get(adduct, {})
                               .get('weight', 0.0) or 0.0)
                except (TypeError, ValueError):
                    wd = 0.0
                if wd <= 0.0 or w0 <= 0.0:
                    continue
                v = cls_intensities.get((c, adduct), 0.0) * w0 / wd
                if v > est[c]:
                    est[c] = v
                    est_from[c] = adduct

        # 衝突ピーク本体の強度(どのクラスから見ても同じピーク)
        peak_intensity = max(
            (cls_intensities.get((c, class_to_adduct.get(c)), 0.0)
             for c in candidate_classes), default=0.0)
        tot_est = float(sum(est.values()))
        fracs = {c: (est[c] / tot_est if tot_est > 0 else 0.0)
                 for c in candidate_classes}
        ratio_ok = bool(
            tot_est > 0.0 and peak_intensity > 0.0
            and ADDUCT_RATIO_CONSIST_LO
            <= tot_est / peak_intensity <= ADDUCT_RATIO_CONSIST_HI)

        # ── 代替: 本数ベースの採点(fix49 まではこれだけ) ──────
        # 「期待 adduct が観測されない = 存在しない」とは限らない(感度
        # 不足で拾えていないだけのことがある)ので、観測なしは減点しない。
        scores: dict[str, float] = {c: 0.0 for c in candidate_classes}
        for adduct in ADDUCT_SET_V54:
            for c in candidate_classes:
                if not fingerprints[c].get(adduct, {}).get('expected', False):
                    continue
                if cls_observations.get((c, adduct), False):
                    scores[c] += 1.0

        winner = None
        method = 'undecided'
        if ratio_ok:
            _top = max(fracs, key=lambda c: fracs[c])
            if fracs[_top] >= ratio_threshold:
                winner = _top
                method = 'ratio'
        if winner is None:
            sorted_scores = sorted(scores.items(), key=lambda x: -x[1])
            top_cls, top_score = sorted_scores[0]
            second_score = (sorted_scores[1][1]
                            if len(sorted_scores) > 1 else 0.0)
            if top_score - second_score >= vote_threshold:
                winner = top_cls
                method = 'count'
        # 勝者が決まったら、採点できなかった候補も含めて棄却する。
        # 1 つの観測ピークに帰属できる分子種は 1 つ、という ⑥/fix49 と
        # 同じ建前を ⑤ でも通す。
        if winner is not None:
            losers = set(c for c in all_classes if c != winner)
        else:
            losers = set()

        # 判別 adduct の集計(誰のみが期待するか)
        discriminator_map: dict[str, list[str]] = {}
        for adduct in ADDUCT_SET_V54:
            expects_classes = [
                c for c in candidate_classes
                if fingerprints[c].get(adduct, {}).get('expected', False)
            ]
            if 0 < len(expects_classes) < len(candidate_classes):
                discriminator_map[adduct] = expects_classes

        # 各 (cls, adduct) の観測フラグも詳細化(class ごとの neutral mass で判定済)
        per_class_obs = {
            c: {a: cls_observations.get((c, a), False)
                for a in ADDUCT_SET_V54}
            for c in candidate_classes
        }

        results[coord] = {
            'winner':              winner,
            'losers':               losers,
            'scores':               scores,
            'method':               method,
            'fracs':                fracs,
            'est':                  est,
            'est_from':             est_from,
            'peak_intensity':       peak_intensity,
            'unjudgeable':          set(unjudgeable),
            'observations':         observations,
            'observation_intensities': observation_intensities,
            'discriminator_map':    discriminator_map,
            'per_class_observations': per_class_obs,
            'class_to_adduct':      class_to_adduct,
            'skipped':              False,
        }

    return results


# ════════════════════════════════════════════════════════════════════
#  Coherence Engine
#  クラス整合性に基づく多変量トレンド解析と衝突帰属エンジン
#
#  設計思想:
#    同じ link 種(ester/ether/plasmalogen)内で、クラスごとに共通の
#    「C×RT、U×RT」依存性を持つと仮定し、クラス別インターセプトだけが
#    異なる pooled 線形モデル RT = α·C + β·U + γ_class を fit する。
#    この fit は以下 3 つの目的を同時に達成する:
#      (1) 非衝突点の残差 → Coherence 外れ値検出
#      (2) U を covariate に含むため SM の多軌跡問題を解決
#      (3) 衝突点で候補クラス仮定下の予測 RT 差 → 自動帰属の根拠
# ════════════════════════════════════════════════════════════════════


def _parse_c_and_u(compound: str) -> tuple[int | None, int | None]:
    """脂質化合物名から総炭素数 C と不飽和度 U を抽出する。

    対応形式:
        "PC 34:1"        → (34, 1)
        "PC O-34:1"      → (34, 1)   # ether 連結(O-)の数値部を取る
        "PC P-36:1"      → (36, 1)   # plasmalogen(P-)
        "SM 34:1;O2"     → (34, 1)   # ;O2 はヒドロキシ数、無関係
        "Cer 18:1;O2/16:0" → (18, 1) # 最初の N:M を採用

    パース不能な場合は (None, None) を返す。
    """
    if not compound:
        return None, None
    # 最初に現れる N:M パターンを採用。O-/P- 接頭辞はオプショナル。
    m = re.search(r'(?:[OP]-)?(\d+):(\d+)', compound.strip())
    if m:
        return int(m.group(1)), int(m.group(2))
    return None, None


def _link_kind(compound: str) -> str:
    """脂質化合物名から側鎖結合種を判定する。

    "ester"       : 通常のエステル結合(デフォルト)
    "ether"       : "クラス名 O-34:1" 形式
    "plasmalogen" : "クラス名 P-36:1" 形式

    SFC での RT 挙動は link 種で系統的に異なるため、Coherence 回帰は
    link 別にモデルを持つ必要がある(実データ検証で ester/ether で
    dRT/dC, dRT/dU が異なることを確認済み)。
    """
    if not compound:
        return 'ester'
    # 化合物名のクラス部直後のスペースに続く "O-" / "P-" を検出
    m = re.search(r'\s([OP])-\d+:\d+', compound.strip())
    if m:
        return 'ether' if m.group(1) == 'O' else 'plasmalogen'
    return 'ester'


def _detect_collision_pairs(
    lib_df: pd.DataFrame,
    ppm_tol: float,
) -> list[tuple[str, str]]:
    """ライブラリから m/z 衝突する可能性のあるクラスペアを自動検出する。

    任意の 2 エントリが |m/z_A - m/z_B| ≤ ppm_tol * m/z / 1e6 を満たす場合、
    それらのクラスをペアとして登録する(同一クラスは除外)。

    戻り値: [(class_A, class_B), ...] (class_A < class_B の辞書順)

    用途: Coherence Engine のペア選択ダイアログ初期化(UI 側で使用)。
    """
    if lib_df.empty:
        return []
    sorted_df = lib_df.sort_values('theoretical_mz').reset_index(drop=True)
    mzs = sorted_df['theoretical_mz'].values
    classes = sorted_df['lipid_class'].values
    pairs: set[tuple[str, str]] = set()
    n = len(sorted_df)
    for i in range(n):
        mz_i = mzs[i]
        tol_i = mz_i * ppm_tol / 1e6
        # 昇順ソート済みなので、j > i で mz_j > mz_i + tol_i になったら打ち切り
        for j in range(i + 1, n):
            if mzs[j] - mz_i > tol_i:
                break
            if classes[i] != classes[j]:
                pair = tuple(sorted([classes[i], classes[j]]))
                pairs.add(pair)
    return sorted(pairs)


def _robust_multivar_fit(
    X: np.ndarray,
    y: np.ndarray,
    mad_k: float = 3.0,
    max_iter: int = 10,
) -> tuple[np.ndarray, np.ndarray, float]:
    """多変量版の反復 MAD 除去ロバスト線形回帰。

    既存 `_robust_linear_fit`(単変数)と同じ戦略を多変数に拡張:
      1. 前処理: y の中央値から極端に離れた点を初期除外(5×MAD)
      2. 内点で最小二乗 fit(X は bias 列を含む設計行列)
      3. 残差の MAD から閾値を作り、閾値を超える点を除外
      4. 収束または最大反復まで繰り返し
      5. 内点残差の標準偏差を rmse として返す

    Parameters
    ----------
    X : ndarray, shape (n, p)
        設計行列。インターセプト用の定数列を自前で含める(本エンジンでは
        クラス別ダミー列がインターセプトの役割を果たす)。
    y : ndarray, shape (n,)
    mad_k : float
        除去閾値の MAD 倍率
    max_iter : int

    戻り値:
      (coef, inlier_mask, rmse)
        coef       : shape (p,) の係数
        inlier_mask: shape (n,) の内点マスク
        rmse       : 内点残差の標準偏差
    """
    n = len(y)
    p = X.shape[1]
    if n < p:
        # 劣決定: fit 不能、ゼロ係数で返す
        return np.zeros(p), np.ones(n, dtype=bool), 0.0

    # 前処理: y の極端な外れ値を初期除外
    y_med = float(np.median(y))
    y_mad = float(np.median(np.abs(y - y_med)))
    if y_mad > 0.0:
        init_threshold = 5.0 * y_mad * 1.4826
        init_mask = np.abs(y - y_med) <= init_threshold
    else:
        init_mask = np.ones(n, dtype=bool)
    # 前処理後でも p 点未満になったらフォールバック
    if init_mask.sum() < p:
        init_mask = np.ones(n, dtype=bool)

    mask = init_mask.copy()
    coef = np.zeros(p)
    rmse = 0.0
    for _ in range(max_iter):
        Xm, ym = X[mask], y[mask]
        if len(Xm) < p:
            break
        try:
            coef, *_ = np.linalg.lstsq(Xm, ym, rcond=None)
        except np.linalg.LinAlgError:
            break
        residuals = y - X @ coef
        inlier_resid = residuals[mask]
        if len(inlier_resid) > 0:
            rmse = float(np.std(inlier_resid, ddof=0))
        mad = float(np.median(np.abs(inlier_resid - np.median(inlier_resid))))
        if mad == 0.0:
            break
        threshold = mad_k * mad * 1.4826
        new_mask = np.abs(residuals) <= threshold
        if new_mask.sum() < p:
            break
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask
    return coef, mask, rmse


# ── ⑥ の予測の不確かさ（外挿の度合いを含む）──────────────────
# rmse を link の pooled rmse 側へ縮める強さ。n/(n+k) が自クラスの重み。
# n=4 なら 4/14 = 0.29 しか自分の rmse を信用しない。
COH_SE_SHRINK_N = 10.0
# 学習点の (C, DB) の分散の下限。1 点しか無い方向で h が発散しないように。
COH_SE_C_VAR_FLOOR = 1.0     # 炭素数 1 個分
COH_SE_U_VAR_FLOOR = 0.25    # 二重結合 0.5 本分
# h の上限。これ以上外挿していても se は 6 倍程度で止める。
COH_SE_LEVERAGE_MAX = 35.0
# link 間 RT 差を信用するのに最低限必要なクラス数。
COH_LINK_DELTA_MIN_CLASSES = 2
# 緩い帰属(片側しか採点できなくてもその候補を採る)の既定。
#   ノンターゲット解析では △ のまま残すより、帰属して信頼度で示す方が
#   使いやすい、という判断。
COH_LOOSE_ATTRIBUTION_DEFAULT = True
# ── 二峰性クラス(位置異性体が分離するクラス)の 2 系列 RT モデル ──
# DG は sn 位置異性体が部分分離するため、1 つの和組成が 2 本のピークに
# なる。1 本のモデルでまとめて fit すると rmse が 10〜12 倍悪化し、
# 判定 B が早い方と遅い方をスポットごとにばらばらに残してしまう。
# α と β は 2 群でほぼ同じで違うのは γ だけなので、γ と γ+Δ の 2 系列
# として扱う。対象クラスは Advanced Setting で足せる。

COH_TWO_SERIES_CLASSES_DEFAULT = ('DG',)
# 各系列に最低これだけ学習点が要る(片側だけに寄っていたら採らない)。
COH_TWO_SERIES_MIN_POINTS = 4
# Δ の下限(min)。これ未満は「2 本に分かれている」とみなさない。
COH_TWO_SERIES_MIN_DELTA = 0.015
# Δ の上限(min)。これを超えるなら別のクラスか外れ値を拾っている。
COH_TWO_SERIES_MAX_DELTA = 0.40
# 2 系列にして rmse がこれだけ良くならなければ採らない。この条件が本体で、
# たまたま残差が 2 つに割れただけのクラスでは rmse はほとんど改善しない。
COH_TWO_SERIES_RMSE_GAIN = 2.5
# 定量エクスポートで同じ compound の複数ピークを合算するか。
QUANT_SUM_TWO_SERIES_DEFAULT = True


def _quant_sum_two_series_rows(rows: list, two_series_classes,
                               sample_cols: list) -> list:
    """同じ分子種の複数ピークを 1 行にまとめ、サンプル値を足す。

    DG は sn 位置異性体が部分分離するので 1 つの和組成が 2 スポットとして
    残る(⑥ の 判定 B を系列単位にしたため)。プリカーサー定量では 2 本の
    合算が求める量なので、ここで足し合わせる。IS も同じ規則でまとまる
    ので、正規化の分母も 2 本の和になる。

    まとめるキーは (Class, Compound, Adduct, Mode)。アダクトが違えば
    別の定量イオンなのでまとめない。対象は `two_series_classes` に
    挙げたクラスだけで、ほかのクラスの行はそのまま通す。

    まとめた行には
      `peaks`      … 合算した本数
      `RT (peaks)` … 合算した各ピークの RT(早い順、';' 区切り)
    を残す。信頼度は最も低いものを採り、`sigma` は最小値、
    `decided_by` と `competitors` は重複を除いてつなぐ。
    """
    if not rows:
        return rows
    _cls_set = set(str(c) for c in (two_series_classes or ()))
    if not _cls_set:
        return rows

    def _num(v):
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return None if f != f else f      # NaN は None にする

    # 低い方を採りたいので「悪いほど大きい」順位を振る
    _RANK = {'high': 0, 'medium': 1, 'no conflict': 2, 'low': 3,
             'none': 4}
    out: list = []
    groups: dict = {}
    order: list = []
    for r in rows:
        if str(r.get('Class', '')) not in _cls_set:
            r = dict(r)
            r['peaks'] = 1
            _rt = _num(r.get('RT (min)'))
            r['RT (peaks)'] = ('' if _rt is None else f'{_rt:.4f}')
            out.append(r)
            continue
        key = (str(r.get('Class', '')), str(r.get('Compound', '')),
               str(r.get('Adduct', '')), str(r.get('Mode', '')))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(r)

    for key in order:
        grp = sorted(groups[key], key=lambda d: (
            _num(d.get('RT (min)')) if _num(d.get('RT (min)')) is not None
            else float('inf')))
        base = dict(grp[0])
        base['peaks'] = len(grp)
        _rts = [_num(r.get('RT (min)')) for r in grp]
        base['RT (peaks)'] = '; '.join(
            f'{v:.4f}' for v in _rts if v is not None)
        if len(grp) > 1:
            for col in sample_cols:
                _s = 0.0
                _any = False
                for r in grp:
                    v = _num(r.get(col))
                    if v is not None:
                        _s += v
                        _any = True
                if _any:
                    base[col] = _s
            _labs = [str(r.get('confidence', '') or '') for r in grp]
            base['confidence'] = max(
                _labs, key=lambda s: _RANK.get(s, 5))
            _sg = [_num(r.get('sigma')) for r in grp]
            _sg = [v for v in _sg if v is not None]
            base['sigma'] = (min(_sg) if _sg else float('nan'))
            base['decided_by'] = ';'.join(dict.fromkeys(
                str(r.get('decided_by', '') or '') for r in grp
                if r.get('decided_by')))
            base['competitors'] = ';'.join(dict.fromkeys(
                _c for r in grp
                for _c in str(r.get('competitors', '') or '').split(';')
                if _c))
        out.append(base)
    return out


def _coherence_series_preds(M: dict, C, U) -> list:
    """モデルの予測 RT を系列ごとに返す。

    ふつうのクラスは 1 要素。二峰性クラス(`series_delta` を持つモデル)は
    [早い系列, 遅い系列] の 2 要素。
    """
    base = (float(M['alpha']) * float(C) + float(M['beta']) * float(U)
            + float(M['gamma']))
    try:
        d = float(M.get('series_delta') or 0.0)
    except (TypeError, ValueError):
        d = 0.0
    if d > 0.0:
        return [base, base + d]
    return [base]


def _coherence_best_pred(M: dict, C, U, obs_rt):
    """観測 RT に最も近い系列を選び (pred, abs_res, series) を返す。

    series は二峰性クラスなら 0(早い)か 1(遅い)、単一系列のクラスなら
    None。判定 B はこの series ごとに「1 compound = 1 ピーク」を課すので、
    DG は 2 本とも残り、ほかのクラスは従来どおり 1 本だけになる。
    """
    preds = _coherence_series_preds(M, C, U)
    try:
        _o = float(obs_rt)
    except (TypeError, ValueError):
        return preds[0], float('nan'), (0 if len(preds) > 1 else None)
    k = min(range(len(preds)), key=lambda j: abs(_o - preds[j]))
    return preds[k], abs(_o - preds[k]), (k if len(preds) > 1 else None)


def _fit_two_series(c_arr, u_arr, y_arr, mad_k, rmse_1,
                    min_points=None, min_delta=None, max_delta=None,
                    gain=None, max_iter=12):
    """2 系列モデル RT = α·C + β·U + γ + Δ·g を当てる(g は系列 0/1)。

    Δ をダミー列 1 本として設計行列に足すだけなので、ふつうの線形回帰で
    済む。系列の所属 g は未知なので初期分割を作り、「最も近い系列に
    付け替えて再 fit」を収束するまで繰り返す(k-means と同じ形)。

    初期分割は**並べた残差の最も大きく空いている所**で切る。中央値で切る
    のは 50 対 50 を仮定しているのと同じで、片方が少ないと見つけられない
    (検証データの neg の DG は 4 対 13 で、中央値分割では無意味な局所解に
    落ちた)。空きの大きい順に 3 つ試し、中央値分割も候補に残して、全点
    RMS が最小になった結果を採る。

    成立条件(どれか 1 つでも外れたら None を返し、1 系列のまま使う):
      - 各系列に `min_points` 以上の学習点がある
      - min_delta <= Δ <= max_delta
      - 1 系列の RMS が 2 系列の RMS の `gain` 倍以上ある

    最後の条件は**全点の二乗平均 (RMS)** で見る。`_robust_multivar_fit` が
    返す rmse は内点だけの広がりなので、少数側の系列が MAD で外れ値として
    除かれ「1 系列で十分」という誤った結論になる(検証データの neg の DG は
    17 点中 4 点が早い系列で、内点 rmse は 0.0045 min しか出なかった)。
    """
    if min_points is None:
        min_points = COH_TWO_SERIES_MIN_POINTS
    if min_delta is None:
        min_delta = COH_TWO_SERIES_MIN_DELTA
    if max_delta is None:
        max_delta = COH_TWO_SERIES_MAX_DELTA
    if gain is None:
        gain = COH_TWO_SERIES_RMSE_GAIN
    min_points = int(min_points)

    c = np.asarray(c_arr, dtype=float).ravel()
    u = np.asarray(u_arr, dtype=float).ravel()
    y = np.asarray(y_arr, dtype=float).ravel()
    n = int(y.size)
    if n != c.size or n != u.size or n < 2 * min_points:
        return None

    # 1 系列 fit。残差は初期分割に、全点 RMS は採否の基準に使う。
    try:
        coef1, _, _ = _robust_multivar_fit(
            np.column_stack([c, u, np.ones(n)]), y, mad_k=mad_k)
        r1 = y - (float(coef1[0]) * c + float(coef1[1]) * u
                  + float(coef1[2]))
    except Exception:
        return None
    rms_1 = float(np.sqrt(np.mean(r1 ** 2)))
    if not (np.isfinite(rms_1) and rms_1 > 0.0):
        return None

    def _em(g0):
        """初期分割 g0 から k-means を回す。(a, b, γ, Δ, 内点 rmse) か None。

        α, β は 2 系列共通とし、γ と Δ は各系列の残差の中央値から取る。
        系列ダミー列を設計行列に入れて頑健回帰すると、少数側の系列が
        まるごと外れ値として落ちてダミー列が定数列と同じになり、γ と
        Δ が縮退して Δ が意味を失う(検証データの neg の DG で、早い
        系列 4 点が遅い系列 13 点より内部ばらつきが大きいために起きた)。
        実測では 2 系列の α と β はほぼ同じなので共通にして困らない。
        """
        g = np.asarray(g0, dtype=float).copy()
        a, b = float(coef1[0]), float(coef1[1])
        gam = delta = float('nan')
        rmse2 = float('nan')
        for _ in range(int(max_iter)):
            n_late = float(g.sum())
            if n_late < min_points or (n - n_late) < min_points:
                return None
            # 系列ごとの中心を引いた y から α, β を取り直す(数回で収束)
            for _k in range(3):
                _r = y - a * c - b * u
                _ge = float(np.median(_r[g < 0.5]))
                _gl = float(np.median(_r[g >= 0.5]))
                try:
                    _cf, _, _rm = _robust_multivar_fit(
                        np.column_stack([c, u, np.ones(n)]),
                        y - np.where(g >= 0.5, _gl, _ge), mad_k=mad_k)
                except Exception:
                    return None
                if not (np.isfinite(_cf[0]) and np.isfinite(_cf[1])):
                    return None
                a, b = float(_cf[0]), float(_cf[1])
                rmse2 = float(_rm)
            _r = y - a * c - b * u
            gam = float(np.median(_r[g < 0.5]))
            _gl = float(np.median(_r[g >= 0.5]))
            delta = _gl - gam
            if not np.isfinite(delta):
                return None
            if delta < 0.0:
                # γ は常に早い系列にする。符号が逆なら系列を入れ替える。
                gam, delta = _gl, -delta
                g = 1.0 - g
            base = a * c + b * u + gam
            g_new = (np.abs(y - (base + delta))
                     < np.abs(y - base)).astype(float)
            if np.array_equal(g_new, g):
                break
            g = g_new
        if not (np.isfinite(rmse2) and rmse2 > 0.0):
            return None
        return a, b, gam, delta, float(rmse2)

    # 初期分割の候補: 並べた残差の空きが大きい所から 3 つ + 中央値分割
    order = np.argsort(r1)
    rs = r1[order]
    gaps = [(float(rs[k] - rs[k - 1]), k)
            for k in range(min_points, n - min_points + 1)]
    if not gaps:
        return None
    gaps.sort(reverse=True)
    inits = []
    for _gv, k in gaps[:3]:
        _g0 = np.zeros(n, dtype=float)
        _g0[order[k:]] = 1.0      # 残差が大きい側 = 遅い系列
        inits.append(_g0)
    inits.append((r1 >= float(np.median(r1))).astype(float))

    best = None   # (rms_2, a, b, gam, delta, rmse2, n_early, n_late)
    for _g0 in inits:
        _out = _em(_g0)
        if _out is None:
            continue
        a, b, gam, delta, rmse2 = _out
        if not (np.isfinite(delta) and min_delta <= delta <= max_delta):
            continue
        base = a * c + b * u + gam
        g_fin = (np.abs(y - (base + delta))
                 < np.abs(y - base)).astype(float)
        n_late = int(g_fin.sum())
        n_early = n - n_late
        if n_early < min_points or n_late < min_points:
            continue
        _r2 = np.where(g_fin > 0.5, y - (base + delta), y - base)
        rms_2 = float(np.sqrt(np.mean(_r2 ** 2)))
        if not (np.isfinite(rms_2) and rms_2 > 0.0):
            continue
        if best is None or rms_2 < best[0]:
            best = (rms_2, a, b, gam, delta, rmse2, n_early, n_late)

    if best is None:
        return None
    rms_2, a, b, gam, delta, rmse2, n_early, n_late = best
    if rms_1 < float(gain) * rms_2:
        return None      # 2 系列にしても当てはまりが良くならない
    return {
        'alpha': a, 'beta': b, 'gamma': gam,
        'series_delta': float(delta), 'rmse': float(rmse2),
        'n_series_early': int(n_early), 'n_series_late': int(n_late),
        'rms_1series': rms_1, 'rms_2series': rms_2,
        'rmse_1series': float(rmse_1),
    }


# df 補正係数の上限。df=1 だと t=12.706 で 6.48 倍になるのでここで止める。
COH_SE_DF_FACTOR_MAX = 6.0
# t(0.975, df)。df が小さいときの予測区間は正規分布より広い。
# scipy に依存させたくないので表引き + 線形補間で済ませる。
_T975 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
    6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
    12: 2.179, 15: 2.131, 20: 2.086, 25: 2.060, 30: 2.042,
}


def _t975(df: float) -> float:
    """t(0.975, df)。df >= 30 は正規分布の 1.96 に寄せる。"""
    try:
        d = float(df)
    except (TypeError, ValueError):
        return 12.706
    if not np.isfinite(d) or d < 1.0:
        return 12.706
    if d >= 30.0:
        # 30 → 2.042、∞ → 1.960 をなだらかにつなぐ
        return max(1.960, 2.042 - (2.042 - 1.960) * min(1.0, (d - 30.0) / 70.0))
    ks = sorted(_T975)
    if d in _T975:
        return _T975[int(d)]
    lo = max(k for k in ks if k <= d)
    hi = min(k for k in ks if k >= d)
    if lo == hi:
        return _T975[lo]
    w = (d - lo) / (hi - lo)
    return _T975[lo] * (1.0 - w) + _T975[hi] * w


def _coherence_df_factor(M: dict) -> tuple[float, float]:
    """df 補正係数と df を返す。

    rmse の自由度が小さいモデルの se を、t 分布ぶん広げる。
      per_class          … df = n − 3 (α, β, γ を自分の点で決めている)
      borrowed / pooled  … df = n − 1 (傾きは借り物、切片だけ自分の点)
    """
    try:
        n = float(M.get('n', 0) or 0)
    except (TypeError, ValueError):
        n = 0.0
    mode = str(M.get('fit_mode', '') or '')
    p = 3.0 if mode == 'per_class' else 1.0
    df = max(n - p, 1.0)
    f = _t975(df) / 1.960
    f = float(min(max(f, 1.0), COH_SE_DF_FACTOR_MAX))
    return f, float(df)


def _coherence_coverage(c_arr, u_arr, rmse_pool) -> dict:
    """モデルの学習点が (C, DB) 空間のどこをどれだけ覆っているか。

    予測の外挿度(leverage)を出すために使う。同じ rmse でも、
    学習点の重心の近くを予測するのと範囲の外へ外挿するのとでは
    予測の確からしさが違う。
    """
    c = np.asarray(c_arr, dtype=float).ravel()
    u = np.asarray(u_arr, dtype=float).ravel()
    out = {
        'c_mean': float(np.mean(c)) if c.size else 0.0,
        'c_var': float(np.var(c)) if c.size > 1 else 0.0,
        'u_mean': float(np.mean(u)) if u.size else 0.0,
        'u_var': float(np.var(u)) if u.size > 1 else 0.0,
        'c_min': (float(np.min(c)) if c.size else float('nan')),
        'c_max': (float(np.max(c)) if c.size else float('nan')),
        'u_min': (float(np.min(u)) if u.size else float('nan')),
        'u_max': (float(np.max(u)) if u.size else float('nan')),
    }
    try:
        out['rmse_pool'] = float(rmse_pool)
    except (TypeError, ValueError):
        out['rmse_pool'] = float('nan')
    return out


def _coherence_pred_se(M: dict, C, U):
    """(C, DB) における予測 RT の標準誤差を返す。

    戻り値: (se, h, rmse_eff)

      se(x0)   = rmse_eff * sqrt(1 + h(x0))
      h(x0)    = 1/n + ((C-C̄)²/var_C + (DB-D̄B)²/var_DB) / (n-1)
      rmse_eff = w * rmse_cls + (1-w) * rmse_pool,  w = n/(n+k)

    h は線形回帰の hat 値で、学習点の重心から離れるほど大きくなる。
    rmse_eff は点数の少ないクラスの rmse を link の pooled rmse 側へ
    縮める。n=4 / パラメータ 3 の rmse は自由度 1 しかなく、scatter の
    推定としてほとんど意味がないため(検証データでは PE[pos] ester が
    rmse=0.00145 min になり、σ が 20 まで膨らんでいた)。

    カバレッジ情報を持たないモデル(旧版の session など)では h = 1/n に
    なり、fix53 までとほぼ同じ挙動になる。
    """
    rmse = float(M.get('rmse', 0.0) or 0.0)
    try:
        n = float(M.get('n', 0) or 0)
    except (TypeError, ValueError):
        n = 0.0
    try:
        rp = float(M.get('rmse_pool', float('nan')))
    except (TypeError, ValueError):
        rp = float('nan')

    if np.isfinite(rp) and rp > 0.0 and n > 0.0:
        w = n / (n + COH_SE_SHRINK_N)
        rmse_eff = w * rmse + (1.0 - w) * rp
    else:
        rmse_eff = rmse

    if n >= 2.0:
        cv = max(float(M.get('c_var', 0.0) or 0.0), COH_SE_C_VAR_FLOOR)
        uv = max(float(M.get('u_var', 0.0) or 0.0), COH_SE_U_VAR_FLOOR)
        try:
            dc = float(C) - float(M.get('c_mean', C))
            du = float(U) - float(M.get('u_mean', U))
        except (TypeError, ValueError):
            dc = du = 0.0
        h = 1.0 / n + (dc * dc / cv + du * du / uv) / max(n - 1.0, 1.0)
    else:
        h = COH_SE_LEVERAGE_MAX
    if not np.isfinite(h) or h < 0.0:
        h = COH_SE_LEVERAGE_MAX
    h = min(h, COH_SE_LEVERAGE_MAX)

    # 自由度の補正。df が小さいモデルの rmse は scatter の推定として
    # 幅が広いので、t 分布ぶん se を広げる。点数が十分あるモデルには
    # ほとんど効かない(df=12 で 1.11 倍、df=3 で 1.62 倍)。
    _f, _df = _coherence_df_factor(M)
    se = float(rmse_eff) * float(np.sqrt(1.0 + h)) * _f
    return float(se), float(h), float(rmse_eff)


def _measure_link_deltas(models: dict, matched: pd.DataFrame) -> dict:
    """link 間の RT 差を実測する。

    同じクラスで 2 つの link の per_class モデルが両方できている
    クラスについて、そのクラスの平均 (C, U) における予測 RT の差を取り、
    中央値を使う。検証データでは ether は ester より 0.10〜0.12 min 早く、
    PC / PE / LPC でよく揃っている。

    戻り値: {(link_from, link_to): delta_min}
            delta は link_from のモデルの γ に足すと link_to になる量。
            `COH_LINK_DELTA_MIN_CLASSES` 未満のクラス数でしか測れない
            組は返さない(当てずっぽうの借用をしないため)。
    """
    by_cls: dict[str, dict[str, dict]] = {}
    for (cls, link_k), M in models.items():
        if M.get('fit_mode') != 'per_class':
            continue
        by_cls.setdefault(str(cls), {})[str(link_k)] = M
    acc: dict[tuple[str, str], list[float]] = {}
    for cls, d in by_cls.items():
        if len(d) < 2:
            continue
        g = matched[matched['lipid_class'].astype(str) == cls]
        if g.empty:
            continue
        try:
            c_bar = float(np.mean(g['C'].values.astype(float)))
            u_bar = float(np.mean(g['U'].values.astype(float)))
        except Exception:
            continue
        for l_from, Mf in d.items():
            for l_to, Mt in d.items():
                if l_from == l_to:
                    continue
                pf = (float(Mf['alpha']) * c_bar + float(Mf['beta']) * u_bar
                      + float(Mf['gamma']))
                pt = (float(Mt['alpha']) * c_bar + float(Mt['beta']) * u_bar
                      + float(Mt['gamma']))
                acc.setdefault((l_from, l_to), []).append(float(pt - pf))
    out: dict[tuple[str, str], float] = {}
    for k, v in acc.items():
        if len(v) >= COH_LINK_DELTA_MIN_CLASSES:
            out[k] = float(np.median(v))
    if out:
        log.info("[Coherence Filter] link deltas: " + ", ".join(
            f"{a}→{b} {d:+.4f} min (n={len(acc[(a, b)])} classes)"
            for (a, b), d in sorted(out.items())))
    elif acc:
        log.info("[Coherence Filter] link deltas: 測れたクラスが "
              f"{max(len(v) for v in acc.values())} 個だけなので使わない "
              f"(最低 {COH_LINK_DELTA_MIN_CLASSES} クラス必要)")
    return out


def _fit_coherence_model(
    match_df: pd.DataFrame,
    conflict_coords: set | None = None,
    rt_outlier_coords: dict[str, set] | None = None,
    min_points: int = 5,
    mad_k: float = 3.0,
    donor_rmse_max: float = 0.02,
    max_model_rmse: float = 0.1,
    min_model_rmse: float = 0.001,
    two_series_classes=None,
) -> dict[tuple[str, str], dict]:
    """link × クラス別の pooled fit を構築する。

    link (ester/ether/plasmalogen) ごとに、その link に属する全クラスの
    マッチ点を使って以下の pooled 線形モデルを fit:

        RT = α · C + β · U + γ_class

    α, β は link 内で全クラス共通、γ_class はクラス別インターセプト。

    残差の広がりが `max_model_rmse` を超えるクラスはモデルを作らない
    (⑥ はそのクラスを採点できず、その衝突は未決になる)。rmse は σ の分母
    なので、成立しないモデルを残すと σ≈0 のまま判定が適用されてしまう。
    逆に縮退した fit で rmse=0 になると σ が発散するため `min_model_rmse`
    を下限として置く。

    クラス単位で点数が足りない場合、link 全体の pooled α/β を
    当てるのではなく、まず「同じクラスの別 link」、次に「中央 RT が最も
    近いクラス」の per_class fit から傾き α, β を借りる(切片 γ は自分の
    点から取り直す)。SFC では溶出時間帯で傾きの符号が変わるため、link
    全体の平均では符号すら合わないことがあるため。

    fit 学習データの除外ルール:
      - conflict_coords に含まれる点(衝突点) → 除外
      - rt_outlier_coords の対応クラスセットに含まれる点(RT 外れ値) → 除外
      - IS 点(is_IS=True)は除外しない(anchor として使う)

    クラス × link の組み合わせで `min_points` 未満のクラスは fit から
    除外する(intercept が不安定になるため)。

    Parameters
    ----------
    match_df : DataFrame
        マッチ結果(matched 列で true の行のみ対象)
    conflict_coords : set of (rt_r6, mz_r6) | None
        除外する衝突座標
    rt_outlier_coords : dict[class, set] | None
        クラスごとの除外 RT 外れ値座標
    min_points : int
        クラス × link 単位での最小学習点数
    mad_k : float
        反復除去の MAD 倍率

    `two_series_classes` に挙げたクラスは 2 系列(γ と γ+Δ)を
    試す。DG は sn 位置異性体が部分分離するため 1 つの和組成が 2 本の
    ピークになり、1 本のモデルでは rmse が 10〜12 倍悪化する。Δ は
    データから測る(`_fit_two_series`)。採れたモデルには
    `series_delta` が入り、以降の予測は近い系列を選ぶ。

    戻り値:
      {(class, link): {'alpha','beta','gamma','rmse','n'}}
      二峰性が成立したクラスは追加で
      {'series_delta','n_series_early','n_series_late','rmse_1series'}

    学習データ抽出を matched=True かつ
    final_status='kept' で絞る。REJ_WINDOW された
    (非 IS-bearing クラスの行が IS-reserved peak に展開された)行は
    `matched=True` のまま `final_status='rejected'` だが、conflict_coords
    から除外されるため従来の filter を素通りしていた。これらが
    ester-link 共通の α, β を歪めると、PC/PE 等の予測 RT がずれて
    Coherence Filter での帰属に失敗する原因となっていた。
    """
    _two_series = set(
        two_series_classes if two_series_classes is not None
        else COH_TWO_SERIES_CLASSES_DEFAULT)
    # 二峰性クラスの 2 系列 fit 用に「判定 B で落ちただけの行」を
    # 拾っておく。判定 B の結果をそのまま 2 系列 fit の入力にすると、前回の
    # 判定 B が片方の系列を落としたクラスでは学習点が片側に寄り、Δ が
    # 測れなくなる(検証データの neg の DG。17 点すべて遅い側で rmse=0.0045)。
    # ③④⑤・手動で落ちている行は入れない。
    _ts_src = None
    if _two_series and not match_df.empty:
        try:
            _m = ((match_df.get('matched', False) == True)
                  & match_df['lipid_class'].astype(str).isin(_two_series)
                  & (match_df.get('coherence_status', STATUS_NA)
                     == STATUS_REJ_RESIDUAL))
            for _col, _bad in (('is_filter_status', STATUS_REJ_WINDOW),
                               ('rt_outlier_status', STATUS_REJ_OUTLIER),
                               ('adduct_filter_status', STATUS_REJ_ADDUCT),
                               ('manual_status', STATUS_REJ_MANUAL)):
                _m &= (match_df.get(_col, STATUS_NA) != _bad)
            _ts_src = match_df[_m].copy()
        except Exception as _e:
            log.warning(f"[Coherence Filter] two-series pool failed: {_e}")
            _ts_src = None

    matched = match_df[match_df['matched']].copy()
    # final_status='kept' に絞ることで上流フィルタ
    # (IS Filter の REJ_WINDOW、Adduct Ion Filter の REJ_ADDUCT 等)で
    # 排除された行が学習データに混入するのを防ぐ。
    if 'final_status' in matched.columns:
        matched = matched[matched['final_status'] == STATUS_KEPT].copy()
    if matched.empty:
        return {}

    # (C, U, link) を付与
    cu = matched['compound'].apply(
        lambda s: pd.Series(_parse_c_and_u(s), index=['C', 'U']))
    matched[['C', 'U']] = cu
    matched['link'] = matched['compound'].apply(_link_kind)
    matched = matched.dropna(subset=['C', 'U']).copy()
    if matched.empty:
        return {}
    matched['C'] = matched['C'].astype(int)
    matched['U'] = matched['U'].astype(int)

    # 学習から除外される点も含めた「存在するクラス×link」。
    # 学習点が 0 でも ⑥ が採点できるよう、あとで link 間の借用で埋める。
    all_pairs = set(zip(matched['lipid_class'].astype(str),
                        matched['link'].astype(str)))

    # 衝突点の除外
    if conflict_coords:
        coords = list(zip(
            matched['obs_rt'].round(6), matched['obs_mz'].round(6)))
        keep = [c not in conflict_coords for c in coords]
        matched = matched[keep].copy()

    # RT 外れ値の除外
    if rt_outlier_coords:
        def _is_rt_outlier(row):
            s = rt_outlier_coords.get(row['lipid_class'], set())
            return (round(row['obs_rt'], 6), round(row['obs_mz'], 6)) in s
        mask = ~matched.apply(_is_rt_outlier, axis=1)
        matched = matched[mask].copy()

    if matched.empty:
        return {}

    # 2 系列 fit 用の点集合。kept 点に「判定 B で落ちただけの点」を
    # 足し、通常 fit と同じ除外(衝突点・RT 外れ値)をかける。
    _ts_pool: dict = {}
    if _two_series:
        _parts = [matched[matched['lipid_class'].astype(str)
                          .isin(_two_series)]]
        if _ts_src is not None and not _ts_src.empty:
            _e2 = _ts_src.copy()
            _cu2 = _e2['compound'].apply(
                lambda s: pd.Series(_parse_c_and_u(s), index=['C', 'U']))
            _e2[['C', 'U']] = _cu2
            _e2['link'] = _e2['compound'].apply(_link_kind)
            _e2 = _e2.dropna(subset=['C', 'U'])
            if not _e2.empty:
                _e2['C'] = _e2['C'].astype(int)
                _e2['U'] = _e2['U'].astype(int)
                if conflict_coords:
                    _c2 = list(zip(_e2['obs_rt'].round(6),
                                   _e2['obs_mz'].round(6)))
                    _e2 = _e2[[c not in conflict_coords for c in _c2]]
                if rt_outlier_coords and len(_e2):
                    _e2 = _e2[~_e2.apply(
                        lambda r: (round(r['obs_rt'], 6),
                                   round(r['obs_mz'], 6))
                        in rt_outlier_coords.get(r['lipid_class'], set()),
                        axis=1)]
                if len(_e2):
                    _parts.append(_e2)
        try:
            _all_ts = pd.concat(_parts, ignore_index=True, sort=False)
        except Exception:
            _all_ts = _parts[0]
        for (_tc, _tl), _gg in _all_ts.groupby(['lipid_class', 'link']):
            _ts_pool[(str(_tc), str(_tl))] = _gg

    # per-class slope 実装
    # クラスごとに独立した α, β, γ を fit する。データ不足クラス
    # (n < min_per_class または C 分散ゼロ等)は pooled α/β + per-class γ
    # にフォールバック。
    MIN_PER_CLASS = 5  # per-class fit に最低必要な点数

    models: dict[tuple[str, str], dict] = {}
    # 傾きを借りる処理は全 link の per_class fit が出そろってから
    # 行うため、フォールバック対象をいったん貯める
    pending: list = []
    for link_k in ['ester', 'ether', 'plasmalogen']:
        d = matched[matched['link'] == link_k]
        if len(d) < min_points:
            # link 全体で点数不足 → この link はスキップ
            continue
        class_counts = d['lipid_class'].value_counts()
        valid_classes = class_counts[class_counts >= 1].index.tolist()
        if len(valid_classes) == 0:
            continue
        # 自由度: n >= 2(α, β) + len(classes) (dummies)
        if len(d) < 2 + len(valid_classes):
            continue
        d_fit = d[d['lipid_class'].isin(valid_classes)].copy()

        # ── Step A: pooled fit(fallback 用 α, β を確保)─────────────
        class_dummies = pd.get_dummies(d_fit['lipid_class'])
        class_names = list(class_dummies.columns)
        X_pool = np.column_stack([
            d_fit['C'].values.astype(float),
            d_fit['U'].values.astype(float),
            class_dummies.values.astype(float),
        ])
        y_pool = d_fit['obs_rt'].values.astype(float)
        try:
            coef_pool, _, rmse_pool = _robust_multivar_fit(
                X_pool, y_pool, mad_k=mad_k)
            pooled_alpha = float(coef_pool[0])
            pooled_beta = float(coef_pool[1])
        except Exception:
            pooled_alpha = 0.0
            pooled_beta = 0.0
            rmse_pool = float('nan')

        # ── Step B: クラスごとに独立 fit を試みる ─────────────────
        for cls in class_names:
            cls_data = d_fit[d_fit['lipid_class'] == cls]
            n_cls = len(cls_data)
            c_arr = cls_data['C'].values.astype(float)
            u_arr = cls_data['U'].values.astype(float)
            y_cls = cls_data['obs_rt'].values.astype(float)

            # per-class fit の前提: 十分な点数 + C に分散がある
            c_var = float(np.var(c_arr)) if n_cls > 1 else 0.0
            can_fit_per_class = (n_cls >= MIN_PER_CLASS and c_var > 0)

            if can_fit_per_class:
                # 設計行列: [C, U, 1] (intercept = γ_cls)
                X_cls = np.column_stack([c_arr, u_arr, np.ones(n_cls)])
                try:
                    coef_cls, _, rmse_cls = _robust_multivar_fit(
                        X_cls, y_cls, mad_k=mad_k)
                    alpha_cls = float(coef_cls[0])
                    beta_cls = float(coef_cls[1])
                    gamma_cls = float(coef_cls[2])
                    # 異常値防御: |α| or |β| が極端(>0.5 min/単位)
                    # または rmse が pooled rmse の 5 倍超過(noisy)
                    # の場合は fit 不安定とみなし fallback
                    if abs(alpha_cls) > 0.5 or abs(beta_cls) > 0.5:
                        raise ValueError(
                            f"unstable slope α={alpha_cls:+.4f} "
                            f"β={beta_cls:+.4f}")
                    if not np.isfinite(rmse_cls):
                        raise ValueError("rmse not finite")
                    if (np.isfinite(rmse_pool) and rmse_pool > 0
                            and rmse_cls > 5.0 * rmse_pool):
                        raise ValueError(
                            f"rmse too large: {rmse_cls:.4f} > 5x pooled "
                            f"{rmse_pool:.4f}")
                    _M = {
                        'alpha': alpha_cls, 'beta': beta_cls,
                        'gamma': gamma_cls, 'rmse': float(rmse_cls),
                        'n': n_cls, 'fit_mode': 'per_class',
                        'rt_med': float(np.median(y_cls)),
                        # 外挿度を出すための学習点の分布
                        **_coherence_coverage(c_arr, u_arr, rmse_pool),
                    }
                    models[(cls, link_k)] = _M
                    continue
                except Exception as exc:
                    # fall through to pooled fallback
                    log.warning(f"[Coherence Filter] per-class fit failed for "
                          f"({cls}, {link_k}): {exc} → fallback to pooled")

            # ── Fallback ────────────────────────────────────────
            # ここで pooled α/β を当てるのをやめ、全 link の
            # per_class fit が出そろってから傾きを借りる(下の第 2 パス)。
            pending.append((cls, link_k, c_arr, u_arr, y_cls, n_cls,
                            pooled_alpha, pooled_beta, rmse_pool))

    # ── fix50 第 2 パス: 傾きを近いクラスから借りる ──────────────
    # SFC では溶出時間帯で傾きの符号が変わる。検証セッション 260828 で
    # fit の良いモデルだけを見ると、クラスの中央 RT と α の相関は
    # pos −0.742 / neg −0.946 で、早く出る中性脂質(TG/DG/MG/FA/Cer)は
    # α≥0、遅く出るリン脂質(PC/PE/PG/PI/PA/SM/LPC)は α≈−0.005〜−0.015。
    # link 全体の pooled α は両者を平均した値になり、符号すら合わない
    # ことがある(pos ester の pooled は +0.0009、PC の実測は −0.0085)。
    # そこで「同じクラスの別 link」を優先しつつ、中央 RT の近いクラスを
    # 候補に並べ、**自分の点への当てはまり(残差の広がり)が最小**になる
    # 相手から傾き α, β を借りる。切片 γ は自分の点から取り直す。
    # 借りても pooled より当てはまらない場合は pooled に戻す。
    donors = [(k, M) for k, M in models.items()
              if M.get('fit_mode') == 'per_class'
              and float(M.get('rmse', 1.0)) <= donor_rmse_max
              and np.isfinite(M.get('rt_med', np.nan))]
    def _fit_with(_a, _b, _c, _u, _y):
        """借りた傾きで γ と残差の広がりを出す。点が少なければ広がりは nan。"""
        _r = _y - _a * _c - _b * _u
        _g = float(np.median(_r)) if len(_r) else 0.0
        if len(_r) < 3:
            return _g, float('nan')
        _s = float(np.median(np.abs(_r - _g))) * 1.4826
        return _g, (_s if np.isfinite(_s) else float('nan'))

    for (cls, link_k, c_arr, u_arr, y_cls, n_cls,
         p_alpha, p_beta, p_rmse) in pending:
        rt_med = float(np.median(y_cls)) if len(y_cls) else float('nan')
        # 同じクラスの別 link を先頭に、あとは中央 RT の近い順
        same = [(k, M) for k, M in donors
                if k[0] == cls and k[1] != link_k]
        other = [(k, M) for k, M in donors
                 if not (k[0] == cls and k[1] != link_k)]
        if np.isfinite(rt_med):
            other.sort(key=lambda kv: abs(float(kv[1]['rt_med']) - rt_med))
        ordered = same + other[:6]

        # 3 点以上あるときは「自分の点への当てはまり」で借用先を選ぶ。
        # RT の近さだけで選ぶと、学習点の分布が kept と違うクラスで
        # 見当外れの傾きを借りてしまう(ST が PI から借りて rmse=1.4 に
        # なった)。2 点以下では広がりを測れないので RT の近さで選ぶ。
        best = None
        for (dc, dl), M in ordered:
            _a, _b = float(M['alpha']), float(M['beta'])
            _g, _s = _fit_with(_a, _b, c_arr, u_arr, y_cls)
            if np.isfinite(_s):
                cost = _s * (0.8 if dc == cls else 1.0)
            else:
                cost = (abs(float(M['rt_med']) - rt_med)
                        if np.isfinite(rt_med) else 1e9)
                cost *= (0.5 if dc == cls else 1.0)
            if best is None or cost < best[0]:
                best = (cost, (dc, dl), _a, _b, _g, _s, float(M['rmse']))

        # pooled との比較。借用が pooled より当てはまらなければ pooled に戻す。
        _pg, _ps = _fit_with(float(p_alpha), float(p_beta),
                             c_arr, u_arr, y_cls)
        use_pooled = best is None
        if (best is not None and np.isfinite(best[5]) and np.isfinite(_ps)
                and _ps < best[5]):
            use_pooled = True
        if use_pooled:
            alpha_v, beta_v = float(p_alpha), float(p_beta)
            gamma_v, spread = _pg, _ps
            rmse_v = float(p_rmse)
            mode_s = 'pooled_fallback'
        else:
            _, (dc, dl), alpha_v, beta_v, gamma_v, spread, drmse = best
            rmse_v = drmse
            mode_s = f'borrowed:{dc}/{dl}'
        # 自分の点で広がりが測れたならそれを rmse にする(下限は借用元の 20%)
        if np.isfinite(spread) and spread > 0:
            rmse_v = max(spread, rmse_v * 0.2)
        models[(cls, link_k)] = {
            'alpha': alpha_v, 'beta': beta_v, 'gamma': gamma_v,
            'rmse': float(rmse_v), 'n': n_cls, 'fit_mode': mode_s,
            'rt_med': rt_med,
            # 外挿度を出すための学習点の分布。傾きを借りていても
            # 「自分の点がどこにあるか」は自分のもので測る。借りた傾きを
            # 学習範囲の外へ延ばすのがまさに危ない場合なので。
            **_coherence_coverage(c_arr, u_arr, p_rmse),
        }

    # ── 二峰性クラスを 2 系列(γ, γ+Δ)に差し替える ──────────
    # 同じ点集合で 1 系列の当てはまりも測り、それを比較の基準にする
    # (per_class fit の rmse は点集合が違うので基準にできない)。
    for (_tc, _tl), _pool in sorted(_ts_pool.items()):
        _M0 = models.get((_tc, _tl))
        if _M0 is None:
            continue
        _c3 = _pool['C'].values.astype(float)
        _u3 = _pool['U'].values.astype(float)
        _y3 = _pool['obs_rt'].values.astype(float)
        try:
            _, _, _r1 = _robust_multivar_fit(
                np.column_stack([_c3, _u3, np.ones(len(_y3))]),
                _y3, mad_k=mad_k)
            _r1 = float(_r1)
        except Exception:
            _r1 = float(_M0.get('rmse', float('nan')))
        _TS = _fit_two_series(_c3, _u3, _y3, mad_k, _r1)
        if not _TS:
            log.info(f"[Coherence Filter] two series not accepted for "
                  f"{_tc}/{_tl} (n={len(_y3)}, 1-series inlier rmse="
                  f"{_r1:.4f})")
            continue
        _rp0 = _M0.get('rmse_pool')
        _M0.update(_TS)
        _M0['fit_mode'] = str(_M0.get('fit_mode', 'per_class')) + ':2series'
        _M0['n'] = int(len(_y3))
        _M0.update(_coherence_coverage(_c3, _u3, _rp0))
        log.info(f"[Coherence Filter] two series for {_tc}/{_tl}: "
              f"\u0394={_TS['series_delta']:+.4f} min, RMS "
              f"{_TS['rms_1series']:.4f} \u2192 {_TS['rms_2series']:.4f} "
              f"({_TS['rms_1series'] / _TS['rms_2series']:.1f}x), "
              f"rmse {_TS['rmse']:.4f} "
              f"(early {_TS['n_series_early']} / "
              f"late {_TS['n_series_late']}, n={len(_y3)})")

    # ── link 間の借用 ──────────────────────────────────────
    # 学習点が 0 のクラス×link はここまでのパスではモデルが作られない。
    # 検証データでは LPS/ester と LPA/ester がこれで、そのクラスの ester 種が
    # 全部「まだ解けていない衝突」の中にしかないため学習点が 0 になる。
    # モデルが無いと ⑥ が採点できず、採点できないから決着せず、決着しない
    # から学習点が増えない、という堂々巡りになっていた。
    # 同じクラスの別 link のモデルに、実測した link 間 RT 差を足して埋める。
    # 借用モデルは n=0 なので _coherence_pred_se の se が最大まで広がり、
    # σ はほぼ 0(= low_confidence)になる。勝ち負けは残差で決まるので、
    # 借用側が勝つことも負けることもある。
    _link_delta = _measure_link_deltas(models, matched)
    _borrowed: list[str] = []
    for (cls, link_k) in sorted(all_pairs):
        if (cls, link_k) in models:
            continue
        donor = None
        for other in ('ester', 'ether', 'plasmalogen'):
            if other != link_k and (cls, other) in models:
                donor = other
                break
        if donor is None:
            continue
        d = _link_delta.get((donor, link_k))
        if d is None:
            continue   # 差が測れないなら作らない
        M = models[(cls, donor)]
        _new = dict(M)
        _new.update({
            'gamma': float(M['gamma']) + float(d),
            'n': 0,
            'fit_mode': f'link_borrowed:{donor}',
            'link_delta': float(d),
        })
        try:
            _new['rt_med'] = float(M.get('rt_med', float('nan'))) + float(d)
        except (TypeError, ValueError):
            pass
        models[(cls, link_k)] = _new
        _borrowed.append(f'{cls}/{link_k}←{donor}({d:+.3f})')
    if _borrowed:
        log.info("[Coherence Filter] link-borrowed models: "
              + ", ".join(_borrowed))

    # ── 成立していないモデルを落とし、rmse に下限を置く ──────
    # rmse は σ(= 残差差 / 平均 rmse)の分母になる。
    #   ・広がりが大きすぎるクラス(ST[pos] は借用先を変えても 1.4 min)を
    #     残すと σ≈0 になるが、⑥ は σ が小さくても判定を適用するため、
    #     意味のない予測で候補を棄却してしまう。モデルごと作らない。
    #   ・逆に縮退した fit で rmse=0 になると σ が発散する
    #     (検証データで PI/ester が rmse=0.00000、σ=318 が出ていた)。
    _dropped = []
    for k in list(models.keys()):
        _r = float(models[k].get('rmse', float('nan')))
        if (not np.isfinite(_r)) or _r > max_model_rmse:
            _dropped.append((k, _r))
            del models[k]
            continue
        if _r < min_model_rmse:
            models[k]['rmse'] = float(min_model_rmse)
            models[k]['rmse_floored'] = True
    if _dropped:
        log.info("[Coherence Filter] model dropped (residual spread too large): "
              + ", ".join(f"{c}/{l} rmse={r:.3f}" for (c, l), r in _dropped))
    return models





def calc_coherence_outliers(
    match_df: pd.DataFrame,
    coherence_models: dict[tuple[str, str], dict],
    mad_k: float = 3.0,
) -> dict[str, set]:
    """Coherence 外れ値判定: 各マッチ点について、自クラス × link の
    pooled fit からの予測残差が大きい点を外れ値としてフラグする。

    単変数の旧 `calc_trend_outliers` の多変数化版。U を回帰の covariate に
    含むため、SM の U 多軌跡問題(U=0 の正常スポットが誤検出される問題)
    が解消される。

    判定基準: 予測残差の絶対値が fit 時の `rmse` × `mad_k` を超える点。

    戻り値: {lipid_class: set of (rt_r6, mz_r6)}
            (既存の `_coherence_outliers` と同じデータ構造)

    判定対象を final_status='kept' の行に限定する。
    REJ_WINDOW 等で既に排除済みの行を再度 outlier と評価しない。
    """
    matched = match_df[match_df['matched']].copy()
    # final_status='kept' のみ判定対象
    if 'final_status' in matched.columns:
        matched = matched[matched['final_status'] == STATUS_KEPT].copy()
    if matched.empty:
        return {}

    cu = matched['compound'].apply(
        lambda s: pd.Series(_parse_c_and_u(s), index=['C', 'U']))
    matched[['C', 'U']] = cu
    matched['link'] = matched['compound'].apply(_link_kind)
    matched = matched.dropna(subset=['C', 'U']).copy()
    if matched.empty:
        return {}
    matched['C'] = matched['C'].astype(int)
    matched['U'] = matched['U'].astype(int)

    outliers: dict[str, set] = {
        cls: set() for cls in matched['lipid_class'].unique()
    }

    for (cls, link_k), M in coherence_models.items():
        if M['rmse'] == 0.0:
            continue
        g = matched[(matched['lipid_class'] == cls) & (matched['link'] == link_k)]
        if g.empty:
            continue
        pred = (M['alpha'] * g['C'].values
                + M['beta'] * g['U'].values
                + M['gamma'])
        residuals = g['obs_rt'].values - pred
        # 二峰性クラスは近い系列からの残差で見る。こうしないと
        # DG の 2 本のうち片方が必ず外れ値として落ちる。
        try:
            _d_ser = float(M.get('series_delta') or 0.0)
        except (TypeError, ValueError):
            _d_ser = 0.0
        if _d_ser > 0.0:
            _r2 = g['obs_rt'].values - (pred + _d_ser)
            residuals = np.where(
                np.abs(_r2) < np.abs(residuals), _r2, residuals)
        # しきい値を点ごとの予測標準誤差にする。mad_k × rmse では
        # 学習範囲の外にある点(外挿)を不当に棄却してしまう。
        _se = np.array([
            _coherence_pred_se(M, _c, _u)[0]
            for _c, _u in zip(g['C'].values, g['U'].values)
        ], dtype=float)
        _se = np.where(np.isfinite(_se) & (_se > 0.0), _se, M['rmse'])
        threshold = mad_k * _se
        mask = np.abs(residuals) > threshold
        out_rows = g[mask]
        outliers[cls].update(
            zip(out_rows['obs_rt'].round(6),
                out_rows['obs_mz'].round(6)))
    return outliers


def resolve_conflicts_by_coherence(
    conflict_coords: set,
    conflict_map: dict[tuple, list[dict]],
    coherence_models: dict[tuple[str, str], dict],
    enabled_pairs: set[frozenset],
    sigma_threshold: float = 3.0,
    loose: bool = COH_LOOSE_ATTRIBUTION_DEFAULT,
) -> dict[tuple, dict]:
    """衝突座標ごとに、候補クラスから Coherence モデルで予測 RT 残差を
    比較し、残差最小のクラスへ帰属する。

    全衝突点は必ず勝者クラスへ帰属される(閾値での「判別不能」は作らない)。
    `sigma_threshold` は「確定/不確定」ではなく「低信頼度フラグ」の
    境界として使う。

    Parameters
    ----------
    conflict_coords : set of (rt_r6, mz_r6)
        calc_conflicts の戻り値
    conflict_map : dict[(rt_r6, mz_r6), list of candidate dict]
        calc_conflicts の戻り値
    coherence_models : dict[(class, link), model dict]
        _fit_coherence_model の戻り値
    enabled_pairs : set of frozenset({class_A, class_B})
        ペア選択 UI で ON のペアのみ処理対象
    sigma_threshold : float
        低信頼度フラグの閾値(σ 単位)
    loose : bool
        採点できる候補が 1 つしかない衝突を、△ のまま残さずに
        その候補へ帰属する。比較していないので σ=0 /
        low_confidence=True / fit_mode に 'single-scored' を残す。
        ノンターゲット解析で △ を減らすための緩和で、既定 ON。

    戻り値:
      {(rt_r6, mz_r6): assignment dict}
      assignment dict の内容:
        'winner_class', 'winner_compound', 'winner_adduct',
        'loser_class',  'loser_compound',  'loser_adduct',
        'losers'  … winner 以外の全候補のリスト。
                    各要素は class / compound / adduct / pred / res /
                    abs_res / rmse / scored。`scored=False` は C が読めない
                    か該当 link のモデルが無かった候補。
        'pred_winner', 'pred_loser',
        'res_winner',  'res_loser',
        'confidence',  'low_confidence' (bool)

      帰属不能(モデル欠如、候補 < 2 等)の衝突座標は戻り値に含まれない。
    """
    from itertools import combinations
    assignments: dict[tuple, dict] = {}

    for coord in conflict_coords:
        candidates = conflict_map.get(coord, [])
        classes_here = {c['lipid_class'] for c in candidates}
        if len(classes_here) < 2:
            continue
        # 有効ペアチェック: この衝突に含まれる全クラス対が enabled_pairs に
        # 含まれていないと処理しない(部分的有効は扱わない)
        valid = True
        for a, b in combinations(classes_here, 2):
            if frozenset({a, b}) not in enabled_pairs:
                valid = False
                break
        if not valid:
            continue

        obs_rt = coord[0]
        scored = []
        for cand in candidates:
            cls = cand['lipid_class']
            compound = cand['compound']
            C, U = _parse_c_and_u(compound)
            if C is None:
                continue
            link = _link_kind(compound)
            if (cls, link) not in coherence_models:
                continue
            M = coherence_models[(cls, link)]
            # 二峰性クラスは観測 RT に近い系列の予測で採点する
            pred, _absres, _ser = _coherence_best_pred(M, C, U, obs_rt)
            res = obs_rt - pred
            # この (C, DB) における予測の標準誤差。学習範囲の外へ
            # 外挿している候補は se が大きくなり、σ に効く。
            _se, _h, _rmse_eff = _coherence_pred_se(M, C, U)
            _dff, _dfv = _coherence_df_factor(M)
            scored.append({
                'class':    cls,
                'compound': compound,
                'adduct':   cand.get('adduct'),
                'pred':     pred,
                'res':      res,
                'abs_res':  abs(res),
                'rmse':     M['rmse'],
                'se':       _se,
                'leverage': _h,
                'rmse_eff': _rmse_eff,
                'df':       _dfv,
                'df_factor': _dff,
                # 計算詳細表示用に保存
                'C':        int(C),
                'U':        int(U),
                'link':     link,
                'alpha':    M['alpha'],
                'beta':     M['beta'],
                'gamma':    M['gamma'],
                'fit_mode': M.get('fit_mode', 'pooled'),
                'series':   _ser,
            })

        if len(scored) < 2:
            # 採点できる候補が 1 つだけのとき。比較はできないが、
            # ノンターゲット用途では △ のまま残すより帰属して信頼度で
            # 示す方が使いやすい。採点できなかった候補は losers に入る。
            if not loose or len(scored) != 1:
                continue
            w = scored[0]
            _unscored = [
                {'class': c['lipid_class'], 'compound': c['compound'],
                 'adduct': c.get('adduct'), 'pred': None, 'res': None,
                 'abs_res': None, 'rmse': None,
                 'se': None, 'leverage': None, 'rmse_eff': None,
                 'scored': False}
                for c in candidates
                if not (c['lipid_class'] == w['class']
                        and c['compound'] == w['compound'])
            ]
            _l0 = _unscored[0] if _unscored else {}
            assignments[coord] = {
                'winner_class':    w['class'],
                'winner_compound': w['compound'],
                'winner_adduct':   w['adduct'],
                'loser_class':     _l0.get('class', ''),
                'loser_compound':  _l0.get('compound', ''),
                'loser_adduct':    _l0.get('adduct', ''),
                'losers':          _unscored,
                'n_candidates':    len(candidates),
                'pred_winner':     w['pred'],
                'pred_loser':      float('nan'),
                'res_winner':      w['res'],
                'res_loser':       float('nan'),
                'confidence':      0.0,
                'low_confidence':  True,
                'single_scored':   True,
                'rmse_ref':        float(w['rmse']),
                'se_winner':       float(w.get('se') or 0.0),
                'se_loser':        float('nan'),
                'se_combined':     float('nan'),
                'leverage_winner': float(w.get('leverage') or 0.0),
                'leverage_loser':  float('nan'),
                'rmse_eff_winner': float(w.get('rmse_eff') or 0.0),
                'rmse_eff_loser':  float('nan'),
                'df_winner':       float(w.get('df') or 0.0),
                'df_loser':        float('nan'),
                'winner_C':        w['C'],
                'winner_U':        w['U'],
                'winner_link':     w['link'],
                'winner_alpha':    w['alpha'],
                'winner_beta':     w['beta'],
                'winner_gamma':    w['gamma'],
                'winner_rmse':     w['rmse'],
                'winner_fit_mode': w['fit_mode'],
                'loser_C':         None, 'loser_U': None,
                'loser_link':      '', 'loser_alpha': float('nan'),
                'loser_beta':      float('nan'),
                'loser_gamma':     float('nan'),
                'loser_rmse':      float('nan'),
                'loser_fit_mode':  'not scored',
            }
            continue
        scored.sort(key=lambda x: x['abs_res'])
        winner = scored[0]
        loser = scored[1]  # 次点。σ の計算と Report の 1 行表示に使う

        # 次点以降も全部 loser として記録する。
        #   fix48 まで assignment に入るのは scored[1] だけで、
        #   _apply_coherence_status も loser_class しか見ていなかったため、
        #   3 候補以上の衝突では 3 番目以降が rejected にならず、
        #   両立しない帰属(PC 31:0 と PE 34:0)が両方とも定量表に
        #   残っていた。
        losers_all = [
            {'class': s['class'], 'compound': s['compound'],
             'adduct': s['adduct'], 'pred': s['pred'], 'res': s['res'],
             'abs_res': s['abs_res'], 'rmse': s['rmse'],
             # 外挿度もレポートに出せるように持たせる
             'se': s.get('se'), 'leverage': s.get('leverage'),
             'rmse_eff': s.get('rmse_eff'), 'scored': True}
            for s in scored[1:]
        ]
        # C が読めない / 該当 link のモデルが無い候補は scored に入らないが、
        # 「winner 以外」であることに変わりはないので棄却対象に含める。
        _scored_keys = {(s['class'], s['compound'], s['adduct'])
                        for s in scored}
        for cand in candidates:
            _k = (cand['lipid_class'], cand['compound'], cand.get('adduct'))
            if _k in _scored_keys:
                continue
            losers_all.append(
                {'class': cand['lipid_class'], 'compound': cand['compound'],
                 'adduct': cand.get('adduct'), 'pred': None, 'res': None,
                 'abs_res': None, 'rmse': None,
                 'se': None, 'leverage': None, 'rmse_eff': None,
                 'scored': False})

        # 信頼度 σ = (敗者絶対残差 - 勝者絶対残差) / 予測差の標準誤差
        #
        # 分母を「平均 rmse」から √(se_w² + se_l²) に替えた。
        #   2 つの予測を比べるのだから、不確かさは二乗和の平方根で合成する。
        #   se は学習範囲の外へ外挿している候補ほど大きくなるので、
        #   「点数の少ないクラスの小さい rmse」で σ が膨らむことがなくなる。
        #   fix53 までは PE 34:0（DB=0、学習範囲は DB 3〜7）を相手に
        #   σ=20 が出て、low_confidence の印もつかないまま通っていた。
        margin = loser['abs_res'] - winner['abs_res']
        rmse_ref = (winner['rmse'] + loser['rmse']) / 2.0   # 参考値として保持
        _se_w = float(winner.get('se') or 0.0)
        _se_l = float(loser.get('se') or 0.0)
        se_comb = float(np.sqrt(_se_w * _se_w + _se_l * _se_l))
        if se_comb <= 0.0:
            sigma = float('inf') if margin > 0 else 0.0
        else:
            sigma = margin / se_comb

        assignments[coord] = {
            'winner_class':    winner['class'],
            'winner_compound': winner['compound'],
            'winner_adduct':   winner['adduct'],
            'loser_class':     loser['class'],
            'loser_compound':  loser['compound'],
            'loser_adduct':    loser['adduct'],
            'losers':          losers_all,   # winner 以外の全候補
            'n_candidates':    len(candidates),
            'pred_winner':     winner['pred'],
            'pred_loser':      loser['pred'],
            'res_winner':      winner['res'],
            'res_loser':       loser['res'],
            'confidence':      sigma,
            'low_confidence':  sigma < sigma_threshold,
            # σ の内訳。detail / report で外挿を見せるために残す。
            'rmse_ref':        rmse_ref,
            'se_winner':       _se_w,
            'se_loser':        _se_l,
            'se_combined':     se_comb,
            'leverage_winner': float(winner.get('leverage') or 0.0),
            'leverage_loser':  float(loser.get('leverage') or 0.0),
            'rmse_eff_winner': float(winner.get('rmse_eff') or 0.0),
            'rmse_eff_loser':  float(loser.get('rmse_eff') or 0.0),
            'df_winner':       float(winner.get('df') or 0.0),
            'df_loser':        float(loser.get('df') or 0.0),
            # 計算詳細(detail panel 用)
            'winner_C':        winner['C'],
            'winner_U':        winner['U'],
            'winner_link':     winner['link'],
            'winner_alpha':    winner['alpha'],
            'winner_beta':     winner['beta'],
            'winner_gamma':    winner['gamma'],
            'winner_rmse':     winner['rmse'],
            'winner_fit_mode': winner['fit_mode'],
            'loser_C':         loser['C'],
            'loser_U':         loser['U'],
            'loser_link':      loser['link'],
            'loser_alpha':     loser['alpha'],
            'loser_beta':      loser['beta'],
            'loser_gamma':     loser['gamma'],
            'loser_rmse':      loser['rmse'],
            'loser_fit_mode':  loser['fit_mode'],
            'rmse_ref':        rmse_ref,
        }
    return assignments


# ════════════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════════════════

def detect_encoding(path: str | Path) -> str:
    det = chardet.detect(Path(path).read_bytes())
    return det.get("encoding", "utf-8") or "utf-8"


def detect_header_row(path: str | Path, encoding: str, max_scan: int = 30) -> int:
    all_kw = set(META_KEYWORDS + HEADER_MARKER_KEYWORDS)
    best_row, best_hits = 0, 0
    with open(path, encoding=encoding, errors="replace") as f:
        for idx, line in enumerate(f):
            if idx > max_scan:
                break
            cells = [c.strip().strip('"').lower() for c in line.split("\t")]
            hits = sum(1 for c in cells
                       if any(c == kw or c.startswith(kw) for kw in all_kw))
            if hits > best_hits:
                best_hits, best_row = hits, idx
    return best_row if best_hits >= 3 else 0


def load_file(path: str | Path) -> tuple[pd.DataFrame, str, int]:
    enc = detect_encoding(path)
    hrow = detect_header_row(path, enc)
    df = pd.read_csv(path, sep="\t", encoding=enc, header=hrow)
    return df, enc, hrow


def _is_meta_col(col_name: str) -> bool:
    cl = str(col_name).strip().lower()
    return any(cl == kw or cl.startswith(kw) for kw in META_KEYWORDS)


def classify_columns(df: pd.DataFrame):
    cols, n = df.columns.tolist(), len(df.columns)
    is_meta  = [_is_meta_col(c) for c in cols]
    is_fixed = [i == 0 and str(cols[0]).strip().lower().startswith("alignment")
                for i in range(n)]
    is_stat  = [False] * n
    for i in range(n - 1, -1, -1):
        cl = str(cols[i]).strip().lower()
        if any(cl.startswith(sk) or sk in cl for sk in STAT_KEYWORDS):
            is_stat[i] = True
            continue
        try:
            float(cl); is_stat[i] = True; continue
        except ValueError:
            pass
        break

    sample_idx = [i for i in range(n)
                  if not is_meta[i] and not is_stat[i] and not is_fixed[i]]
    s_start = sample_idx[0]  if sample_idx else 0
    s_end   = sample_idx[-1] if sample_idx else n - 1
    return (s_start, s_end,
            [i for i in range(n) if is_meta[i]],
            [i for i in range(n) if is_stat[i]],
            [i for i in range(n) if is_fixed[i]])


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lmap = {str(c).strip().lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lmap:
            return lmap[cand.lower()]
    return None


def make_unique_mz(mz_array: np.ndarray) -> np.ndarray:
    out = np.round(mz_array.astype(float).copy(), 6)
    seen: dict[float, int] = {}
    for i in range(len(out)):
        v = out[i]
        while v in seen:
            v = round(v + v * 1e-6, 6)
        seen[v] = i
        out[i] = v
    return out


# ════════════════════════════════════════════════════════════════════
#  Per-sample peak list support (v1.1.0)
#  ----
#  MS-DIAL の Peak list result エクスポート(per-sample TXT)を読み込み、
#  LipidZoner 内部で simple alignment + isotope grouping を実施する
#  ユーティリティ群。v1.1.0-alpha.1 では UI 結合なし、CLI / programmatic
#  使用のみ。v1.0.0 既存パスの挙動には影響しない。
# ════════════════════════════════════════════════════════════════════

# ¹³C natural mass difference. Used for isotope spacing detection.
DELTA_13C: float = 1.00336

# MS-DIAL per-sample peak list で必須となる最小カラム
_PER_SAMPLE_REQUIRED_COLS = (
    'Peak ID',
    'Precursor m/z',
    'RT (min)',
    'Height',
    'Isotope',
)


def load_per_sample_peaklist(
    file_path: str | Path,
    encoding: str = 'utf-8',
) -> pd.DataFrame:
    """
    MS-DIAL Peak list result(per-sample export)の 1 ファイル読込。

    Parameters
    ----------
    file_path : str | Path
        per-sample TXT のパス
    encoding : str
        既定 utf-8(MS-DIAL の英数字主の TXT は通常 OK)

    Returns
    -------
    pd.DataFrame
        最小カラム: Peak ID, Precursor m/z, RT (min), Height, Isotope, Adduct, S/N
        その他 MS-DIAL の全カラムも保持する
    """
    df = pd.read_csv(file_path, sep='\t', encoding=encoding)
    missing = [c for c in _PER_SAMPLE_REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Per-sample peak list missing columns: {missing}. "
            f"File: {file_path}"
        )
    # 数値カラムを float / int に揃える(欠損や型不一致対策)
    for col in ('Precursor m/z', 'RT (min)', 'Height'):
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['Isotope'] = pd.to_numeric(df['Isotope'], errors='coerce').fillna(0).astype(int)
    # NaN 行(数値変換失敗)を除去
    df = df.dropna(subset=['Precursor m/z', 'RT (min)', 'Height']).reset_index(drop=True)
    return df


def load_per_sample_peaklists(
    folder: str | Path,
    pattern: str = '*.txt',
    exclude_patterns: tuple[str, ...] = ('Height_', 'Mz_', 'Rt_', 'SN_', 'PeakID_', 'IdentificationMethod_'),
) -> list[tuple[str, pd.DataFrame]]:
    """
    フォルダ内の per-sample peak list TXT を全てロード。

    MS-DIAL の Alignment result export(Height_*, Mz_* など)を除外し、
    per-sample export ファイル名(サンプル名.txt)のみを対象とする。

    Returns
    -------
    list[tuple[str, pd.DataFrame]]
        [(sample_name, peak_df), ...] の順序付きリスト
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise ValueError(f"Not a directory: {folder}")
    candidates = sorted(folder.glob(pattern))
    out: list[tuple[str, pd.DataFrame]] = []
    for f in candidates:
        if any(f.name.startswith(p) for p in exclude_patterns):
            continue
        sample_name = f.stem
        try:
            df = load_per_sample_peaklist(f)
            out.append((sample_name, df))
        except ValueError:
            # 非 per-sample TXT は静かにスキップ
            continue
    return out


# ════════════════════════════════════════════════════════════════════
#  Alignment core  (fix18 / S2 で Alignment v0.5.27 から移植)
#  ----
#  以下は standalone ツール Alignment_v0.5.27.py のコア関数をそのまま
#  取り込んだもの。LipidZoner が独自に持っていた旧実装
#    - simple_align_per_sample  (single-linkage)
#    - group_isotopes_in_aligned(単調減少 ISO_MONO_TOL=1.2)
#  は fix18 で削除した。取り違え防止のため、旧実装は残していない。
#
#  置換の根拠:
#   * alignment  — v0.5.20 の centroid 精緻化クラスタリングは、到着順に
#     依存する旧 greedy 法と違って決定的で、許容内にある 2 ピークを
#     取りこぼさない。6 サンプル Pos で features 4771 -> 4640、
#     全サンプル検出 2850 -> 2910、境界付近ペア 204 -> 66。
#     さらに avg_mz / avg_rt は強度加重平均になり、area_s* /
#     fwhm_warning / intra_sample_dup の QC 列も得られる。
#   * isotope    — S0 実測比較(2026-08-16、検証データ Pos/Neg 各 6 サンプル)。
#     IS 消失ゼロ(pos/HexCer は 1 -> 2 に増加)。判定が割れた feature の
#     観測同位体比 ÷ 理論同位体比は、v0523=実在/fix9=同位体 の 374 件で
#     中央値 2.02(同位体では説明できない量 = 実在種)、v0523=同位体/
#     fix9=実在 の 132 件で中央値 1.02(理論どおり = 純同位体)。
#     両方向とも v0.5.23 の判定が正しく、旧 fix9 は実在種を過剰除去して
#     いた。ライブラリ実在ヒットは pos +25 / neg +2。
#     レポート: Alignment/S0_isotope_comparison/
# ════════════════════════════════════════════════════════════════════

# ── A) クラスタリング補助 ─────────────────────────────────────────
def _mz_tol_array(mz, mz_tol, mz_ppm, mz_min_da):
    """Per-peak absolute m/z tolerance.

    If mz_ppm is None -> constant mz_tol (backward compatible).
    Else -> max(mz_ppm*1e-6*mz, floor) where floor = mz_min_da or mz_tol.
    """
    mz = np.asarray(mz, dtype=float)
    if mz_ppm is None:
        return np.full(mz.shape, float(mz_tol))
    floor = float(mz_min_da) if mz_min_da is not None else float(mz_tol)
    return np.maximum(mz_ppm * 1e-6 * mz, floor)


def _assign_nearest(mz, rt, cen_mz, cen_rt, mz_tol_arr, rt_tol):
    """Assign each peak to nearest centroid within (tol, rt_tol); -1 if none.
    `cen_mz` must be sorted ascending. `mz_tol_arr` is per-peak (len == len(mz))."""
    n = len(mz)
    assign = np.full(n, -1, dtype=int)
    for i in range(n):
        mt = mz_tol_arr[i]
        lo = np.searchsorted(cen_mz, mz[i] - mt)
        hi = np.searchsorted(cen_mz, mz[i] + mt)
        best = -1
        best_d = np.inf
        for c in range(lo, hi):
            dmz = abs(cen_mz[c] - mz[i])
            drt = abs(cen_rt[c] - rt[i])
            if dmz < mt and drt < rt_tol:
                d = (dmz / mt) ** 2 + (drt / rt_tol) ** 2
                if d < best_d:
                    best_d = d
                    best = c
        assign[i] = best
    return assign


def _renucleate_orphans(assign, mz, rt, k):
    assign = assign.copy()
    for j, idx in enumerate(np.where(assign < 0)[0]):
        assign[idx] = k + j
    return assign


def _recompute_centroids(assign, mz, rt, h):
    new_mz, new_rt = [], []
    for c in np.unique(assign):
        m = assign == c
        w = h[m]
        if w.sum() > 0:
            new_mz.append(float(np.average(mz[m], weights=w)))
            new_rt.append(float(np.average(rt[m], weights=w)))
        else:
            new_mz.append(float(mz[m].mean()))
            new_rt.append(float(rt[m].mean()))
    return np.asarray(new_mz), np.asarray(new_rt)


# ── B) サンプル間 alignment (v0.5.20 centroid 精緻化 + v0.5.22 ppm) ──
def simple_align_per_sample(peak_dfs, mz_tol=0.005, rt_tol=0.05, n_refine=3,
                            mz_ppm=None, mz_min_da=None,
                            height_col='Height'):
    rows = []
    # height_col='Area' のときは area_s* が height_s* と重複するので出さない
    has_area = (height_col != 'Area') and any('Area' in df.columns for df in peak_dfs)
    has_fwhm = any('fwhm_sec' in df.columns for df in peak_dfs)
    for i, df in enumerate(peak_dfs):
        for _, r in df.iterrows():
            rec = {
                'sample_idx': i, 'peak_id': r.get('Peak ID', -1),
                'mz': float(r['Precursor m/z']), 'rt': float(r['RT (min)']),
                'height': (float(r[height_col])
                           if pd.notna(r.get(height_col, None)) else 0.0),
                'iso_msd': int(r.get('Isotope', 0)) if pd.notna(r.get('Isotope', 0)) else 0,
                'add_msd': str(r.get('Adduct', '')),
            }
            if has_area:
                rec['area'] = float(r['Area']) if pd.notna(r.get('Area', None)) else 0.0
            if has_fwhm:
                rec['fwhm_sec'] = float(r['fwhm_sec']) if pd.notna(r.get('fwhm_sec', None)) else 0.0
            rows.append(rec)
    base_cols = ['feature_id', 'avg_mz', 'avg_rt', 'n_samples',
                 'majority_iso_msdial', 'majority_adduct_msdial',
                 'fwhm_warning', 'intra_sample_dup',
                 # ギャップフィリングの窓幅に使う平均ピーク幅
                 'mean_fwhm_sec']
    if not rows:
        cols = base_cols + [f'height_s{i}' for i in range(len(peak_dfs))]
        if has_area:
            cols += [f'area_s{i}' for i in range(len(peak_dfs))]
        cols += [f'filled_s{i}' for i in range(len(peak_dfs))]
        return pd.DataFrame(columns=cols)
    pk = pd.DataFrame(rows).sort_values('mz').reset_index(drop=True)
    if pk.empty:
        cols = base_cols + [f'height_s{i}' for i in range(len(peak_dfs))]
        if has_area:
            cols += [f'area_s{i}' for i in range(len(peak_dfs))]
        cols += [f'filled_s{i}' for i in range(len(peak_dfs))]
        return pd.DataFrame(columns=cols)

    mz = pk['mz'].to_numpy(); rt = pk['rt'].to_numpy(); h = pk['height'].to_numpy()
    n = len(pk)
    mz_tol_arr = _mz_tol_array(mz, mz_tol, mz_ppm, mz_min_da)

    # seeds (intensity-anchored)
    order = np.lexsort((rt, mz, -h))
    seed_mz, seed_rt = [], []
    for si in order:
        smz = np.asarray(seed_mz)
        if smz.size:
            near = (np.abs(smz - mz[si]) < mz_tol_arr[si]) & \
                   (np.abs(np.asarray(seed_rt) - rt[si]) < rt_tol)
            if near.any():
                continue
        seed_mz.append(float(mz[si])); seed_rt.append(float(rt[si]))
    cen_mz = np.asarray(seed_mz, float); cen_rt = np.asarray(seed_rt, float)

    assign = np.zeros(n, dtype=int)
    for _ in range(max(1, n_refine)):
        o = np.argsort(cen_mz); cen_mz = cen_mz[o]; cen_rt = cen_rt[o]
        # per-peak tol must follow the peak, not reorder -> pass mz_tol_arr as is
        assign = _assign_nearest(mz, rt, cen_mz, cen_rt, mz_tol_arr, rt_tol)
        assign = _renucleate_orphans(assign, mz, rt, len(cen_mz))
        cen_mz, cen_rt = _recompute_centroids(assign, mz, rt, h)
    o = np.argsort(cen_mz); cen_mz = cen_mz[o]; cen_rt = cen_rt[o]
    assign = _assign_nearest(mz, rt, cen_mz, cen_rt, mz_tol_arr, rt_tol)
    assign = _renucleate_orphans(assign, mz, rt, len(cen_mz))

    uniq = np.unique(assign)
    remap = {old: newi for newi, old in enumerate(uniq)}
    pk['feature_id'] = np.array([remap[a] for a in assign], dtype=int)

    n_samples = len(peak_dfs)
    out_rows = []
    for fid, grp in pk.groupby('feature_id'):
        fwhm_warn = False
        mean_fwhm = 0.0
        if has_fwhm and 'fwhm_sec' in grp.columns:
            fwhm_warn = bool((grp['fwhm_sec'] > 30.0).any())
            # 検出できた検体の平均ピーク幅。ギャップフィリングの窓幅に
            # 使う（MS-DIAL と同じ規則）。0 や欠損は平均から外す。
            _fw = pd.to_numeric(grp['fwhm_sec'], errors='coerce')
            _fw = _fw[_fw > 0]
            if len(_fw):
                mean_fwhm = float(_fw.mean())
        w = grp['height'].to_numpy(); gm = grp['mz'].to_numpy(); gr = grp['rt'].to_numpy()
        if w.sum() > 0:
            avg_mz = float(np.average(gm, weights=w)); avg_rt = float(np.average(gr, weights=w))
        else:
            avg_mz = float(gm.mean()); avg_rt = float(gr.mean())
        dup = bool((grp['sample_idx'].value_counts() > 1).any())
        row = {'feature_id': int(fid), 'avg_mz': avg_mz, 'avg_rt': avg_rt,
               'n_samples': int(grp['sample_idx'].nunique()),
               'majority_iso_msdial': int(grp['iso_msd'].mode().iloc[0]) if len(grp) else 0,
               'majority_adduct_msdial': str(grp['add_msd'].mode().iloc[0]) if len(grp) else '',
               'fwhm_warning': fwhm_warn, 'intra_sample_dup': dup,
               'mean_fwhm_sec': mean_fwhm}
        for s in range(n_samples):
            sub = grp[grp['sample_idx'] == s]
            row[f'height_s{s}'] = float(sub['height'].max()) if len(sub) else 0.0
            if has_area:
                row[f'area_s{s}'] = float(sub['area'].max()) if len(sub) else 0.0
            # 検出由来か埋めた値かの旗。ここでは全部 False。
            row[f'filled_s{s}'] = False
        out_rows.append(row)
    return pd.DataFrame(out_rows).sort_values('avg_mz').reset_index(drop=True)


# ── C) サンプル間 RT ドリフト補正 (v0.5.22) ───────────────────────
def estimate_rt_offsets(peak_dfs, mz_tol=0.01, rt_gen=0.30,
                        mz_ppm=None, mz_min_da=None, min_intensity_quantile=0.5):
    """Per-sample RT offsets that center every sample on the cross-sample
    consensus RT, estimated from peaks detected in ALL samples.

    Returns (offsets, stats). offsets[s] is subtracted from sample s RT.
    A sample's corrected RT = rt - offsets[s]; sum(offsets) ~ 0.
    """
    N = len(peak_dfs)
    if N < 2:
        return [0.0] * N, {'n_common': 0}
    parts = []
    for i, df in enumerate(peak_dfs):
        parts.append(pd.DataFrame({
            'mz': pd.to_numeric(df['Precursor m/z'], errors='coerce'),
            'rt': pd.to_numeric(df['RT (min)'], errors='coerce'),
            'h':  pd.to_numeric(df['Height'], errors='coerce'), 's': i}))
    pk = pd.concat(parts, ignore_index=True).dropna(subset=['mz', 'rt', 'h'])
    pk = pk.sort_values('mz').reset_index(drop=True)
    if pk.empty:
        return [0.0] * N, {'n_common': 0}
    mz = pk['mz'].to_numpy(); rt = pk['rt'].to_numpy()
    n = len(pk)
    mtol = _mz_tol_array(mz, mz_tol, mz_ppm, mz_min_da)
    # m/z anchor clustering (generous), then RT anchor clustering (generous)
    mzc = np.zeros(n, dtype=np.int64); cur = 0; a = mz[0]
    for i in range(1, n):
        if mz[i] - a < mtol[i]:
            mzc[i] = cur
        else:
            cur += 1; mzc[i] = cur; a = mz[i]
    order = np.lexsort((rt, mzc)); mzc_s = mzc[order]; rt_s = rt[order]
    fs = np.zeros(n, dtype=np.int64); fc = 0; i = 0
    while i < n:
        j = i
        while j < n and mzc_s[j] == mzc_s[i]:
            j += 1
        a_rt = rt_s[i]; fs[i] = fc
        for k in range(i + 1, j):
            if rt_s[k] - a_rt < rt_gen:
                fs[k] = fc
            else:
                fc += 1; fs[k] = fc; a_rt = rt_s[k]
        fc += 1; i = j
    fid = np.zeros(n, dtype=np.int64); fid[order] = fs
    pk['fid'] = fid
    # tallest peak per (feature, sample)
    idx = pk.groupby(['fid', 's'])['h'].idxmax()
    tall = pk.loc[idx]
    # common features: present in all N samples
    cov = tall.groupby('fid')['s'].nunique()
    common = cov[cov == N].index
    if len(common) == 0:
        return [0.0] * N, {'n_common': 0}
    sub = tall[tall['fid'].isin(common)]
    # optional: keep only the more intense common features (robustness)
    feat_med_h = sub.groupby('fid')['h'].median()
    thr = feat_med_h.quantile(min_intensity_quantile) if len(feat_med_h) else 0
    keep = feat_med_h[feat_med_h >= thr].index
    sub = sub[sub['fid'].isin(keep)]
    # consensus RT per feature = mean across samples
    consensus = sub.groupby('fid')['rt'].mean()
    sub = sub.merge(consensus.rename('crt'), on='fid')
    sub['dev'] = sub['rt'] - sub['crt']
    offsets = [0.0] * N
    for s in range(N):
        d = sub[sub['s'] == s]['dev']
        offsets[s] = float(np.median(d)) if len(d) else 0.0
    # re-center so mean offset is exactly 0
    m = float(np.mean(offsets))
    offsets = [o - m for o in offsets]
    stats = {'n_common': int(len(keep)),
             'max_abs_offset': float(np.max(np.abs(offsets))) if offsets else 0.0}
    return offsets, stats


def correct_rt_drift(peak_dfs, mz_tol=0.01, rt_gen=0.30,
                     mz_ppm=None, mz_min_da=None):
    """Return (corrected_dfs, offsets, stats). Each sample's 'RT (min)' is
    shifted by -offset_s so all samples share a consensus RT axis."""
    offsets, stats = estimate_rt_offsets(
        peak_dfs, mz_tol=mz_tol, rt_gen=rt_gen, mz_ppm=mz_ppm, mz_min_da=mz_min_da)
    out = []
    for i, df in enumerate(peak_dfs):
        d = df.copy()
        d['RT (min)'] = pd.to_numeric(d['RT (min)'], errors='coerce') - offsets[i]
        out.append(d)
    return out, offsets, stats


# ── D) 同位体グルーピング (v0.5.23 脂質特化 averagine) ────────────
def group_isotopes(aligned_df, n_samples, intensity_prefix="height_s",
                   do_grouping=True, kC=0.050, kO=0.017, tol=0.5,
                   mz_tol=0.006, rt_tol=0.03, max_iso=3):
    """v0.5.19: post-alignment isotope grouping by averagine theoretical ratios.

    For each feature, look for a more-intense "parent" at avg_mz - k*1.003355
    (k=1..max_iso, z=1) at the same RT. The expected isotope intensity is the
    parent intensity * theoretical ratio R_k from an averagine model:
        nC = kC*mz, nO = kO*mz
        R1 = nC*0.0107 + nO*0.00038   (13C + 17O)
        R2 = C(nC,2)*0.0107^2 + nO*0.00205   (13C2 + 18O)
        R3 = C(nC,3)*0.0107^3
    If observed <= expected*(1+tol) the feature is a PURE isotope (iso_weight=k);
    otherwise it is a real co-eluting compound (iso_weight=0) and the part above
    the predicted isotope is recorded as per-sample excess (for NNLS).

    v0.5.23 lipid-specific averagine: kC/kO default to values derived from
    the project's own curated lipid library (260427_LipidQuant_lipidlibrary_
    curated.xlsm, 795 lipids across 11 classes, 343-1241 Da). Mean C/mass
    = 0.049 (~ generic 0.050, kept 0.050) and mean O/mass = 0.017 (generic
    averagine used 0.011; lipids carry ~50% more oxygen per Da from phosphate
    and ester groups). The change lifts the predicted M+2 ratio ~8-12% (via
    18O) while M+1 is essentially unchanged (13C-dominated), so a few borderline
    M+2 features are now correctly called pure isotopes.

    Validated 2026-06-04 on real plasma data: separates pure-isotope
    mis-annotations from real species (e.g. PI 38:3 kept, PC 36:0 = M+2 of
    PC 36:1 flagged). Adds columns iso_weight, iso_parent_fid, and (per sample)
    excess_<intensity col>. Does not drop rows; collapse is done at export.
    """
    df = aligned_df.copy().sort_values("avg_mz").reset_index(drop=True)
    if "feature_id" not in df.columns:
        df["feature_id"] = np.arange(len(df))
    icols = [f"{intensity_prefix}{i}" for i in range(n_samples)
             if f"{intensity_prefix}{i}" in df.columns]
    n = len(df)
    df["iso_weight"] = 0
    df["iso_parent_fid"] = -1
    if (not do_grouping) or n == 0 or not icols:
        return df
    I = df[icols].to_numpy(dtype=float)
    # aggregate feature intensity = median of nonzero per-sample values
    _Inan = np.where(I > 0, I, np.nan)
    with np.errstate(all="ignore"):
        agg = np.nanmedian(_Inan, axis=1)
    agg = np.where(np.isfinite(agg), agg, 0.0)
    mz = df["avg_mz"].to_numpy(); rt = df["avg_rt"].to_numpy()
    fid = df["feature_id"].to_numpy()
    C13 = 1.003355

    def Rk(m, k):
        nC = kC * m; nO = kO * m
        if k == 1:
            return nC * 0.0107 + nO * 0.00038
        if k == 2:
            return (nC * (nC - 1) / 2.0) * 0.0107 ** 2 + nO * 0.00205
        if k == 3:
            return (nC * (nC - 1) * (nC - 2) / 6.0) * 0.0107 ** 3
        return 0.0

    weight = np.zeros(n, int); parent_fid = np.full(n, -1, dtype=np.int64)
    excess = np.zeros((n, len(icols)))
    # mz is already ascending (sorted above) -> binary search on mz itself
    for i in range(n):
        for k in range(1, max_iso + 1):
            tgt = mz[i] - k * C13
            lo = np.searchsorted(mz, tgt - mz_tol)
            hi = np.searchsorted(mz, tgt + mz_tol)
            if hi <= lo:
                continue
            cand = np.arange(lo, hi)
            cand = cand[(np.abs(rt[cand] - rt[i]) < rt_tol) & (mz[cand] < mz[i] - 0.5)]
            if len(cand) == 0:
                continue
            j = cand[np.argmax(agg[cand])]
            if agg[j] <= 0:
                continue
            # per-sample ratio where BOTH parent and feature are detected
            both = (I[i] > 0) & (I[j] > 0)
            if not both.any():
                continue
            med_ratio = float(np.median(I[i][both] / I[j][both]))
            Rkj = Rk(mz[j], k)
            if med_ratio <= Rkj * (1 + tol):
                weight[i] = k; parent_fid[i] = fid[j]
                break
            else:
                excess[i] = np.maximum(excess[i], I[i] - I[j] * Rkj)
                parent_fid[i] = fid[j]
                break
    df["iso_weight"] = weight
    df["iso_parent_fid"] = parent_fid
    for sidx, c in enumerate(icols):
        df[f"excess_{c}"] = np.clip(excess[:, sidx], 0, None)
    return df


# ── E) iso_weight -> iso_position アダプタ (fix18 新規) ────────────
def derive_iso_labels(aligned_df: pd.DataFrame) -> pd.DataFrame:
    """group_isotopes の出力を PerSampleFileEntry が要求する形に変換する。

    group_isotopes は
      iso_weight     : 0 = 実在種、k = 直前の親に対する M+k
      iso_parent_fid : 親 feature の feature_id(-1 = なし)
    を返す。PerSampleFileEntry は
      iso_position   : 0 = M+0(monoisotopic)、k = M+k
      iso_group_id   : 同一 isotope envelope に付く ID
    を要求するので、親の連鎖を root まで辿って変換する。

    重要な点が 2 つある。

    1. iso_weight は「直前の親に対する相対位置」であって絶対位置ではない。
       group_isotopes は k=1,2,3 の順に親を探して最初に見つかった時点で
       break するため、M+2 は(M+1 が存在すれば)M+1 を親として
       iso_weight=1 を持つ。したがって iso_position は連鎖に沿った
       iso_weight の累積和になる。

    2. iso_weight=0 だが iso_parent_fid != -1 の行が存在する。これは
       「同位体と共溶出しているが、同位体だけでは説明できない量を持つ
       実在種」(超過分は excess_* 列に記録済み)であり、親の envelope に
       吸収してはいけない。よって連鎖を辿るのは iso_weight > 0 の間だけ
       とし、この行は自分自身を root とする。

    壊れた連鎖(親が見つからない)や循環参照はその時点で打ち切る。
    """
    df = aligned_df.copy()
    n = len(df)
    if n == 0:
        df['iso_position'] = pd.Series(dtype=int)
        df['iso_group_id'] = pd.Series(dtype=int)
        return df
    if 'iso_weight' not in df.columns or 'iso_parent_fid' not in df.columns:
        raise ValueError(
            "derive_iso_labels: aligned_df must come from group_isotopes "
            "(iso_weight / iso_parent_fid columns are required)")

    fid = df['feature_id'].to_numpy()
    weight = df['iso_weight'].to_numpy()
    parent = df['iso_parent_fid'].to_numpy()
    row_of_fid = {int(f): i for i, f in enumerate(fid)}

    iso_position = np.zeros(n, dtype=int)
    root_of = np.arange(n, dtype=int)
    n_broken = 0
    for i in range(n):
        cur = i
        acc = 0
        seen = {i}
        while weight[cur] > 0:
            j = row_of_fid.get(int(parent[cur]))
            if j is None or j in seen:
                n_broken += 1
                break
            acc += int(weight[cur])
            cur = j
            seen.add(j)
        iso_position[i] = acc
        root_of[i] = cur

    gid_of_root: dict[int, int] = {}
    iso_group_id = np.empty(n, dtype=int)
    for i in range(n):
        r = int(root_of[i])
        if r not in gid_of_root:
            gid_of_root[r] = len(gid_of_root)
        iso_group_id[i] = gid_of_root[r]

    if n_broken:
        log.info(f"[derive_iso_labels] rows whose parent chain broke: {n_broken} "
              f"(truncated there; each is treated as its own root)")

    df['iso_position'] = iso_position
    df['iso_group_id'] = iso_group_id
    return df


# ── F) alignment パラメータの単一の出所 (fix18 新規) ───────────────
#  fix19 で「Peak Detection & Alignment」タブから編集できるようにする。
#  それまでは create_per_sample_file_entry / 読込経路がこの既定値を使う。
#  値の根拠は project memory の Alignment 一貫性方針(2026-05-31)および
#  S0 実測比較(2026-08-16)。
# ── ギャップフィリング ────────────────────────────────────
# MS-DIAL の "Gap filling by compulsion" に相当。既定 ON も本家に合わせる。
#
# 窓幅は「その feature を検出できた検体の平均ピーク幅」× 倍率。MS-DIAL の
# 規則そのまま。FILL_RT_HALF_FALLBACK は fwhm が取れない入力（MS-DIAL TXT
# など）のための固定値で、検証データの実測 FWHM 2.2 秒 ≒ 0.037 分に由来する。
#
# align_rt_tol（既定 0.05 分）を流用してはいけない。DG の sn 位置異性体は
# RT 差 0.052 分で並ぶため、±0.05 分の窓では隣のピークを拾う（実測 24 組中
# 7 組が汚染。±0.035 分なら 0 組）。
GAP_FILL_ENABLED_DEFAULT = True
GAP_FILL_WIDTH_FACTOR = 1.0      # 平均 FWHM に掛ける倍率（片側窓 = FWHM×倍率）
GAP_FILL_RT_HALF_FALLBACK = 0.035   # 分。fwhm 不明時の片側窓
GAP_FILL_MZ_PPM = 10.0           # アノテーションと同じ許容幅
GAP_FILL_SANITY_RATIO = 5.0      # 検出検体の中央値の何倍を超えたら要注意にするか


DEFAULT_ALIGN_PARAMS: dict = {
    # --- サンプル間 alignment ---
    'align_mz_tol':   0.005,   # Da。mz_ppm を使う場合は下限として働く
    'align_rt_tol':   0.05,    # min
    'align_mz_ppm':   None,    # None = 絶対値のみ。数値を入れると ppm 併用
    'align_n_refine': 3,       # centroid 精緻化の反復回数
    'height_col':     'Height',  # 'Height' または 'Area'
    'rt_drift_correct': False,   # サンプル間 RT ドリフト補正(既定 off)
    # --- 同位体グルーピング ---
    'iso_mz_tol':     0.006,   # Da
    'iso_rt_tol':     0.03,    # min
    'iso_max_iso':    3,       # M+3 まで探索
    'iso_ratio_tol':  0.5,     # 観測 <= 理論 x (1 + これ) なら純同位体
    'iso_kC':         0.050,   # averagine 炭素係数(脂質ライブラリ由来)
    'iso_kO':         0.017,   # averagine 酸素係数(脂質ライブラリ由来)
    'iso_enabled':    True,
    # --- ギャップフィリング ---
    #  MS-DIAL の "Gap filling by compulsion" と同じで既定 ON。
    #  raw mzML 経路でのみ働く(TXT 経路には生データが無い)。
    'gap_fill':              GAP_FILL_ENABLED_DEFAULT,
    'gap_fill_width_factor': GAP_FILL_WIDTH_FACTOR,
    'gap_fill_mz_ppm':       GAP_FILL_MZ_PPM,
}


# ── F-2) raw mzML 検出パラメータの単一の出所 (fix21 新規) ────────
#  Step 2 の spinbox 既定値と session 保存の両方がこの dict を見る。
#  値の根拠は Alignment v0.5 の pyOpenMS 一貫性チューニング
#  (noise=200 / snr=3)と min_trace_length 検出バグ修正(2026-06-03)。
DEFAULT_DETECT_PARAMS: dict = {
    'mass_error_ppm':      20.0,   # MassTraceDetection の m/z 許容(ppm)
    'noise_threshold_int': 200.0,  # 強度の下限
    'chrom_peak_snr':      3.0,    # クロマトピークの S/N 下限
    'chrom_fwhm':          5.0,    # 想定ピーク幅(秒)
    'min_trace_length':    3.0,    # これより短い trace を捨てる(秒)
    'centroid':            'auto',  # 'auto' | 'on' | 'off'
}


def align_per_sample_folder(
    folder: str | Path,
    params: dict | None = None,
    sample_pattern: str = '*.txt',
) -> tuple[pd.DataFrame, list[str], dict]:
    """per-sample フォルダを読み、alignment + 同位体グルーピングまで通す。

    create_per_sample_file_entry の中身をここに切り出した。
    fix19 の「Peak Detection & Alignment」タブからも同じ関数を呼ぶ。

    Returns
    -------
    (aligned_df, sample_names, info)
        aligned_df は iso_position / iso_group_id 付き。
        info は実行時の統計(RT オフセット等)。
    """
    p = dict(DEFAULT_ALIGN_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if k in p})

    loaded = load_per_sample_peaklists(folder, pattern=sample_pattern)
    if not loaded:
        raise ValueError(f"No per-sample peak list TXT found in: {folder}")
    sample_names = [name for name, _ in loaded]
    peak_dfs = [df for _, df in loaded]
    info: dict = {'n_samples': len(peak_dfs),
                  'n_peaks_in': int(sum(len(d) for d in peak_dfs))}

    # RT ドリフト補正(任意)
    if p['rt_drift_correct'] and len(peak_dfs) > 1:
        peak_dfs, offsets, stats = correct_rt_drift(
            peak_dfs, mz_tol=max(p['align_mz_tol'], 0.01),
            mz_ppm=p['align_mz_ppm'], mz_min_da=p['align_mz_tol'])
        info['rt_offsets'] = offsets
        info['rt_drift_stats'] = stats

    aligned = simple_align_per_sample(
        peak_dfs,
        mz_tol=p['align_mz_tol'], rt_tol=p['align_rt_tol'],
        n_refine=p['align_n_refine'], mz_ppm=p['align_mz_ppm'],
        mz_min_da=p['align_mz_tol'], height_col=p['height_col'])
    info['n_features'] = int(len(aligned))

    aligned = group_isotopes(
        aligned, len(sample_names), intensity_prefix='height_s',
        do_grouping=bool(p['iso_enabled']),
        kC=p['iso_kC'], kO=p['iso_kO'], tol=p['iso_ratio_tol'],
        mz_tol=p['iso_mz_tol'], rt_tol=p['iso_rt_tol'],
        max_iso=p['iso_max_iso'])
    aligned = derive_iso_labels(aligned)
    info['n_monoisotopic'] = int((aligned['iso_position'] == 0).sum())
    info['params'] = p
    return aligned, sample_names, info


# ── G) アライメント結果のエクスポート (fix19 で Alignment から移植) ──
def write_msdial_compatible_txt(
    aligned_df: pd.DataFrame,
    sample_names: list[str],
    out_path: str | Path,
    ion_mode: str = 'neg',
    intensity_col: str = 'height',
    monoisotopic_only: bool = False,
) -> int:
    """Write all aligned features in MS-DIAL Alignment result-compatible format.

    No isotope filtering is done here (v0.2 design); downstream tools handle
    isotopes. LipidZoner v1.0.0-beta.13's existing 'Alignment result' input
    mode reads this format. 4 header rows + 1 column header row, then data.

    Parameters
    ----------
    intensity_col : 'height' | 'area'
        Which intensity to put in the sample columns (default 'height').
        Downstream NNLS can be configured to use either.

    Returns the number of features written.
    """
    m0 = aligned_df.copy().sort_values('avg_mz').reset_index(drop=True)
    if monoisotopic_only and 'iso_weight' in m0.columns:
        m0 = m0[m0['iso_weight'] == 0].reset_index(drop=True)  # collapse: M+0 only
    fid_to_alnid = ({int(f): k + 1 for k, f in enumerate(m0['feature_id'])}
                    if 'feature_id' in m0.columns else {})
    if m0.empty:
        with open(out_path, 'w', encoding='utf-8', newline='') as f:
            f.write("# No features\n")
        return 0
    has_area = any(c.startswith('area_s') for c in m0.columns)
    if intensity_col == 'area' and not has_area:
        intensity_col = 'height'   # fallback

    # MS-DIAL header (4 rows of metadata + 1 row of column names)
    n_samples = len(sample_names)
    leading_blanks = 35  # number of blank meta columns before sample columns

    lines = []
    # Row 1: Class
    lines.append('\t'.join([''] * leading_blanks + ['Class'] + ['1'] * n_samples + ['NA', 'NA']))
    # Row 2: File type
    lines.append('\t'.join([''] * leading_blanks + ['File type'] + ['Sample'] * n_samples + ['NA', 'NA']))
    # Row 3: Injection order
    lines.append('\t'.join([''] * leading_blanks + ['Injection order'] + [str(i + 1) for i in range(n_samples)] + ['NA', 'NA']))
    # Row 4: Batch ID
    lines.append('\t'.join([''] * leading_blanks + ['Batch ID'] + ['1'] * n_samples + ['Average', 'Stdev']))

    # Row 5: column header
    meta_cols = [
        'Alignment ID', 'Average Rt(min)', 'Average Mz', 'Metabolite name',
        'Adduct type', 'Post curation result', 'Fill %', 'MS/MS assigned',
        'Reference RT', 'Reference m/z', 'Formula', 'Ontology', 'INCHIKEY',
        'SMILES', 'Annotation tag (VS1.0)', 'RT matched', 'm/z matched',
        'MS/MS matched', 'Comment', 'Manually modified for quantification',
        'Manually modified for annotation', 'Isotope tracking parent ID',
        'Isotope tracking weight number', 'RT similarity', 'm/z similarity',
        'Simple dot product', 'Weighted dot product', 'Reverse dot product',
        'Matched peaks count', 'Matched peaks percentage', 'Total score',
        'S/N average', 'Spectrum reference file name', 'MS1 isotopic spectrum',
        'MS/MS spectrum',
    ]
    assert len(meta_cols) == leading_blanks, f"meta cols len {len(meta_cols)} != {leading_blanks}"
    lines.append('\t'.join(meta_cols + sample_names + ['Average', 'Stdev']))

    # Data rows
    intensity_prefix = 'area_s' if intensity_col == 'area' else 'height_s'
    for i, row in m0.iterrows():
        meta = [''] * leading_blanks
        meta[0] = str(i + 1)                       # Alignment ID
        meta[1] = f"{row['avg_rt']:.4f}"           # Average Rt(min)
        meta[2] = f"{row['avg_mz']:.5f}"           # Average Mz
        meta[3] = 'Unknown'                        # Metabolite name
        adduct = row.get('majority_adduct_msdial', '')
        if not adduct:
            adduct = '[M-H]-' if ion_mode == 'neg' else '[M+H]+'
        meta[4] = str(adduct)                      # Adduct type
        _w = int(row.get('iso_weight', 0) or 0)
        _pf = int(row.get('iso_parent_fid', -1) or -1)
        # 21 = Isotope tracking parent ID, 22 = Isotope tracking weight number
        meta[21] = str(fid_to_alnid.get(_pf, i + 1)) if _w > 0 else str(i + 1)
        meta[22] = str(_w)
        vals = meta
        # Per-sample intensities (height or area, user choice)
        h = np.array([row[f'{intensity_prefix}{s}'] for s in range(n_samples)], dtype=float)
        # v0.5.24: Fill % = fraction of samples where the feature was detected
        # (>0). MS-DIAL reports this 0-100 with one decimal (e.g. 5/6 -> 83.3).
        n_det = int((h > 0).sum())
        meta[6] = f"{n_det / n_samples * 100:.1f}" if n_samples else "0.0"
        for s in range(n_samples):
            vals.append(f"{h[s]:.0f}")
        # v0.5.24: Average / Stdev over DETECTED samples only. Including the
        # zeros of undetected samples deflated both; MS-DIAL's own export
        # averages detected samples, so this matches it.
        hd = h[h > 0]
        if hd.size:
            vals.append(f"{hd.mean():.0f}")
            vals.append(f"{hd.std(ddof=0):.0f}")
        else:
            vals.append("0")
            vals.append("0")
        lines.append('\t'.join(vals))

    with open(out_path, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(lines) + '\n')
    return len(m0)


def write_extended_csv(
    aligned_df: pd.DataFrame,
    sample_names: list[str],
    out_path: str | Path,
) -> int:
    """Write all isotope features (M+0..M+max_iso) as extended CSV.

    For downstream NNLS isotope deconvolution. Columns:
      feature_id, iso_group_id, iso_position, avg_mz, avg_rt, n_samples,
      majority_iso_msdial, majority_adduct_msdial,
      <sample_name_0>, <sample_name_1>, ..., <sample_name_N-1>

    Returns the number of rows written.
    """
    df = aligned_df.copy()
    # Rename height columns to actual sample names
    rename_map = {f'height_s{i}': sample_names[i] for i in range(len(sample_names))}
    df = df.rename(columns=rename_map)
    df.to_csv(out_path, index=False)
    return len(df)


# ════════════════════════════════════════════════════════════════════
#  H) Raw mzML ピーク検出 (fix20 / S4 で Alignment v0.5.27 から移植)
#  ----
#  pyOpenMS は **別プロセス**(Alignment_raw_worker.py)で呼ぶ。理由は
#  pyOpenMS が Qt5 を同梱しており、PySide6(Qt6)と同一プロセスに同居すると
#  Windows で DLL がぶつかるため(症状: "DLL load failed while importing
#  QtCore")。この制約により worker を単一 .py に取り込むことはできない。
#  worker は LipidZoner 本体と同じフォルダに置いておくこと。
# ════════════════════════════════════════════════════════════════════

ALIGN_WORKER_NAME = "Alignment_raw_worker.py"


def resolve_align_worker() -> "Path | None":
    """Alignment_raw_worker.py の場所を解決する。見つからなければ None。

    探索順: 本体と同じフォルダ → 起動スクリプトのフォルダ →
    PyInstaller の展開先(sys._MEIPASS)→ カレントディレクトリ。
    """
    cands: list[Path] = []
    try:
        cands.append(Path(__file__).resolve().parent / ALIGN_WORKER_NAME)
    except Exception:
        pass
    try:
        cands.append(Path(sys.argv[0]).resolve().parent / ALIGN_WORKER_NAME)
    except Exception:
        pass
    _mei = getattr(sys, '_MEIPASS', None)
    if _mei:
        cands.append(Path(_mei) / ALIGN_WORKER_NAME)
    try:
        cands.append(Path.cwd() / ALIGN_WORKER_NAME)
    except Exception:
        pass
    for c in cands:
        try:
            if c.is_file():
                return c
        except Exception:
            continue
    return None


# pyOpenMS 可用性は **遅延チェック**(結果はキャッシュ)。
# Alignment 本体は import 時に worker を --check で起動していたが、
# LipidZoner は per-sample TXT だけでも使えるので、起動のたびに
# サブプロセスを立ち上げて数秒待たせるのは割に合わない。raw mzML
# モードを選んだ時点で初めて確認する。
_PYOPENMS_STATE: "tuple[bool, str | None, str | None] | None" = None


def check_pyopenms(force: bool = False) -> "tuple[bool, str | None, str | None]":
    """(available, version, error) を返す。初回のみ worker を --check で起動。"""
    global _PYOPENMS_STATE
    if _PYOPENMS_STATE is not None and not force:
        return _PYOPENMS_STATE
    worker = resolve_align_worker()
    if worker is None:
        _PYOPENMS_STATE = (
            False, None,
            f"{ALIGN_WORKER_NAME} not found. "
            f"Place it in the same folder as LipidZoner.")
        return _PYOPENMS_STATE
    try:
        result = _subprocess.run(
            [sys.executable, str(worker), "--check"],
            capture_output=True, text=True, timeout=60)
        if result.returncode == 0:
            _PYOPENMS_STATE = (True, result.stdout.strip(), None)
        else:
            err = (result.stderr or result.stdout or "Unknown error").strip()
            _PYOPENMS_STATE = (False, None, err[:500])
    except Exception as e:
        _PYOPENMS_STATE = (False, None, f"{type(e).__name__}: {e}")
    return _PYOPENMS_STATE


def detect_mzml_polarity(path) -> "str | None":
    """mzML の先頭にある polarity cvParam を読んで 'pos'/'neg' を返す。

    positive=MS:1000130 / negative=MS:1000129。バイト列を先頭から走査して
    最初のタグで打ち切るので、巨大な mzML でも速く、全体を読み込まない。
    gzip 圧縮 mzML にも対応。判定できなければ None。
    """
    import gzip
    POS = b"MS:1000130"
    NEG = b"MS:1000129"
    try:
        with open(path, "rb") as fh:
            magic = fh.read(2)
        opener = gzip.open if magic == b"\x1f\x8b" else open
        with opener(path, "rb") as fh:
            tail = b""
            while True:
                chunk = fh.read(65536)
                if not chunk:
                    break
                buf = tail + chunk
                ip = buf.find(POS)
                ineg = buf.find(NEG)
                if ip != -1 or ineg != -1:
                    if ineg == -1:
                        return "pos"
                    if ip == -1:
                        return "neg"
                    return "pos" if ip < ineg else "neg"
                tail = buf[-16:]   # タグがチャンク境界で割れるのを防ぐ
    except Exception:
        return None
    return None


def load_raw_mzml_peaks(
    file_path: str | Path,
    mass_error_ppm: float = 20.0,
    noise_threshold_int: float = 1000.0,
    chrom_peak_snr: float = 3.0,
    chrom_fwhm: float = 5.0,
    min_trace_length: float = 3.0,
    centroid: str = "auto",
    progress_cb=None,
) -> 'pd.DataFrame':
    """Read .mzML and run pyOpenMS Stage 1+2 peak detection via subprocess.

    pyOpenMS is invoked via a separate Python process (Alignment_raw_worker.py)
    to avoid DLL conflicts with PySide6's Qt6 on Windows. The worker writes
    the result to a temporary CSV which is then loaded back here.

    Parameters
    ----------
    file_path : str | Path
        .mzML file path
    mass_error_ppm : float
        Mass trace clustering tolerance (ppm).
    noise_threshold_int : float
        Minimum intensity threshold.
    chrom_peak_snr : float
        Minimum chromatographic peak S/N.
    chrom_fwhm : float
        Expected chromatographic peak FWHM (sec).
    progress_cb : callable | None
        Optional progress callback(stage_name: str, percent: float)

    Returns
    -------
    pd.DataFrame with columns:
      Peak ID, Precursor m/z, RT (min), Height, Area, fwhm_sec, Isotope, Adduct
    """
    ok, _ver, err = check_pyopenms()
    if not ok:
        raise RuntimeError(
            "pyOpenMS not available. " + (err or "Install via `pip install pyopenms`."))
    import subprocess as _sub
    worker = resolve_align_worker()
    if worker is None:
        raise RuntimeError(f"Worker script not found: {ALIGN_WORKER_NAME}")

    if progress_cb:
        progress_cb("Spawning pyOpenMS worker", 5.0)

    # Temporary output CSV
    with _tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
        out_csv = tmp.name

    try:
        cmd = [
            sys.executable, str(worker),
            "--mzml", str(file_path),
            "--out",  out_csv,
            "--mass-error-ppm", f"{mass_error_ppm}",
            "--noise", f"{noise_threshold_int}",
            "--snr",   f"{chrom_peak_snr}",
            "--fwhm",  f"{chrom_fwhm}",
            "--min-trace-length", f"{min_trace_length}",
            "--centroid", f"{centroid}",
        ]
        if progress_cb:
            progress_cb("Running peak detection", 20.0)
        result = _sub.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            raise RuntimeError(
                f"Worker subprocess failed (exit {result.returncode}):\n"
                f"--- stderr ---\n{result.stderr[-2000:]}\n"
                f"--- stdout ---\n{result.stdout[-500:]}")
        if progress_cb:
            progress_cb("Loading worker output", 90.0)
        df = pd.read_csv(out_csv)
    finally:
        try:
            os.unlink(out_csv)
        except Exception:
            pass

    if progress_cb:
        progress_cb("Done", 100.0)
    return df


_DETECT_MAX_WORKERS = 4


_WORKER_TIMEOUT_SEC = 600.0


def _build_worker_cmd(worker, file_path, out_csv, mass_error_ppm,
                      noise_threshold_int, chrom_peak_snr, chrom_fwhm,
                      min_trace_length, centroid):
    """Build the Alignment_raw_worker.py command line for one mzML file."""
    return [
        sys.executable, str(worker),
        "--mzml", str(file_path),
        "--out",  str(out_csv),
        "--mass-error-ppm", f"{mass_error_ppm}",
        "--noise", f"{noise_threshold_int}",
        "--snr",   f"{chrom_peak_snr}",
        "--fwhm",  f"{chrom_fwhm}",
        "--min-trace-length", f"{min_trace_length}",
        "--centroid", f"{centroid}",
    ]


def _finalize_worker(returncode, stdout, stderr, out_csv):
    """Turn a finished worker process into {'df': DataFrame|None, 'error': str|None}."""
    if returncode != 0:
        return {'df': None, 'error': (
            f"worker exit {returncode}: "
            f"{(stderr or '')[-800:].strip() or (stdout or '')[-400:].strip()}")}
    try:
        df = pd.read_csv(out_csv)
    except Exception as e:
        return {'df': None, 'error': f"could not read worker output: {e}"}
    return {'df': df, 'error': None}


def _run_worker_pool(jobs, max_workers, should_cancel=None, on_tick=None,
                     timeout=_WORKER_TIMEOUT_SEC, poll_interval=0.05):
    """Run detection worker subprocesses with bounded concurrency.

    Parameters
    ----------
    jobs : list of dict, each {'key': hashable, 'cmd': list[str], 'out_csv': str}
    max_workers : int   maximum simultaneous subprocesses
    should_cancel : callable() -> bool   polled each tick; True aborts
    on_tick : callable(completed:int, total:int, running_keys:list)  UI pump
    timeout : float   per-job wall-clock limit (seconds)

    Returns
    -------
    (results, canceled) where results maps key -> {'df', 'error'} and canceled
    is True if aborted via should_cancel.
    """
    import subprocess as _sub
    import time as _time

    results = {}
    pending = list(jobs)
    running = {}          # key -> (proc, start_time, job)
    total = len(jobs)
    completed = 0
    canceled = False
    max_workers = max(1, int(max_workers))

    while pending or running:
        if should_cancel is not None and should_cancel():
            canceled = True
            for _key, (proc, _s, _j) in running.items():
                try:
                    proc.terminate()
                except Exception:
                    pass
            break

        # fill up to max_workers
        while pending and len(running) < max_workers:
            job = pending.pop(0)
            try:
                proc = _sub.Popen(job['cmd'], stdout=_sub.PIPE,
                                  stderr=_sub.PIPE, text=True)
            except Exception as e:
                results[job['key']] = {'df': None, 'error': f"spawn failed: {e}"}
                completed += 1
                continue
            running[job['key']] = (proc, _time.time(), job)

        # poll running processes
        done_keys = []
        for key, (proc, start, job) in running.items():
            rc = proc.poll()
            if rc is not None:
                out, err = proc.communicate()
                results[key] = _finalize_worker(rc, out, err, job['out_csv'])
                done_keys.append(key)
            elif (_time.time() - start) > timeout:
                try:
                    proc.kill()
                except Exception:
                    pass
                try:
                    out, err = proc.communicate(timeout=5)
                except Exception:
                    out, err = '', ''
                results[key] = {'df': None,
                                'error': f"timeout after {timeout:.0f}s"}
                done_keys.append(key)
        for key in done_keys:
            running.pop(key, None)
            completed += 1

        if on_tick is not None:
            on_tick(completed, total, list(running.keys()))

        if pending or running:
            _time.sleep(poll_interval)

    return results, canceled


# ── H2) ギャップフィリング ────────────────────────────────
#  MS-DIAL の "Gap filling by compulsion"(既定 ON)の移植。
#
#  アラインメント後、ある検体で 0 になっているセルを raw mzML から積分し
#  直す。窓幅は **その feature を検出できた検体の平均ピーク幅**（
#  `mean_fwhm_sec`）で、MS-DIAL と同じ規則。極大が無くても積分する。
#
#  なぜ align_rt_tol を流用しないか: DG の sn 位置異性体は RT 差 0.052 分で
#  並ぶ。±0.05 分の窓では隣のピークを拾い、実測で 24 組中 7 組が汚染した。
#  平均 FWHM 由来の ±0.035 分では 0 組。
#
#  面積について: 高さは検出値をよく再現する（検出済みセルで検算して中央比
#  0.99、600 件中 594 件が ±10% 以内）が、面積は固定窓の台形積分なので
#  pyOpenMS の computePeakArea（マストレース全体を積分）とは一致しない
#  （中央比 0.84、四分位 0.61〜0.93）。Height 基準の定量では問題にならない
#  が、Area 基準を選んでいるときは埋めた値の精度が落ちることを UI で伝える。

def build_gap_fill_targets(aligned, n_samples, *,
                           width_factor=GAP_FILL_WIDTH_FACTOR,
                           rt_half_fallback=GAP_FILL_RT_HALF_FALLBACK,
                           height_prefix='height_s'):
    """検体ごとの埋め対象を返す。

    Returns
    -------
    dict  sample_idx -> DataFrame[feature_id, mz, rt, rt_half]
        空の検体は含めない。
    """
    hcols = [f'{height_prefix}{i}' for i in range(n_samples)]
    hcols = [c for c in hcols if c in aligned.columns]
    if not hcols or aligned.empty:
        return {}
    H = aligned[hcols].to_numpy(dtype=float)
    H = np.nan_to_num(H, nan=0.0)
    miss = H <= 0

    # 窓幅: 平均 FWHM(秒) → 分。取れていない feature は固定値へ落とす。
    if 'mean_fwhm_sec' in aligned.columns:
        fw = pd.to_numeric(aligned['mean_fwhm_sec'], errors='coerce')
        fw = fw.fillna(0.0).to_numpy(dtype=float)
    else:
        fw = np.zeros(len(aligned), dtype=float)
    half = fw / 60.0 * float(width_factor)
    half = np.where(half > 0, half, float(rt_half_fallback))

    fid = aligned['feature_id'].to_numpy()
    mz = pd.to_numeric(aligned['avg_mz'], errors='coerce').to_numpy(float)
    rt = pd.to_numeric(aligned['avg_rt'], errors='coerce').to_numpy(float)
    ok = np.isfinite(mz) & np.isfinite(rt)

    out = {}
    for s in range(len(hcols)):
        # その検体で 0、かつ他のどこかで検出されている行
        sel = miss[:, s] & (~miss).any(axis=1) & ok
        if not sel.any():
            continue
        out[s] = pd.DataFrame({'feature_id': fid[sel], 'mz': mz[sel],
                               'rt': rt[sel], 'rt_half': half[sel]})
    return out


def _build_fill_worker_cmd(worker, mzml_path, targets_csv, out_csv, *,
                           mz_ppm=GAP_FILL_MZ_PPM,
                           rt_half_fallback=GAP_FILL_RT_HALF_FALLBACK):
    """worker の --fill コマンドラインを 1 ファイル分組み立てる。"""
    return [
        sys.executable, str(worker), "--fill",
        "--mzml", str(mzml_path),
        "--targets", str(targets_csv),
        "--out", str(out_csv),
        "--fill-ppm", f"{mz_ppm}",
        "--fill-rt-half", f"{rt_half_fallback}",
    ]


def gap_fill_aligned(aligned, raw_paths, *, worker=None,
                     width_factor=GAP_FILL_WIDTH_FACTOR,
                     rt_half_fallback=GAP_FILL_RT_HALF_FALLBACK,
                     mz_ppm=GAP_FILL_MZ_PPM,
                     sanity_ratio=GAP_FILL_SANITY_RATIO,
                     height_col='Height',
                     max_workers=None, should_cancel=None, on_tick=None):
    """欠測セルを raw mzML から埋める。(aligned_filled, info) を返す。

    aligned は変更しない（コピーを返す）。失敗しても元の表はそのまま返し、
    info['error'] に理由を入れる。つまり呼び出し側は戻り値をそのまま使える。

    Parameters
    ----------
    raw_paths : list[str]   検体 0..n-1 に対応する mzML のパス。
    height_col : 'Height' or 'Area'
        アラインメントで主強度として使った列。'Area' のときは埋める値も
        worker の area 列から採る。

    Notes
    -----
    worker は窓内で RT が最も近い極大を height にする。極大が 1 つも無い
    ときだけ窓内の最大に落とし、source 列に 'forced' と書く。その件数は
    info['n_forced'] に入る（検証データでは 2,251 件中 14 件）。

    面積は固定窓の台形積分なので、pyOpenMS の computePeakArea（マストレース
    全体を積分）とは一致しない（検出済みセルで検算して中央比 0.84、四分位
    0.61〜0.93）。高さは中央比 0.99（500 件中 500 件が ±25% 以内）。
    Area 基準の定量では埋めた値の精度が落ちる。
    """
    info = {'enabled': True, 'n_targets': 0, 'n_filled': 0, 'n_suspect': 0,
            'n_forced': 0,
            'per_sample': {}, 'error': None, 'canceled': False,
            'width_factor': float(width_factor),
            'rt_half_fallback': float(rt_half_fallback),
            'mz_ppm': float(mz_ppm)}
    if aligned is None or aligned.empty:
        info['error'] = 'aligned table is empty'
        return aligned, info

    raw_paths = [str(p) for p in (raw_paths or [])]
    hcols = [c for c in aligned.columns if c.startswith('height_s')]
    n_samples = len(hcols)
    if len(raw_paths) != n_samples:
        info['error'] = (f'raw mzML の数({len(raw_paths)})が検体数'
                         f'({n_samples})と合わない')
        return aligned, info
    missing_files = [p for p in raw_paths if not Path(p).is_file()]
    if missing_files:
        info['error'] = ('raw mzML が見つからない: '
                         + ', '.join(Path(p).name for p in missing_files[:3]))
        return aligned, info

    if worker is None:
        worker = resolve_align_worker()
    if worker is None:
        info['error'] = f'{ALIGN_WORKER_NAME} が見つからない'
        return aligned, info

    targets = build_gap_fill_targets(
        aligned, n_samples, width_factor=width_factor,
        rt_half_fallback=rt_half_fallback)
    if not targets:
        info['error'] = None
        return aligned.copy(), info
    info['n_targets'] = int(sum(len(t) for t in targets.values()))

    import tempfile as _tempfile
    import os as _os

    jobs, tmp_files = [], []
    for s, tdf in targets.items():
        with _tempfile.NamedTemporaryFile(suffix='_tg.csv', delete=False,
                                          mode='w', newline='',
                                          encoding='utf-8') as tf:
            tdf.to_csv(tf, index=False)
            tg_csv = tf.name
        with _tempfile.NamedTemporaryFile(suffix='_fill.csv',
                                          delete=False) as of:
            out_csv = of.name
        tmp_files += [tg_csv, out_csv]
        jobs.append({'key': s, 'out_csv': out_csv,
                     'cmd': _build_fill_worker_cmd(
                         worker, raw_paths[s], tg_csv, out_csv,
                         mz_ppm=mz_ppm,
                         rt_half_fallback=rt_half_fallback)})

    if max_workers is None:
        max_workers = min(len(jobs), _DETECT_MAX_WORKERS,
                          max(1, os.cpu_count() or _DETECT_MAX_WORKERS))
    try:
        results, canceled = _run_worker_pool(
            jobs, max_workers, should_cancel=should_cancel, on_tick=on_tick)
    finally:
        for p in tmp_files:
            try:
                _os.unlink(p)
            except Exception:
                pass
    info['canceled'] = bool(canceled)

    out = aligned.copy()
    for s in range(n_samples):
        fc = f'filled_s{s}'
        if fc not in out.columns:
            out[fc] = False
        else:
            out[fc] = out[fc].fillna(False).astype(bool)

    # 検出された検体の中央値（怪しい埋め値の判定に使う）
    H0 = out[[f'height_s{i}' for i in range(n_samples)]].to_numpy(float)
    H0 = np.nan_to_num(H0, nan=0.0)
    _pos = np.where(H0 > 0, H0, np.nan)
    with np.errstate(all='ignore'):
        med_det = np.nanmedian(_pos, axis=1)
    fid_pos = {int(f): i for i, f in enumerate(out['feature_id'].to_numpy())}

    use_area_as_height = (str(height_col) == 'Area')
    errors = []
    for s in range(n_samples):
        r = results.get(s)
        if r is None:
            continue
        if r.get('df') is None:
            errors.append(f'{Path(raw_paths[s]).name}: {r.get("error")}')
            continue
        fdf = r['df']
        if fdf.empty or 'feature_id' not in fdf.columns:
            continue
        rows = [fid_pos.get(int(f), -1) for f in fdf['feature_id'].to_numpy()]
        hv = pd.to_numeric(fdf.get('height', 0), errors='coerce')
        hv = hv.fillna(0.0).to_numpy(float)
        av = pd.to_numeric(fdf.get('area', 0), errors='coerce')
        av = av.fillna(0.0).to_numpy(float)
        main = av if use_area_as_height else hv
        # 'forced' = 窓に極大が無く、窓内の最大で埋めた行。
        #   MS-DIAL の compulsion そのものだが、根拠はいちばん弱い。
        _srcv = (fdf['source'].astype(str).to_numpy()
                 if 'source' in fdf.columns
                 else np.full(len(fdf), '', dtype=object))
        n_f = 0
        n_s = 0
        n_forced = 0
        hcol = f'height_s{s}'
        acol = f'area_s{s}'
        hnew = out[hcol].to_numpy(float).copy()
        anew = (out[acol].to_numpy(float).copy()
                if acol in out.columns else None)
        fnew = out[f'filled_s{s}'].to_numpy(bool).copy()
        for q, ri in enumerate(rows):
            if ri < 0 or main[q] <= 0:
                continue
            hnew[ri] = float(main[q])
            if anew is not None:
                anew[ri] = float(av[q])
            fnew[ri] = True
            n_f += 1
            if str(_srcv[q]) == 'forced':
                n_forced += 1
            m = med_det[ri]
            if np.isfinite(m) and m > 0 and main[q] > sanity_ratio * m:
                n_s += 1
        out[hcol] = hnew
        if anew is not None:
            out[acol] = anew
        out[f'filled_s{s}'] = fnew
        info['per_sample'][s] = {'n_targets': int(len(targets.get(s, []))),
                                 'n_filled': n_f, 'n_suspect': n_s,
                                 'n_forced': n_forced}
        info['n_filled'] += n_f
        info['n_suspect'] += n_s
        info['n_forced'] += n_forced

    # 埋めた後は検出検体数が増えるので n_samples 列を引き直す
    if 'n_samples' in out.columns:
        Hf = out[[f'height_s{i}' for i in range(n_samples)]].to_numpy(float)
        out['n_samples'] = (np.nan_to_num(Hf, nan=0.0) > 0).sum(axis=1).astype(int)

    if errors:
        info['error'] = '; '.join(errors[:3])
    return out, info


# ── I) パラメータの自動推定 (fix20 で移植) ────────────────────────
def consistency_noise_sweep(
    base_dfs: list,
    mz_tol: float = 0.010,
    rt_tol: float = 0.100,
    candidates=(50, 100, 150, 200, 300, 500, 800, 1200),
    cv_max: float = 0.20,
):
    """v0.5.10: Tune the noise threshold purely from cross-sample consistency
    (no external reference). Strategy: cluster the peaks ONCE (anchor-based,
    same logic as simple_align_per_sample), then for each candidate noise level
    simulate detection by keeping only peaks with Height >= candidate and score
    the resulting feature table.

    Objective (maximised):  score = n_robust * sqrt(full_fraction)
      - n_robust      = # features present in ALL samples with height CV <= cv_max
      - full_fraction = (# features in all samples) / (# total features)
    n_robust rewards reproducible, fully-detected features; sqrt(full_fraction)
    mildly favours purity but with diminishing returns, so a noisy polarity
    (e.g. Neg) is not pushed to a high threshold that discards real features
    for no CV improvement. (v0.5.13: was full_fraction, which over-weighted
    purity and over-tuned Neg.)

    Returns (best_noise, rows) where each row is a dict with keys:
      noise, n_total, n_full, full_frac, median_cv, n_robust, score
    """
    N = len(base_dfs)
    parts = []
    for i, df in enumerate(base_dfs):
        parts.append(pd.DataFrame({
            'mz': pd.to_numeric(df['Precursor m/z'], errors='coerce'),
            'rt': pd.to_numeric(df['RT (min)'], errors='coerce'),
            'height': pd.to_numeric(df['Height'], errors='coerce'),
            'sidx': i,
        }))
    pk = pd.concat(parts, ignore_index=True).dropna(subset=['mz', 'rt', 'height'])
    pk = pk.sort_values(['mz', 'rt']).reset_index(drop=True)
    n = len(pk)
    if n == 0:
        return None, []
    mz = pk['mz'].to_numpy(); rt = pk['rt'].to_numpy()
    # m/z anchor clustering
    mzc = np.zeros(n, dtype=np.int64); cur = 0; a_mz = mz[0]
    for i in range(1, n):
        if mz[i] - a_mz < mz_tol:
            mzc[i] = cur
        else:
            cur += 1; mzc[i] = cur; a_mz = mz[i]
    # RT anchor clustering within each m/z cluster -> feature id
    order = np.lexsort((rt, mzc))
    mzc_s = mzc[order]; rt_s = rt[order]
    feat_sorted = np.zeros(n, dtype=np.int64); fc = 0
    i = 0
    while i < n:
        j = i
        while j < n and mzc_s[j] == mzc_s[i]:
            j += 1
        a_rt = rt_s[i]; feat_sorted[i] = fc
        for k in range(i + 1, j):
            if rt_s[k] - a_rt < rt_tol:
                feat_sorted[k] = fc
            else:
                fc += 1; feat_sorted[k] = fc; a_rt = rt_s[k]
        fc += 1
        i = j
    fid = np.zeros(n, dtype=np.int64); fid[order] = feat_sorted
    sidx = pk['sidx'].to_numpy(); height = pk['height'].to_numpy()

    rows = []
    for C in candidates:
        m = height >= C
        if not m.any():
            rows.append(dict(noise=C, n_total=0, n_full=0, full_frac=0.0,
                             median_cv=float('nan'), n_robust=0, score=0.0))
            continue
        sub = pd.DataFrame({'fid': fid[m], 'sidx': sidx[m], 'h': height[m]})
        g = sub.groupby(['fid', 'sidx'])['h'].max().reset_index()
        piv = g.pivot(index='fid', columns='sidx', values='h')
        n_total = piv.shape[0]
        full = piv.dropna(axis=0, how='any')
        n_full = full.shape[0]
        if n_full and N > 1:
            H = full.to_numpy()
            mean = H.mean(axis=1); std = H.std(axis=1, ddof=1)
            cv = np.divide(std, mean, out=np.full_like(std, np.nan), where=mean > 0)
            n_robust = int(np.nansum(cv <= cv_max))
            median_cv = float(np.nanmedian(cv))
        else:
            n_robust = 0; median_cv = float('nan')
        full_frac = (n_full / n_total) if n_total else 0.0
        score = n_robust * (full_frac ** 0.5)   # v0.5.13: sqrt-damped purity
        rows.append(dict(noise=C, n_total=n_total, n_full=n_full,
                         full_frac=full_frac, median_cv=median_cv,
                         n_robust=n_robust, score=score))
    valid = [r for r in rows if r['n_total'] > 0]
    best = max(valid, key=lambda r: r['score'])['noise'] if valid else None
    return best, rows


def estimate_alignment_tolerances(base_dfs, mz_gen=0.02, rt_gen=0.2):
    """v0.5.14: Suggest alignment m/z & RT tolerances from cross-sample scatter.

    Cluster peaks generously, keep full-coverage features (one peak per sample,
    max height), drop gross mis-merges (spread > 0.5 * generous window), and use
    the 99th percentile of the per-feature m/z and RT spread (max-min across
    samples). This is robust to nearby-compound mis-merges while still covering
    real run-to-run drift. Returns (mz_tol, rt_tol, stats) or (None, None, {}).
    """
    N = len(base_dfs)
    parts = []
    for i, df in enumerate(base_dfs):
        parts.append(pd.DataFrame({
            'mz': pd.to_numeric(df['Precursor m/z'], errors='coerce'),
            'rt': pd.to_numeric(df['RT (min)'], errors='coerce'),
            'h':  pd.to_numeric(df['Height'], errors='coerce'),
            's':  i}))
    pk = pd.concat(parts, ignore_index=True).dropna(subset=['mz', 'rt', 'h'])
    pk = pk.sort_values(['mz', 'rt']).reset_index(drop=True)
    n = len(pk)
    if n == 0:
        return None, None, {}
    mz = pk['mz'].to_numpy(); rt = pk['rt'].to_numpy()
    mzc = np.zeros(n, dtype=np.int64); cur = 0; a_mz = mz[0]
    for i in range(1, n):
        if mz[i] - a_mz < mz_gen:
            mzc[i] = cur
        else:
            cur += 1; mzc[i] = cur; a_mz = mz[i]
    order = np.lexsort((rt, mzc)); mzc_s = mzc[order]; rt_s = rt[order]
    fs = np.zeros(n, dtype=np.int64); fc = 0; i = 0
    while i < n:
        j = i
        while j < n and mzc_s[j] == mzc_s[i]:
            j += 1
        a_rt = rt_s[i]; fs[i] = fc
        for k in range(i + 1, j):
            if rt_s[k] - a_rt < rt_gen:
                fs[k] = fc
            else:
                fc += 1; fs[k] = fc; a_rt = rt_s[k]
        fc += 1
        i = j
    fid = np.zeros(n, dtype=np.int64); fid[order] = fs; pk['fid'] = fid
    nss = pk.groupby('fid')['s'].nunique()
    full_ids = nss[nss == N].index
    mzsp = []; rtsp = []
    for fid_ in full_ids:
        sub = pk[pk['fid'] == fid_].sort_values('h').groupby('s').tail(1)
        mzv = sub['mz'].to_numpy(); rtv = sub['rt'].to_numpy()
        ms = float(mzv.max() - mzv.min()); rs = float(rtv.max() - rtv.min())
        if ms > 0.5 * mz_gen or rs > 0.5 * rt_gen:
            continue
        mzsp.append(ms); rtsp.append(rs)
    if len(mzsp) < 20:
        return None, None, {'n_full': len(full_ids), 'n_clean': len(mzsp)}
    mzsp = np.array(mzsp); rtsp = np.array(rtsp)
    mz_tol = round(float(np.percentile(mzsp, 99)), 4)
    rt_tol = round(float(np.percentile(rtsp, 99)), 3)
    stats = {
        'n_full': int(len(full_ids)), 'n_clean': int(len(mzsp)),
        'mz_med': float(np.median(mzsp)), 'mz_p95': float(np.percentile(mzsp, 95)),
        'mz_p99': float(np.percentile(mzsp, 99)),
        'rt_med': float(np.median(rtsp)), 'rt_p95': float(np.percentile(rtsp, 95)),
        'rt_p99': float(np.percentile(rtsp, 99)),
    }
    return mz_tol, rt_tol, stats




# ════════════════════════════════════════════════════════════════════
#  FileEntry
# ════════════════════════════════════════════════════════════════════

class FileEntry:
    def __init__(self, path: str | Path, tag: str, ion_mode: str = "pos"):
        self.path     = Path(path)
        self.tag      = tag
        self.ion_mode = ion_mode   # "pos" or "neg"
        self.df, self.encoding, self.header_row = load_file(path)
        (self.sample_start, self.sample_end,
         self.meta_indices, self.stat_indices,
         self.fixed_indices) = classify_columns(self.df)
        self.rt_col = find_column(
            self.df, ["Average Rt(min)", "Average Rt", "RT", "Retention time"])
        self.mz_col = find_column(
            self.df, ["Average Mz", "Average m/z", "m/z", "Precursor m/z"])

    @property
    def label(self) -> str:
        return f"{self.tag} [{self.ion_mode}] {self.path.name}"

    def sample_columns(self) -> list[str]:
        return self.df.columns[self.sample_start: self.sample_end + 1].tolist()

    def rt_array(self) -> np.ndarray | None:
        return (pd.to_numeric(self.df[self.rt_col], errors="coerce").values
                if self.rt_col else None)

    def mz_array(self) -> np.ndarray | None:
        return (pd.to_numeric(self.df[self.mz_col], errors="coerce").values
                if self.mz_col else None)

# ════════════════════════════════════════════════════════════════════
#  StubFileEntry
# ════════════════════════════════════════════════════════════════════

class StubFileEntry:
    """FileEntry-compatible placeholder with empty data.

    Used by main() to construct an empty AP that the user can populate
    via the Manage Files dialog. All array-returning methods yield
    empty numpy arrays, and sample_columns() returns []. The real
    FileEntry replaces this once data is loaded via _load_fe().
    """

    def __init__(self):
        self.path = Path("<no file>")
        self.tag = "#0"
        self.ion_mode = "pos"
        self.df = pd.DataFrame()
        self.encoding = "utf-8"
        self.header_row = 0
        self.sample_start = 0
        self.sample_end = -1
        self.meta_indices = []
        self.stat_indices = []
        self.fixed_indices = []
        self.rt_col = None
        self.mz_col = None

    @property
    def label(self) -> str:
        return "(no data loaded)"

    def sample_columns(self) -> list:
        return []

    def rt_array(self):
        import numpy as _np
        return _np.array([], dtype=float)

    def mz_array(self):
        import numpy as _np
        return _np.array([], dtype=float)


# ════════════════════════════════════════════════════════════════════
#  PerSampleFileEntry (v1.1.0)
#  ----
#  FileEntry と同じ interface を持ち、per-sample peak list を裏で保持する
#  クラス。rt_array() / mz_array() / sample_columns() は M+0 features の
#  みを返すので、既存の match_library_to_peaks 以降のパイプラインに透過的
#  に乗せられる。Isotope (M+1, M+2, M+3) は self._aligned_df 内に保持し、
#  isotope_envelope_for_feature() などのアクセサで個別に取得する。
# ════════════════════════════════════════════════════════════════════

class PerSampleFileEntry:
    """FileEntry-compatible wrapper for per-sample peak list workflow.

    Attributes
    ----------
    path : Path
        合成されたパス(`PerSample:<folder name>` 形式、表示用のみ)
    tag : str
        FE タグ(例: "#1")
    ion_mode : str
        "pos" または "neg"
    df : pd.DataFrame
        M+0 features のみを含む MS-DIAL alignment-like テーブル
        (Alignment ID, Average Rt(min), Average Mz, <sample names>)
    encoding : str
        "utf-8"(固定、互換性のため)
    header_row : int
        0(固定、互換性のため)
    sample_start, sample_end : int
        df 内のサンプル列範囲(classify_columns 互換)
    meta_indices, stat_indices, fixed_indices : list[int]
        df のカラム分類(classify_columns 互換)
    rt_col, mz_col : str
        df 上の RT, m/z カラム名

    Internal
    --------
    _aligned_df : pd.DataFrame
        全 features(M+0, M+1, M+2, M+3)を含む aligned テーブル
        cols: feature_id, avg_mz, avg_rt, n_samples, iso_group_id,
              iso_position, majority_iso_msdial, majority_adduct_msdial,
              height_s0, ..., height_s{N-1}
    _sample_names : list[str]
        サンプル名(per-sample TXT の stem)
    """

    def __init__(
        self,
        aligned_df: 'pd.DataFrame',
        sample_names: list[str],
        ion_mode: str = "pos",
        tag: str = "#0",
        source_folder: str | Path | None = None,
    ):
        if 'iso_position' not in aligned_df.columns:
            raise ValueError(
                "aligned_df must include 'iso_position' column "
                "(call group_isotopes() then derive_iso_labels() first)")
        n_samples_in_df = sum(1 for c in aligned_df.columns
                              if str(c).startswith('height_s'))
        if n_samples_in_df != len(sample_names):
            raise ValueError(
                f"sample_names length ({len(sample_names)}) != "
                f"height_s* columns count ({n_samples_in_df})")

        self._aligned_df = aligned_df.copy()
        self._sample_names = list(sample_names)
        self.ion_mode = ion_mode
        self.tag = tag
        # 実フォルダパスを保持(Session 保存/復元で per-sample データを
        # 自動再ロードするために使う)。
        self.source_folder = str(source_folder) if source_folder else None
        # サンプル順の生データ(mzML)実パス。Candidate Picker が
        # EIC を引くのに使う。① タブ経由なら実際に読んだファイル、
        # session 復元なら保存された値が入る。無ければ None。
        self.raw_paths: "list[str] | None" = None
        folder_label = Path(source_folder).name if source_folder else "memory"
        self.path = Path(f"<PerSample:{folder_label}>")

        # M+0 features のみ抽出して MS-DIAL-like df を構築
        m0 = aligned_df[aligned_df['iso_position'] == 0].copy()
        m0 = m0.sort_values('avg_mz').reset_index(drop=True)

        synth = pd.DataFrame({
            'Alignment ID'   : np.arange(len(m0), dtype=int),
            'Average Rt(min)': m0['avg_rt'].values,
            'Average Mz'     : m0['avg_mz'].values,
        })
        # M+0 内部の feature_id を保持(isotope 取得用)
        synth['_feature_id_internal'] = m0['feature_id'].values
        synth['_iso_group_id_internal'] = m0['iso_group_id'].values

        # サンプル列(height_s0, ..., を実サンプル名にリネーム)
        for i, name in enumerate(sample_names):
            synth[name] = m0[f'height_s{i}'].values

        self.df = synth
        self.encoding = "utf-8"
        self.header_row = 0

        # classify_columns 互換のカラム分類を手動構築
        cols = self.df.columns.tolist()
        meta_cols = ['Alignment ID', 'Average Rt(min)', 'Average Mz',
                     '_feature_id_internal', '_iso_group_id_internal']
        self.meta_indices = [i for i, c in enumerate(cols) if c in meta_cols]
        self.stat_indices = []
        self.fixed_indices = [0]  # Alignment ID
        sample_idx = [i for i, c in enumerate(cols) if c in sample_names]
        self.sample_start = sample_idx[0] if sample_idx else 0
        self.sample_end = sample_idx[-1] if sample_idx else -1

        self.rt_col = 'Average Rt(min)'
        self.mz_col = 'Average Mz'

    # ───── FileEntry compatible API ─────
    @property
    def label(self) -> str:
        return f"{self.tag} [{self.ion_mode}] {self.path.name}"

    def sample_columns(self) -> list[str]:
        return list(self._sample_names)

    def rt_array(self) -> 'np.ndarray':
        return pd.to_numeric(self.df[self.rt_col], errors='coerce').values

    def mz_array(self) -> 'np.ndarray':
        return pd.to_numeric(self.df[self.mz_col], errors='coerce').values

    # ───── New accessors (v1.1.0) ─────
    def all_features_df(self) -> 'pd.DataFrame':
        """全 features(M+0, M+1, M+2, M+3)を含む aligned table を返す。"""
        return self._aligned_df.copy()

    def m0_features_df(self) -> 'pd.DataFrame':
        """M+0 features のみを返す(self.df と概ね同じだが内部列付き)。"""
        return self._aligned_df[self._aligned_df['iso_position'] == 0].copy()

    def iso_group_id_array(self) -> 'np.ndarray':
        """M+0 features の iso_group_id 配列(rt_array() / mz_array() と同じ長さ)。"""
        return self.df['_iso_group_id_internal'].values

    def isotope_envelope_for_feature(
        self,
        m0_row_idx: int,
        max_iso: int = 3,
    ) -> dict:
        """
        M+0 features 内の row index(self.df の行)から、その iso_group の
        全 isotope (M+0, M+1, ..., M+max_iso) を per-sample height で返す。

        Returns
        -------
        dict
            {
              'iso_group_id': int,
              'mz_M0': float, 'rt_M0': float,
              'heights': {
                0: np.ndarray (n_samples,),  # M+0
                1: np.ndarray (n_samples,),  # M+1 (or zeros if not detected)
                2: np.ndarray (n_samples,),  # M+2
                ...
              },
              'detected_iso_positions': list[int],  # 実際に linked された positions
            }
        """
        if m0_row_idx < 0 or m0_row_idx >= len(self.df):
            raise IndexError(f"m0_row_idx out of range: {m0_row_idx}")
        gid = int(self.df.iloc[m0_row_idx]['_iso_group_id_internal'])
        grp = self._aligned_df[self._aligned_df['iso_group_id'] == gid]
        n_samples = len(self._sample_names)
        heights: dict[int, 'np.ndarray'] = {}
        detected: list[int] = []
        for k in range(0, max_iso + 1):
            row = grp[grp['iso_position'] == k]
            if len(row) == 0:
                heights[k] = np.zeros(n_samples, dtype=float)
            else:
                r = row.iloc[0]
                heights[k] = np.array(
                    [float(r[f'height_s{i}']) for i in range(n_samples)],
                    dtype=float)
                detected.append(k)
        m0_row = grp[grp['iso_position'] == 0].iloc[0] if (grp['iso_position'] == 0).any() else None
        return {
            'iso_group_id': gid,
            'mz_M0': float(m0_row['avg_mz']) if m0_row is not None else np.nan,
            'rt_M0': float(m0_row['avg_rt']) if m0_row is not None else np.nan,
            'heights': heights,
            'detected_iso_positions': detected,
        }


def resolve_entry_raw_paths(fe) -> "list[str]":
    """エントリに対応する生データ(mzML)のパスをサンプル順で返す。

    優先順:
      1. fe.raw_paths に入っていて、実在するもの
      2. fe.source_folder の中から、サンプル名 + '.mzML' を探す
         (session を別マシンで開いた / raw_paths が無い旧 session 用)

    1 件も見つからなければ空リストを返す。呼び出し側は
    「EIC は出せない」と表示すること(TXT 取り込みのデータなど、
    そもそも生データが存在しない経路がある)。
    """
    names = list(getattr(fe, '_sample_names', None)
                 or (fe.sample_columns() if hasattr(fe, 'sample_columns') else []))
    stored = list(getattr(fe, 'raw_paths', None) or [])
    if stored and len(stored) == len(names) and all(Path(p).exists() for p in stored):
        return stored

    folder = getattr(fe, 'source_folder', None)
    if not folder:
        return [p for p in stored if Path(p).exists()]
    d = Path(folder)
    if not d.is_dir():
        return [p for p in stored if Path(p).exists()]
    out: list[str] = []
    for nm in names:
        hit = None
        for ext in ('.mzML', '.mzml', '.MZML'):
            cand = d / f"{nm}{ext}"
            if cand.exists():
                hit = cand
                break
        if hit is None:
            return []          # 1 つでも欠けたら並び順が崩れるので諦める
        out.append(str(hit))
    return out


def create_per_sample_file_entry(
    folder: str | Path,
    ion_mode: str = "pos",
    tag: str = "#0",
    params: dict | None = None,
    sample_pattern: str = '*.txt',
) -> PerSampleFileEntry:
    """
    Per-sample peak list フォルダから PerSampleFileEntry を一発生成。

    alignment / 同位体グルーピングの実体は align_per_sample_folder に
    移した。個別の tol 引数は params 辞書(DEFAULT_ALIGN_PARAMS と同じキー)
    にまとめてある。params=None なら既定値。

    Parameters
    ----------
    folder : str | Path
        per-sample TXT が並んだフォルダ
    ion_mode : str
        "pos" or "neg"
    tag : str
        FileEntry タグ
    params : dict | None
        DEFAULT_ALIGN_PARAMS のキーの一部または全部。未指定キーは既定値。
    sample_pattern : str
        ファイル glob パターン

    Returns
    -------
    PerSampleFileEntry
    """
    aligned, sample_names, info = align_per_sample_folder(
        folder, params=params, sample_pattern=sample_pattern)
    fe = PerSampleFileEntry(
        aligned_df=aligned,
        sample_names=sample_names,
        ion_mode=ion_mode,
        tag=tag,
        source_folder=folder,
    )
    # fix19 の UI 表示・再アライメント用に、実行時の情報を持たせておく
    fe.align_info = info
    fe.align_params = info.get('params', dict(DEFAULT_ALIGN_PARAMS))
    return fe




# ════════════════════════════════════════════════════════════════════
#  InspectDialog
# ════════════════════════════════════════════════════════════════════

# ════════════════════════════════════════════════════════════════════
#  ① タブの可視化 (fix20 / S4 で Alignment v0.5.27 から移植)
#  ----
#  ChromatogramWidget / ChromatogramDialog : TIC・EIC ビューア(raw モード)
#  AlignmentPlotWindow                     : RT×m/z 散布図 + サンプル別強度 + EIC
#  AdvancedSettingsDialog                  : 検出/アライメントパラメータの自動推定
#
#  AdvancedSettingsDialog は既存の ParametersDialog("Advanced Setting")とは
#  対象が異なる(こちらは検出パラメータ、あちらはアノテーションのフィルタ)
#  ため、統合せず別ダイアログのまま残している。
# ════════════════════════════════════════════════════════════════════

# matplotlib は LipidZoner 本体が冒頭で無条件に import しているので、
# ここに到達している時点で利用可能。移植元の分岐をそのまま活かすため
# フラグだけ定義しておく。
MATPLOTLIB_AVAILABLE = True
MATPLOTLIB_ERROR = ""
# 移植元は NavigationToolbar2QT を NavToolbar という名前で import している。
# LipidZoner 側は NavigationToolbar という名前なので別名を張る。
NavToolbar = NavigationToolbar


class ChromatogramWidget(QWidget):
    """Embedded chromatogram viewer for Step 2.

    Shows TIC (Total Ion Chromatogram) or EIC for selected sample.
    Uses Alignment_raw_worker.py via subprocess to extract data.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._sample_paths: list[Path] = []
        self._sample_names: list[str] = []
        self._tic_cache = {}  # name -> df

        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)

        # Controls row
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Sample:"))
        self._cb_sample = QComboBox()
        self._cb_sample.setMinimumWidth(200)
        self._cb_sample.currentIndexChanged.connect(self._redraw)
        ctrl.addWidget(self._cb_sample)
        ctrl.addWidget(QLabel("View:"))
        self._cb_mode = QComboBox()
        self._cb_mode.addItems(["TIC", "EIC"])
        self._cb_mode.currentIndexChanged.connect(self._mode_changed)
        ctrl.addWidget(self._cb_mode)
        ctrl.addWidget(QLabel("m/z:"))
        self._sp_eic_mz = QDoubleSpinBox()
        self._sp_eic_mz.setDecimals(4)
        self._sp_eic_mz.setRange(50.0, 5000.0)
        self._sp_eic_mz.setValue(760.5851)
        self._sp_eic_mz.setEnabled(False)
        ctrl.addWidget(self._sp_eic_mz)
        ctrl.addWidget(QLabel("± Da:"))
        self._sp_eic_tol = QDoubleSpinBox()
        self._sp_eic_tol.setDecimals(4)
        self._sp_eic_tol.setRange(0.0005, 0.1)
        self._sp_eic_tol.setValue(0.005)
        self._sp_eic_tol.setEnabled(False)
        ctrl.addWidget(self._sp_eic_tol)
        self._btn_refresh = QPushButton("View")
        self._btn_refresh.clicked.connect(self._redraw)
        ctrl.addWidget(self._btn_refresh)
        ctrl.addStretch()
        lay.addLayout(ctrl)

        # Plot canvas
        if MATPLOTLIB_AVAILABLE:
            self._fig = Figure(figsize=(8, 2.5), tight_layout=True)
            self._canvas = FigureCanvas(self._fig)
            self._ax = self._fig.add_subplot(111)
            lay.addWidget(self._canvas)
        else:
            self._canvas = None
            lay.addWidget(QLabel(f"matplotlib unavailable: {MATPLOTLIB_ERROR}"))

    def set_samples(self, sample_names: list[str], sample_paths: list):
        """Update the sample selector."""
        self._sample_names = list(sample_names)
        self._sample_paths = [Path(p) for p in sample_paths]
        self._tic_cache.clear()
        self._cb_sample.blockSignals(True)
        self._cb_sample.clear()
        for name in self._sample_names:
            self._cb_sample.addItem(name)
        self._cb_sample.blockSignals(False)
        if self._canvas is not None and self._sample_names:
            self._redraw()

    def _mode_changed(self):
        is_eic = self._cb_mode.currentText() == "EIC"
        self._sp_eic_mz.setEnabled(is_eic)
        self._sp_eic_tol.setEnabled(is_eic)
        self._redraw()

    def _redraw(self):
        if self._canvas is None or not self._sample_names:
            return
        idx = self._cb_sample.currentIndex()
        if idx < 0 or idx >= len(self._sample_paths):
            return
        sample_name = self._sample_names[idx]
        path = self._sample_paths[idx]
        mode = self._cb_mode.currentText()
        ax = self._ax
        ax.clear()
        ax.text(0.5, 0.5, f"Loading {mode}...", ha='center', va='center',
                transform=ax.transAxes, color='gray')
        ax.set_axis_off()
        self._canvas.draw_idle()
        QApplication.processEvents()

        worker = resolve_align_worker()
        import subprocess as _sub, tempfile as _tempfile, os as _os
        try:
            with _tempfile.NamedTemporaryFile(suffix='.csv', delete=False, mode='w') as tmp:
                out_csv = tmp.name
            if mode == "TIC":
                cmd = [sys.executable, str(worker), "--tic",
                       "--mzml", str(path), "--out", out_csv]
            else:
                cmd = [sys.executable, str(worker), "--eic",
                       "--mzml", str(path), "--out", out_csv,
                       "--eic-mz", f"{self._sp_eic_mz.value():.6f}",
                       "--eic-mz-tol", f"{self._sp_eic_tol.value():.6f}"]
            r = _sub.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode != 0:
                raise RuntimeError(r.stderr[-500:])
            df = pd.read_csv(out_csv)
        except Exception as e:
            ax.clear()
            ax.text(0.5, 0.5, f"Failed: {e}", ha='center', va='center',
                    transform=ax.transAxes, color='red', fontsize=8)
            ax.set_axis_off()
            self._canvas.draw_idle()
            return
        finally:
            try: _os.unlink(out_csv)
            except Exception: pass

        ax.clear()
        ax.set_axis_on()
        x = df['rt_min'].values
        y_col = 'tic' if mode == "TIC" else 'intensity'
        y = df[y_col].values
        ax.plot(x, y, lw=1.0, color='steelblue')
        ax.fill_between(x, 0, y, alpha=0.2, color='steelblue')
        if mode == "EIC":
            ax.set_title(f"EIC m/z = {self._sp_eic_mz.value():.4f} ± {self._sp_eic_tol.value():.4f} Da   |   {sample_name}",
                         fontsize=9)
        else:
            ax.set_title(f"TIC   |   {sample_name}", fontsize=9)
        ax.set_xlabel('RT (min)')
        ax.set_ylabel('Intensity')
        ax.grid(alpha=0.3)
        self._canvas.draw_idle()


class ChromatogramDialog(QDialog):
    """Resizable dialog wrapping ChromatogramWidget (v0.4.1).

    Provides a larger, more comfortable viewing area than the embedded
    widget in Step 2. The widget itself is reused inside.
    """
    def __init__(self, sample_names, sample_paths, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Chromatogram Viewer")
        self.resize(1100, 600)
        self.setSizeGripEnabled(True)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        self._viewer = ChromatogramWidget(self)
        lay.addWidget(self._viewer)
        self._viewer.set_samples(sample_names, sample_paths)
        # Add Navigation toolbar for matplotlib (zoom/pan)
        if MATPLOTLIB_AVAILABLE and self._viewer._canvas is not None:
            from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as _NT
            tb = _NT(self._viewer._canvas, self)
            lay.addWidget(tb)
        # Close button
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_row.addWidget(btn_close)
        lay.addLayout(btn_row)


class AlignmentPlotWindow(QMainWindow):
    """Interactive visualization for aligned features.

    Three plots:
      1. Scatter (RT vs m/z, color by n_samples) — click to select feature
      2. Per-sample heights bar chart — updates on feature selection
      3. EIC (extracted ion chromatogram) — extracted via worker subprocess
         (only available when source data was Raw mode .mzML)
    """

    def __init__(self, aligned_df, sample_names, sample_paths,
                 raw_mode=False, intensity_col='height',
                 mz_tol_eic=0.005, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Alignment Visualization")
        self.resize(1200, 800)
        self._aligned_df = aligned_df.reset_index(drop=True).copy()
        self._sample_names = sample_names
        self._sample_paths = sample_paths
        self._raw_mode = raw_mode
        self._intensity_col = intensity_col  # 'height' or 'area'
        self._mz_tol_eic = mz_tol_eic
        self._selected_idx = None
        self._eic_cache = {}  # (feat_idx) -> dict[sample_idx, df]

        central = QSplitter(Qt.Vertical)
        self.setCentralWidget(central)

        # ── Top: scatter ────────────────────────────────────
        scatter_w = QWidget()
        scatter_lay = QVBoxLayout(scatter_w)
        scatter_lay.setContentsMargins(2, 2, 2, 2)
        ctrl_row = QHBoxLayout()
        ctrl_row.addWidget(QLabel("Color by:"))
        self._cb_color = QComboBox()
        self._cb_color.addItems(["n_samples", "log10(mean intensity)",
                                  "fwhm_warning"])
        self._cb_color.currentIndexChanged.connect(self._redraw_scatter)
        ctrl_row.addWidget(self._cb_color)
        ctrl_row.addStretch()
        ctrl_row.addWidget(QLabel(f"Features: {len(self._aligned_df)}"))
        scatter_lay.addLayout(ctrl_row)
        self._scatter_fig = Figure(figsize=(10, 4), tight_layout=True)
        self._scatter_canvas = FigureCanvas(self._scatter_fig)
        self._scatter_ax = self._scatter_fig.add_subplot(111)
        self._scatter_cbar = None   # v0.5.17: reuse one colorbar (avoid stacking)
        scatter_lay.addWidget(NavToolbar(self._scatter_canvas, self))
        scatter_lay.addWidget(self._scatter_canvas)
        self._scatter_canvas.mpl_connect('button_press_event',
                                         self._on_scatter_click)
        central.addWidget(scatter_w)

        # ── Bottom: bar + EIC side by side ──────────────────
        bottom_split = QSplitter(Qt.Horizontal)
        # Bar
        bar_w = QWidget()
        bar_lay = QVBoxLayout(bar_w)
        bar_lay.setContentsMargins(2, 2, 2, 2)
        self._bar_label = QLabel("Click a point in the scatter to view")
        self._bar_label.setStyleSheet("color:#555; padding:2px;")
        bar_lay.addWidget(self._bar_label)
        self._bar_fig = Figure(figsize=(5, 3), tight_layout=True)
        self._bar_canvas = FigureCanvas(self._bar_fig)
        self._bar_ax = self._bar_fig.add_subplot(111)
        bar_lay.addWidget(self._bar_canvas)
        bottom_split.addWidget(bar_w)
        # EIC
        eic_w = QWidget()
        eic_lay = QVBoxLayout(eic_w)
        eic_lay.setContentsMargins(2, 2, 2, 2)
        self._eic_label = QLabel("EIC")
        self._eic_label.setStyleSheet("color:#555; padding:2px;")
        eic_lay.addWidget(self._eic_label)
        self._eic_fig = Figure(figsize=(5, 3), tight_layout=True)
        self._eic_canvas = FigureCanvas(self._eic_fig)
        self._eic_ax = self._eic_fig.add_subplot(111)
        eic_lay.addWidget(self._eic_canvas)
        bottom_split.addWidget(eic_w)
        central.addWidget(bottom_split)

        if not self._raw_mode:
            self._eic_label.setText(
                "EIC: not available (input was TXT mode, no raw mzML)")

        self._redraw_scatter()

    # ─── Scatter ────────────────────────────────────────────
    def _redraw_scatter(self):
        ax = self._scatter_ax
        ax.clear()
        df = self._aligned_df
        mode = self._cb_color.currentText()
        if mode == "n_samples":
            c = df['n_samples'].values
            cmap = 'viridis'
            cbar_label = 'n_samples detected'
        elif mode.startswith("log10"):
            h_cols = [f'{self._intensity_col}_s{i}' for i in range(len(self._sample_names))]
            h_cols = [c for c in h_cols if c in df.columns]
            mean_h = df[h_cols].mean(axis=1).clip(lower=1.0)
            import numpy as _np
            c = _np.log10(mean_h.values)
            cmap = 'plasma'
            cbar_label = f'log10(mean {self._intensity_col})'
        else:  # fwhm_warning
            if 'fwhm_warning' in df.columns:
                c = df['fwhm_warning'].astype(int).values
            else:
                c = [0] * len(df)
            cmap = 'coolwarm'
            cbar_label = 'fwhm_warning (1=warning)'
        sc = ax.scatter(df['avg_rt'].values, df['avg_mz'].values,
                        c=c, cmap=cmap, s=8, alpha=0.7, edgecolors='none',
                        picker=True)
        ax.set_xlabel('RT (min)')
        ax.set_ylabel('m/z')
        ax.set_title(f'Aligned features ({len(df)})')
        # v0.5.17: create the colorbar once and update it on later redraws;
        #   previously a new colorbar was added every click -> bars stacked up
        #   and the plot kept shrinking.
        if self._scatter_cbar is None:
            self._scatter_cbar = self._scatter_fig.colorbar(sc, ax=ax,
                                                            label=cbar_label)
        else:
            self._scatter_cbar.update_normal(sc)
            self._scatter_cbar.set_label(cbar_label)
        # Mark selected
        if self._selected_idx is not None and 0 <= self._selected_idx < len(df):
            sel = df.iloc[self._selected_idx]
            ax.scatter([sel['avg_rt']], [sel['avg_mz']], s=120,
                       facecolors='none', edgecolors='red', linewidths=2,
                       zorder=10)
        self._scatter_canvas.draw_idle()

    def _on_scatter_click(self, event):
        if event.inaxes != self._scatter_ax:
            return
        if event.xdata is None or event.ydata is None:
            return
        # Find nearest feature in display coordinates
        import numpy as _np
        df = self._aligned_df
        # Normalize axes ranges so distance comparison is scale-invariant
        x_rng = max(df['avg_rt'].max() - df['avg_rt'].min(), 1e-6)
        y_rng = max(df['avg_mz'].max() - df['avg_mz'].min(), 1e-6)
        dx = (df['avg_rt'].values - event.xdata) / x_rng
        dy = (df['avg_mz'].values - event.ydata) / y_rng
        d2 = dx*dx + dy*dy
        idx = int(_np.argmin(d2))
        self._selected_idx = idx
        self._show_feature(idx)

    # ─── Bar (per-sample heights) ───────────────────────────
    def _show_feature(self, idx):
        row = self._aligned_df.iloc[idx]
        title = (f"Feature #{idx}  m/z = {row['avg_mz']:.4f}  "
                 f"RT = {row['avg_rt']:.3f} min  "
                 f"({int(row['n_samples'])}/{len(self._sample_names)} samples)")
        self._bar_label.setText(title)
        ax = self._bar_ax
        ax.clear()
        n = len(self._sample_names)
        cols = [f'{self._intensity_col}_s{i}' for i in range(n)]
        vals = [float(row.get(c, 0.0)) for c in cols]
        x = list(range(n))
        ax.bar(x, vals, color='steelblue', edgecolor='navy')
        ax.set_xticks(x)
        ax.set_xticklabels([s[-12:] for s in self._sample_names],
                           rotation=30, ha='right', fontsize=7)
        ax.set_ylabel(self._intensity_col)
        ax.grid(axis='y', alpha=0.3)
        self._bar_canvas.draw_idle()
        self._redraw_scatter()  # update red highlight
        # Trigger EIC if possible
        if self._raw_mode and self._sample_paths:
            self._update_eic(idx, row)
        else:
            ax = self._eic_ax
            ax.clear()
            ax.text(0.5, 0.5,
                    'EIC not available\n(no raw mzML loaded)',
                    ha='center', va='center', transform=ax.transAxes,
                    color='gray')
            ax.set_axis_off()
            self._eic_canvas.draw_idle()

    # ─── EIC (subprocess) ───────────────────────────────────
    def _update_eic(self, feat_idx, row):
        ax = self._eic_ax
        ax.clear()
        # Show "loading" message
        ax.text(0.5, 0.5, 'Loading EIC...', ha='center', va='center',
                transform=ax.transAxes, color='gray')
        ax.set_axis_off()
        self._eic_canvas.draw_idle()
        QApplication.processEvents()

        # Run subprocess for each sample (cached)
        import subprocess as _sub
        import tempfile as _tempfile
        worker = resolve_align_worker()
        eic_data = self._eic_cache.get(feat_idx)
        if eic_data is None:
            eic_data = {}
            for s_idx, mzml_path in enumerate(self._sample_paths):
                if not mzml_path or not Path(mzml_path).exists():
                    continue
                with _tempfile.NamedTemporaryFile(suffix='.csv',
                                                  delete=False, mode='w') as tmp:
                    out_csv = tmp.name
                try:
                    cmd = [
                        sys.executable, str(worker), "--eic",
                        "--mzml", str(mzml_path),
                        "--out", out_csv,
                        "--eic-mz", f"{row['avg_mz']:.6f}",
                        "--eic-mz-tol", f"{self._mz_tol_eic}",
                    ]
                    r = _sub.run(cmd, capture_output=True, text=True, timeout=120)
                    if r.returncode != 0:
                        continue
                    eic_data[s_idx] = pd.read_csv(out_csv)
                except Exception:
                    continue
                finally:
                    try: os.unlink(out_csv)
                    except Exception: pass
            self._eic_cache[feat_idx] = eic_data

        # Plot
        ax = self._eic_ax
        ax.clear()
        ax.set_axis_on()
        if not eic_data:
            ax.text(0.5, 0.5, 'EIC extraction failed', ha='center', va='center',
                    transform=ax.transAxes, color='red')
            ax.set_axis_off()
        else:
            for s_idx, df in eic_data.items():
                label = self._sample_names[s_idx][-12:]
                ax.plot(df['rt_min'], df['intensity'], label=label, lw=1.0, alpha=0.85)
            ax.axvline(float(row['avg_rt']), color='red', ls='--', lw=0.8,
                       alpha=0.5, label='avg_rt')
            ax.set_xlabel('RT (min)')
            ax.set_ylabel('Intensity')
            ax.set_title(f"EIC m/z = {row['avg_mz']:.4f} ± {self._mz_tol_eic} Da")
            ax.legend(fontsize=6, loc='best')
            ax.grid(alpha=0.3)
        self._eic_canvas.draw_idle()


class AdvancedSettingsDialog(QDialog):
    """v0.5.9: Hidden home for automatic parameter estimation, so the main
    Step 2 form stays simple. Each estimator measures a property of the
    detected peaks, writes the value into Step 2, and the user re-runs
    Detect Peaks to apply it."""

    def __init__(self, main_window, section="detection", parent=None):
        super().__init__(parent)
        self._mw = main_window
        self._section = section
        lay = QVBoxLayout(self)

        # fix34 で 2 列になったので、どちらの極性を触っているか明示する
        _pol = str(getattr(main_window, '_al_cur_pol', 'pos')).upper()
        if section == "alignment":
            self.setWindowTitle(
                f"Advanced settings - Alignment (Step 3)   [{_pol}]")
            self.resize(620, 460)
            intro = QLabel(
                f"Estimate alignment tolerances for {_pol} from the data. The "
                "button reads the cross-sample scatter of full-coverage "
                "features and <b>proposes</b> values. Nothing changes until you "
                "press [Apply proposed values]. If peaks are not detected yet, "
                "you will be offered to run detection first.")
            intro.setWordWrap(True)
            lay.addWidget(intro)
            b5 = QPushButton("Estimate  m/z & RT tolerance  from data")
            b5.setMinimumHeight(30)
            b5.setToolTip(
                "Suggest alignment tolerances from the cross-sample scatter of "
                "full-coverage features (99th percentile of m/z and RT spread).")
            b5.clicked.connect(self._do_align)
            lay.addWidget(b5)
        else:
            self.setWindowTitle(
                f"Advanced settings - Peak detection (Step 2)   [{_pol}]")
            self.resize(660, 560)
            intro = QLabel(
                f"Estimate detection parameters for {_pol} from the loaded data. "
                "Each button reads the measured properties of the detected "
                "peaks and <b>proposes</b> a value. Nothing changes until you "
                "press [Apply proposed values]. If peaks are not detected yet, "
                "you will be offered to run detection first.")
            intro.setWordWrap(True)
            lay.addWidget(intro)
            # Order matches the Step 2 form: Mass accuracy, Noise, S/N, FWHM
            b1 = QPushButton("Estimate  Mass trace m/z accuracy (ppm)  from data")
            b1.setMinimumHeight(30)
            b1.clicked.connect(self._do_mass)
            lay.addWidget(b1)
            b2 = QPushButton("Auto-tune  Noise threshold  by cross-sample consistency")
            b2.setMinimumHeight(30)
            b2.setToolTip(
                "Optimise the noise threshold using only the data: maximise the "
                "number of reproducible features present in all samples "
                "(CV<=20%).\nRuns one detection pass at low noise, then sweeps. "
                "Needs >=3 samples.")
            b2.clicked.connect(self._do_noise)
            lay.addWidget(b2)
            b3 = QPushButton("Auto-tune  Chromatographic S/N  by cross-sample consistency")
            b3.setMinimumHeight(30)
            b3.setToolTip(
                "Optimise S/N by the same consistency objective. Unlike noise, "
                "each S/N value needs a real detection pass, so this is "
                "slower.\nNeeds >=3 samples.")
            b3.clicked.connect(self._do_snr)
            lay.addWidget(b3)
            b4 = QPushButton("Estimate  Expected peak FWHM (sec)  from data")
            b4.setMinimumHeight(30)
            b4.clicked.connect(self._do_fwhm)
            lay.addWidget(b4)
            note = QLabel(
                "S/N = 3 is a standard detection convention, so manual is "
                "usually fine; an optional consistency tuner is provided. Noise "
                "has no single optimum in the abstract, but is tuned here "
                "against a chosen target - cross-sample consistency.")
            note.setWordWrap(True)
            note.setStyleSheet("color:#666;")
            lay.addWidget(note)

        self._log = QTextEdit()
        self._log.setReadOnly(True)
        self._log.setPlaceholderText("Estimation results appear here.")
        _f = QFont("Courier New")
        _f.setStyleHint(QFont.Monospace)
        _f.setPointSize(max(9, UI_FONT_POINT_SIZE - 1))
        self._log.setFont(_f)
        lay.addWidget(self._log, 1)

        # 推定した値は「提案」として並べ、Apply を押して初めて反映する
        self._prop_box = QGroupBox("Proposed values  (current → proposed)")
        _pv = QVBoxLayout(self._prop_box)
        self._prop_lbl = QLabel("No estimate yet.")
        self._prop_lbl.setStyleSheet("color:#555;")
        _pf = QFont("Courier New")
        _pf.setStyleHint(QFont.Monospace)
        _pf.setPointSize(max(9, UI_FONT_POINT_SIZE - 1))
        self._prop_lbl.setFont(_pf)
        _pv.addWidget(self._prop_lbl)
        lay.addWidget(self._prop_box)

        row = QHBoxLayout()
        self._btn_apply = QPushButton("Apply proposed values")
        self._btn_apply.setMinimumHeight(30)
        _fa = self._btn_apply.font()
        _fa.setBold(True)
        self._btn_apply.setFont(_fa)
        self._btn_apply.setToolTip(
            "Write the proposed values into the form. Nothing is changed "
            "before you press this.")
        self._btn_apply.setEnabled(False)
        self._btn_apply.clicked.connect(self._on_apply)
        row.addWidget(self._btn_apply)
        self._btn_discard = QPushButton("Discard")
        self._btn_discard.setToolTip("Throw the proposal away and keep the "
                                     "values you already have.")
        self._btn_discard.setEnabled(False)
        self._btn_discard.clicked.connect(self._on_discard)
        row.addWidget(self._btn_discard)
        row.addStretch()
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        row.addWidget(close)
        lay.addLayout(row)

    # ── 提案の表示と反映 ──────────────────────────────────
    def _refresh_proposal(self):
        try:
            txt = self._mw._al_proposal_summary()
        except Exception:
            txt = ""
        has = bool(txt)
        self._prop_lbl.setText(txt if has else "No estimate yet.")
        self._prop_lbl.setStyleSheet("color:#036;" if has else "color:#555;")
        self._btn_apply.setEnabled(has)
        self._btn_discard.setEnabled(has)

    def _on_apply(self):
        msg = self._mw._al_apply_proposal()
        self._append(msg)
        self._refresh_proposal()
        _re_run = ("Run alignment" if self._section == "alignment"
                   else "Detect Peaks")
        self._append(f"Now press [{_re_run}] in the main window to use them.")

    def _on_discard(self):
        self._mw._al_clear_proposal()
        self._refresh_proposal()
        self._append("Proposal discarded. The form is unchanged.")

    def _append(self, msg):
        if msg:
            self._log.append(msg + "\n")
        self._refresh_proposal()

    def _do_mass(self):
        self._append(self._mw._al_estimate_mass_accuracy_from_data(parent=self))

    def _do_fwhm(self):
        self._append(self._mw._al_estimate_fwhm_from_data(parent=self))

    def _do_noise(self):
        self._append(self._mw._al_estimate_noise_by_consistency(parent=self))

    def _do_snr(self):
        self._append(self._mw._al_estimate_snr_by_consistency(parent=self))

    def _do_align(self):
        self._append(self._mw._al_estimate_alignment_tolerances_from_data(parent=self))


class InspectDialog(QDialog):
    def __init__(self, fe: FileEntry, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Inspect: {fe.label}")
        self.resize(850, 600)
        lay = QVBoxLayout(self)

        info = QGroupBox("File Information")
        form = QFormLayout(info)
        for lbl, val in [
            ("Path:", str(fe.path)), ("Encoding:", fe.encoding),
            ("Header row:", str(fe.header_row)), ("Rows:", str(len(fe.df))),
            ("Columns:", str(len(fe.df.columns))),
            ("m/z column:", fe.mz_col or "(not found)"),
            ("RT column:",  fe.rt_col or "(not found)"),
        ]:
            form.addRow(lbl, QLabel(val))
        lay.addWidget(info)

        col_box = QGroupBox("Column Classification")
        col_lay = QVBoxLayout(col_box)
        tbl = QTableWidget(len(fe.df.columns), 4)
        tbl.setHorizontalHeaderLabels(
            ["Index", "Column Name", "Classification", "Action"])
        tbl.horizontalHeader().setStretchLastSection(True)
        tbl.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        meta_s  = set(fe.meta_indices)
        stat_s  = set(fe.stat_indices)
        fixed_s = set(fe.fixed_indices)
        for i, col in enumerate(fe.df.columns):
            if i in meta_s:
                cls_t, bg, act = "META",   QColor(255, 220, 220), "→ Excluded"
            elif i in stat_s:
                cls_t, bg, act = "STAT",   QColor(255, 235, 200), "→ Excluded"
            elif i in fixed_s:
                cls_t, bg, act = "FIXED",  QColor(220, 220, 255), "→ Excluded"
            elif fe.sample_start <= i <= fe.sample_end:
                cls_t, bg, act = "SAMPLE", QColor(220, 255, 220), "→ Output"
            else:
                cls_t, bg, act = "OTHER",  QColor(240, 240, 240), "→ Excluded"
            for c, v in enumerate([str(i), str(col), cls_t, act]):
                item = QTableWidgetItem(v)
                item.setBackground(bg)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                tbl.setItem(i, c, item)
        col_lay.addWidget(tbl)
        n_s = (fe.sample_end - fe.sample_start + 1
               if fe.sample_end >= fe.sample_start else 0)
        col_lay.addWidget(QLabel(
            f"Sample columns: {n_s} (index {fe.sample_start}–{fe.sample_end})"
            f"  |  Meta excluded: {len(fe.meta_indices)}"
            f"  |  Stat excluded: {len(fe.stat_indices)}"
        ))
        lay.addWidget(col_box)

        ov = QGroupBox("Manual Override (optional)")
        ov_lay = QHBoxLayout(ov)
        ov_lay.addWidget(QLabel("Sample start col:"))
        self._spin_s = QSpinBox()
        self._spin_s.setRange(0, len(fe.df.columns) - 1)
        self._spin_s.setValue(fe.sample_start)
        ov_lay.addWidget(self._spin_s)
        ov_lay.addWidget(QLabel("Sample end col:"))
        self._spin_e = QSpinBox()
        self._spin_e.setRange(0, len(fe.df.columns) - 1)
        self._spin_e.setValue(fe.sample_end)
        ov_lay.addWidget(self._spin_e)
        self._fe = fe
        btn = QPushButton("Apply Override")
        btn.clicked.connect(self._apply)
        ov_lay.addWidget(btn)
        lay.addWidget(ov)

        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        lay.addWidget(close)

    def _apply(self):
        self._fe.sample_start = self._spin_s.value()
        self._fe.sample_end   = self._spin_e.value()
        QMessageBox.information(
            self, "Override Applied",
            f"Sample range set to columns "
            f"{self._fe.sample_start}–{self._fe.sample_end}.")


# ════════════════════════════════════════════════════════════════════
#  MatchTableDialog
# ════════════════════════════════════════════════════════════════════

class MatchTableDialog(QDialog):
    """PreviewRTDialog から切り出した化合物テーブル専用ダイアログ。

    - 親(PreviewRTDialog)が保持する QTableWidget (_match_table) と
      QLabel (_lbl_match_count) を自分の layout に格納してホストする。
    - Modeless で開き、Annotation Pipeline 本体と同時に操作可能。
    - close 時は deleteLater せず hide のみ。再表示で同じインスタンスを使い回す。
    """

    def __init__(self, parent: QDialog,
                 table_widget: QTableWidget,
                 count_label: QLabel,
                 title: str = "Compounds Table"):
        super().__init__(parent)
        self.setWindowTitle(title)
        # Window フラグ: 親より独立したサブウィンドウ。閉じる/最小化/最大化可能。
        self.setWindowFlags(
            Qt.Window
            | Qt.WindowTitleHint
            | Qt.WindowSystemMenuHint
            | Qt.WindowCloseButtonHint
            | Qt.WindowMinMaxButtonsHint)
        self.setModal(False)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)
        # SizePolicy を再調整: ダイアログ内ではテーブルを伸縮させる
        table_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        lay.addWidget(table_widget, stretch=1)
        lay.addWidget(count_label, stretch=0)

        self.resize(720, 560)

    def closeEvent(self, event):
        """× ボタン: deleteLater せず hide のみ(子ウィジェットを保持)。"""
        event.ignore()
        self.hide()


# ════════════════════════════════════════════════════════════════════
#  PreviewRTDialog
# ════════════════════════════════════════════════════════════════════

class PreviewRTDialog(QDialog):
    """RT vs m/z プレビュー + ライブラリオーバーレイ + タスク登録"""

    task_ready     = Signal(dict)
    settings_saved = Signal(dict)   # ダイアログ終了時に現在のPreview設定を通知
    reserved_updated = Signal(str, set)  # (file_tag, reserved_coords) を通知
    # (file_tag, {loser_class: set of (rt, mz)}) を通知
    coherence_updated = Signal(str, dict)
    # match_df を MainWindow に通知(_export_one で final_status='kept' を使う)
    match_df_updated = Signal(str, object)

    # チェックボックスの特殊タグ
    _TAG_NOT_MATCHED            = "__not_matched__"
    _TAG_RT_OUTLIERS            = "__rt_outliers__"
    _TAG_COHERENCE_OUTLIERS     = "__coherence_outliers__"
    _TAG_CONFLICTS              = "__conflicts__"
    _TAG_LOWCONF                = "__lowconf__"
    _TAG_SHOW_IS_REJECTED       = "__show_is_rejected__"
    _TAG_SHOW_ADDUCT_REJECTED   = "__show_adduct_rejected__"
    _TAG_CURATED_VIEW           = "__curated_view__"
    _TAG_IS_ONLY                = "__is_only__"

    def __init__(
        self,
        fe: FileEntry,
        all_entries: list[FileEntry],
        initial_settings: dict | None = None,
        all_initial_settings: dict | None = None,
        parent=None,
    ):
        # parent=None で top-level Qt window として登録。Windows での
        # タスクバーボタンを安定して出すため、QDialog の owner-chain を
        # 切る。host (MainWindow) への delegate は self._host で行う。
        super().__init__(None)
        self._host       = parent   # MainWindow への参照を明示保持
        self.fe           = fe
        self._all_entries = all_entries
        # 全モード分の初期設定を保持
        # (key: ion_mode, value: preview_settings dict)
        self._all_initial_settings = dict(all_initial_settings or {})
        # "Preview RT" → "Annotation Pipeline" にリネーム
        # 実行ファイル名を固定名にしたので、版番号はここに出す。
        # 見えているのはこのウィンドウで、MainWindow のタイトル
        # ("LipidZoner v…")は隠れていて読めないため。
        self.setWindowTitle(ap_window_title(fe.label))
        # 独立した top-level window として表示。Qt.Window を明示しないと
        # Qt.Dialog 扱いになり、Windows のタスクバーに固有ボタンが出ない。
        # 親 (MainWindow) が hidden のため、最小化すると行き場を失って
        # 復元不能になる現象が起きる。
        self.setWindowFlags(
            Qt.Window
            | Qt.WindowTitleHint
            | Qt.WindowSystemMenuHint
            | Qt.WindowCloseButtonHint
            | Qt.WindowMinMaxButtonsHint)
        self.resize(1100, 750)

        # 状態変数
        self._selected_range: tuple[float, float] | None = None
        # 範囲選択の永続ハイライト矩形 (mode -> Artist)
        self._persistent_range_artists: dict[str, object] = {}
        self._lib_path: Path | None = None
        self._match_df: pd.DataFrame | None = None
        self._class_colors: dict[str, str] = {}
        self._matched_coords: set          = set()
        self._outlier_coords: dict[str, set] = {}
        # Coherence 外れ値（link × クラス別 pooled 回帰からの残差ベース、
        # 旧 _coherence_outliers の置換）
        self._coherence_outliers: dict[str, set] = {}
        # Coherence Engine の pooled fit 結果 {(class, link): {alpha,beta,gamma,rmse,n}}
        self._coherence_models: dict[tuple[str, str], dict] = {}
        # 衝突帰属結果 {(rt_r6, mz_r6): {winner_*, loser_*, confidence, low_confidence, ...}}
        self._coherence_assignments: dict[tuple, dict] = {}
        # 有効な衝突ペア(ペア選択ダイアログで ON のもの)
        self._coherence_pairs: list[tuple[str, str]] = []
        # ライブラリから自動検出された全衝突ペア(UI 表示用)
        self._coherence_all_pairs: list[tuple[str, str]] = []
        self._conflict_coords: set           = set()
        # 競合詳細: {(rt, mz): [{lipid_class, compound, adduct}, ...]}
        self._conflict_map: dict[tuple, list[dict]] = {}
        # Filter ごとのクラス別 exempt 集合。
        # チェックを外したクラスはその Filter で判定対象外となり、
        # 既存の rejected フラグが付かない(既に付いていたら剥がす)。
        # キー: 'is' / 'adduct' / 'rt_outlier' / 'coherence'
        # 値: そのフィルタを適用しないクラス名の set
        # Curated view mode フラグ(ボタンの checked 状態を
        # 持たない代わりに内部で管理)
        self._curated_mode: bool = False
        self._filter_class_exemptions: dict[str, set[str]] = {
            'is': set(),
            'adduct': set(),
            'rt_outlier': set(),
            'coherence': set(),
        }
        # auto-replay 中は IS Filter のプレビューダイアログ等を
        # スキップして保存済み choices で直接処理するためのフラグ
        self._auto_replay_in_progress: bool = False
        # Run All 実行中フラグ。③ のプレビューダイアログを抑制する。
        self._run_all_in_progress: bool = False
        # ③ が 2 パスマッチングまで到達したか。Run All の中止判定用。
        self._is_filter_applied_ok: bool = False
        self._overlay_artists: list        = []
        self._annot                        = None
        # チェックボックス辞書: {クラス名 or タグ: QCheckBox}
        self._class_checks: dict[str, QCheckBox] = {}
        # モード切替で _class_checks が再構築される際に、
        # ユーザーが明示的に切り替えたチェック状態を保持するための永続辞書。
        # 他モード固有のクラス(モード切替で _class_checks から消える)も
        # この辞書には残るため、再登場時に復元できる。
        self._user_filter_prefs: dict[str, bool] = {}
        # Fix RT 設定: {lipid_class: (reference_rt, tolerance_min, source)}
        # source = "manual" | "auto_is"
        self._class_ref_rt: dict[str, tuple[float, float, str]] = {}
        # reserved_coords: IS優先マッチングで帰属された観測ピークの座標
        # （エクスポート時の除外に使用）
        self._reserved_coords: set = set()
        # IS Filter プレビューダイアログでの選択状態
        # {cls: [{compound, adduct, use, manual_rt}, ...]}
        self._is_filter_choices: dict[str, list[dict]] = {}
        # ランタイム IS 由来アダクト fingerprint
        # ② IS Filter 完了時に _compute_is_fingerprints() で生成。
        # 形式は CURATED_ADDUCT_FINGERPRINTS と同じ:
        #   {class_name: {adduct: {'expected': bool, 'weight': float,
        #                          'ratio_target': None}}}
        # Adduct Ion Filter (Step ⑤) が衝突解決の参照として使う。
        # これは **表示用のマージ済みビュー**。正は下の
        # _runtime_is_fingerprints_by_mode（(クラス, 極性) キー）。
        # クラス名だけをキーにすると、両極性に IS があるクラス
        # (検証データでは HexCer / Cer) が後に走った極性で上書きされる。
        self._runtime_is_fingerprints: dict[
            str, dict[str, dict[str, object]]] = {}
        # (クラス, 極性) → {adduct: spec}
        self._runtime_is_fingerprints_by_mode: dict[
            tuple[str, str], dict[str, dict[str, object]]] = {}
        # Adduct Ion Filter の判定結果
        #   {coord: {'winner': str|None, 'losers': set, 'scores': dict,
        #            'observations': dict, 'skipped': bool}}
        self._adduct_attribution: dict = {}
        # 手動で winner 指定された座標(プロット右クリックで設定)
        self._manual_winner_coords: set = set()
        # 各脂質クラスの定量イオンモード選択
        # {class_name: 'pos' | 'neg' | 'skip'}
        # フィルタ後・タスク登録前に「PCはposで定量」「LPSはskip」等を決定する。
        self._quant_ion_choices: dict[str, str] = {}
        # redraw 抑制フラグ(両モード処理中の中間描画を回避)
        self._suppress_overlay_draw: bool = False
        # per-mode 状態スナップショット
        # キー = ion_mode ('pos' | 'neg')、値 = mode ごとの状態 dict
        # _on_mode_radio_changed で snapshot 保存/復元する。
        self._mode_state: dict[str, dict] = {}
        self._active_mode: str = fe.ion_mode

        rt = fe.rt_array()
        mz = fe.mz_array()

        # 上位レイアウトは QTabWidget(Analysis / Utility)
        main_lay = QVBoxLayout(self)
        main_lay.setContentsMargins(0, 0, 0, 0)
        main_lay.setSpacing(0)

        # Session の Save/Load は QTabWidget の Session タブに移動。
        # ここでは Ctrl+S / Ctrl+O ショートカット用に QAction を AP 自身に
        # 付与しておく(タブ表示なしでも shortcut が効く)。
        _act_save = QAction("Save LipidZoner Session File…", self)
        _act_save.setShortcut(QKeySequence("Ctrl+S"))
        _act_save.triggered.connect(self._ap_save_session)
        self.addAction(_act_save)
        _act_load = QAction("Load LipidZoner Session File…", self)
        _act_load.setShortcut(QKeySequence("Ctrl+O"))
        _act_load.triggered.connect(self._ap_load_session)
        self.addAction(_act_load)

        self._tabs = QTabWidget(self)
        main_lay.addWidget(self._tabs)

        # ── Session タブ(最初に配置。起動時デフォルトは Analysis に切替) ──
        session_tab = QWidget()
        session_lay = QVBoxLayout(session_tab)
        session_lay.setContentsMargins(12, 12, 12, 12)
        session_lay.setSpacing(8)
        sess_box = QGroupBox("LipidZoner Session File")
        sess_box_lay = QVBoxLayout(sess_box)
        sess_box_lay.setSpacing(6)
        _sess_desc = QLabel(
            "<b>Save</b>: Save the unified LipidZoner Session File (JSON). "
            "Includes input files, parameters, ion selection, and filter exemptions.<br>"
            "<b>Load (Full)</b>: Load everything from a session file "
            "(input files, parameters, ion selection, filter exemptions, tasks).<br>"
            "<b>Load (Settings only)</b>: Apply only <i>parameters</i> and "
            "<i>filter exemptions</i> from a session file to the currently loaded data. "
            "Useful when reusing pipeline settings across different datasets. "
            "Run All is invoked automatically after loading."
            "<br><b>Note:</b> Before clicking <b>Load Settings Only</b>, "
            "load your data and lipid library in the <b>Analysis</b> tab first.")
        _sess_desc.setStyleSheet("color:#444;")
        _sess_desc.setWordWrap(True)
        _sess_desc.setTextFormat(Qt.RichText)
        sess_box_lay.addWidget(_sess_desc)
        _sess_btn_row = QHBoxLayout()
        _sess_btn_save = QPushButton("Save LipidZoner Session File…")
        _sess_btn_save.setToolTip(
            "Save the unified LipidZoner Session File (JSON).  [Ctrl+S]")
        _sess_btn_save.clicked.connect(self._ap_save_session)
        _sess_btn_row.addWidget(_sess_btn_save)
        _sess_btn_load = QPushButton("Load LipidZoner Session File…")
        _sess_btn_load.setToolTip(
            "Load all sections (input files, parameters, ion selection,\n"
            "filter exemptions, tasks) from a session file.  [Ctrl+O]")
        _sess_btn_load.clicked.connect(self._ap_load_session)
        _sess_btn_row.addWidget(_sess_btn_load)
        _sess_btn_load_settings = QPushButton("Load Settings Only…")
        _sess_btn_load_settings.setToolTip(
            "Apply only parameters and filter exemptions from a session file\n"
            "to the currently loaded data, then run the full pipeline.\n"
            "Input files, tasks, and ion selection are NOT changed.")
        _sess_btn_load_settings.clicked.connect(self._ap_load_session_settings_only)
        _sess_btn_row.addWidget(_sess_btn_load_settings)
        _sess_btn_row.addStretch()
        sess_box_lay.addLayout(_sess_btn_row)
        session_lay.addWidget(sess_box)
        session_lay.addStretch()
        # Session はトップレベルタブではなく Utility 内サブタブへ移す。
        # ここでは widget を保持するだけで、実際の addTab は
        # _build_utility_tab() の末尾で行う。
        self._session_tab_widget = session_tab

        # ── Analysis タブ ────────────────────────────────────────────
        analysis_tab = QWidget()
        analysis_lay = QVBoxLayout(analysis_tab)
        analysis_lay.setContentsMargins(6, 6, 6, 6)
        analysis_lay.setSpacing(4)
        self._tabs.addTab(analysis_tab, TAB_LABELS['annot'])

        # ── ライブラリ操作パネル──────────
        lib_box = QGroupBox("Lipid Annotation")
        lib_box.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        lib_lay = QVBoxLayout(lib_box)
        lib_lay.setContentsMargins(8, 4, 8, 4)
        lib_lay.setSpacing(2)

        # ── Data / Library / Clear / Run All を 1 行に統合(左詰め) ──
        top_row = QHBoxLayout()
        top_row.setSpacing(6)

        # Data セクション
        top_row.addWidget(QLabel("Data:"))
        self._lbl_data_summary = QLabel("(no data loaded)")
        self._lbl_data_summary.setStyleSheet("color:#666;")
        self._lbl_data_summary.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        top_row.addWidget(self._lbl_data_summary)
        self._btn_load_data = QPushButton("Load data…")
        self._btn_load_data.setToolTip(
            "Open the Manage Data Files dialog to add / remove / inspect\n"
            "MS-DIAL data files (pos + neg pair).")
        self._btn_load_data.clicked.connect(self._open_manage_files_dialog)
        top_row.addWidget(self._btn_load_data)

        # Library セクション
        top_row.addWidget(QLabel("Library:"))
        self._lbl_lib = QLabel("No library loaded")
        self._lbl_lib.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self._lbl_lib.setStyleSheet("color: #555; padding: 0 4px;")
        top_row.addWidget(self._lbl_lib)
        btn_lib = QPushButton("Load library…")
        btn_lib.clicked.connect(self._load_library)
        top_row.addWidget(btn_lib)

        # Clear / Run All
        btn_clear = QPushButton("Clear")
        btn_clear.clicked.connect(self._clear_overlay)
        top_row.addWidget(btn_clear)

        self._btn_run_all = QPushButton("▶ Run All")
        self._btn_run_all.setToolTip(
            "Run all steps from ① Match Overlay to ⑥ Coherence Filter\n"
            "in sequence (full restart).\n"
            "Each call starts from scratch using current parameter values.")
        self._btn_run_all.clicked.connect(self._run_all)
        top_row.addWidget(self._btn_run_all)

        # 左詰め: 右側を stretch で埋める
        top_row.addStretch()

        lib_lay.addLayout(top_row)

        # ステップ行共通のラベル幅（統一感のため）
        STEP_LABEL_W = 150
        INPUT_W      = 110

        # ── 横一列ステップ行 + Parameters ダイアログ ──────────
        # spinbox 群は不可視 holder に保持し、Parameters ダイアログが UI で
        # 表示・編集する。各ステップメソッドは self._spin_xxx.value() を読む。

        # ── 不可視 spinbox holder(各 spinbox は self.* で参照可能)──────
        self._params_holder = QWidget()
        self._params_holder.setVisible(False)
        # ① ppm tol
        self._spin_ppm = QDoubleSpinBox(self._params_holder)
        self._spin_ppm.setRange(0.1, 200.0)
        self._spin_ppm.setDecimals(1)
        self._spin_ppm.setSingleStep(1.0)
        self._spin_ppm.setValue(10.0)
        self._spin_ppm.setSuffix(" ppm")
        self._spin_ppm.setFixedWidth(110)
        # ボイドボリューム(非保持成分・キャリーオーバー)の RT 上限。
        # メソッド固有の値なのでユーザーが指定できるようにする。
        self._spin_void_rt = QDoubleSpinBox(self._params_holder)
        self._spin_void_rt.setRange(0.0, 10.0)
        self._spin_void_rt.setDecimals(3)
        self._spin_void_rt.setSingleStep(0.05)
        self._spin_void_rt.setValue(VOID_RT_CUTOFF_DEFAULT)
        self._spin_void_rt.setSuffix(" min")
        self._spin_void_rt.setSpecialValueText("(off)")
        self._spin_void_rt.setFixedWidth(110)
        self._spin_void_rt.setToolTip(
            "Observed peaks eluting earlier than this RT are not taken as\n"
            "library match candidates. This removes unretained compounds\n"
            "and carry-over that appear in the void volume region.\n\n"
            "The value is specific to the column, flow rate and method,\n"
            "so decide it from the RT distribution of your own data.\n"
            "Set 0 to disable.\n\n"
            f"The default {VOID_RT_CUTOFF_DEFAULT} min comes from the SFC verification data:\n"
            "  RT 0.35-0.50 holds a dense band of unretained compounds,\n"
            "  0.55-0.75 is a valley, and the real signal starts at 0.80\n"
            "  (earliest IS: SE at 0.818).")
        # ② IS tol
        self._spin_is_tol = QDoubleSpinBox(self._params_holder)
        self._spin_is_tol.setRange(0.001, 5.0)
        self._spin_is_tol.setDecimals(3)
        self._spin_is_tol.setSingleStep(0.05)
        self._spin_is_tol.setValue(0.10)
        self._spin_is_tol.setSuffix(" min")
        self._spin_is_tol.setFixedWidth(110)
        # IS 自動選択に強度を加味するか(既定 OFF)。
        #   手動選択(Candidate Picker → Manual RT)が最優先なのは変わらない。
        #   これはピン留めが無い IS のフォールバックの挙動を変えるだけ。
        self._cb_is_auto_intensity = QCheckBox(
            "Use intensity in IS auto-pick", self._params_holder)
        self._cb_is_auto_intensity.setChecked(False)
        self._cb_is_auto_intensity.setToolTip(
            "Changes how a non-pinned IS is chosen when several\n"
            "candidates fall inside the ppm window.\n\n"
            "OFF (default): pick the smallest |Δm/z| (legacy behaviour)\n"
            f"ON: drop candidates below {int(IS_AUTO_MIN_REL_INT*100)}% of the strongest\n"
            "    candidate, then pick the smallest |Δm/z|\n\n"
            "On the verification data one IS improved and one got worse.\n"
            "For an IS with several candidates, choosing it explicitly in\n"
            "the Candidate Picker is the reliable option.")
        # ⑤ Adduct Ion Filter — ppm / RT / int(将来 vote)
        self._spin_adduct_ppm = QDoubleSpinBox(self._params_holder)
        self._spin_adduct_ppm.setRange(0.1, 200.0)
        self._spin_adduct_ppm.setDecimals(1)
        self._spin_adduct_ppm.setSingleStep(1.0)
        self._spin_adduct_ppm.setValue(10.0)
        self._spin_adduct_ppm.setFixedWidth(80)
        self._spin_adduct_rt = QDoubleSpinBox(self._params_holder)
        self._spin_adduct_rt.setRange(0.001, 5.0)
        self._spin_adduct_rt.setDecimals(3)
        self._spin_adduct_rt.setSingleStep(0.01)
        self._spin_adduct_rt.setValue(0.05)
        self._spin_adduct_rt.setFixedWidth(80)
        self._spin_adduct_int = QDoubleSpinBox(self._params_holder)
        self._spin_adduct_int.setRange(0.0, 1e9)
        self._spin_adduct_int.setDecimals(0)
        self._spin_adduct_int.setSingleStep(100.0)
        self._spin_adduct_int.setValue(0.0)
        self._spin_adduct_int.setFixedWidth(90)
        # vote_threshold を可視化。
        # binary discriminator スコアになり、score は整数値
        # (-2, -1, 0, 1, 2 ...)。default=1.0 で「勝者は次点より少なくとも
        # 1 つ多い discriminator が必要」を意味する直感的な閾値になる。
        self._spin_adduct_vote = QDoubleSpinBox(self._params_holder)
        self._spin_adduct_vote.setRange(0.0, 10.0)
        self._spin_adduct_vote.setDecimals(1)
        self._spin_adduct_vote.setSingleStep(0.5)
        self._spin_adduct_vote.setValue(1.0)
        self._spin_adduct_vote.setFixedWidth(80)
        # ④ MAD×
        self._spin_iqr = QDoubleSpinBox(self._params_holder)
        self._spin_iqr.setRange(0.1, 10.0)
        self._spin_iqr.setDecimals(1)
        self._spin_iqr.setSingleStep(0.5)
        self._spin_iqr.setValue(3.0)
        self._spin_iqr.setFixedWidth(110)
        # ⑥ σ
        self._spin_sigma = QDoubleSpinBox(self._params_holder)
        self._spin_sigma.setDecimals(1)
        self._spin_sigma.setRange(1.0, 10.0)
        self._spin_sigma.setSingleStep(0.5)
        self._spin_sigma.setValue(3.0)
        self._spin_sigma.setFixedWidth(70)
        # 緩い帰属。⑥ が採点できる候補を 1 つしか作れない衝突を
        #   △ のまま残さず、その候補へ帰属する。ノンターゲット解析で
        #   「ある程度の間違いは許容するので △ を減らしたい」という要望。
        self._cb_loose_coherence = QCheckBox(
            "Loose attribution for unscorable conflicts",
            self._params_holder)
        self._cb_loose_coherence.setChecked(COH_LOOSE_ATTRIBUTION_DEFAULT)
        self._cb_loose_coherence.setToolTip(
            "What to do with a conflict where only one candidate can be\n"
            "scored (the other class has no RT model for its link).\n\n"
            "ON (default): attribute the spot to the candidate that could\n"
            "    be scored. No comparison was made, so it is reported as\n"
            "    low confidence (sigma = 0, decided_by = coherence:single).\n"
            "OFF: leave the spot undecided.\n\n"
            "Intended for non-target work, where an attribution with a\n"
            "confidence label is more useful than no attribution.")
        # 二峰性クラス。DG は sn 位置異性体(1,3- と 1,2-)が部分
        #   分離するので、1 つの和組成が 2 本のピークになる。ここに
        #   挙げたクラスは ⑥ が 2 系列(γ と γ+Δ)でモデル化し、
        #   判定 B の「1 compound = 1 ピーク」を系列単位に緩める。
        self._ed_two_series = QLineEdit(
            ', '.join(COH_TWO_SERIES_CLASSES_DEFAULT), self._params_holder)
        self._ed_two_series.setFixedWidth(170)
        self._ed_two_series.setToolTip(
            "Lipid classes whose regioisomers partly separate, so one sum\n"
            "composition appears as two peaks (comma separated).\n\n"
            "DG is the default: its sn-1,3 and sn-1,2 isomers elute about\n"
            "0.05 min apart in SFC. For a class listed here the RT model\n"
            "gets two parallel series (gamma and gamma+delta), delta is\n"
            "measured from the data, and both peaks are kept instead of\n"
            "only the one closer to a single-series prediction.\n\n"
            "The two series are only used when the data supports them:\n"
            "each series needs at least a few points, the spacing must be\n"
            "within range, and the fit has to improve several fold.\n\n"
            "Leave empty to disable.")
        # 定量エクスポートで 2 本を合算するか。プリカーサー定量
        #   では位置異性体の合計が求める量なので既定で ON。IS も同じ
        #   規則で合算される。
        self._cb_sum_two_series = QCheckBox(
            "Sum peaks in quant export", self._params_holder)
        self._cb_sum_two_series.setChecked(QUANT_SUM_TWO_SERIES_DEFAULT)
        self._cb_sum_two_series.setToolTip(
            "For the bimodal classes above, add up the sample values of\n"
            "the peaks that share one compound when exporting the quant\n"
            "ion table.\n\n"
            "ON (default): one row per compound, sample values summed.\n"
            "    The internal standard is summed the same way, so the\n"
            "    normalisation denominator is also the sum.\n"
            "    Columns 'peaks' and 'RT (peaks)' show what was added.\n"
            "OFF: one row per peak, nothing summed.\n\n"
            "This only affects the quant ion table. The non-target table\n"
            "always lists one row per feature.")

        # ── 横一列ステップ行(ステップ名そのものをボタン化) ──
        # 「① Match & Overlay」「② IS Filter」… のようにステップ名 = 押下で
        # 該当ステップを両モード実行するボタンとして配置する。Apply 文字は廃止。
        step_row = QHBoxLayout()
        step_row.setSpacing(6)
        step_row.setContentsMargins(0, 2, 0, 2)

        # 各ステップグループ間に縦区切りを挟むヘルパ
        def _add_separator():
            sep = QFrame(self)
            sep.setFrameShape(QFrame.VLine)
            sep.setFrameShadow(QFrame.Sunken)
            sep.setStyleSheet("color: #C8C6BD;")
            step_row.addWidget(sep)

        # ① Match & Overlay
        # 他ステップと同様 checkable 化(実行済みで青表示)
        self._btn_match = QPushButton("① Match & Overlay")
        self._btn_match.setCheckable(True)
        self._btn_match.setToolTip(
            "Step ① Match Overlay: match every library entry to the\n"
            "closest observed peak by m/z error.\n"
            "Click to run for all loaded modes (pos + neg).\n"
            "Button turns blue once matched. Click again re-runs match.")
        self._btn_match.clicked.connect(self._on_match_clicked)
        step_row.addWidget(self._btn_match)
        _add_separator()

        # ② Conflict Search
        # トグル化(両モード一括 ON/OFF、青で適用中表示)
        self._btn_conflict = QPushButton("② Conflict Search")
        self._btn_conflict.setCheckable(True)
        self._btn_conflict.setToolTip(
            "Step ② Conflict Search: detect spots where multiple lipid\n"
            "classes are assigned to the same observed peak.\n"
            "Toggle: click to apply (button stays pressed), click again\n"
            "to revert. Reverting also clears Step ⑤'s attribution.")
        self._btn_conflict.clicked.connect(self._toggle_conflict_search)
        self._btn_conflict.setEnabled(False)
        step_row.addWidget(self._btn_conflict)
        _add_separator()

        # ③ IS Filter
        # トグル化
        self._btn_is_filter = QPushButton("③ IS Filter")
        self._btn_is_filter.setCheckable(True)
        self._btn_is_filter.setToolTip(
            "Step ③ IS Filter: switch to IS-based two-pass matching.\n"
            "Also generates IS-derived adduct fingerprints used by Step ⑤.\n"
            "Toggle: click to apply (opens IS preview dialog),\n"
            "click again to revert to non-IS-filtered match.")
        self._btn_is_filter.clicked.connect(self._toggle_is_filter)
        self._btn_is_filter.setEnabled(False)
        step_row.addWidget(self._btn_is_filter)
        _add_separator()

        # ④ RT Outlier Filter
        # トグル化(既存 _rt_outlier_search が内部で toggle)
        self._btn_rt_outlier = QPushButton("④ RT Outlier Filter")
        self._btn_rt_outlier.setCheckable(True)
        self._btn_rt_outlier.setToolTip(
            "Step ④ RT Outlier Filter: detect within-class RT outliers.\n"
            "Toggle: click to apply, click again to revert.")
        self._btn_rt_outlier.clicked.connect(
            lambda: self._apply_to_all_modes(self._rt_outlier_search))
        self._btn_rt_outlier.setEnabled(False)
        step_row.addWidget(self._btn_rt_outlier)
        _add_separator()

        # ⑤ Adduct Ion Filter (+ サブボタン Patterns / Report は同グループ)
        # トグル化(両モード一括 ON/OFF)
        self._btn_adduct_filter = QPushButton("⑤ Adduct Ion Filter")
        self._btn_adduct_filter.setCheckable(True)
        self._btn_adduct_filter.setToolTip(
            "Step ⑤ Adduct Ion Filter: resolve conflicts using observed\n"
            "adduct ion patterns (requires Step ③ IS Filter first).\n"
            "Toggle: click to apply (button stays pressed), click again\n"
            "to revert. Applies to both pos and neg modes.")
        self._btn_adduct_filter.clicked.connect(self._toggle_adduct_filter)
        self._btn_adduct_filter.setEnabled(False)
        step_row.addWidget(self._btn_adduct_filter)
        # Patterns / Report ボタンは Advanced Dialog に移動
        self._btn_adduct_patterns = QPushButton("Patterns…")
        self._btn_adduct_patterns.clicked.connect(self._open_adduct_patterns)
        self._btn_adduct_patterns.setEnabled(False)
        self._btn_adduct_patterns.setVisible(False)   # 非表示(widget は保持)
        self._btn_adduct_report = QPushButton("Report…")
        self._btn_adduct_report.clicked.connect(self._open_adduct_report)
        self._btn_adduct_report.setEnabled(False)
        self._btn_adduct_report.setVisible(False)
        _add_separator()

        # ⑥ Coherence Filter (+ サブボタン Pairs / Report は同グループ)
        # トグル化(既存 _coherence_search が内部で toggle)
        self._btn_coherence = QPushButton("⑥ Coherence Filter")
        self._btn_coherence.setCheckable(True)
        self._btn_coherence.setToolTip(
            "Step ⑥ Coherence Filter: pooled multivariate regression\n"
            "for cross-class RT-based attribution.\n"
            "Toggle: click to apply, click again to revert.")
        self._btn_coherence.clicked.connect(
            lambda: self._apply_to_all_modes(self._coherence_search))
        self._btn_coherence.setEnabled(False)
        step_row.addWidget(self._btn_coherence)
        # Pairs / Report ボタンは Advanced Dialog に移動
        self._btn_coherence_pairs = QPushButton("Pairs…")
        self._btn_coherence_pairs.clicked.connect(self._open_coherence_pairs)
        self._btn_coherence_pairs.setEnabled(False)
        self._btn_coherence_pairs.setVisible(False)
        self._btn_coherence_report = QPushButton("Report…")
        self._btn_coherence_report.clicked.connect(self._open_coherence_report)
        self._btn_coherence_report.setEnabled(False)
        self._btn_coherence_report.setVisible(False)

        # ⑦ Fix RT 機能撤廃
        # (IS Filter の手動 RT 入力で代替可能なため不要)

        # ⑦ Select Quant Ion (フィルタ後の定量モード選択)
        # 他のステップボタンと同じく完了で青く表示するため
        # checkable 化。click で常にダイアログを開き、ダイアログ後に
        # _quant_ion_choices の有無で checked 状態を再同期する。
        _add_separator()
        self._btn_quant_ion = QPushButton("⑦ Select Quant Ion…")
        self._btn_quant_ion.setCheckable(True)
        self._btn_quant_ion.setToolTip(
            "After filtering and manual attribution, decide which ion mode\n"
            "(pos / neg) to use for quantification per lipid class.\n"
            "Uncertain classes can be excluded (skip). Save / Load JSON.")
        self._btn_quant_ion.clicked.connect(self._open_quant_ion_selector)
        step_row.addWidget(self._btn_quant_ion)

        # 右側に Advanced ボタン配置
        step_row.addStretch()
        btn_params = QPushButton("Advanced…")
        btn_params.setToolTip(
            "Open the Advanced dialog: edit step parameters and access\n"
            "step-specific reports (③ IS rejected, ⑤ Patterns/Report,\n"
            "④ RT outlier list, ⑥ Pairs/Report).")
        btn_params.clicked.connect(self._open_parameters_dialog)
        step_row.addWidget(btn_params)

        lib_lay.addLayout(step_row)

        analysis_lay.addWidget(lib_box, stretch=0)

        # RT/mz が取得できない場合
        if rt is None or mz is None:
            analysis_lay.addWidget(
                QLabel("RT or m/z column not found – cannot preview."))
            self._restore_settings(initial_settings)
            return

        mask = ~(np.isnan(rt) | np.isnan(mz))
        self._rt = rt[mask]
        self._mz = mz[mask]
        # サンプル強度行列(同じ mask を適用)。
        # 各候補の mean_intensity / sample_intensities を表示するため。
        # shape: (n_valid_peaks, n_samples)
        try:
            sample_cols = fe.sample_columns()
            self._sample_intensities = (
                fe.df[sample_cols].apply(pd.to_numeric, errors="coerce")
                  .values[mask].astype(float))
            self._sample_col_names = list(sample_cols)
        except Exception:
            self._sample_intensities = None
            self._sample_col_names = []

        # ── 左右分割レイアウト（QSplitter でドラッグ調整可能）────────
        h_split = QSplitter(Qt.Horizontal)
        h_split.setChildrenCollapsible(False)
        analysis_lay.addWidget(h_split, stretch=1)

        # 左: 散布図 + ヒストグラム
        left_w = QWidget()
        left_lay = QVBoxLayout(left_w)
        left_lay.setContentsMargins(0, 0, 0, 0)
        h_split.addWidget(left_w)

        self._fig = Figure(constrained_layout=True)
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._toolbar = NavigationToolbar(self._canvas, left_w)
        self._toolbar.setMaximumHeight(32)
        left_lay.addWidget(self._toolbar)
        left_lay.addWidget(self._canvas, stretch=1)

        # + 修正: 散布図 2 段(物理位置固定: pos=top、neg=bottom)
        # _ax_sc_pos / _ax_sc_neg は物理位置に紐付き、active モード切替後も
        # 不変。_ax_sc / _ax_sc_other はエイリアスとして active モードに従う。
        gs = self._fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.20)
        self._ax_sc_pos = self._fig.add_subplot(gs[0])  # 上 = pos(常時)
        self._ax_sc_neg = self._fig.add_subplot(
            gs[1], sharex=self._ax_sc_pos)               # 下 = neg(常時)

        # active モードに従って _ax_sc / _ax_sc_other を割り当てる
        if self.fe.ion_mode == 'pos':
            self._ax_sc       = self._ax_sc_pos  # primary = top
            self._ax_sc_other = self._ax_sc_neg  # secondary = bottom
        else:
            self._ax_sc       = self._ax_sc_neg  # primary = bottom
            self._ax_sc_other = self._ax_sc_pos  # secondary = top

        # 後方互換: _ax_hist 属性(廃止だが None で残置)
        self._ax_hist = None

        # 各物理軸専用の base_scatter を作成
        # primary 軸の base_scatter は self._rt/self._mz(active mode のデータ)
        if self.fe.ion_mode == 'pos':
            self._base_scatter_pos = self._ax_sc_pos.scatter(
                self._rt, self._mz, s=4, alpha=0.4,
                color='steelblue', zorder=2)
            # neg 用は他モードからの取得
            other_rt, other_mz, _ = self._get_other_mode_data()
            if other_rt is not None and len(other_rt) > 0:
                self._base_scatter_neg = self._ax_sc_neg.scatter(
                    other_rt, other_mz, s=4, alpha=0.4,
                    color='steelblue', zorder=2)
            else:
                self._base_scatter_neg = self._ax_sc_neg.scatter(
                    [], [], s=4, alpha=0.4, color='steelblue', zorder=2)
        else:
            self._base_scatter_neg = self._ax_sc_neg.scatter(
                self._rt, self._mz, s=4, alpha=0.4,
                color='steelblue', zorder=2)
            other_rt, other_mz, _ = self._get_other_mode_data()
            if other_rt is not None and len(other_rt) > 0:
                self._base_scatter_pos = self._ax_sc_pos.scatter(
                    other_rt, other_mz, s=4, alpha=0.4,
                    color='steelblue', zorder=2)
            else:
                self._base_scatter_pos = self._ax_sc_pos.scatter(
                    [], [], s=4, alpha=0.4, color='steelblue', zorder=2)

        # 後方互換: _base_scatter は active mode 側の base
        self._base_scatter = (self._base_scatter_pos
                              if self.fe.ion_mode == 'pos'
                              else self._base_scatter_neg)

        # タイトル / ラベル設定
        self._ax_sc_pos.set_xlabel("RT (min)")
        self._ax_sc_pos.set_ylabel("m/z")
        self._ax_sc_pos.set_title("RT vs m/z  (pos)")
        self._ax_sc_neg.set_xlabel("RT (min)")
        self._ax_sc_neg.set_ylabel("m/z")
        self._ax_sc_neg.set_title("RT vs m/z  (neg)")

        # Tooltip annotation
        # constrained_layout 由来の図サイズ変動を抑止
        # 両軸それぞれに annotation を持つ。各軸の data
        # coords で正しく解釈されるようにし、モード切替で位置がズレない。
        def _make_annot(ax):
            a = ax.annotate(
                "", xy=(0, 0), xytext=(10, 10),
                textcoords="offset points",
                bbox=dict(boxstyle="round,pad=0.3", fc="#EEEDFE",
                          ec="#AFA9EC", alpha=0.95),
                fontsize=8, zorder=10, annotation_clip=False)
            a.set_visible(False)
            a.set_in_layout(False)
            return a

        self._annot_pos = _make_annot(self._ax_sc_pos)
        self._annot_neg = _make_annot(self._ax_sc_neg)
        # 後方互換: 既存コードは self._annot を参照する
        self._annot = (self._annot_pos
                       if self.fe.ion_mode == 'pos'
                       else self._annot_neg)

        # イベント接続
        self._canvas.mpl_connect("motion_notify_event", self._on_hover)
        self._canvas.mpl_connect("scroll_event",        self._on_scroll)
        self._canvas.mpl_connect("button_press_event",  self._on_mouse_press)

        # 初期範囲を保存（ホームリセット用）
        # pos/neg 別々に保存(両軸の natural な m/z 範囲を保つ)
        self._canvas.draw()
        self._home_xlim_pos = self._ax_sc_pos.get_xlim()
        self._home_ylim_pos = self._ax_sc_pos.get_ylim()
        self._home_xlim_neg = self._ax_sc_neg.get_xlim()
        self._home_ylim_neg = self._ax_sc_neg.get_ylim()
        # 旧 attribute は active 軸のものとして互換性維持
        self._home_xlim = self._ax_sc.get_xlim()
        self._home_ylim = self._ax_sc.get_ylim()

        # SpanSelector pos / neg 両軸に取り付け
        # interactive=False に変更
        #   each drag fully replaces previous selection. No resize handles,
        #   no "drag-to-move" trap inside the existing red span.
        #   The most predictable UX for "select a new range every time".
        # button=1 (左クリックのみ) に制限。
        #   default では全ボタンに反応するため、右クリック(curation メニュー)
        #   後に SpanSelector の state が drag-in-progress のまま残り、
        #   続くマウス移動が誤って範囲選択として処理されてしまっていた。
        self._span_pos = SpanSelector(
            self._ax_sc_pos,
            lambda lo, hi: self._on_select(lo, hi, 'pos'),
            "horizontal",
            useblit=True,
            props=dict(alpha=0.3, facecolor="red"),
            interactive=False,
            button=1)
        self._span_neg = SpanSelector(
            self._ax_sc_neg,
            lambda lo, hi: self._on_select(lo, hi, 'neg'),
            "horizontal",
            useblit=True,
            props=dict(alpha=0.3, facecolor="red"),
            interactive=False,
            button=1)

        h_bar = QHBoxLayout()
        self._lbl_range = QLabel("Select a range on the scatter plot")
        h_bar.addWidget(self._lbl_range)
        # 範囲選択を消すボタン
        btn_clear_range = QPushButton("Clear range")
        btn_clear_range.setToolTip(
            "Clear the current RT range selection on both Pos and Neg plots.")
        btn_clear_range.clicked.connect(self._clear_range_selection)
        h_bar.addWidget(btn_clear_range)
        btn_add = QPushButton("Add Task …")
        btn_add.setToolTip(
            "Add a single task for the currently selected RT range.")
        btn_add.clicked.connect(self._open_add_task)
        h_bar.addWidget(btn_add)
        # パイプライン結果から Task 自動生成
        btn_auto = QPushButton("Auto-add Tasks…")
        btn_auto.setToolTip(
            "Automatically add one task per confirmed lipid class.\n"
            "RT range is taken from the matched & kept spots' obs_rt\n"
            "(min - padding to max + padding). Runs on both pos and neg.")
        btn_auto.clicked.connect(self._auto_add_tasks_from_pipeline)
        h_bar.addWidget(btn_auto)
        # "Set as Fix RT" 撤廃
        # (IS Filter の手動 RT 入力で代替可能)
        left_lay.addLayout(h_bar)

        # 右: マッチリスト
        right_w = QGroupBox("Matched features")
        # テーブルの内容で右パネル全体がダイアログを押し広げないようにする
        right_w.setSizePolicy(QSizePolicy.Ignored, QSizePolicy.Expanding)
        right_w.setMinimumWidth(0)
        right_lay = QVBoxLayout(right_w)
        h_split.addWidget(right_w)

        # スプリッター初期比率: 左60% / 右40%
        h_split.setSizes([660, 440])

        # ── Mode セレクタを撤廃し、両モード統合表示に変更 ──
        # Pos/Neg ラジオは UI から外したが、_apply_to_all_modes 内のラジオ
        # 同期処理など既存コードが参照するため widget は残す(layout には
        # 追加しないので非表示)。
        self._mode_radio_pos = QRadioButton("Pos")
        self._mode_radio_neg = QRadioButton("Neg")
        if self.fe.ion_mode == 'pos':
            self._mode_radio_pos.setChecked(True)
        else:
            self._mode_radio_neg.setChecked(True)
        self._mode_radio_pos.setVisible(False)
        self._mode_radio_neg.setVisible(False)
        # 既存ハンドラは残す(他コードからの programmatic setChecked にも反応)
        self._mode_radio_pos.toggled.connect(self._on_mode_radio_changed)
        self._mode_radio_neg.toggled.connect(self._on_mode_radio_changed)

        # Show classes ヘッダ行 (All / None: bulk マクロ)
        # Curated view は View filters > Restrict display
        # セクションに移動。ここでは All / None のみ。
        cls_header = QHBoxLayout()
        cls_header.addWidget(QLabel("Show classes:"))
        self._btn_all_classes  = QPushButton("All")
        self._btn_none_classes = QPushButton("None")
        for b in (self._btn_all_classes, self._btn_none_classes):
            b.setFixedWidth(44)
            b.setFixedHeight(20)
            b.setStyleSheet("font-size:11px; padding:1px 4px;")
        self._btn_all_classes.clicked.connect(self._select_all_classes)
        self._btn_none_classes.clicked.connect(self._select_no_classes)
        cls_header.addWidget(self._btn_all_classes)
        cls_header.addWidget(self._btn_none_classes)
        cls_header.addStretch()
        # 化合物テーブルを別ダイアログで開くボタン
        self._btn_show_compounds = QPushButton("Compounds Table…")
        self._btn_show_compounds.setToolTip(
            "Open the matched compounds table in a separate window\n"
            "(modeless: can be used alongside the Annotation Pipeline window).")
        self._btn_show_compounds.clicked.connect(
            self._open_match_table_dialog)
        cls_header.addWidget(self._btn_show_compounds)
        right_lay.addLayout(cls_header)

        self._filter_container = QWidget()
        self._filter_grid = QGridLayout(self._filter_container)
        self._filter_grid.setContentsMargins(4, 2, 4, 2)
        self._filter_grid.setSpacing(3)
        self._filter_grid.setHorizontalSpacing(8)
        self._filter_container.setSizePolicy(
            QSizePolicy.Preferred, QSizePolicy.Minimum)

        # QScrollArea に入れて高さを固定する
        # チェックボックス数が増えてもウィンドウが押し広がらないようにする
        self._filter_scroll = QScrollArea()
        self._filter_scroll.setWidget(self._filter_container)
        self._filter_scroll.setWidgetResizable(True)
        self._filter_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarAlwaysOff)
        self._filter_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarAsNeeded)
        self._filter_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        # 最小高さのみ保証（比率は right_lay の stretch で制御）
        self._filter_scroll.setMinimumHeight(100)
        right_lay.addWidget(self._filter_scroll, stretch=1)

        # テーブル + count label は別ダイアログ
        # (MatchTableDialog) で表示する。ここでは widget だけ生成し、
        # right_lay には追加しない。MatchTableDialog 初回オープン時に
        # ダイアログ側 layout に格納する。
        # Mode 列を追加(両モード統合表示)
        self._match_table = QTableWidget(0, 7)
        self._match_table.setHorizontalHeaderLabels(
            ["Mode", "Class", "Compound", "Adduct", "Δppm", "RT(min)",
             "Conflicts with"])
        self._match_table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.Stretch)
        self._match_table.horizontalHeader().setSectionResizeMode(
            6, QHeaderView.Stretch)
        # 別ダイアログで表示するため、size policy はダイアログ側で再調整
        self._match_table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._match_table.setMinimumWidth(0)
        self._match_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self._match_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._match_table.itemSelectionChanged.connect(self._on_table_select)
        # 右クリックで手動キュレーション
        self._match_table.setContextMenuPolicy(Qt.CustomContextMenu)
        self._match_table.customContextMenuRequested.connect(
            self._on_match_table_context_menu)

        self._lbl_match_count = QLabel("No library loaded")

        # MatchTableDialog インスタンス(初回 _open_match_table_dialog で生成)
        self._match_table_dialog = None

        self._canvas.draw_idle()

        # ── Utility タブを構築 ────────────────────────────────
        self._build_utility_tab()

        # ── 前回設定の復元 ────────────────────────────────────────────
        self._restore_settings(initial_settings)

        # ── 終了時にPreview設定とreserved_coordsをMainWindowへ通知 ─────
        self.finished.connect(self._on_finished)

        # ── 修正: ダイアログサイズ + 中央配置 ───────────────────
        try:
            screen = QApplication.primaryScreen()
            if screen is not None:
                geom = screen.availableGeometry()
                target_w = min(1700, max(1200, int(geom.width() * 0.92)))
                target_h = min(1100, max(800, int(geom.height() * 0.92)))
                self.resize(target_w, target_h)
                self.move(
                    geom.x() + (geom.width() - target_w) // 2,
                    geom.y() + (geom.height() - target_h) // 2,
                )
        except Exception as _e:
            log.warning(f"[PreviewRTDialog position] failed: {_e}")

    # ── Utility タブ構築 ────────────────────────────────────


    def _build_external_import_tab(self) -> QWidget:
        """Utility > External data import を構築する。

        本ソフトの既定は「① タブで raw mzML からピーク検出 + アライメント」。
        外部で作ったテーブルを読むのは拡張的な使い方なので、こちらに集めた。

          - MS-DIAL Per-sample peak list フォルダ
                サンプルごとのピークリスト。アライメントは本ソフトが行う。
                ① タブ Step 3 のパラメータ(host._align_params)がそのまま使われる。
          - MS-DIAL Alignment result(1 ファイル)
                検出もアライメントも外部で済んだ表。そのままアノテーションへ。
        """
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        lay = QVBoxLayout(content)
        lay.setContentsMargins(10, 10, 10, 10)
        lay.setSpacing(10)

        intro = QLabel(
            "Import a peak table produced elsewhere (extension).\n\n"
            # タブ名をベタ書きすると改名時に食い違うので参照にする
            f"Normally you run detection and alignment from raw mzML in the\n"
            f"'{TAB_LABELS['align']}' tab. Use this when you bring in data that\n"
            "was already processed in MS-DIAL or a similar tool.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color:#444; padding:4px;")
        lay.addWidget(intro)

        box = QGroupBox("Import")
        blay = QVBoxLayout(box)
        _b1 = QPushButton("Manage data files…")
        _b1.setMinimumHeight(30)
        _b1.setToolTip(
            "Add / remove / inspect input files.\n"
            "Input mode offers:\n"
            "  - Per-sample peak list folder (MS-DIAL)\n"
            "      One TXT per sample; this app does the alignment\n"
            "  - Alignment result (1 file)\n"
            "      Detection and alignment already done elsewhere")
        _b1.clicked.connect(self._open_manage_files_dialog)
        blay.addWidget(_b1)

        note = QLabel(
            "When you load a per-sample peak list, the alignment and isotope\n"
            "grouping parameters set in Step 3 of the ① tab are used as-is.\n"
            "A confirmation dialog shows those values before it runs.")
        note.setWordWrap(True)
        note.setStyleSheet("color:#666; padding:2px;")
        blay.addWidget(note)
        lay.addWidget(box)

        self._lbl_import_summary = QLabel("(no data loaded)")
        self._lbl_import_summary.setStyleSheet("color:#555; padding:4px;")
        self._lbl_import_summary.setWordWrap(True)
        lay.addWidget(self._lbl_import_summary)

        lay.addStretch()
        return tab

    def _refresh_import_summary(self):
        """External data import タブの読込状況ラベルを更新する。"""
        host = self._host
        if host is None or not hasattr(host, '_file_entries'):
            return
        entries = list(host._file_entries or [])
        if not entries:
            text = "(no data loaded)"
        else:
            text = "Loaded:\n" + "\n".join(
                f"  {fe.tag} [{fe.ion_mode}] {fe.path.name}"
                f"  ({len(fe.df)} features)" for fe in entries)
        try:
            self._lbl_import_summary.setText(text)
        except Exception:
            pass

    # ══════════════════════════════════════════════════════════════
    #  ① Peak Detection & Alignment タブ (fix19 / S3)
    #  ----
    #  Alignment_v0.5.27.py の Step 1〜4 構成を LipidZoner のタブとして
    #  作り直したもの。standalone 版は QMainWindow だったので、レイアウト
    #  だけ QWidget 用に組み替え、パラメータの意味と既定値は踏襲した。
    #
    #  fix19 時点の対応:
    #    Step 1  per-sample peak list TXT / raw mzML の 2 モード
    #            (Alignment は v0.5.2 で TXT モードを削除していたが、
    #             検証データセットが TXT なので復活させた)
    #    Step 2  raw mzML のピーク検出 → fix20 で実装(今は無効表示)
    #    Step 3  アライメント + 同位体グルーピング。パラメータは
    #            DEFAULT_ALIGN_PARAMS のキーと 1:1 で対応する
    #    Step 4  サマリ表示 / TXT・CSV 保存 / Annotation タブへ引き渡し
    #
    #  これにより「フォルダを選んだ瞬間に既定値でアライメントが走る」
    #  という fix16 以前の挙動は不要になった。
    # ══════════════════════════════════════════════════════════════

    # ── ① タブ: 極性ごとの状態 ───────────────────────────
    #  pos と neg を両方読み込むのが通常フローなので、Ion mode を
    #  切り替えても互いの状態が消えないようにする。
    #
    #  設計: 状態は self._al_state = {'pos': {...}, 'neg': {...}} に持ち、
    #  _al_loaded_files などは「今表示している極性」を指す property に
    #  する。こうすると ① タブの既存メソッド(37 個)は 1 行も書き換えず
    #  にそのまま動く。fix24 までの実装は属性への直接代入と
    #  item 代入(self._al_loaded_files[k] = ...)の両方をするので、
    #  getter が実体の list を返し、setter が差し替える形にしてある。
    #
    #  UI のパラメータ(Step 2 / Step 3)はウィジェットが 1 組しかないので、
    #  極性を切り替える瞬間に _al_collect_* で退避し _al_apply_* で
    #  復元する。この収集・適用の仕組みは fix21 の session 保存で作った
    #  ものをそのまま使い回している。

    # 極性ごとに 1 つずつ持つウィジェット。
    # self._al_<name> は「今対象にしている列」の widget を返す
    # プロパティになる(クラス定義の直後で生成)。
    _AL_PER_POL_WIDGETS = (
        # Step 1
        'file_table',
        # Step 2
        'sp_mass_ppm', 'btn_adv_det', 'sp_noise', 'sp_snr', 'sp_fwhm',
        'sp_min_trace', 'cmb_centroid', 'btn_detect', 'btn_chrom',
        'lbl_detect',
        # Step 3
        'sp_mz_tol', 'cb_ppm', 'sp_ppm', 'sp_rt_tol', 'sp_refine',
        'cb_rt_drift', 'rb_height', 'rb_area', 'int_group', 'cb_iso',
        'sp_iso_mz', 'sp_iso_rt', 'sp_max_iso', 'sp_iso_ratio',
        'sp_kC', 'sp_kO', 'btn_adv_align', 'btn_run', 'lbl_status',
        # Step 4
        'txt_summary', 'btn_txt', 'btn_csv', 'btn_plots',
    )

    def _al_with_pol(self, pol: str, fn, *args, **kw):
        """列 pol を対象にして fn を呼ぶ。

        2 列になったので「今どの列を操作しているか」を呼び出しごとに
        決める必要がある。_al_cur_pol を一時的に差し替えて実行し、
        必ず元に戻す。
        """
        old = getattr(self, '_al_cur_pol', 'pos')
        self._al_cur_pol = pol
        try:
            return fn(*args, **kw)
        finally:
            self._al_cur_pol = old

    @staticmethod
    def _al_blank_state() -> dict:
        return {'loaded_files': [], 'aligned_df': None, 'sample_names': [],
                'info': {}, 'detect_params_used': None,
                'params': None, 'detect': None, 'sent': False}

    def _al_init_state(self):
        self._al_state = {'pos': PreviewRTDialog._al_blank_state(),
                          'neg': PreviewRTDialog._al_blank_state()}
        # 極性ごとのウィジェット置き場。プロパティがここを引く。
        self._al_w = {'pos': {}, 'neg': {}}
        self._al_cur_pol = 'pos'
        self._al_switching = False

    def _al_ensure_state(self) -> dict:
        if not hasattr(self, '_al_state'):
            self._al_init_state()
        return self._al_state

    @property
    def _al_cur_state(self) -> dict:
        st = self._al_ensure_state()
        return st.get(getattr(self, '_al_cur_pol', 'pos'), st['pos'])

    @property
    def _al_loaded_files(self):
        return self._al_cur_state['loaded_files']

    @_al_loaded_files.setter
    def _al_loaded_files(self, v):
        self._al_cur_state['loaded_files'] = v

    @property
    def _al_aligned_df(self):
        return self._al_cur_state['aligned_df']

    @_al_aligned_df.setter
    def _al_aligned_df(self, v):
        self._al_cur_state['aligned_df'] = v

    @property
    def _al_sample_names(self):
        return self._al_cur_state['sample_names']

    @_al_sample_names.setter
    def _al_sample_names(self, v):
        self._al_cur_state['sample_names'] = v

    @property
    def _al_info(self):
        return self._al_cur_state['info']

    @_al_info.setter
    def _al_info(self, v):
        self._al_cur_state['info'] = v

    @property
    def _al_detect_params_used(self):
        return self._al_cur_state['detect_params_used']

    @_al_detect_params_used.setter
    def _al_detect_params_used(self, v):
        self._al_cur_state['detect_params_used'] = v

    # ── ① タブ: 極性の切り替え ───────────────────────────
    def _al_current_polarity(self) -> str:
        """ラジオボタンが指している極性。"""
        try:
            return 'neg' if self._al_rb_neg.isChecked() else 'pos'
        except Exception:
            return 'pos'

    def _al_stash_params(self, pol: str):
        """pol の列に出ている Step2 / Step3 の値を状態へ退避する。

        2 列になったので、**その列の widget** を読む。
        """
        try:
            st = self._al_ensure_state()[pol]
        except KeyError:
            return
        if not (getattr(self, '_al_w', None) or {}).get(pol):
            return
        try:
            st['params'] = self._al_with_pol(pol, self._al_collect_params)
            st['detect'] = self._al_with_pol(
                pol, self._al_collect_detect_params)
        except Exception as e:
            log.warning(f"stash params failed ({pol}): {e}")

    def _al_restore_params(self, pol: str) -> bool:
        """pol の状態に入っている Step2 / Step3 の値を pol の列へ戻す。

        まだ一度も設定していない極性は、その列に今出ている値を
        そのままその極性の値として確定させる。

        2 列になったので読み書きとも **その列の widget** に対して
        行う(fix33 までは「表示中の列」だった)。

        Returns
        -------
        bool : 引き継ぎが起きたら True
        """
        try:
            st = self._al_ensure_state()[pol]
        except KeyError:
            return False
        if not (getattr(self, '_al_w', None) or {}).get(pol):
            return False
        inherited = False
        try:
            if st.get('params') is None and st.get('detect') is None:
                st['params'] = self._al_with_pol(pol, self._al_collect_params)
                st['detect'] = self._al_with_pol(
                    pol, self._al_collect_detect_params)
                inherited = True
            else:
                if st.get('params'):
                    self._al_with_pol(
                        pol, self._al_apply_align_params, st['params'])
                if st.get('detect'):
                    self._al_with_pol(
                        pol, self._al_apply_detect_params, st['detect'])
        except Exception as e:
            log.warning(f"restore params failed ({pol}): {e}")
        return inherited

    def _al_set_polarity(self, pol: str):
        """これから操作する極性を決める。

        fix33 までは「表示を切り替える」操作だったが、2 列になったので
        表示は動かない。_al_cur_pol と隠しラジオを合わせるだけにする。
        _al_run_both がここを使って pos → neg と回す。
        """
        if pol not in ('pos', 'neg'):
            return
        old = getattr(self, '_al_cur_pol', 'pos')
        if old != pol:
            self._al_stash_params(old)
        self._al_switching = True
        try:
            (self._al_rb_neg if pol == 'neg'
             else self._al_rb_pos).setChecked(True)
        except Exception:
            pass
        finally:
            self._al_switching = False
        self._al_cur_pol = pol

    def _al_on_ion_mode_changed(self, *_a):
        """Ion mode ラジオが動いたとき。

        fix34 で 2 列になり、画面上の切り替えは無くなった。ラジオは
        _al_cur_pol を映す内部状態として残っているだけなので、ここでは
        UI を組み直さない(組み直すと両列とも同じ極性の値で埋まる)。
        """
        return

    def _al_sync_ui_to_state(self):
        """状態に合わせて画面を作り直す(fix34: 両方の列に効かせる)。

        fix33 までは「表示中の極性」だけを更新していたが、2 列になって
        両方が常に見えているので、片方だけ更新すると古い数字が残る。
        """
        for pol in ('pos', 'neg'):
            self._al_with_pol(pol, self._al_sync_one_column, pol)
        self._al_update_pol_summary()

    def _al_sync_one_column(self, pol: str):
        """1 つの列ぶんの Step1 / 2 / 3 / 4 を状態に合わせ直す。"""
        if not (getattr(self, '_al_w', None) or {}).get(pol):
            return
        try:
            self._al_refresh_file_table()
        except Exception:
            pass
        try:
            self._al_update_status()
        except Exception:
            pass
        has_align = self._al_cur_state['aligned_df'] is not None
        for name in ('_al_btn_txt', '_al_btn_csv', '_al_btn_plots'):
            b = getattr(self, name, None)
            if b is not None:
                b.setEnabled(has_align)
        try:
            if has_align:
                self._al_show_summary()
            else:
                self._al_txt_summary.clear()
        except Exception:
            pass
        try:
            dets = [df for _, _, df in self._al_loaded_files if df is not None]
            if dets:
                self._al_lbl_detect.setText(
                    f"✓ {len(dets)} samples, "
                    f"{sum(len(d) for d in dets)} peaks detected")
            elif self._al_loaded_files:
                self._al_lbl_detect.setText(
                    "Press [Detect Peaks] to run detection for this polarity.")
            else:
                self._al_lbl_detect.setText(
                    "Add mzML files, then press [Detect Peaks].")
        except Exception:
            pass

    def _al_update_pol_summary(self):
        """両極性の進捗を 1 行にまとめて出す(fix25 / fix34 で文言整理)。

        fix34 で 2 列になったので、これは「もう一方の様子」を伝える
        ためではなく、**Send to Annotation を押してよい状態か**を
        1 行で確かめるためのもの。
        """
        state = self._al_ensure_state()
        parts = []
        ready = 0
        for pol in ('pos', 'neg'):
            st = state[pol]
            nfile = len(st['loaded_files'])
            if not nfile:
                parts.append(f"{pol.upper()}: no files")
            elif st['aligned_df'] is not None:
                ready += 1
                parts.append(f"{pol.upper()}: {nfile} files → "
                             f"{len(st['aligned_df'])} features"
                             + ("  [sent]" if st.get('sent') else ""))
            elif any(df is not None for _, _, df in st['loaded_files']):
                parts.append(f"{pol.upper()}: {nfile} files, detected, "
                             f"not aligned")
            else:
                parts.append(f"{pol.upper()}: {nfile} files, not detected")
        # fix31 で Annotation は極性ごとに 1 エントリになったので、
        # 「両方そろってから Send」が通常フローになる。
        try:
            self._al_btn_send.setEnabled(ready > 0)
        except Exception:
            pass
        try:
            self._al_lbl_pol_summary.setText("        ".join(parts))
        except Exception:
            pass

    @staticmethod
    def _al_two_col(left, right):
        """左右 2 列に並べる。真ん中に区切り線を入れる。"""
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(14)
        h.addWidget(left, 1)
        sep = QFrame()
        sep.setFrameShape(QFrame.VLine)
        sep.setFrameShadow(QFrame.Sunken)
        h.addWidget(sep)
        h.addWidget(right, 1)
        return w

    def _al_col_header(self, pol: str) -> QLabel:
        lab = QLabel("Pos" if pol == 'pos' else "Neg")
        f = lab.font()
        f.setBold(True)
        lab.setFont(f)
        lab.setStyleSheet("color:#036;")
        return lab

    def _build_alignment_tab(self) -> QWidget:
        """① Detection タブを構築して返す。

        pos / neg を左右 2 列で同時に見せる。Ion mode の切り替えは
        廃止した。各 Step の枠の中が Pos | Neg の 2 列になっているので、
        同じパラメータが真横に並んで比べられる。
        """
        self._al_init_state()

        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        outer.addWidget(scroll)
        content = QWidget()
        scroll.setWidget(content)
        root = QVBoxLayout(content)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # ── 共有ウィジェット ──────────────────────────────────────
        # 入力モードのラジオは撤去した。このタブは raw mzML 専用で、
        #   外部で作った表(per-sample peak list / Alignment result)の
        #   取り込みは Utility > External data import に一本化した。

        # Ion mode ラジオは画面から消えたが、_al_cur_pol を映す
        #   内部状態として残す。_al_current_polarity() など 20 箇所が
        #   参照しているので、消すより隠すほうが安全(fix24 と同じ判断)。
        self._al_rb_pos = QRadioButton("Pos")
        self._al_rb_neg = QRadioButton("Neg")
        self._al_rb_pos.setChecked(True)
        self._al_ion_group = QButtonGroup(self)
        self._al_ion_group.addButton(self._al_rb_pos)
        self._al_ion_group.addButton(self._al_rb_neg)
        self._al_rb_pos.setVisible(False)
        self._al_rb_neg.setVisible(False)

        # "Input: Raw data (.mzML) — To load a table made elsewhere…"
        #   の行は撤去した。このタブが raw mzML 専用なのは
        #   Step 2 の見出しに書いてあるので、重ねて説明する必要がない。

        # ── 各極性の列を作る ──────────────────────────────────────
        cols = {}
        _keep = getattr(self, '_al_cur_pol', 'pos')
        for _pol in ('pos', 'neg'):
            self._al_cur_pol = _pol
            cols[_pol] = self._al_build_pol_columns(_pol)
        self._al_cur_pol = _keep

        # Step 1
        step1 = QGroupBox("Step 1: Data input")
        s1 = QVBoxLayout(step1)
        s1.addWidget(self._al_two_col(cols['pos']['c1'], cols['neg']['c1']))
        # 常時出ていた説明文は撤去。ファイルを追加したときの
        #   結果や警告を出すときだけ姿を見せる。
        self._al_lbl_ion_auto = QLabel("")
        self._al_lbl_ion_auto.setStyleSheet("color:#666;")
        self._al_lbl_ion_auto.setVisible(False)
        s1.addWidget(self._al_lbl_ion_auto)
        root.addWidget(step1)

        # Step 2
        self._al_step2 = QGroupBox("Step 2: Peak detection")
        s2 = QVBoxLayout(self._al_step2)
        s2.addWidget(self._al_two_col(cols['pos']['c2'], cols['neg']['c2']))
        self._al_lbl_pyopenms = QLabel("")
        self._al_lbl_pyopenms.setWordWrap(True)
        self._al_lbl_pyopenms.setStyleSheet("color:#555; padding:2px;")
        s2.addWidget(self._al_lbl_pyopenms)
        root.addWidget(self._al_step2)

        # Step 3
        step3 = QGroupBox("Step 3: Alignment + isotope grouping")
        s3 = QVBoxLayout(step3)
        s3.addWidget(self._al_two_col(cols['pos']['c3'], cols['neg']['c3']))
        root.addWidget(step3)

        # Step 4
        step4 = QGroupBox("Step 4: Results / export")
        s4 = QVBoxLayout(step4)
        s4.addWidget(self._al_two_col(cols['pos']['c4'], cols['neg']['c4']))
        root.addWidget(step4)

        # ── 両極性まとめて動かすボタン(共有)──────────────────────
        both_row = QHBoxLayout()
        both_row.addStretch()
        self._al_btn_run_both = QPushButton(
            "▶▶  Run both  (detect + align, pos and neg)")
        self._al_btn_run_both.setMinimumHeight(32)
        self._al_btn_run_both.setMinimumWidth(320)
        self._al_btn_run_both.setToolTip(
            "For every polarity that has files, run peak detection (Step 2)\n"
            "if it has not been done yet, then continue through\n"
            "alignment (Step 3).")
        _fb = self._al_btn_run_both.font()
        _fb.setBold(True)
        self._al_btn_run_both.setFont(_fb)
        self._al_btn_run_both.clicked.connect(self._al_run_both)
        both_row.addWidget(self._al_btn_run_both)
        self._al_btn_send = QPushButton("→  Send to Annotation")
        self._al_btn_send.setMinimumHeight(32)
        self._al_btn_send.setMinimumWidth(220)
        self._al_btn_send.setToolTip(
            "Register the alignment result as a PerSampleFileEntry in the\n"
            "Annotation tab. Only M+0 is annotated; isotopes are kept\n"
            "internally for quantification.")
        self._al_btn_send.setFont(_fb)
        self._al_btn_send.setEnabled(False)
        self._al_btn_send.clicked.connect(self._al_send_to_annotation)
        both_row.addWidget(self._al_btn_send)
        both_row.addStretch()
        root.addLayout(both_row)

        # 両極性の進捗を 1 行で(Send できる状態かの確認用)
        self._al_lbl_pol_summary = QLabel("")
        self._al_lbl_pol_summary.setStyleSheet(
            "color:#046; background:#eef4f8; padding:4px; border-radius:3px;")
        root.addWidget(self._al_lbl_pol_summary)
        root.addStretch()

        self._al_step2.setVisible(True)
        for _pol in ('pos', 'neg'):
            self._al_with_pol(
                _pol, self._al_on_iso_toggled,
                self._al_w[_pol]['cb_iso'].isChecked())
        return tab

    def _al_build_pol_columns(self, pol: str) -> dict:
        """1 つの極性ぶんの Step 1〜4 の中身を作る。

        呼び出し前に self._al_cur_pol = pol にしておくこと。
        self._al_xxx = ... の代入は、プロパティ経由で
        self._al_w[pol]['xxx'] に入る。
        """
        D = DEFAULT_DETECT_PARAMS
        P = DEFAULT_ALIGN_PARAMS
        W = lambda f, *a, **k: self._al_with_pol(pol, f, *a, **k)

        # ══ Step 1 列 ═══════════════════════════════════════════
        c1 = QWidget()
        v1 = QVBoxLayout(c1)
        v1.setContentsMargins(0, 0, 0, 0)
        v1.addWidget(self._al_col_header(pol))
        btn_row = QHBoxLayout()
        for txt, slot, tip in (
            ("Add file(s)…", self._al_add_files,
             "Select .mzML files to add (multiple allowed)"),
            ("Add folder…", self._al_add_folder,
             "Add every .mzML inside a folder"),
            ("Remove selected", self._al_remove_selected,
             "Remove the selected rows"),
        ):
            b = QPushButton(txt)
            b.setToolTip(tip)
            b.clicked.connect(lambda _c=False, f=slot: W(f))
            btn_row.addWidget(b)
        b_clear = QPushButton("Clear all")
        b_clear.setToolTip("Discard every loaded file of this polarity")
        b_clear.clicked.connect(
            lambda _c=False: W(self._al_clear_all, polarity='current'))
        btn_row.addWidget(b_clear)
        btn_row.addStretch()
        v1.addLayout(btn_row)

        self._al_file_table = QTableWidget(0, 4)
        t = self._al_file_table
        t.setHorizontalHeaderLabels(["#", "Sample name", "Peaks", "Path"])
        _h = t.horizontalHeader()
        # 列を圧縮する。Path は末尾だけ見えれば pos/neg の取り違えは
        # 分かるので狭くし、全文はツールチップに出す。
        _h.setSectionResizeMode(0, QHeaderView.Fixed)
        t.setColumnWidth(0, 30)
        _h.setSectionResizeMode(1, QHeaderView.Stretch)
        # 検出後は "4820 / 780 / 166 (M+0/M+1/M+2)" が入るので
        # 66 px では読めなかった。実測 219-244 px 必要。
        _h.setSectionResizeMode(2, QHeaderView.Fixed)
        t.setColumnWidth(2, 230)
        _h.setSectionResizeMode(3, QHeaderView.Fixed)
        t.setColumnWidth(3, 130)
        t.setWordWrap(False)
        t.setTextElideMode(Qt.ElideMiddle)
        t.setSelectionBehavior(QTableWidget.SelectRows)
        t.setEditTriggers(QTableWidget.NoEditTriggers)
        t.setMinimumHeight(150)
        t.verticalHeader().setVisible(False)
        v1.addWidget(t)

        # ══ Step 2 列 ═══════════════════════════════════════════
        c2 = QWidget()
        v2 = QVBoxLayout(c2)
        v2.setContentsMargins(0, 0, 0, 0)
        v2.addWidget(self._al_col_header(pol))
        s2f = QFormLayout()
        s2f.setSpacing(6)
        v2.addLayout(s2f)

        self._al_sp_mass_ppm = QDoubleSpinBox()
        self._al_sp_mass_ppm.setDecimals(1)
        self._al_sp_mass_ppm.setRange(1.0, 200.0)
        self._al_sp_mass_ppm.setValue(D['mass_error_ppm'])
        self._al_sp_mass_ppm.setSuffix(" ppm")
        self._al_sp_mass_ppm.setMaximumWidth(180)
        self._al_sp_mass_ppm.setToolTip(
            "m/z clustering tolerance for MassTraceDetection (ppm).\n"
            "[Advanced settings...] can estimate it from the data.")
        self._al_btn_adv_det = QPushButton("Advanced settings...")
        self._al_btn_adv_det.setToolTip(
            "Estimate mass accuracy / noise / S/N / peak width from the data.")
        self._al_btn_adv_det.clicked.connect(
            lambda _c=False: W(self._al_advanced_detection))
        _mrow = QHBoxLayout()
        _mrow.setContentsMargins(0, 0, 0, 0)
        _mrow.addWidget(self._al_sp_mass_ppm)
        _mrow.addStretch()
        _mrow.addWidget(self._al_btn_adv_det)
        _mw_ = QWidget()
        _mw_.setLayout(_mrow)
        s2f.addRow("Mass trace m/z accuracy:", _mw_)

        self._al_sp_noise = QDoubleSpinBox()
        self._al_sp_noise.setDecimals(0)
        self._al_sp_noise.setRange(0.0, 1e9)
        self._al_sp_noise.setSingleStep(100.0)
        self._al_sp_noise.setValue(D['noise_threshold_int'])
        self._al_sp_noise.setMaximumWidth(180)
        s2f.addRow("Noise threshold (intensity):", self._al_sp_noise)

        self._al_sp_snr = QDoubleSpinBox()
        self._al_sp_snr.setDecimals(1)
        self._al_sp_snr.setRange(0.0, 100.0)
        self._al_sp_snr.setValue(D['chrom_peak_snr'])
        self._al_sp_snr.setMaximumWidth(180)
        s2f.addRow("Chromatographic peak S/N:", self._al_sp_snr)

        self._al_sp_fwhm = QDoubleSpinBox()
        self._al_sp_fwhm.setDecimals(1)
        self._al_sp_fwhm.setRange(0.5, 60.0)
        self._al_sp_fwhm.setValue(D['chrom_fwhm'])
        self._al_sp_fwhm.setSuffix(" sec")
        self._al_sp_fwhm.setMaximumWidth(180)
        self._al_sp_fwhm.setToolTip(
            "Expected chromatographic peak width (FWHM, seconds).\n"
            "Min trace length follows this at 0.6x automatically.")
        s2f.addRow("Expected peak FWHM:", self._al_sp_fwhm)

        self._al_sp_min_trace = QDoubleSpinBox()
        self._al_sp_min_trace.setDecimals(1)
        self._al_sp_min_trace.setRange(0.5, 30.0)
        self._al_sp_min_trace.setValue(3.0)
        self._al_sp_min_trace.setSuffix(" sec")
        self._al_sp_min_trace.setMaximumWidth(180)
        self._al_sp_min_trace.setToolTip(
            "Discard mass traces shorter than this. Always keep it below the\n"
            "real peak width. The pyOpenMS default of 5.0 sec targets LC and\n"
            "discards the ~4 sec peaks of SFC (on 2026-06-03 this dropped\n"
            "coverage to 63%).")
        s2f.addRow("Min trace length:", self._al_sp_min_trace)
        self._al_sp_fwhm.valueChanged.connect(
            lambda v: W(self._al_sync_min_trace_length, v))
        W(self._al_sync_min_trace_length, self._al_sp_fwhm.value())

        self._al_cmb_centroid = QComboBox()
        self._al_cmb_centroid.addItems(["auto", "on", "off"])
        self._al_cmb_centroid.setCurrentText(D['centroid'])
        self._al_cmb_centroid.setMaximumWidth(180)
        self._al_cmb_centroid.setToolTip(
            "Whether to centroid profile-mode spectra before detection\n"
            "(PeakPickerHiRes). Feeding profile data to MassTrace detection\n"
            "without centroiding splits one peak into several.\n"
            "auto: only for profile data / on: always / off: never")
        s2f.addRow("Centroid mode:", self._al_cmb_centroid)

        det_row = QHBoxLayout()
        self._al_btn_detect = QPushButton("Detect Peaks")
        self._al_btn_detect.setMinimumHeight(28)
        _fd = self._al_btn_detect.font()
        _fd.setBold(True)
        self._al_btn_detect.setFont(_fd)
        self._al_btn_detect.clicked.connect(
            lambda _c=False: W(self._al_detect_peaks))
        det_row.addWidget(self._al_btn_detect)
        self._al_btn_chrom = QPushButton("View Chromatogram…")
        self._al_btn_chrom.setToolTip(
            "Open the TIC / EIC viewer in a separate window.")
        self._al_btn_chrom.clicked.connect(
            lambda _c=False: W(self._al_view_chromatogram))
        det_row.addWidget(self._al_btn_chrom)
        det_row.addStretch()
        v2.addLayout(det_row)
        self._al_lbl_detect = QLabel(
            "Add mzML files, then press [Detect Peaks].")
        self._al_lbl_detect.setStyleSheet("color:#555; padding:2px;")
        self._al_lbl_detect.setWordWrap(True)
        v2.addWidget(self._al_lbl_detect)

        # ══ Step 3 列 ═══════════════════════════════════════════
        c3 = QWidget()
        v3 = QVBoxLayout(c3)
        v3.setContentsMargins(0, 0, 0, 0)
        v3.addWidget(self._al_col_header(pol))
        form = QFormLayout()
        form.setSpacing(6)
        v3.addLayout(form)

        self._al_sp_mz_tol = QDoubleSpinBox()
        self._al_sp_mz_tol.setDecimals(4)
        self._al_sp_mz_tol.setRange(0.0001, 0.1)
        self._al_sp_mz_tol.setSingleStep(0.001)
        self._al_sp_mz_tol.setValue(P['align_mz_tol'])
        self._al_sp_mz_tol.setSuffix(" Da")
        self._al_sp_mz_tol.setMaximumWidth(180)
        self._al_sp_mz_tol.setToolTip(
            "m/z tolerance for treating peaks from different samples as the\n"
            "same feature. With ppm enabled this acts as the floor on the\n"
            "low-mass side.")
        form.addRow("m/z tolerance:", self._al_sp_mz_tol)

        self._al_cb_ppm = QCheckBox("Use ppm (floor = m/z tolerance above)")
        self._al_sp_ppm = QDoubleSpinBox()
        self._al_sp_ppm.setDecimals(1)
        self._al_sp_ppm.setRange(1.0, 100.0)
        self._al_sp_ppm.setValue(15.0)
        self._al_sp_ppm.setSuffix(" ppm")
        self._al_sp_ppm.setMaximumWidth(130)
        self._al_sp_ppm.setEnabled(bool(P['align_mz_ppm']))
        self._al_cb_ppm.setChecked(P['align_mz_ppm'] is not None)
        if P['align_mz_ppm']:
            self._al_sp_ppm.setValue(float(P['align_mz_ppm']))
        self._al_cb_ppm.toggled.connect(self._al_sp_ppm.setEnabled)
        _ppm_row = QHBoxLayout()
        _ppm_row.setContentsMargins(0, 0, 0, 0)
        _ppm_row.addWidget(self._al_cb_ppm)
        _ppm_row.addWidget(self._al_sp_ppm)
        _ppm_row.addStretch()
        _ppm_w = QWidget()
        _ppm_w.setLayout(_ppm_row)
        form.addRow("m/z tolerance mode:", _ppm_w)

        self._al_sp_rt_tol = QDoubleSpinBox()
        self._al_sp_rt_tol.setDecimals(3)
        self._al_sp_rt_tol.setRange(0.001, 5.0)
        self._al_sp_rt_tol.setSingleStep(0.01)
        self._al_sp_rt_tol.setValue(P['align_rt_tol'])
        self._al_sp_rt_tol.setSuffix(" min")
        self._al_sp_rt_tol.setMaximumWidth(180)
        form.addRow("RT tolerance:", self._al_sp_rt_tol)

        self._al_sp_refine = QSpinBox()
        self._al_sp_refine.setRange(1, 10)
        self._al_sp_refine.setValue(P['align_n_refine'])
        self._al_sp_refine.setMaximumWidth(180)
        self._al_sp_refine.setToolTip(
            "Number of centroid refinement iterations. Strong peaks seed the\n"
            "centroids, then every peak is reassigned to its nearest one.")
        form.addRow("Centroid refine iterations:", self._al_sp_refine)

        self._al_cb_rt_drift = QCheckBox("Correct cross-sample RT drift")
        self._al_cb_rt_drift.setChecked(bool(P['rt_drift_correct']))
        self._al_cb_rt_drift.setToolTip(
            "Use high-intensity peaks found in every sample to median-shift\n"
            "each sample's RT onto the consensus. Needs 2 or more samples.")
        form.addRow("RT drift correction:", self._al_cb_rt_drift)

        self._al_rb_height = QRadioButton("Height")
        self._al_rb_area = QRadioButton("Area")
        self._al_rb_height.setChecked(P['height_col'] == 'Height')
        self._al_rb_area.setChecked(P['height_col'] == 'Area')
        # QButtonGroup は必ず極性ごとに作る。1 個を共有すると
        # pos の Height と neg の Height が排他になってしまう(fix26 で
        # Candidate Picker が踏んだのと同じ罠)。
        self._al_int_group = QButtonGroup(self)
        self._al_int_group.addButton(self._al_rb_height)
        self._al_int_group.addButton(self._al_rb_area)
        _int_row = QHBoxLayout()
        _int_row.setContentsMargins(0, 0, 0, 0)
        _int_row.addWidget(self._al_rb_height)
        _int_row.addWidget(self._al_rb_area)
        _int_row.addStretch()
        _int_w = QWidget()
        _int_w.setLayout(_int_row)
        form.addRow("Intensity column:", _int_w)

        _sep_gf = QFrame()
        _sep_gf.setFrameShape(QFrame.HLine)
        _sep_gf.setStyleSheet("color:#ddd;")
        form.addRow(_sep_gf)

        # ギャップフィリング。MS-DIAL と同じく既定 ON。
        self._al_cb_gapfill = QCheckBox("Fill missing values from raw data")
        self._al_cb_gapfill.setChecked(bool(P['gap_fill']))
        self._al_cb_gapfill.setToolTip(
            "After alignment, re-integrate the raw mzML for every sample where\n"
            "a feature was not detected, so the cell is a measured value\n"
            "instead of a zero. Integration is forced: the window is used even\n"
            "when there is no local maximum.\n"
            "\n"
            "The window is the average peak width of the samples that DID\n"
            "detect the feature, times the factor below. This follows MS-DIAL's\n"
            "'Gap filling by compulsion', which is also on by default there.\n"
            "\n"
            "Needs the raw .mzML files. Filled cells are flagged, and the\n"
            "quantification table marks species and IS that used one.")
        self._al_cb_gapfill.toggled.connect(
            lambda on: W(self._al_on_gapfill_toggled, on))
        form.addRow("Gap filling:", self._al_cb_gapfill)

        self._al_sp_gf_factor = QDoubleSpinBox()
        self._al_sp_gf_factor.setDecimals(2)
        self._al_sp_gf_factor.setRange(0.10, 5.00)
        self._al_sp_gf_factor.setSingleStep(0.1)
        self._al_sp_gf_factor.setValue(float(P['gap_fill_width_factor']))
        self._al_sp_gf_factor.setSuffix(" x peak width")
        self._al_sp_gf_factor.setMaximumWidth(180)
        self._al_sp_gf_factor.setToolTip(
            "Half width of the integration window, as a multiple of the\n"
            "average peak width (FWHM) of the samples that detected the\n"
            "feature. 1.0 means +/- one peak width.\n"
            "\n"
            "Do not raise this much. On the SFC method the peak width is\n"
            "about 2.2 s (0.037 min), while DG sn-regioisomers sit 0.052 min\n"
            "apart: at +/- 0.05 min the fill starts picking up the neighbour\n"
            "(7 of 24 doublets were contaminated in testing, versus 0 at\n"
            "+/- 0.035 min).")
        form.addRow("  Window half width:", self._al_sp_gf_factor)

        self._al_sp_gf_ppm = QDoubleSpinBox()
        self._al_sp_gf_ppm.setDecimals(1)
        self._al_sp_gf_ppm.setRange(1.0, 100.0)
        self._al_sp_gf_ppm.setSingleStep(1.0)
        self._al_sp_gf_ppm.setValue(float(P['gap_fill_mz_ppm']))
        self._al_sp_gf_ppm.setSuffix(" ppm")
        self._al_sp_gf_ppm.setMaximumWidth(180)
        self._al_sp_gf_ppm.setToolTip(
            "m/z window for the fill, in ppm. Keep it equal to the annotation\n"
            "tolerance so a filled cell means the same thing as a matched one.")
        form.addRow("  Fill m/z tolerance:", self._al_sp_gf_ppm)

        _sep = QFrame()
        _sep.setFrameShape(QFrame.HLine)
        _sep.setStyleSheet("color:#ddd;")
        form.addRow(_sep)

        # 画面には手法名を出さない。中身の説明はツールチップに残す。
        self._al_cb_iso = QCheckBox("Enable isotope grouping")
        self._al_cb_iso.setChecked(bool(P['iso_enabled']))
        self._al_cb_iso.setToolTip(
            "Use the theoretical isotope ratios of a lipid-specific averagine\n"
            "to tell pure isotopes from real species (Alignment v0.5.23).\n"
            "Turning it off makes every feature M+0.")
        self._al_cb_iso.toggled.connect(
            lambda on: W(self._al_on_iso_toggled, on))
        form.addRow("Isotope grouping:", self._al_cb_iso)

        self._al_sp_iso_mz = QDoubleSpinBox()
        self._al_sp_iso_mz.setDecimals(4)
        self._al_sp_iso_mz.setRange(0.0001, 0.1)
        self._al_sp_iso_mz.setSingleStep(0.001)
        self._al_sp_iso_mz.setValue(P['iso_mz_tol'])
        self._al_sp_iso_mz.setSuffix(" Da")
        self._al_sp_iso_mz.setMaximumWidth(180)
        self._al_sp_iso_mz.setToolTip(
            "m/z matching for isotopes. Can differ from the alignment value\n"
            "(isotope search looks around the theoretical Δm/z = k × 1.003355).")
        form.addRow("  Isotope m/z tolerance:", self._al_sp_iso_mz)

        self._al_sp_iso_rt = QDoubleSpinBox()
        self._al_sp_iso_rt.setDecimals(3)
        self._al_sp_iso_rt.setRange(0.001, 5.0)
        self._al_sp_iso_rt.setSingleStep(0.01)
        self._al_sp_iso_rt.setValue(P['iso_rt_tol'])
        self._al_sp_iso_rt.setSuffix(" min")
        self._al_sp_iso_rt.setMaximumWidth(180)
        form.addRow("  Isotope RT tolerance:", self._al_sp_iso_rt)

        self._al_sp_max_iso = QSpinBox()
        self._al_sp_max_iso.setRange(1, 6)
        self._al_sp_max_iso.setValue(P['iso_max_iso'])
        self._al_sp_max_iso.setMaximumWidth(180)
        self._al_sp_max_iso.setToolTip(
            "Highest isotope position to search (3 = up to M+3)")
        form.addRow("  Max isotope position:", self._al_sp_max_iso)

        self._al_sp_iso_ratio = QDoubleSpinBox()
        self._al_sp_iso_ratio.setDecimals(2)
        self._al_sp_iso_ratio.setRange(0.05, 2.0)
        self._al_sp_iso_ratio.setSingleStep(0.05)
        self._al_sp_iso_ratio.setValue(P['iso_ratio_tol'])
        self._al_sp_iso_ratio.setMaximumWidth(180)
        self._al_sp_iso_ratio.setToolTip(
            "Treat a peak as a pure isotope when the observed ratio <=\n"
            "theoretical ratio × (1 + this value). Lowering it keeps\n"
            "borderline peaks as real species; raising it removes isotopes\n"
            "more aggressively. Default 0.5.")
        form.addRow("  Isotope ratio tolerance:", self._al_sp_iso_ratio)

        _kk_row = QHBoxLayout()
        _kk_row.setContentsMargins(0, 0, 0, 0)
        self._al_sp_kC = QDoubleSpinBox()
        self._al_sp_kC.setDecimals(4)
        self._al_sp_kC.setRange(0.001, 0.2)
        self._al_sp_kC.setSingleStep(0.001)
        self._al_sp_kC.setValue(P['iso_kC'])
        self._al_sp_kC.setMaximumWidth(120)
        self._al_sp_kC.setPrefix("kC ")
        self._al_sp_kO = QDoubleSpinBox()
        self._al_sp_kO.setDecimals(4)
        self._al_sp_kO.setRange(0.001, 0.2)
        self._al_sp_kO.setSingleStep(0.001)
        self._al_sp_kO.setValue(P['iso_kO'])
        self._al_sp_kO.setMaximumWidth(120)
        self._al_sp_kO.setPrefix("kO ")
        for _w in (self._al_sp_kC, self._al_sp_kO):
            _w.setToolTip(
                "Averagine element coefficients (C / O atoms per unit mass).\n"
                "The defaults were derived from this project's curated lipid\n"
                "library (795 lipids, 11 classes). Normally left unchanged.")
            _kk_row.addWidget(_w)
        _kk_row.addStretch()
        _kk_w = QWidget()
        _kk_w.setLayout(_kk_row)
        form.addRow("  Averagine coefficients:", _kk_w)

        run_row = QHBoxLayout()
        self._al_btn_run = QPushButton("▶  Run alignment")
        self._al_btn_run.setMinimumHeight(30)
        self._al_btn_run.setToolTip(
            "Run alignment for this polarity.")
        _f = self._al_btn_run.font()
        _f.setBold(True)
        self._al_btn_run.setFont(_f)
        self._al_btn_run.clicked.connect(
            lambda _c=False: W(self._al_run_alignment))
        run_row.addWidget(self._al_btn_run)
        _btn_reset = QPushButton("Reset to defaults")
        _btn_reset.clicked.connect(lambda _c=False: W(self._al_reset_params))
        run_row.addWidget(_btn_reset)
        self._al_btn_adv_align = QPushButton("Advanced settings...")
        self._al_btn_adv_align.setToolTip(
            "Estimate the m/z / RT tolerances from the cross-sample scatter\n"
            "of the loaded data.")
        self._al_btn_adv_align.clicked.connect(
            lambda _c=False: W(self._al_advanced_alignment))
        run_row.addWidget(self._al_btn_adv_align)
        run_row.addStretch()
        v3.addLayout(run_row)
        self._al_lbl_status = QLabel("Status: ready (no files loaded)")
        self._al_lbl_status.setStyleSheet("color:#555; padding:2px;")
        self._al_lbl_status.setWordWrap(True)
        v3.addWidget(self._al_lbl_status)

        # ══ Step 4 列 ═══════════════════════════════════════════
        c4 = QWidget()
        v4 = QVBoxLayout(c4)
        v4.setContentsMargins(0, 0, 0, 0)
        v4.addWidget(self._al_col_header(pol))
        self._al_txt_summary = QTextEdit()
        self._al_txt_summary.setReadOnly(True)
        self._al_txt_summary.setMinimumHeight(170)
        self._al_txt_summary.setPlaceholderText(
            "Press Run alignment and the result summary appears here.")
        v4.addWidget(self._al_txt_summary)
        save_row = QHBoxLayout()
        self._al_btn_txt = QPushButton("Save aligned TXT…")
        self._al_btn_txt.setToolTip(
            "Save M+0 features as a TXT compatible with the MS-DIAL Alignment\n"
            "result format. To load one back in, use\n"
            "Utility > External data import.")
        self._al_btn_txt.clicked.connect(lambda _c=False: W(self._al_save_txt))
        self._al_btn_txt.setEnabled(False)
        save_row.addWidget(self._al_btn_txt)
        self._al_btn_csv = QPushButton("Save extended CSV…")
        self._al_btn_csv.setToolTip(
            "Save every feature, including M+1 / M+2 / M+3 and the excess_*\n"
            "columns, as CSV (isotope-resolved feature table).")
        self._al_btn_csv.clicked.connect(lambda _c=False: W(self._al_save_csv))
        self._al_btn_csv.setEnabled(False)
        save_row.addWidget(self._al_btn_csv)
        self._al_btn_plots = QPushButton("Show Plots…")
        self._al_btn_plots.setToolTip(
            "Open the RT×m/z scatter, per-sample intensities and EIC (raw mode\n"
            "only) in a separate window. Clicking a point selects that feature.")
        self._al_btn_plots.clicked.connect(
            lambda _c=False: W(self._al_show_plots))
        self._al_btn_plots.setEnabled(False)
        save_row.addWidget(self._al_btn_plots)
        save_row.addStretch()
        v4.addLayout(save_row)

        return {'c1': c1, 'c2': c2, 'c3': c3, 'c4': c4}

    # ── ① タブ Step 2: raw mzML ピーク検出 (fix20 / S4) ──────────
    def _al_sync_min_trace_length(self, fwhm):
        """Min trace length を FWHM の 0.6 倍に追随させる。

        MassTraceDetection の min_trace_length は実ピーク幅より必ず下に
        なければならず、上回ると実ピークごと捨てられる。FWHM を変えたときに
        古い値が取り残されないよう連動させる(Alignment v0.5.18 と同じ)。
        """
        try:
            self._al_sp_min_trace.setValue(round(0.6 * float(fwhm), 1))
        except Exception:
            pass

    def _al_check_pyopenms_ui(self) -> bool:
        """pyOpenMS の可用性を確認し、Step 2 のラベルに状態を出す。"""
        ok, ver, err = check_pyopenms()
        if ok:
            self._al_lbl_pyopenms.setText(
                f"Using pyOpenMS {ver} through the worker process.")
            self._al_lbl_pyopenms.setStyleSheet("color:#070; padding:2px;")
        else:
            self._al_lbl_pyopenms.setText(
                "⚠ pyOpenMS is unavailable, so raw mzML peak detection is off.\n"
                f"{(err or '')[:300]}\n"
                f"(Place {ALIGN_WORKER_NAME} in the same folder as LipidZoner and "
                "run `pip install pyopenms`.)")
            self._al_lbl_pyopenms.setStyleSheet("color:#b00; padding:2px;")
        try:
            self._al_btn_detect.setEnabled(bool(ok))
        except Exception:
            pass
        return bool(ok)

    def _al_auto_set_ion_mode(self):
        """表示中の極性のファイル群と Ion mode の整合をラベルに出す。

        fix24 まではここで polarity を判定して Ion mode を強制設定し、
        混在なら警告していた。fix25 では _al_add_paths がファイル単位で
        pos / neg に振り分けるので、混在はもう異常ではない。
        この関数は表示の更新だけを行う(呼び出し互換のため残してある)。
        """
        if not hasattr(self, '_al_polarity_cache'):
            self._al_polarity_cache = {}
        cur = getattr(self, '_al_cur_pol', 'pos')
        pols = []
        has_mzml = False
        for _name, p, _df in self._al_loaded_files:
            if not str(p).lower().endswith(".mzml"):
                continue
            has_mzml = True
            key = str(p)
            pol = self._al_polarity_cache.get(key)
            if pol is None:
                try:
                    pol = detect_mzml_polarity(p)
                except Exception:
                    pol = None
                self._al_polarity_cache[key] = pol
            if pol:
                pols.append(pol)
        uniq = set(pols)
        if not has_mzml:
            # 説明文は出さない。何も言うことが無いときは隠す。
            self._al_lbl_ion_auto.setText("")
            self._al_lbl_ion_auto.setVisible(False)
        elif not uniq:
            self._al_lbl_ion_auto.setText("⚠ Could not read the polarity from the mzML")
            self._al_lbl_ion_auto.setStyleSheet("color:#b00;")
        elif uniq == {cur}:
            self._al_lbl_ion_auto.setText(
                f"{cur.upper()}  ({len(pols)} files / read from the mzML)")
            self._al_lbl_ion_auto.setStyleSheet("color:#070;")
        else:
            self._al_lbl_ion_auto.setText(
                f"⚠ Showing {cur.upper()}, but the loaded files are "
                f"{'/'.join(sorted(uniq)).upper()}")
            self._al_lbl_ion_auto.setStyleSheet("color:#b00;")
        # 中身があるときだけ見せる
        self._al_lbl_ion_auto.setVisible(
            bool(self._al_lbl_ion_auto.text()))
        self._al_update_pol_summary()

    # ── 推定値は提案として溜め、[Apply] で初めて反映する ──
    #   fix36 までは推定ボタンが即座に入力欄を書き換えていたので、
    #   結果を見て「採らない」と判断しても元の値が失われていた。
    _AL_PROPOSAL_LABELS = {
        'sp_mass_ppm':  ('Mass trace m/z accuracy', ' ppm'),
        'sp_noise':     ('Noise threshold', ''),
        'sp_snr':       ('Chromatographic peak S/N', ''),
        'sp_fwhm':      ('Expected peak FWHM', ' sec'),
        'sp_mz_tol':    ('m/z tolerance', ' Da'),
        'sp_rt_tol':    ('RT tolerance', ' min'),
    }

    def _al_clear_proposal(self):
        self._al_proposal = {}

    def _al_propose(self, attr: str, value: float) -> float:
        """推定値を「提案」として記録する(まだ入力欄には入れない)。

        戻り値は、その入力欄に実際に入る値(range / decimals で丸めた後)。
        報告文に「提案値」として出すために、丸めた結果を返している。
        """
        if not hasattr(self, '_al_proposal'):
            self._al_proposal = {}
        w = getattr(self, '_al_' + attr)
        cur = w.value()
        # 実際に入る値を知るため一度入れて、すぐ元に戻す
        w.blockSignals(True)
        try:
            w.setValue(value)
            eff = w.value()
        finally:
            w.setValue(cur)
            w.blockSignals(False)
        self._al_proposal[attr] = (cur, eff)
        return eff

    def _al_proposal_summary(self) -> str:
        """溜まっている提案を「現在 → 提案」で 1 行ずつ返す。"""
        p = getattr(self, '_al_proposal', None) or {}
        if not p:
            return ""
        out = []
        for attr, (cur, new) in p.items():
            name, unit = self._AL_PROPOSAL_LABELS.get(attr, (attr, ''))
            mark = '' if cur == new else '  ← changes'
            out.append(f"  {name:28s} {cur:g}{unit}  →  {new:g}{unit}{mark}")
        return "\n".join(out)

    def _al_apply_proposal(self) -> str:
        """溜まっている提案を実際に入力欄へ入れる。"""
        p = getattr(self, '_al_proposal', None) or {}
        if not p:
            return "Nothing to apply."
        done = []
        for attr, (cur, new) in p.items():
            try:
                getattr(self, '_al_' + attr).setValue(new)
                name, unit = self._AL_PROPOSAL_LABELS.get(attr, (attr, ''))
                done.append(f"  {name}: {cur:g}{unit} -> {new:g}{unit}")
            except Exception as e:
                done.append(f"  {attr}: failed ({e})")
        self._al_clear_proposal()
        return "Applied:\n" + "\n".join(done)

    def _al_advanced_detection(self):
        self._al_clear_proposal()
        AdvancedSettingsDialog(self, section="detection", parent=self).exec()

    def _al_advanced_alignment(self):
        self._al_clear_proposal()
        AdvancedSettingsDialog(self, section="alignment", parent=self).exec()

    def _al_detect_peaks(self, _checked=False, *, quiet: bool = False):
        """読み込んだ mzML 全部にピーク検出をかける(worker を並列実行)。

        quiet=True なら完了ダイアログを出さない([Run both] 用)。
        失敗・中断の通知は quiet でも出す。
        """
        if not self._al_loaded_files:
            QMessageBox.information(self, "No files",
                                    "Add .mzML files in Step 1.")
            return
        ok, _ver, err = check_pyopenms()
        if not ok:
            QMessageBox.warning(
                self, "pyOpenMS unavailable",
                "Raw mzML detection needs pyOpenMS (through the worker).\n\n"
                + (err or ""))
            self._al_check_pyopenms_ui()
            return
        worker = resolve_align_worker()
        if worker is None:
            QMessageBox.critical(
                self, "Worker missing",
                f"{ALIGN_WORKER_NAME} not found.\n"
                "Place it in the same folder as LipidZoner.")
            return

        det_params = dict(
            mass_error_ppm=float(self._al_sp_mass_ppm.value()),
            noise_threshold_int=float(self._al_sp_noise.value()),
            chrom_peak_snr=float(self._al_sp_snr.value()),
            chrom_fwhm=float(self._al_sp_fwhm.value()),
            min_trace_length=float(self._al_sp_min_trace.value()),
            centroid=self._al_cmb_centroid.currentText(),
        )
        # 実際に使った検出条件を保持し、session とサマリに載せる
        self._al_detect_params_used = dict(det_params)
        jobs, tmp_csvs = [], {}
        for k, (sample_name, p, _) in enumerate(self._al_loaded_files):
            with _tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
                out_csv = tmp.name
            tmp_csvs[k] = out_csv
            jobs.append({"key": k, "out_csv": out_csv,
                         "cmd": _build_worker_cmd(worker, p, out_csv, **det_params)})

        max_workers = min(len(jobs), _DETECT_MAX_WORKERS,
                          max(1, os.cpu_count() or _DETECT_MAX_WORKERS))
        progress = QProgressDialog(
            f"Detecting peaks ({max_workers} parallel)…", "Cancel",
            0, len(jobs), self)
        progress.setWindowTitle("Peak detection")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)

        def _on_tick(completed, total, running_keys):
            progress.setValue(completed)
            running = ", ".join(
                Path(self._al_loaded_files[rk][1]).name for rk in running_keys[:3])
            more = "..." if len(running_keys) > 3 else ""
            progress.setLabelText(
                f"Detected {completed}/{total}   |   running "
                f"({max_workers} parallel): {running}{more}")
            QApplication.processEvents()

        results, canceled = _run_worker_pool(
            jobs, max_workers,
            should_cancel=progress.wasCanceled, on_tick=_on_tick)

        errors = []
        for k, (sample_name, p, _) in enumerate(self._al_loaded_files):
            r = results.get(k)
            if r is None:
                continue          # canceled で未着手
            if r.get("df") is not None:
                self._al_loaded_files[k] = (sample_name, p, r["df"])
            elif r.get("error"):
                errors.append(f"{Path(p).name}: {r['error']}")
        for out_csv in tmp_csvs.values():
            try:
                os.unlink(out_csv)
            except Exception:
                pass
        progress.setValue(len(self._al_loaded_files))

        self._al_refresh_file_table()
        self._al_update_status()
        n_total = sum(len(df) for _, _, df in self._al_loaded_files
                      if df is not None)
        self._al_lbl_detect.setText(
            f"✓ {len(self._al_loaded_files)} samples, {n_total} peaks detected")
        # クロマトグラムを開いていたら更新
        dlg = getattr(self, '_al_chrom_dialog', None)
        if dlg is not None:
            try:
                dlg._viewer.set_samples(
                    [n for n, _, _ in self._al_loaded_files],
                    [p for _, p, _ in self._al_loaded_files])
            except Exception:
                pass

        if canceled:
            QMessageBox.information(
                self, "Peak detection canceled",
                "Cancelled. Results for the files already finished are kept.")
        elif errors:
            QMessageBox.warning(
                self, "Some files failed",
                f"Detection failed for {len(errors)} file(s):\n\n" + "\n".join(errors[:10]))
        elif not quiet:
            QMessageBox.information(
                self, "Peak detection complete",
                f"Detected {n_total} peaks across {len(self._al_loaded_files)} samples.")

    def _al_view_chromatogram(self):
        """TIC / EIC ビューアを開く(raw モードのみ)。"""
        if not self._al_loaded_files:
            QMessageBox.information(
                self, "No files", "Add .mzML files in Step 1.")
            return
        names = [n for n, _, _ in self._al_loaded_files]
        paths = [p for _, p, _ in self._al_loaded_files]
        try:
            dlg = ChromatogramDialog(names, paths, parent=self)
            self._al_chrom_dialog = dlg
            dlg.show()
            dlg.raise_()
            dlg.activateWindow()
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(
                self, "Chromatogram failed", f"{type(e).__name__}: {e}")

    def _al_show_plots(self):
        """アライメント結果の可視化ウィンドウを開く。"""
        if self._al_aligned_df is None or self._al_aligned_df.empty:
            QMessageBox.information(
                self, "No data", "Run alignment in Step 3 first.")
            return
        names = [n for n, _, _ in self._al_loaded_files]
        paths = [str(p) for _, p, _ in self._al_loaded_files]
        params = self._al_info.get('params', {})
        try:
            dlg = AlignmentPlotWindow(
                self._al_aligned_df, names, paths,
                raw_mode=True,          # raw 専用
                intensity_col=('area' if params.get('height_col') == 'Area'
                               else 'height'),
                mz_tol_eic=float(params.get('align_mz_tol', 0.005)),
                parent=self)
            dlg.show()
            self._al_plot_window = dlg   # 参照を保持しないと即座に破棄される
        except Exception as e:
            import traceback
            traceback.print_exc()
            QMessageBox.critical(
                self, "Plot failed", f"{type(e).__name__}: {e}")

    # ── ① タブ: パラメータ自動推定 (fix20 で Alignment から移植) ──
    def _al_detected_dfs_with(self, col):
        return [df for _n, _p, df in self._al_loaded_files
                if df is not None and col in df.columns]

    def _al_ensure_detected(self, col, parent=None):
        """Return detected per-sample DataFrames that contain `col`. If none
        are available yet, offer to run peak detection first."""
        have = self._al_detected_dfs_with(col)
        if have:
            return have
        if not self._al_loaded_files:
            QMessageBox.information(parent or self, "No files",
                "Add .mzML files first (Step 1).")
            return []
        r = QMessageBox.question(parent or self, "Detect peaks first",
            "Peaks have not been detected yet (or were detected with an older\n"
            "version that lacks this measurement). Run peak detection now with\n"
            "the current parameters?",
            QMessageBox.Yes | QMessageBox.No)
        if r != QMessageBox.Yes:
            return []
        self._al_detect_peaks()
        return self._al_detected_dfs_with(col)

    def _al_estimate_fwhm_from_data(self, parent=None):
        """Set 'Expected peak FWHM' to the median measured peak width.
        Returns a human-readable result string, or None if it could not run."""
        have = self._al_ensure_detected('fwhm_sec', parent)
        if not have:
            return None
        vals = []
        for df in have:
            v = pd.to_numeric(df['fwhm_sec'], errors='coerce')
            v = v[(v > 0.0) & (v < 120.0)].dropna()
            if len(v):
                vals.append(v.to_numpy())
        if not vals:
            return None
        allv = np.concatenate(vals)
        med = round(float(np.median(allv)), 1)
        q1, q3 = (float(x) for x in np.percentile(allv, [25, 75]))
        # 提案するだけ
        prop = self._al_propose('sp_fwhm', med)
        return (f"Expected peak FWHM\n"
                f"  median measured FWHM = {med} s  "
                f"(N={allv.size}, IQR {q1:.1f}-{q3:.1f} s)\n"
                f"  -> proposed FWHM = {prop} s (not applied yet).\n"
                f"     Min trace length follows FWHM at 0.6x when applied.")

    def _al_estimate_mass_accuracy_from_data(self, parent=None):
        """Set 'Mass trace m/z accuracy' to ~3x the median per-trace m/z
        scatter (the instrument's measured mass precision, in ppm).
        Returns a result string, or None if it could not run."""
        have = self._al_ensure_detected('mz_sd_ppm', parent)
        if not have:
            return None
        vals = []
        for df in have:
            v = pd.to_numeric(df['mz_sd_ppm'], errors='coerce')
            v = v[(v > 0.0) & (v < 100.0)].dropna()
            if len(v):
                vals.append(v.to_numpy())
        if not vals:
            return None
        allv = np.concatenate(vals)
        med = float(np.median(allv))
        p95 = float(np.percentile(allv, 95))
        tol = round(3.0 * med, 1)
        # 提案するだけ。[Apply] を押すまで入力欄は変えない。
        prop = self._al_propose('sp_mass_ppm', tol)
        return (f"Mass trace m/z accuracy\n"
                f"  median per-trace m/z scatter = {med:.2f} ppm  "
                f"(N={allv.size}, p95 {p95:.2f} ppm)\n"
                f"  -> proposed tolerance = 3x median = {prop} ppm "
                f"(not applied yet)")

    def _al_estimate_noise_by_consistency(self, parent=None):
        """v0.5.10: Auto-tune the Noise threshold from cross-sample consistency.
        Runs ONE detection pass at a low base noise on all samples (using the
        current m/z accuracy / S/N / FWHM / centroid settings), then sweeps
        candidate noise levels by Height-thresholding and scores each by
        n_robust * full_fraction. Sets the Noise field to the best value.
        Returns a result string, or None if it could not run."""
        owner = parent or self
        if not self._al_loaded_files:
            QMessageBox.information(owner, "No files",
                "Add .mzML files first (Step 1).")
            return None
        N = len(self._al_loaded_files)
        if N < 3:
            QMessageBox.information(owner, "Need more samples",
                "Cross-sample consistency tuning needs at least 3 samples.")
            return None
        if not check_pyopenms()[0]:
            QMessageBox.warning(owner, "pyOpenMS unavailable",
                "Raw peak detection is unavailable.")
            return None
        base_noise = 50.0
        r = QMessageBox.question(owner, "Consistency-based noise tuning",
            f"This runs ONE peak-detection pass on all {N} samples at a low "
            f"base noise ({int(base_noise)}), then sweeps candidate noise "
            f"levels to maximise cross-sample consistency.\n\n"
            f"It may take tens of seconds. Continue?",
            QMessageBox.Yes | QMessageBox.No)
        if r != QMessageBox.Yes:
            return None
        progress = QProgressDialog(
            f"Detecting at base noise {int(base_noise)} ...",
            "Cancel", 0, N, owner)
        progress.setWindowTitle("Consistency tuning")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        base_dfs = []
        errors = []
        for k, (name, pth, _df) in enumerate(self._al_loaded_files):
            progress.setValue(k)
            progress.setLabelText(f"Detecting (base noise): {Path(pth).name}")
            QApplication.processEvents()
            if progress.wasCanceled():
                return None
            try:
                df = load_raw_mzml_peaks(
                    pth,
                    mass_error_ppm=float(self._al_sp_mass_ppm.value()),
                    noise_threshold_int=base_noise,
                    chrom_peak_snr=float(self._al_sp_snr.value()),
                    chrom_fwhm=float(self._al_sp_fwhm.value()),
                    min_trace_length=float(self._al_sp_min_trace.value()),
                    centroid=self._al_cmb_centroid.currentText(),
                )
                base_dfs.append(df)
            except Exception as e:
                errors.append(f"{Path(pth).name}: {e}")
        progress.setValue(N)
        if errors or len(base_dfs) < 3:
            QMessageBox.warning(owner, "Detection failed",
                "Base detection failed:\n\n" + "\n".join(errors[:8]))
            return None
        mz_tol = float(self._al_sp_mz_tol.value())
        rt_tol = float(self._al_sp_rt_tol.value())
        cv_max = 0.20
        best, rows = consistency_noise_sweep(
            base_dfs, mz_tol=mz_tol, rt_tol=rt_tol, cv_max=cv_max)
        if best is None:
            QMessageBox.warning(owner, "No result",
                "Consistency sweep produced no usable features.")
            return None
        prop = self._al_propose('sp_noise', float(best))
        lines = [
            "Noise threshold - cross-sample consistency tuning",
            f"  (N={N} samples, mz_tol={mz_tol:.4f}, rt_tol={rt_tol:.3f}, "
            f"CV<={int(cv_max*100)}%)",
            f"  {'noise':>6} {'total':>6} {'full':>6} {'full%':>6} "
            f"{'medCV':>6} {'robust':>6} {'score':>8}",
        ]
        for r_ in rows:
            mc = (f"{r_['median_cv']*100:.1f}%"
                  if r_['median_cv'] == r_['median_cv'] else "  -  ")
            mark = "  <= best" if r_['noise'] == best else ""
            lines.append(
                f"  {r_['noise']:>6} {r_['n_total']:>6} {r_['n_full']:>6} "
                f"{r_['full_frac']*100:>5.1f}% {mc:>6} {r_['n_robust']:>6} "
                f"{r_['score']:>8.1f}{mark}")
        lines.append(f"  -> best noise = {best}; proposed value = {prop:g} "
                     f"(not applied yet).")
        lines.append("  (approximate sweep: one base detection + Height "
                     "thresholding; apply via real re-detection.)")
        return "\n".join(lines)

    def _al_estimate_snr_by_consistency(self, parent=None):
        """v0.5.11: Auto-tune Chromatographic S/N by cross-sample consistency.
        Unlike noise, S/N filtering happens inside ElutionPeakDetection and
        cannot be simulated by post-hoc thresholding, so each candidate needs a
        real detection pass on all samples (slower). Scores each candidate with
        the same objective (n_robust * full_fraction) evaluated at the current
        noise threshold, and sets the S/N field to the best value."""
        owner = parent or self
        if not self._al_loaded_files:
            QMessageBox.information(owner, "No files",
                "Add .mzML files first (Step 1).")
            return None
        N = len(self._al_loaded_files)
        if N < 3:
            QMessageBox.information(owner, "Need more samples",
                "Cross-sample consistency tuning needs at least 3 samples.")
            return None
        if not check_pyopenms()[0]:
            QMessageBox.warning(owner, "pyOpenMS unavailable",
                "Raw peak detection is unavailable.")
            return None
        snr_cands = [1.0, 2.0, 3.0, 5.0]
        noise = float(self._al_sp_noise.value())
        mz_tol = float(self._al_sp_mz_tol.value())
        rt_tol = float(self._al_sp_rt_tol.value())
        cv_max = 0.20
        r = QMessageBox.question(owner, "Consistency-based S/N tuning",
            f"This re-detects all {N} samples for each candidate S/N "
            f"({', '.join(str(int(x)) for x in snr_cands)}) at the current "
            f"noise ({int(noise)}) -- {len(snr_cands)} x {N} = "
            f"{len(snr_cands)*N} detection passes.\n\n"
            f"This can take a minute or more. Continue?",
            QMessageBox.Yes | QMessageBox.No)
        if r != QMessageBox.Yes:
            return None
        total = len(snr_cands) * N
        progress = QProgressDialog("Consistency S/N tuning ...",
                                   "Cancel", 0, total, owner)
        progress.setWindowTitle("Consistency tuning (S/N)")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)
        step = 0
        results = []
        for snr in snr_cands:
            dfs = []
            failed = False
            for (name, pth, _df) in self._al_loaded_files:
                progress.setValue(step); step += 1
                progress.setLabelText(f"S/N={int(snr)}: {Path(pth).name}")
                QApplication.processEvents()
                if progress.wasCanceled():
                    return None
                try:
                    dfs.append(load_raw_mzml_peaks(
                        pth,
                        mass_error_ppm=float(self._al_sp_mass_ppm.value()),
                        noise_threshold_int=noise,
                        chrom_peak_snr=snr,
                        chrom_fwhm=float(self._al_sp_fwhm.value()),
                    min_trace_length=float(self._al_sp_min_trace.value()),
                        centroid=self._al_cmb_centroid.currentText()))
                except Exception:
                    failed = True
                    break
            if failed or len(dfs) < 3:
                continue
            _b, rows = consistency_noise_sweep(
                dfs, mz_tol=mz_tol, rt_tol=rt_tol,
                candidates=[noise], cv_max=cv_max)
            if rows:
                row = rows[0]
                results.append((snr, row['n_total'], row['n_full'],
                                row['full_frac'], row['median_cv'],
                                row['n_robust'], row['score']))
        progress.setValue(total)
        if not results:
            QMessageBox.warning(owner, "No result",
                "S/N sweep produced no usable result.")
            return None
        best = max(results, key=lambda t: t[6])[0]
        prop = self._al_propose('sp_snr', float(best))
        lines = [
            "Chromatographic S/N - cross-sample consistency tuning",
            f"  (N={N} samples, noise={int(noise)}, CV<={int(cv_max*100)}%)",
            f"  {'S/N':>5} {'total':>6} {'full':>6} {'full%':>6} "
            f"{'medCV':>6} {'robust':>6} {'score':>8}",
        ]
        for (snr, nt, nf, ff, mcv, nr, sc) in results:
            mc = f"{mcv*100:.1f}%" if mcv == mcv else "  -  "
            mark = "  <= best" if snr == best else ""
            lines.append(f"  {int(snr):>5} {nt:>6} {nf:>6} {ff*100:>5.1f}% "
                         f"{mc:>6} {nr:>6} {sc:>8.1f}{mark}")
        lines.append(f"  -> best S/N = {int(best)}; proposed value = {prop:g} "
                     f"(not applied yet).")
        return "\n".join(lines)

    def _al_estimate_alignment_tolerances_from_data(self, parent=None):
        """v0.5.14: Suggest Step 3 m/z & RT tolerances from the cross-sample
        scatter of full-coverage features in the current detection."""
        owner = parent or self
        have = self._al_ensure_detected('Height', owner)
        if not have:
            return None
        if len(have) < 3:
            QMessageBox.information(owner, "Need more samples",
                "Tolerance estimation needs at least 3 detected samples.")
            return None
        mz_tol, rt_tol, st = estimate_alignment_tolerances(have)
        if mz_tol is None:
            QMessageBox.warning(owner, "Not enough data",
                "Could not find enough full-coverage features to estimate "
                "tolerances. Detect peaks on >=3 samples first.")
            return None
        # 提案するだけ
        a_mz = self._al_propose('sp_mz_tol', mz_tol)
        a_rt = self._al_propose('sp_rt_tol', rt_tol)
        return (
            "Alignment tolerances (cross-sample scatter of full-coverage "
            "features)\n"
            f"  clean features used: {st['n_clean']} / {st['n_full']} full\n"
            f"  m/z spread mDa: median {st['mz_med']*1e3:.2f}, "
            f"p95 {st['mz_p95']*1e3:.2f}, p99 {st['mz_p99']*1e3:.2f}\n"
            f"  RT  spread sec: median {st['rt_med']*60:.2f}, "
            f"p95 {st['rt_p95']*60:.2f}, p99 {st['rt_p99']*60:.2f}\n"
            f"  -> proposed m/z tolerance = p99 = {a_mz} Da; "
            f"RT tolerance = p99 = {a_rt} min (not applied yet).")


    # ── ① タブ: 検出パラメータの収集 / 適用 ──────────────
    def _al_collect_detect_params(self) -> dict:
        """Step 2 の UI 値を DEFAULT_DETECT_PARAMS と同じキーで返す。"""
        return {
            'mass_error_ppm':      float(self._al_sp_mass_ppm.value()),
            'noise_threshold_int': float(self._al_sp_noise.value()),
            'chrom_peak_snr':      float(self._al_sp_snr.value()),
            'chrom_fwhm':          float(self._al_sp_fwhm.value()),
            'min_trace_length':    float(self._al_sp_min_trace.value()),
            'centroid':            str(self._al_cmb_centroid.currentText()),
        }

    def _al_on_gapfill_toggled(self, on: bool):
        """ギャップフィリングを切ると下の 2 つを無効にする。"""
        for w in (getattr(self, '_al_sp_gf_factor', None),
                  getattr(self, '_al_sp_gf_ppm', None)):
            if w is not None:
                w.setEnabled(bool(on))

    def _al_apply_align_params(self, p: dict):
        """dict(DEFAULT_ALIGN_PARAMS のキー)を Step 3 の UI に反映する。"""
        if not p:
            return
        g = lambda k: p.get(k, DEFAULT_ALIGN_PARAMS[k])
        self._al_sp_mz_tol.setValue(float(g('align_mz_tol')))
        self._al_sp_rt_tol.setValue(float(g('align_rt_tol')))
        _ppm = g('align_mz_ppm')
        self._al_cb_ppm.setChecked(_ppm is not None)
        if _ppm is not None:
            self._al_sp_ppm.setValue(float(_ppm))
        self._al_sp_refine.setValue(int(g('align_n_refine')))
        _hc = str(g('height_col'))
        self._al_rb_area.setChecked(_hc == 'Area')
        self._al_rb_height.setChecked(_hc != 'Area')
        self._al_cb_rt_drift.setChecked(bool(g('rt_drift_correct')))
        self._al_cb_iso.setChecked(bool(g('iso_enabled')))
        self._al_sp_iso_mz.setValue(float(g('iso_mz_tol')))
        self._al_sp_iso_rt.setValue(float(g('iso_rt_tol')))
        self._al_sp_max_iso.setValue(int(g('iso_max_iso')))
        self._al_sp_iso_ratio.setValue(float(g('iso_ratio_tol')))
        self._al_sp_kC.setValue(float(g('iso_kC')))
        self._al_sp_kO.setValue(float(g('iso_kO')))

        self._al_cb_gapfill.setChecked(bool(g('gap_fill')))
        self._al_sp_gf_factor.setValue(float(g('gap_fill_width_factor')))
        self._al_sp_gf_ppm.setValue(float(g('gap_fill_mz_ppm')))
        self._al_on_gapfill_toggled(bool(g('gap_fill')))

    def _al_apply_detect_params(self, p: dict):
        """dict(DEFAULT_DETECT_PARAMS のキー)を Step 2 の UI に反映する。"""
        if not p:
            return
        g = lambda k: p.get(k, DEFAULT_DETECT_PARAMS[k])
        # FWHM は valueChanged で min_trace_length を上書きするので、
        # FWHM を先に入れてから min_trace_length を最後に入れ直す。
        self._al_sp_mass_ppm.setValue(float(g('mass_error_ppm')))
        self._al_sp_noise.setValue(float(g('noise_threshold_int')))
        self._al_sp_snr.setValue(float(g('chrom_peak_snr')))
        self._al_sp_fwhm.setValue(float(g('chrom_fwhm')))
        self._al_sp_min_trace.setValue(float(g('min_trace_length')))
        self._al_cmb_centroid.setCurrentText(str(g('centroid')))

    def _al_collect_session_section(self) -> dict:
        """session に保存する alignment セクションを組み立てる。

        pos / neg のパラメータをそれぞれ 'by_polarity' に入れる。
        旧版が読めるよう、表示中の極性の値は従来どおり 'params' /
        'detect' にも入れておく。

        2 列になったので **両方の列を退避してから** 読む。
        fix33 までは片方しか退避しておらず、未設定の極性に表示中の列の
        値が入ってしまっていた(neg の noise に pos の 200 が入る等)。
        """
        cur = getattr(self, '_al_cur_pol', 'pos')
        for pol in ('pos', 'neg'):
            self._al_stash_params(pol)
        state = self._al_ensure_state()
        by_pol = {}
        for pol in ('pos', 'neg'):
            st = state[pol]
            by_pol[pol] = {
                'params': dict(st['params']) if st.get('params')
                          else self._al_with_pol(pol, self._al_collect_params),
                'detect': dict(st['detect']) if st.get('detect')
                          else self._al_with_pol(
                              pol, self._al_collect_detect_params),
            }
        return {
            'params': dict(by_pol[cur]['params']),
            'detect': dict(by_pol[cur]['detect']),
            'input_mode': 'raw',          # raw 専用
            'ion_mode': cur,
            'by_polarity': by_pol,
        }

    def _al_apply_session_section(self, sec: dict):
        """session の alignment セクションを ① タブに復元する。

        by_polarity があれば **列ごとに** 反映する。無い古い
        session は 'params' / 'detect' を両列に入れる(fix33 までの
        単一 UI で保存されたものなので、それが両極性の値だった)。
        """
        if not isinstance(sec, dict):
            return
        by_pol = sec.get('by_polarity')
        # 入力モードのラジオは撤去した。古い session に
        #   input_mode: 'txt' が入っていても raw として扱う。
        state = self._al_ensure_state()
        for pol in ('pos', 'neg'):
            d = (by_pol or {}).get(pol) if isinstance(by_pol, dict) else None
            if not isinstance(d, dict):
                d = {'params': sec.get('params') or {},
                     'detect': sec.get('detect') or {}}
            try:
                if isinstance(d.get('params'), dict) and d['params']:
                    state[pol]['params'] = dict(d['params'])
                    self._al_with_pol(
                        pol, self._al_apply_align_params, d['params'])
                if isinstance(d.get('detect'), dict) and d['detect']:
                    state[pol]['detect'] = dict(d['detect'])
                    self._al_with_pol(
                        pol, self._al_apply_detect_params, d['detect'])
            except Exception as e:
                log.warning(f"apply session section failed ({pol}): {e}")
        try:
            want = sec.get('ion_mode')
            if want in ('pos', 'neg'):
                self._al_cur_pol = want
                self._al_switching = True
                try:
                    (self._al_rb_neg if want == 'neg'
                     else self._al_rb_pos).setChecked(True)
                finally:
                    self._al_switching = False
        except Exception as e:
            log.warning(f"ion mode restore failed: {e}")
        try:
            self._al_sync_ui_to_state()
        except Exception:
            pass

    # ── ① タブ: パラメータ収集 ────────────────────────────────────
    def _al_collect_params(self) -> dict:
        """UI の値を DEFAULT_ALIGN_PARAMS と同じキーの dict にまとめる。"""
        return {
            'align_mz_tol':     float(self._al_sp_mz_tol.value()),
            'align_rt_tol':     float(self._al_sp_rt_tol.value()),
            'align_mz_ppm':     (float(self._al_sp_ppm.value())
                                 if self._al_cb_ppm.isChecked() else None),
            'align_n_refine':   int(self._al_sp_refine.value()),
            'height_col':       ('Area' if self._al_rb_area.isChecked()
                                 else 'Height'),
            'rt_drift_correct': bool(self._al_cb_rt_drift.isChecked()),
            'iso_mz_tol':       float(self._al_sp_iso_mz.value()),
            'iso_rt_tol':       float(self._al_sp_iso_rt.value()),
            'iso_max_iso':      int(self._al_sp_max_iso.value()),
            'iso_ratio_tol':    float(self._al_sp_iso_ratio.value()),
            'iso_kC':           float(self._al_sp_kC.value()),
            'iso_kO':           float(self._al_sp_kO.value()),
            'iso_enabled':      bool(self._al_cb_iso.isChecked()),

            'gap_fill':              bool(self._al_cb_gapfill.isChecked()),
            'gap_fill_width_factor': float(self._al_sp_gf_factor.value()),
            'gap_fill_mz_ppm':       float(self._al_sp_gf_ppm.value()),
        }

    def _al_reset_params(self):
        """パラメータを DEFAULT_ALIGN_PARAMS に戻す。"""
        P = DEFAULT_ALIGN_PARAMS
        self._al_sp_mz_tol.setValue(P['align_mz_tol'])
        self._al_sp_rt_tol.setValue(P['align_rt_tol'])
        self._al_cb_ppm.setChecked(P['align_mz_ppm'] is not None)
        if P['align_mz_ppm']:
            self._al_sp_ppm.setValue(float(P['align_mz_ppm']))
        self._al_sp_refine.setValue(P['align_n_refine'])
        self._al_rb_height.setChecked(P['height_col'] == 'Height')
        self._al_rb_area.setChecked(P['height_col'] == 'Area')
        self._al_cb_rt_drift.setChecked(bool(P['rt_drift_correct']))
        self._al_cb_iso.setChecked(bool(P['iso_enabled']))

        self._al_cb_gapfill.setChecked(bool(P['gap_fill']))
        self._al_sp_gf_factor.setValue(float(P['gap_fill_width_factor']))
        self._al_sp_gf_ppm.setValue(float(P['gap_fill_mz_ppm']))
        self._al_on_gapfill_toggled(bool(P['gap_fill']))
        self._al_sp_iso_mz.setValue(P['iso_mz_tol'])
        self._al_sp_iso_rt.setValue(P['iso_rt_tol'])
        self._al_sp_max_iso.setValue(P['iso_max_iso'])
        self._al_sp_iso_ratio.setValue(P['iso_ratio_tol'])
        self._al_sp_kC.setValue(P['iso_kC'])
        self._al_sp_kO.setValue(P['iso_kO'])
        # Step 2 の検出パラメータも既定に戻す
        try:
            self._al_apply_detect_params(dict(DEFAULT_DETECT_PARAMS))
        except Exception:
            pass

    def _al_on_mode_changed(self):
        """入力モード切替(fix35 で空実装)。

        fix24 で raw mzML 固定になり、fix35 でラジオそのものを撤去した。
        session の復元経路など外から呼ばれる可能性が残っているので、
        名前だけ残して no-op にしてある。
        """
        try:
            self._al_step2.setVisible(True)
        except Exception:
            pass

    def _al_ensure_pyopenms_checked(self):
        """① タブが開かれたときに一度だけ pyOpenMS を確認する。"""
        if getattr(self, '_al_pyopenms_checked', False):
            return
        self._al_pyopenms_checked = True
        try:
            self._al_check_pyopenms_ui()
        except Exception as e:
            log.warning(f"pyopenms check failed: {e}")

    def _al_on_iso_toggled(self, on: bool):
        for w in (self._al_sp_iso_mz, self._al_sp_iso_rt, self._al_sp_max_iso,
                  self._al_sp_iso_ratio, self._al_sp_kC, self._al_sp_kO):
            w.setEnabled(bool(on))

    # ── ① タブ: ファイル操作 ─────────────────────────────────────
    def _al_add_files(self):
        # raw mzML 専用。TXT を選ばせる分岐は撤去した。
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add raw data .mzML file(s)", "",
            "mzML files (*.mzML);;All (*)")
        if paths:
            self._al_add_paths(paths)

    def _al_add_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Add folder", "")
        if not folder:
            return
        f = Path(folder)
        # raw mzML 専用。TXT を拾う分岐は撤去した。
        paths = sorted(str(x) for x in f.glob('*.mzML'))
        if not paths:
            QMessageBox.warning(
                self, "No files",
                f"No mzML file was found in:\n{folder}")
            return
        self._al_add_paths(paths)

    def _al_add_paths(self, paths: list):
        """ファイルを極性ごとに振り分けて登録する。

        fix24 までは「読み込んだ mzML の polarity が混在していたら警告」
        だったが、pos と neg を両方読むのが通常フローなので、
        ファイル単位で polarity を読んで pos / neg のリストへ自動で
        振り分けるように変えた。polarity を読めなかったものは
        今対象にしている極性へ入れる。
        """
        raw = True          # このタブは raw mzML 専用
        state = self._al_ensure_state()
        cur = getattr(self, '_al_cur_pol', 'pos')
        added = {'pos': 0, 'neg': 0}
        errors: list[str] = []
        unknown: list[str] = []
        if not hasattr(self, '_al_polarity_cache'):
            self._al_polarity_cache = {}
        for p in paths:
            p = Path(p)
            name = p.stem
            pol = None
            if str(p).lower().endswith('.mzml'):
                pol = self._al_polarity_cache.get(str(p))
                if pol is None:
                    try:
                        pol = detect_mzml_polarity(p)
                    except Exception as e:
                        log.warning(f"polarity detect failed {p.name}: {e}")
                        pol = None
                    self._al_polarity_cache[str(p)] = pol
            if pol not in ('pos', 'neg'):
                unknown.append(p.name)
                pol = cur
            if any(nm == name for nm, _, _ in state[pol]['loaded_files']):
                continue
            try:
                # raw mzML 専用。検出は Step 2 で行うので
                # ここでは中身を読まない。
                df = None
            except Exception as e:
                errors.append(f"{p.name}: {e}")
                continue
            state[pol]['loaded_files'].append((name, p, df))
            added[pol] += 1

        # 今の極性に 1 件も入らず、もう一方に入ったならそちらへ寄せる。
        # (mzML/neg だけを追加したのに pos の空の表が出る、を防ぐ)
        if added[cur] == 0:
            other = 'neg' if cur == 'pos' else 'pos'
            if added[other]:
                self._al_set_polarity(other)

        msgs = [f"{pol.upper()} {n}" for pol, n in added.items() if n]
        if msgs:
            txt = "Added: " + " / ".join(msgs)
            if unknown:
                txt += (f"  ({len(unknown)} with unknown polarity went to "
                        f"{cur.upper()})")
            self._al_lbl_ion_auto.setText(txt)
            self._al_lbl_ion_auto.setStyleSheet("color:#070;")
            self._al_lbl_ion_auto.setVisible(True)
        elif errors or unknown:
            self._al_lbl_ion_auto.setText("No file could be added")
            self._al_lbl_ion_auto.setStyleSheet("color:#b00;")
            self._al_lbl_ion_auto.setVisible(True)

        self._al_sync_ui_to_state()
        if errors:
            QMessageBox.warning(
                self, "Some files failed",
                f"Failed to load {len(errors)} file(s):\n\n" + '\n'.join(errors[:10]))

    def _al_remove_selected(self):
        rows = sorted((i.row() for i in
                       self._al_file_table.selectionModel().selectedRows()),
                      reverse=True)
        if not rows:
            QMessageBox.information(
                self, "No selection", "Select a row first.")
            return
        for r in rows:
            if 0 <= r < len(self._al_loaded_files):
                self._al_loaded_files.pop(r)
        self._al_refresh_file_table()
        self._al_update_status()

    def _al_clear_all(self, _checked=False, *, polarity: str = 'ask'):
        """読み込んだファイルと結果を破棄する。

        polarity='current' なら表示中の極性だけ、'both' なら両方、
        'ask'(ボタンから押されたとき)なら、もう一方にもデータが
        あるときだけどちらを消すか聞く。
        呼び出し側が極性を明示した場合は確認しない。
        """
        state = self._al_ensure_state()
        cur = getattr(self, '_al_cur_pol', 'pos')
        other = 'neg' if cur == 'pos' else 'pos'
        if polarity == 'ask' and not state[other]['loaded_files']:
            polarity = 'current'
        if polarity == 'ask':
            box = QMessageBox(self)
            box.setWindowTitle("Clear all")
            box.setIcon(QMessageBox.Question)
            box.setText(
                f"{other.upper()} also holds {len(state[other]['loaded_files'])}"
                f" file(s). Which do you want to discard?")
            b_cur = box.addButton(f"{cur.upper()} only", QMessageBox.AcceptRole)
            b_both = box.addButton("Both", QMessageBox.DestructiveRole)
            box.addButton("Cancel", QMessageBox.RejectRole)
            box.exec()
            clicked = box.clickedButton()
            if clicked is b_both:
                polarity = 'both'
            elif clicked is b_cur:
                polarity = 'current'
            else:
                return
        targets = ('pos', 'neg') if polarity == 'both' else (cur,)
        for pol in targets:
            # 消すのはファイルと結果だけ。調整したパラメータは残す。
            # (fix24 の _al_clear_all もパラメータには触っていなかった。
            #  ここで state を丸ごと作り直すと Step2/Step3 のチューニングが
            #  黙って既定値へ戻ってしまう)
            keep_p = state[pol].get('params')
            keep_d = state[pol].get('detect')
            state[pol] = PreviewRTDialog._al_blank_state()
            state[pol]['params'] = keep_p
            state[pol]['detect'] = keep_d
        self._al_sync_ui_to_state()

    def _al_refresh_file_table(self):
        t = self._al_file_table
        t.setRowCount(len(self._al_loaded_files))
        for r, (name, path, df) in enumerate(self._al_loaded_files):
            if df is None:
                summary = "(not yet detected)"
            elif 'Isotope' in df.columns and (df['Isotope'] != 0).any():
                n0 = int((df['Isotope'] == 0).sum())
                n1 = int((df['Isotope'] == 1).sum())
                n2 = int((df['Isotope'] == 2).sum())
                summary = f"{n0} / {n1} / {n2} (M+0/M+1/M+2)"
            else:
                summary = f"{len(df)} peaks"
            t.setItem(r, 0, QTableWidgetItem(str(r + 1)))
            it = QTableWidgetItem(name); it.setToolTip(name)
            t.setItem(r, 1, it)
            t.setItem(r, 2, QTableWidgetItem(summary))
            ip = QTableWidgetItem(str(path)); ip.setToolTip(str(path))
            t.setItem(r, 3, ip)

    def _al_update_status(self):
        n = len(self._al_loaded_files)
        if n == 0:
            self._al_lbl_status.setText("Status: ready (no files loaded)")
        elif self._al_aligned_df is None:
            self._al_lbl_status.setText(
                f"Status: {n} file(s) loaded, not yet aligned")
        else:
            self._al_lbl_status.setText(
                f"Status: {n} file(s), "
                f"{len(self._al_aligned_df)} aligned features")

    # ── ① タブ: 実行 ─────────────────────────────────────────────
    def _al_run_alignment(self, _checked=False, *, quiet: bool = False):
        """Step 3 を実行する。表示中の極性に対して走る。

        quiet は現状ダイアログを出さないが、将来の通知追加に
        備えて _al_detect_peaks と signature を揃えてある。
        """
        if not self._al_loaded_files:
            QMessageBox.information(
                self, "No data", "Add files in Step 1.")
            return
        missing = [n for n, _, df in self._al_loaded_files if df is None]
        if missing:
            QMessageBox.warning(
                self, "Peaks not detected",
                "Some files have not been detected yet:\n  "
                + "\n  ".join(missing[:10])
                + "\n\nRun [Detect Peaks] in Step 2 first.")
            return

        params = self._al_collect_params()
        sample_names = [n for n, _, _ in self._al_loaded_files]
        peak_dfs = [df for _, _, df in self._al_loaded_files]
        info = {'n_samples': len(peak_dfs),
                'n_peaks_in': int(sum(len(d) for d in peak_dfs))}
        # 検出条件も一緒に記録する(fix35: raw 専用になった)
        info['input_mode'] = 'raw'
        info['detect_params'] = dict(
            getattr(self, '_al_detect_params_used', None)
            or self._al_collect_detect_params())

        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            if params['rt_drift_correct'] and len(peak_dfs) > 1:
                peak_dfs, offsets, stats = correct_rt_drift(
                    peak_dfs, mz_tol=max(params['align_mz_tol'], 0.01),
                    mz_ppm=params['align_mz_ppm'],
                    mz_min_da=params['align_mz_tol'])
                info['rt_offsets'] = offsets
                info['rt_drift_stats'] = stats
            aligned = simple_align_per_sample(
                peak_dfs,
                mz_tol=params['align_mz_tol'], rt_tol=params['align_rt_tol'],
                n_refine=params['align_n_refine'],
                mz_ppm=params['align_mz_ppm'],
                mz_min_da=params['align_mz_tol'],
                height_col=params['height_col'])
            info['n_features'] = int(len(aligned))
            # 同位体グルーピングの**前**に欠測セルを埋める。
            #   group_isotopes は検体間の強度中央値で親子を決めるので、
            #   埋めたあとの表を渡した方が判定が安定する。
            if params.get('gap_fill', False):
                QApplication.restoreOverrideCursor()
                aligned, gf_info = self._al_do_gap_fill(
                    aligned, params, sample_names)
                QApplication.setOverrideCursor(Qt.WaitCursor)
                info['gap_fill'] = gf_info
            aligned = group_isotopes(
                aligned, len(sample_names), intensity_prefix='height_s',
                do_grouping=bool(params['iso_enabled']),
                kC=params['iso_kC'], kO=params['iso_kO'],
                tol=params['iso_ratio_tol'],
                mz_tol=params['iso_mz_tol'], rt_tol=params['iso_rt_tol'],
                max_iso=params['iso_max_iso'])
            aligned = derive_iso_labels(aligned)
            info['n_monoisotopic'] = int((aligned['iso_position'] == 0).sum())
            info['params'] = params
        except Exception as e:
            QApplication.restoreOverrideCursor()
            import traceback
            traceback.print_exc()
            QMessageBox.critical(
                self, "Alignment failed", f"{type(e).__name__}: {e}")
            return
        QApplication.restoreOverrideCursor()

        self._al_aligned_df = aligned
        self._al_sample_names = sample_names
        self._al_info = info
        self._al_cur_state['sent'] = False   # 再アライメント後は未登録
        # 以後の読込経路(Manage Files / session)も同じパラメータを使うよう共有
        try:
            if self._host is not None:
                self._host._align_params = params
        except Exception:
            pass
        self._al_show_summary()
        for b in (self._al_btn_send, self._al_btn_txt, self._al_btn_csv,
                  self._al_btn_plots):
            b.setEnabled(True)
        self._al_update_status()
        self._al_update_pol_summary()

    def _al_do_gap_fill(self, aligned, params, sample_names):
        """欠測セルを raw mzML から埋める。(aligned, info) を返す。

        raw mzML が揃っていなければ何もせず、理由を info に入れて返す。
        アラインメント自体は成功しているので、ここで失敗しても落とさない。
        """
        raw_paths = [p for _, p, _ in self._al_loaded_files]
        worker = resolve_align_worker()
        if worker is None:
            return aligned, {'enabled': True,
                             'error': f'{ALIGN_WORKER_NAME} が見つからない'}

        progress = QProgressDialog(
            "Filling missing values from raw data…", "Cancel",
            0, len(raw_paths), self)
        progress.setWindowTitle("Gap filling")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setValue(0)

        def _on_tick(completed, total, running_keys):
            progress.setMaximum(total)
            progress.setValue(completed)
            running = ", ".join(Path(raw_paths[rk]).name
                                for rk in running_keys[:2]
                                if 0 <= rk < len(raw_paths))
            more = "..." if len(running_keys) > 2 else ""
            progress.setLabelText(
                f"Filled {completed}/{total} file(s)   |   running: "
                f"{running}{more}")
            QApplication.processEvents()

        try:
            out, gf_info = gap_fill_aligned(
                aligned, raw_paths, worker=worker,
                width_factor=float(params.get(
                    'gap_fill_width_factor', GAP_FILL_WIDTH_FACTOR)),
                mz_ppm=float(params.get('gap_fill_mz_ppm', GAP_FILL_MZ_PPM)),
                height_col=str(params.get('height_col', 'Height')),
                should_cancel=progress.wasCanceled, on_tick=_on_tick)
        except Exception as e:
            import traceback
            traceback.print_exc()
            return aligned, {'enabled': True,
                             'error': f'{type(e).__name__}: {e}'}
        finally:
            progress.close()

        if gf_info.get('error'):
            QMessageBox.warning(
                self, "Gap filling incomplete",
                "Alignment finished, but filling missing values did not:\n\n"
                f"{gf_info['error']}\n\n"
                "The table keeps the zeros for those cells.")
        return out, gf_info

    def _al_run_both(self, _checked=False):
        """ファイルのある極性を pos → neg の順に通しで実行する。

        各極性について「未検出なら Step 2 → Step 3」を回す。
        どこかで失敗したらそこで止め、済んだところまでは残す。
        """
        state = self._al_ensure_state()
        todo = [pol for pol in ('pos', 'neg') if state[pol]['loaded_files']]
        if not todo:
            QMessageBox.information(
                self, "No files", "Add .mzML files in Step 1.")
            return
        if not check_pyopenms()[0]:          # raw 専用
            self._al_check_pyopenms_ui()
            QMessageBox.warning(
                self, "pyOpenMS unavailable",
                "Raw mzML detection needs pyOpenMS.")
            return
        done = []
        for pol in todo:
            self._al_set_polarity(pol)
            self._al_sync_ui_to_state()
            QApplication.processEvents()
            if any(df is None for _, _, df in self._al_loaded_files):
                self._al_detect_peaks(quiet=True)
                if any(df is None for _, _, df in self._al_loaded_files):
                    QMessageBox.warning(
                        self, "Stopped",
                        f"Stopped: peak detection for {pol.upper()} did not finish.\n"
                        + ("\n".join(done) if done else ""))
                    return
            self._al_run_alignment(quiet=True)
            if self._al_cur_state['aligned_df'] is None:
                QMessageBox.warning(
                    self, "Stopped",
                    f"Stopped: alignment for {pol.upper()} failed.\n"
                    + ("\n".join(done) if done else ""))
                return
            done.append(f"  {pol.upper()}: {len(self._al_loaded_files)} samples "
                        f"→ {len(self._al_cur_state['aligned_df'])} features")
        self._al_sync_ui_to_state()
        QMessageBox.information(
            self, "Run both complete",
            "Both polarities finished:\n\n" + "\n".join(done)
            + "\n\nNext: press [Send to Annotation] to register both at once.")

    def _al_show_summary(self):
        a, info = self._al_aligned_df, self._al_info
        p = info.get('params', {})
        n = len(a)
        n_m0 = int((a['iso_position'] == 0).sum())
        dist = a['iso_position'].value_counts().sort_index()
        dist_txt = ', '.join(f"M+{int(k)}: {int(v)}" for k, v in dist.items())
        lines = [
            f"Samples            : {info.get('n_samples')}",
            f"Input peaks        : {info.get('n_peaks_in')}",
            f"Aligned features   : {n}",
            f"Monoisotopic (M+0) : {n_m0}  ({n_m0/n*100:.1f}%)" if n else "",
            f"Isotope positions  : {dist_txt}",
            "",
            "── Parameters ─────────────────────────────",
            f"  align m/z tol    : {p.get('align_mz_tol')} Da"
            + (f"  + {p.get('align_mz_ppm')} ppm"
               if p.get('align_mz_ppm') else ""),
            f"  align RT tol     : {p.get('align_rt_tol')} min",
            f"  refine iters     : {p.get('align_n_refine')}",
            f"  intensity column : {p.get('height_col')}",
            f"  RT drift correct : {p.get('rt_drift_correct')}",
            f"  isotope grouping : {p.get('iso_enabled')}"
            + (f"  (mz {p.get('iso_mz_tol')} / rt {p.get('iso_rt_tol')} / "
               f"max M+{p.get('iso_max_iso')} / ratio tol "
               f"{p.get('iso_ratio_tol')} / kC {p.get('iso_kC')} "
               f"kO {p.get('iso_kO')})" if p.get('iso_enabled') else ""),
        ]
        if info.get('input_mode') == 'raw' and info.get('detect_params'):
            d = info['detect_params']
            lines += ["",
                      "── Peak detection (raw mzML) ──────────────",
                      f"  mass accuracy    : {d.get('mass_error_ppm')} ppm",
                      f"  noise threshold  : {d.get('noise_threshold_int')}",
                      f"  peak S/N         : {d.get('chrom_peak_snr')}",
                      f"  expected FWHM    : {d.get('chrom_fwhm')} s",
                      f"  min trace length : {d.get('min_trace_length')} s",
                      f"  centroid mode    : {d.get('centroid')}"]
        # ギャップフィリング
        gf = info.get('gap_fill') or {}
        if gf:
            if gf.get('error'):
                lines += ["",
                          "── Gap filling ────────────────────────────",
                          f"  NOT done: {gf['error']}"]
            else:
                _nt = int(gf.get('n_targets', 0))
                _nf = int(gf.get('n_filled', 0))
                _ns = int(gf.get('n_suspect', 0))
                _pct = (100.0 * _nf / _nt) if _nt else 0.0
                lines += ["",
                          "── Gap filling ────────────────────────────",
                          f"  missing cells    : {_nt}",
                          f"  filled from raw  : {_nf}  ({_pct:.1f}%)",
                          f"  window half width: "
                          f"{gf.get('width_factor')} x peak width",
                          f"  m/z window       : {gf.get('mz_ppm')} ppm"]
                _nfo = int(gf.get('n_forced', 0))
                if _nfo:
                    lines += [f"  no local maximum : {_nfo}  (filled with the "
                              f"window maximum anyway)"]
                if _ns:
                    lines += [f"  CHECK: {_ns} filled value(s) are more than "
                              f"5x the median of the samples that",
                              "         detected the feature — possibly a "
                              "different peak in the window."]
                if gf.get('canceled'):
                    lines += ["  (canceled part way — some cells are still "
                              "zero)"]
        if 'rt_offsets' in info:
            offs = ', '.join(f"{o:+.4f}" for o in info['rt_offsets'])
            st = info.get('rt_drift_stats', {})
            lines += ["",
                      "── RT drift correction ────────────────────",
                      f"  common peaks used : {st.get('n_common')}",
                      f"  offsets (min)     : {offs}"]
        qc = []
        if 'fwhm_warning' in a.columns:
            qc.append(f"FWHM warning: {int(a['fwhm_warning'].sum())} features")
        if 'intra_sample_dup' in a.columns:
            qc.append(
                f"intra-sample duplicate: {int(a['intra_sample_dup'].sum())} features")
        # 何セルが埋めた値かを QC にも出す
        _fcols = [c for c in a.columns if str(c).startswith('filled_s')]
        if _fcols:
            try:
                _nfill = int(a[_fcols].fillna(False).to_numpy(dtype=bool).sum())
            except Exception:
                _nfill = 0
            _ncell = len(a) * len(_fcols)
            if _ncell:
                qc.append(f"filled cells: {_nfill} / {_ncell} "
                          f"({100.0 * _nfill / _ncell:.1f}%)")
        if qc:
            lines += ["", "── QC ─────────────────────────────────────"] \
                     + [f"  {x}" for x in qc]
        lines += ["",
                  "Next: press [Send to Annotation] to register this in the",
                  "Annotation tab (only M+0 is annotated; isotopes are",
                  "kept internally for quantification)."]
        self._al_txt_summary.setPlainText(
            '\n'.join(x for x in lines if x is not None))

    # ── ① タブ: 出力 ─────────────────────────────────────────────
    @staticmethod
    def _al_entry_desc(fe) -> str:
        """エントリの規模と検出条件を 1 行で表す。"""
        ai = getattr(fe, 'align_info', None) or {}
        dp = ai.get('detect_params') or {}
        bits = []
        try:
            bits.append(f"{len(fe._aligned_df):,} features")
        except Exception:
            pass
        if ai.get('n_samples'):
            bits.append(f"{ai['n_samples']} samples")
        if dp.get('noise_threshold_int') is not None:
            bits.append(f"noise {dp['noise_threshold_int']:g}")
        if dp.get('mass_error_ppm') is not None:
            bits.append(f"{dp['mass_error_ppm']:g} ppm")
        return "  ".join(bits) if bits else "(no information)"

    def _al_confirm_replace(self, existing, st: dict, pol: str) -> bool:
        """同じ極性を置き換えてよいか尋ねる。

        新旧の規模と検出条件を並べ、破棄されるものを列挙する。
        Advanced Setting の設定値は破棄しないこともここで明示する。
        """
        new_bits = []
        try:
            new_bits.append(f"{len(st['aligned_df']):,} features")
        except Exception:
            pass
        info = st.get('info') or {}
        dp = (info.get('detect_params')
              or self._al_collect_detect_params() or {})
        if st.get('loaded_files'):
            new_bits.append(f"{len(st['loaded_files'])} samples")
        if dp.get('noise_threshold_int') is not None:
            new_bits.append(f"noise {dp['noise_threshold_int']:g}")
        if dp.get('mass_error_ppm') is not None:
            new_bits.append(f"{dp['mass_error_ppm']:g} ppm")
        lost = self._al_describe_purge(pol)
        msg = (
            f"Annotation can hold only one per-sample dataset per polarity.\n"
            f"[{pol}] is already registered.\n\n"
            f"  current  {existing.tag}:  {self._al_entry_desc(existing)}\n"
            f"  new      →         {'  '.join(new_bits) or '(no information)'}\n\n"
            f"Replacing it discards the Annotation results and the manual\n"
            f"curation for [{pol}]:\n  " + ("\n  ".join(lost) if lost else "(none)")
            + "\n\nThe Advanced Setting values (ppm / Void volume cutoff /\n"
              "MADx / IS tol / Adduct / σ / Filter Classes) are kept, so\n"
              "running Run All again reproduces the same conditions.\n\n"
              "Replace it?")
        mb = QMessageBox(self)
        mb.setIcon(QMessageBox.Warning)
        mb.setWindowTitle(f"Replace [{pol}]")
        mb.setText(msg)
        b_yes = mb.addButton("Replace", QMessageBox.AcceptRole)
        b_no = mb.addButton("Cancel", QMessageBox.RejectRole)
        mb.setDefaultButton(b_no)
        mb.exec()
        return mb.clickedButton() is b_yes

    def _al_describe_purge(self, pol: str) -> list:
        """置き換えで失われる Annotation の成果を数える。"""
        out = []
        snap = (getattr(self, '_mode_state', None) or {}).get(pol)
        cur = (getattr(self, '_active_mode', None) == pol)
        mdf = (getattr(self, '_match_df', None) if cur
               else (snap or {}).get('_match_df'))
        try:
            if mdf is not None and not mdf.empty:
                _m = mdf[mdf['matched']]
                out.append(f"match results: {len(_m)} rows / "
                           f"{_m['lipid_class'].nunique()} classes")
        except Exception:
            pass
        ch = (getattr(self, '_is_filter_choices', None) if cur
              else (snap or {}).get('_is_filter_choices')) or {}
        try:
            n_pin = sum(1 for lst in ch.values() for e in lst
                        if e.get('selected_peak_idx') is not None)
            n_rt = sum(1 for lst in ch.values() for e in lst
                       if e.get('manual_rt') is not None)
            n_off = sum(1 for lst in ch.values() for e in lst
                        if not e.get('use'))
            if ch:
                out.append(f"③ IS choices for {len(ch)} classes"
                           f" (pinned {n_pin} / Manual RT {n_rt} / "
                           f"Use=off {n_off})")
        except Exception:
            pass
        win = (getattr(self, '_manual_winner_coords', None) if cur
               else (snap or {}).get('_manual_winner_coords')) or set()
        if win:
            out.append(f"{len(win)} manually chosen winners")
        qi = getattr(self, '_quant_ion_choices', None) or {}
        n_q = sum(1 for v in qi.values() if v == pol)
        if n_q:
            out.append(f"{n_q} classes whose ⑦ Quant Ion pointed at {pol}")
        return out

    def _al_purge_annotation_state(self, pol: str, tag: str) -> str:
        """置き換えた極性の Annotation の成果を捨てる。

        Advanced Setting の設定値(_spin_ppm / _spin_void_rt / _spin_is_tol /
        _spin_iqr / _spin_sigma / _spin_adduct_* / _cb_is_auto_intensity /
        _filter_class_exemptions)には **触らない**。
        チューニング済みの条件は残して解析だけやり直せるようにする。

        戻り値: 何を捨てたかの短い説明。
        """
        lost = self._al_describe_purge(pol)
        # AP 側: その極性の per-mode state
        try:
            (getattr(self, '_mode_state', None) or {}).pop(pol, None)
        except Exception:
            pass
        if getattr(self, '_active_mode', None) == pol:
            self._match_df = None
            for _a, _v in (('_matched_coords', set()),
                           ('_reserved_coords', set()),
                           ('_conflict_coords', set()),
                           ('_conflict_map', {}),
                           ('_outlier_coords', {}),
                           ('_coherence_outliers', {}),
                           ('_coherence_assignments', {}),
                           ('_coherence_models', {}),
                           ('_coherence_pairs', []),
                           ('_coherence_all_pairs', []),
                           ('_adduct_attribution', {}),
                           ('_manual_winner_coords', set()),
                           ('_class_ref_rt', {}),
                           ('_is_filter_choices', {})):
                try:
                    setattr(self, _a, type(_v)(_v))
                except Exception:
                    pass
        # ⑦ でこの極性を指していたクラスの選択も外す
        try:
            qi = getattr(self, '_quant_ion_choices', None)
            if isinstance(qi, dict):
                for _c in [c for c, v in qi.items() if v == pol]:
                    qi.pop(_c, None)
        except Exception:
            pass
        # host 側: tag / mode をキーにした残骸
        host = getattr(self, '_host', None)
        if host is not None:
            for _name, _key in (('_preview_settings_by_mode', pol),
                                ('_reserved_coords_by_tag', tag),
                                ('_coherence_loser_by_tag', tag),
                                ('_preview_match_df_by_tag', tag)):
                try:
                    d = getattr(host, _name, None)
                    if isinstance(d, dict):
                        d.pop(_key, None)
                except Exception:
                    pass
        log.info(f"discarded the Annotation state of {tag} [{pol}]: "
              + ("; ".join(lost) if lost else "(none)"))
        return "; ".join(lost)

    def _al_send_to_annotation(self, _checked=False):
        """アライメント済みの極性をまとめて host に登録する。

        fix24 までは表示中の極性 1 つだけを登録していた。pos と neg を
        両方使うのが通常フローなので、済んでいる方を全部登録する。
        既に登録した極性は、再アライメントするまで二重登録しない
        (state[pol]['sent'] で管理)。
        """
        host = self._host
        if host is None or not hasattr(host, '_file_entries'):
            QMessageBox.warning(
                self, "Cannot register", "Host main window not available.")
            return
        state = self._al_ensure_state()
        ready = [pol for pol in ('pos', 'neg')
                 if state[pol]['aligned_df'] is not None]
        if not ready:
            QMessageBox.information(
                self, "Nothing to register",
                "Run the Step 3 alignment first.")
            return
        fresh = [pol for pol in ready if not state[pol].get('sent')]
        if not fresh:
            QMessageBox.information(
                self, "Already registered",
                "Every aligned polarity is already registered.\n"
                "Run Step 3 again to register them anew.")
            return

        registered, failed, replaced = [], [], []
        for pol in fresh:
            st = state[pol]
            folders = {str(p.parent) for _, p, _ in st['loaded_files']}
            source_folder = folders.pop() if len(folders) == 1 else None
            # Annotation は極性ごとに 1 エントリだけ。既に同じ
            # 極性が登録されていれば、確認したうえで置き換える。
            existing = next(
                (f for f in (host._file_entries or [])
                 if getattr(f, 'ion_mode', None) == pol), None)
            if existing is not None and not self._al_confirm_replace(
                    existing, st, pol):
                continue
            _bumped = False
            try:
                if existing is not None:
                    tag = existing.tag
                else:
                    host._file_counter += 1
                    _bumped = True
                    tag = f"#{host._file_counter}"
                fe = PerSampleFileEntry(
                    aligned_df=st['aligned_df'],
                    sample_names=st['sample_names'],
                    ion_mode=pol,
                    tag=tag,
                    source_folder=source_folder,
                )
                fe.align_info = st['info']
                fe.align_params = (st['info'] or {}).get(
                    'params', dict(DEFAULT_ALIGN_PARAMS))
                # 実際に読んだ生データのパスをサンプル順で持たせる
                fe.raw_paths = [str(p) for _, p, _ in st['loaded_files']]
                if existing is not None:
                    _i = host._file_entries.index(existing)
                    host._file_entries[_i] = fe
                    for _j in range(host._cmb_files.count()):
                        if host._cmb_files.itemData(_j) == tag:
                            host._cmb_files.setItemText(_j, fe.label)
                            break
                    _purged = self._al_purge_annotation_state(pol, tag)
                    replaced.append(f"  replaced {tag} [{pol}]"
                                    + (f" / {_purged}" if _purged else ""))
                else:
                    host._file_entries.append(fe)
                    host._cmb_files.addItem(fe.label, tag)
                st['sent'] = True
                registered.append(
                    f"  {tag} [{pol}]  M+0 {len(fe.df)} / "
                    f"{len(st['aligned_df'])} features total")
            except Exception as e:
                if _bumped:
                    try:
                        host._file_counter -= 1
                    except Exception:
                        pass
                import traceback
                traceback.print_exc()
                failed.append(f"  {pol}: {type(e).__name__}: {e}")
        if registered:
            try:
                host._cmb_files.setCurrentIndex(host._cmb_files.count() - 1)
            except Exception:
                pass
        try:
            self._refresh_data_summary()
        except Exception:
            pass
        self._al_update_pol_summary()
        if replaced:
            log.info("replaced:\n" + "\n".join(replaced))
        if failed and not registered:
            QMessageBox.critical(
                self, "Register failed", "\n".join(failed))
        elif failed:
            QMessageBox.warning(
                self, "Partially registered",
                "Registered:\n" + "\n".join(registered)
                + "\n\nFailed:\n" + "\n".join(failed))
        else:
            QMessageBox.information(
                self, "Registered",
                "Registered in the Annotation tab:\n\n" + "\n".join(registered)
                + "\n\nContinue the analysis in the Annotation tab.")

    def _al_save_txt(self):
        if self._al_aligned_df is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save aligned TXT",
            timestamped_filename("aligned_peaks.txt"),
            "Text files (*.txt);;All (*)")
        if not path:
            return
        ion_mode = 'pos' if self._al_rb_pos.isChecked() else 'neg'
        icol = ('area' if self._al_rb_area.isChecked() else 'height')
        try:
            n = write_msdial_compatible_txt(
                self._al_aligned_df, self._al_sample_names, path,
                ion_mode=ion_mode, intensity_col=icol,
                monoisotopic_only=True)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", f"{type(e).__name__}: {e}")
            return
        QMessageBox.information(
            self, "Saved", f"Wrote {n} M+0 features to:\n{path}")

    def _al_save_csv(self):
        if self._al_aligned_df is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save extended CSV",
            timestamped_filename("extended_peaks.csv"),
            "CSV files (*.csv);;All (*)")
        if not path:
            return
        try:
            n = write_extended_csv(
                self._al_aligned_df, self._al_sample_names, path)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", f"{type(e).__name__}: {e}")
            return
        QMessageBox.information(
            self, "Saved", f"Wrote {n} rows (all features, isotopes included) to:\n{path}")

    def _build_quantification_tab(self):
        """Placeholder tab. Absolute quantification will be added in a
        future release; use the LipidQuant relay export (Utility tab) for now."""
        tab = QWidget()
        lay = QVBoxLayout(tab)
        lay.setContentsMargins(16, 16, 16, 16)
        lay.setSpacing(8)
        lay.addWidget(QLabel("<h2>Quantification</h2>"))
        lbl = QLabel(
            "To be added in a future release.<br><br>"
            "For now, export the annotated peak table with "
            "<b>Utility &gt; Export &gt; Quant ion table</b> or the "
            "<b>LipidQuant relay</b> sub-tab and perform quantification "
            "externally.")
        lbl.setStyleSheet("color:#666;")
        lbl.setWordWrap(True)
        lay.addWidget(lbl)
        lay.addStretch()
        return tab
    def _per_sample_entry(self):
        # 今表示中の per-sample データと match_df は PreviewRTDialog
        # 自身が self.fe / self._match_df として保持している。まずそれを使い、
        # 無ければ host(MainWindow)側を探す。
        fe = getattr(self, "fe", None)
        if isinstance(fe, PerSampleFileEntry):
            return fe
        host = getattr(self, "_host", None)
        for f in (getattr(host, "_file_entries", None) or []):
            if isinstance(f, PerSampleFileEntry):
                return f
        return None

    def _entry_for_mode(self, mode: str):
        """指定極性の PerSampleFileEntry と match_df を返す。

        戻り値: (fe, mdf) 。見つからなければ (None, None)。
        """
        fe = None
        for f in (self._all_entries or []):
            if (isinstance(f, PerSampleFileEntry)
                    and getattr(f, 'ion_mode', None) == mode):
                fe = f
                break
        if fe is None:
            host = getattr(self, '_host', None)
            for f in (getattr(host, '_file_entries', None) or []):
                if (isinstance(f, PerSampleFileEntry)
                        and getattr(f, 'ion_mode', None) == mode):
                    fe = f
                    break
        if fe is None:
            return None, None
        mdf = self._get_match_df_for_mode(mode)
        if mdf is None or getattr(mdf, 'empty', True):
            mdf = (getattr(getattr(self, '_host', None),
                           '_preview_match_df_by_tag', {}) or {}).get(fe.tag)
        return fe, mdf

    def _nontarget_table(self, fe, mdf, values=None, attribution=None):
        """1 極性ぶんのノンターゲット表を組み立てる。

        アライメント後の M+0 feature を全部並べ、アノテーションが
        付いたものには化合物名・クラス・アダクトと final_status を、
        付かなかったものには 'not annotated' を入れる。

        1 つの feature に複数の候補が当たることがある(multi-label)ので、
        kept のものを代表にし、残りは other_candidates に並べる。

        values に fe_intensity_tables() が返す DataFrame を渡すと、
        サンプル列の値をそれで置き換える(height / area の出し分け)。
        None なら従来どおり fe.df のサンプル列を使う。

        attribution に {(obs_rt, obs_mz): merged_entry} を渡すと、
        confidence / sigma / decided_by / competitors の 4 列を足す。
        ノンターゲット解析では「帰属はしておいて信頼度で絞る」運用を
        するため。
        """
        import pandas as _pd
        base = getattr(fe, 'df', None)
        if base is None or getattr(base, 'empty', True):
            return None
        try:
            smp = list(fe.sample_columns())
        except Exception:
            smp = []
        _rt_col = getattr(fe, 'rt_col', None)
        _mz_col = getattr(fe, 'mz_col', None)
        out = _pd.DataFrame()
        out['Feature ID'] = (base['Alignment ID'].values
                             if 'Alignment ID' in base.columns
                             else np.arange(len(base)))
        out['RT (min)'] = (base[_rt_col].values if _rt_col in base.columns
                           else np.full(len(base), np.nan))
        out['m/z'] = (base[_mz_col].values if _mz_col in base.columns
                      else np.full(len(base), np.nan))
        # RT × m/z のユニーク ID。多変量解析で pos/neg を積んでも
        # 衝突しないよう極性を先頭に付ける。丸めの結果ぶつかった場合は
        # #2, #3 … を足して、結合キーとして必ず一意になるようにする。
        _mode = str(getattr(fe, 'ion_mode', '') or 'na')
        _uid, _seen = [], {}
        for _rt_v, _mz_v in zip(out['RT (min)'], out['m/z']):
            _u = feature_uid(_mode, _rt_v, _mz_v)
            _c = _seen.get(_u, 0) + 1
            _seen[_u] = _c
            _uid.append(_u if _c == 1 else f"{_u}#{_c}")
        out.insert(0, 'Feature UID', _uid)
        # feature index -> 当たったライブラリ行
        ann: dict = {}
        if (mdf is not None and not getattr(mdf, 'empty', True)
                and 'matched_idx' in mdf.columns):
            try:
                _hit = mdf[mdf['matched_idx'] >= 0]
            except Exception:
                _hit = mdf.iloc[0:0]
            for _, r in _hit.iterrows():
                try:
                    i = int(r['matched_idx'])
                except (TypeError, ValueError):
                    continue
                if 0 <= i < len(base):
                    ann.setdefault(i, []).append(r)
        cmpd, cls_, add, sts, why, ncand, alt = (
            [], [], [], [], [], [], [])
        # 信頼度 4 列
        conf_l, conf_s, conf_by, conf_vs = [], [], [], []
        try:
            _sig_thr = float(self._spin_sigma.value())
        except Exception:
            _sig_thr = 3.0
        for i in range(len(base)):
            hits = ann.get(i) or []
            if not hits:
                cmpd.append(''); cls_.append(''); add.append('')
                sts.append(NONTARGET_UNANNOTATED); why.append('')
                ncand.append(0); alt.append('')
                conf_l.append(''); conf_s.append(float('nan'))
                conf_by.append(''); conf_vs.append('')
                continue
            # Compound / Class / Adduct は **kept のときだけ** 入れる。
            # フィルタで落ちた帰属が「その feature の正体」として
            # 残ってしまうのを防ぐ。落ちた候補は candidates に残す。
            kept = next((h for h in hits
                         if str(h.get('final_status', '')) == STATUS_KEPT),
                        None)
            if kept is not None:
                cmpd.append(str(kept.get('compound', '') or ''))
                cls_.append(str(kept.get('lipid_class', '') or ''))
                add.append(str(kept.get('adduct', '') or ''))
                sts.append(STATUS_KEPT)
                why.append('')
            else:
                cmpd.append(''); cls_.append(''); add.append('')
                sts.append(str(hits[0].get('final_status', '') or ''))
                # 代表(m/z 誤差が最小のもの)がどのフィルタで落ちたか。
                # 候補ごとの顛末は candidates 側で追える。
                why.append(rejection_reasons(hits[0]))
            # その座標の帰属から信頼度をまとめる
            _entry = None
            if attribution:
                try:
                    _entry = attribution.get(
                        (round(float(hits[0]['obs_rt']), 6),
                         round(float(hits[0]['obs_mz']), 6)))
                except Exception:
                    _entry = None
            _lab, _sg, _by, _vs = annotation_confidence(
                _entry, kept if kept is not None else hits[0], _sig_thr)
            conf_l.append(_lab)
            try:
                conf_s.append(float(_sg))
            except (TypeError, ValueError):
                conf_s.append(float('nan'))
            conf_by.append(_by); conf_vs.append(_vs)
            ncand.append(len(hits))
            # 候補は「残り」ではなく **全部** 並べる。
            # 化合物名の列から外したぶん、ここに情報を残しておく。
            alt.append('; '.join(
                f"{h.get('compound', '')}[{h.get('lipid_class', '')}]"
                f":{h.get('final_status', '')}" for h in hits[:8])
                + (' ...' if len(hits) > 8 else ''))
        out['Compound'] = cmpd
        out['Class'] = cls_
        out['Adduct'] = add
        out['final_status'] = sts
        out['rejected_by'] = why
        out['n_candidates'] = ncand
        out['candidates'] = alt
        # 信頼度(ノンターゲット解析で絞り込むため)
        out['confidence'] = conf_l
        out['sigma'] = conf_s
        out['decided_by'] = conf_by
        out['competitors'] = conf_vs
        # どのセルがギャップフィリング由来かを出す。
        _fl = fe_filled_table(fe)
        if _fl is not None:
            _fa = _fl[smp].to_numpy(dtype=bool)
            out['n_filled'] = _fa.sum(axis=1).astype(int)
            out['filled_samples'] = [
                '; '.join(s for s, f in zip(smp, _fa[i]) if f)
                for i in range(len(base))]
        _vals = values if values is not None else base
        for s in smp:
            if s in _vals.columns:
                out[s] = np.asarray(_vals[s].values)
        return out

    def _export_nontarget_features(self):
        """アノテーションに通らなかった feature も含めて出力する。

        ノンターゲット解析の入口。定量とは独立なので、
        [Apply Quantification] を押していなくても使える。

        height と area の両方を書き、シート名に指標を明記する
        (pos_height / pos_area / neg_height / neg_area)。area が出るのは
        Detection で Height を選んで走らせたときだけ — Area を選ぶと
        area_s* は書かれず height 側も残らないため、その場合は 1 指標のみ。
        """
        import pandas as _pd
        sheets, info, notes = {}, [], []
        # 信頼度列のもとになる帰属（両モードぶん）を 1 回だけ作る
        _attr_by_mode: dict = {}
        try:
            for (_m, _rt, _mz), _e in (
                    self._build_merged_attribution() or {}).items():
                _attr_by_mode.setdefault(str(_m), {})[
                    (round(float(_rt), 6), round(float(_mz), 6))] = _e
        except Exception as _e:
            log.warning(f"attribution build failed: {_e}")
        for mode in ('pos', 'neg'):
            fe, mdf = self._entry_for_mode(mode)
            if fe is None:
                continue
            try:
                tabs = fe_intensity_tables(fe)
            except Exception as e:
                log.warning(f"intensity tables failed ({mode}): {e}")
                tabs = []
            if not tabs:
                continue
            for metric, vals in tabs:
                try:
                    tbl = self._nontarget_table(
                        fe, mdf, values=vals,
                        attribution=_attr_by_mode.get(mode))
                except Exception as e:
                    import traceback
                    QMessageBox.critical(
                        self, "Export failed",
                        f"{mode}/{metric}: {type(e).__name__}: {e}\n\n"
                        + traceback.format_exc()[-1200:])
                    return
                if tbl is None or tbl.empty:
                    continue
                sheets[f"{mode}_{metric}"] = tbl
                if metric == tabs[0][0]:
                    _ann = int((tbl['final_status']
                                != NONTARGET_UNANNOTATED).sum())
                    _kept = int((tbl['final_status'] == STATUS_KEPT).sum())
                    info.append(f"{mode}: {len(tbl)} features, "
                                f"{_ann} annotated ({_kept} kept)")
            if len(tabs) == 1:
                notes.append(
                    f"{mode}: only '{tabs[0][0]}' is available — the "
                    f"alignment was run with that intensity column, so the "
                    f"other one was never stored.")
        if not sheets:
            QMessageBox.information(
                self, "No features",
                "No aligned features were found. Load per-sample data "
                "(Session) or send an alignment from the Detection tab "
                "first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export all features (non-target)",
            timestamped_filename("nontarget_features.xlsx"),
            "Excel (*.xlsx);;CSV (*.csv)")
        if not path:
            return
        try:
            if path.lower().endswith(".csv"):
                # 極性でサンプル列が違い得るので、混ぜずに分ける
                written = []
                for name, tbl in sheets.items():
                    _p = path[:-4] + f"_{name}.csv"
                    tbl.to_csv(_p, index=False)
                    written.append(_p)
            else:
                with _pd.ExcelWriter(path) as _xw:
                    for name, tbl in sheets.items():
                        tbl.to_excel(_xw, sheet_name=name, index=False)
                written = [path]
        except Exception as e:
            QMessageBox.critical(self, "Export failed",
                                 f"{type(e).__name__}: {e}")
            return
        QMessageBox.information(
            self, "Exported",
            "\n".join(info)
            + "\n\nSheets: " + ", ".join(sheets)
            + ("\n\n" + "\n".join(notes) if notes else "")
            + "\n\nWrote:\n" + "\n".join(written))

    def _build_utility_tab(self):
        """Utility タブを構築する。図エクスポート、Joint Plots、テーブル、
        パラメータ JSON 保存・読込などの二次機能を集約する。
        """
        utility_tab = QWidget()
        outer = QVBoxLayout(utility_tab)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # スクロール可能なコンテナ
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer.addWidget(scroll)

        content = QWidget()
        scroll.setWidget(content)
        lay = QVBoxLayout(content)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(10)

        # ── 階層構造化 ──
        # ■ Reports
        #   ➤ Filtering Results
        #       Figures (PNG):  ③ | ④ | ⑤ | ⑥ | Summary All
        #       Tables (xlsx):  ③ | ④ | ⑤ | ⑥ | Summary All
        #   ➤ Scatter Plots
        #       □ Include histograms
        #       Current view:  [Pos] [Neg] [Merged]
        #       Curated view:  [Pos] [Neg] [Merged]   (TBD)
        #   ➤ Heatmap
        #       [IS adduct ion patterns] [Observed adduct ion patterns]
        # ■ Parameters
        #       [Save JSON] [Load JSON] [Reset]

        SECTION_STYLE = (
            "font-size: 14px; font-weight: bold;"
            " color: #1A1A1A; padding: 8px 0 2px 0;")
        SUB_STYLE = (
            "font-size: 11px; color: #555; padding-top: 2px;")

        # ============ ■ Reports ============================
        hdr_reports = QLabel("■  Reports")
        hdr_reports.setStyleSheet(SECTION_STYLE)
        lay.addWidget(hdr_reports)

        # ── ➤ Filtering Results ─────────────────────────────
        # Figures 行は削除(全ステップに自然な図がないため
        # Tables (xlsx) のみに集約。⑤ Adduct は dialog 内で PNG/xlsx 両方
        # 保存可能なため Tables 側からアクセス可能)。
        gb_filter = QGroupBox("➤  Filtering Results")
        gb_filter_lay = QVBoxLayout(gb_filter)

        tbl_hdr = QLabel("Tables (xlsx):  "
                         "(each opens a preview dialog with Export…)")
        tbl_hdr.setStyleSheet(SUB_STYLE)
        gb_filter_lay.addWidget(tbl_hdr)
        tbl_row = QHBoxLayout()
        btn_tbl_is = QPushButton("③ IS Filter")
        btn_tbl_is.setToolTip(
            "Open ③ IS-rejected spots table dialog (xlsx export inside).")
        btn_tbl_is.clicked.connect(self._open_is_filter_report)
        tbl_row.addWidget(btn_tbl_is)
        btn_tbl_rt = QPushButton("④ RT outlier")
        btn_tbl_rt.setToolTip(
            "Open ④ RT outlier table dialog (xlsx export inside).")
        btn_tbl_rt.clicked.connect(self._open_rt_outlier_report)
        tbl_row.addWidget(btn_tbl_rt)
        btn_tbl_ad = QPushButton("⑤ Adduct Ion Filter")
        btn_tbl_ad.setToolTip(
            "Open ⑤ Adduct attribution report dialog (xlsx export inside).")
        btn_tbl_ad.clicked.connect(self._open_adduct_report)
        tbl_row.addWidget(btn_tbl_ad)
        btn_tbl_co = QPushButton("⑥ Coherence Filter")
        btn_tbl_co.setToolTip(
            "Open ⑥ Coherence attribution dialog (xlsx export inside).")
        btn_tbl_co.clicked.connect(self._open_coherence_report)
        tbl_row.addWidget(btn_tbl_co)
        btn_tbl_sum = QPushButton("Summary All")
        btn_tbl_sum.setToolTip(
            "Open combined summary dialog (xlsx export inside).")
        btn_tbl_sum.clicked.connect(self._open_full_report)
        tbl_row.addWidget(btn_tbl_sum)
        tbl_row.addStretch()
        gb_filter_lay.addLayout(tbl_row)
        lay.addWidget(gb_filter)

        # ── ➤ Scatter Plots ─────────────────────────────────
        gb_scatter = QGroupBox("➤  Scatter Plots")
        gb_scatter_lay = QVBoxLayout(gb_scatter)

        hist_row = QHBoxLayout()
        self._chk_include_hist = QCheckBox(
            "Include histograms (marginal RT / m/z)")
        self._chk_include_hist.setToolTip(
            "If checked, scatter previews include marginal histograms\n"
            "(joint plot style). Affects subsequent Pos/Neg/Merged buttons.")
        hist_row.addWidget(self._chk_include_hist)
        hist_row.addStretch()
        gb_scatter_lay.addLayout(hist_row)

        cv_hdr = QLabel("Current view:")
        cv_hdr.setStyleSheet(SUB_STYLE)
        gb_scatter_lay.addWidget(cv_hdr)
        cv_row = QHBoxLayout()
        btn_cv_pos = QPushButton("Pos")
        btn_cv_pos.setToolTip(
            "Open Pos scatter (with histograms if checkbox is on).")
        btn_cv_pos.clicked.connect(
            lambda: self._utility_scatter_action('current', 'pos'))
        cv_row.addWidget(btn_cv_pos)
        btn_cv_neg = QPushButton("Neg")
        btn_cv_neg.setToolTip(
            "Open Neg scatter (with histograms if checkbox is on).")
        btn_cv_neg.clicked.connect(
            lambda: self._utility_scatter_action('current', 'neg'))
        cv_row.addWidget(btn_cv_neg)
        btn_cv_merged = QPushButton("Merged")
        btn_cv_merged.setToolTip(
            "Save the combined Pos+Neg scatter.")
        btn_cv_merged.clicked.connect(
            lambda: self._utility_scatter_action('current', 'merged'))
        cv_row.addWidget(btn_cv_merged)
        cv_row.addStretch()
        gb_scatter_lay.addLayout(cv_row)

        # Curated view 実装済み。⑦ Select Quant Ion で
        # 選択されたクラス×モードのスポットのみを描画する。
        curated_hdr = QLabel(
            "Curated view:  (uses ⑦ Select Quant Ion choices)")
        curated_hdr.setStyleSheet(SUB_STYLE)
        gb_scatter_lay.addWidget(curated_hdr)
        curated_row = QHBoxLayout()
        btn_curv_pos = QPushButton("Pos")
        btn_curv_pos.setToolTip(
            "Pos scatter showing only classes selected as 'pos'\n"
            "in ⑦ Select Quant Ion.")
        btn_curv_pos.clicked.connect(
            lambda: self._utility_scatter_action('curated', 'pos'))
        curated_row.addWidget(btn_curv_pos)
        btn_curv_neg = QPushButton("Neg")
        btn_curv_neg.setToolTip(
            "Neg scatter showing only classes selected as 'neg'\n"
            "in ⑦ Select Quant Ion.")
        btn_curv_neg.clicked.connect(
            lambda: self._utility_scatter_action('curated', 'neg'))
        curated_row.addWidget(btn_curv_neg)
        btn_curv_merged = QPushButton("Merged")
        btn_curv_merged.setToolTip(
            "Merged scatter combining pos- and neg-selected classes\n"
            "(per ⑦ Select Quant Ion choice). Skip classes excluded.")
        btn_curv_merged.clicked.connect(
            lambda: self._utility_scatter_action('curated', 'merged'))
        curated_row.addWidget(btn_curv_merged)
        curated_row.addStretch()
        gb_scatter_lay.addLayout(curated_row)
        lay.addWidget(gb_scatter)

        # ── ➤ Heatmap ───────────────────────────────────────
        gb_heatmap = QGroupBox("➤  Heatmap")
        gb_heatmap_lay = QHBoxLayout(gb_heatmap)
        btn_hm_is = QPushButton("IS adduct ion patterns")
        btn_hm_is.setToolTip(
            "Save the heatmap of CURATED + runtime IS-derived adduct "
            "patterns (reference).")
        btn_hm_is.clicked.connect(self._utility_save_patterns_png)
        gb_heatmap_lay.addWidget(btn_hm_is)
        btn_hm_obs = QPushButton("Observed adduct ion patterns")
        btn_hm_obs.setToolTip(
            "Save the heatmap of observed adduct patterns at non-IS "
            "conflict spots.")
        btn_hm_obs.clicked.connect(self._utility_save_report_png)
        gb_heatmap_lay.addWidget(btn_hm_obs)
        gb_heatmap_lay.addStretch()
        lay.addWidget(gb_heatmap)

        # ── ➤ Quant Ion Table ──────────────────────────────
        # ⑦ Select Quant Ion で選択されたクラス×モードの kept スポットの
        # サンプル別の値を 1 つの xlsx に出力する。
        # 出ていたのは「エリア」ではなく Detection で選んだ指標
        # だった。名前を直し、height と area を別シートで書く。
        gb_qarea = QGroupBox("➤  Quant Ion Table")
        gb_qarea_lay = QVBoxLayout(gb_qarea)
        qarea_desc = QLabel(
            "Export sample-wise intensities for spots selected in "
            "⑦ Select Quant Ion (final_status = kept only). "
            "Height and area are written to separate sheets.")
        qarea_desc.setStyleSheet(SUB_STYLE)
        qarea_desc.setWordWrap(True)
        gb_qarea_lay.addWidget(qarea_desc)
        qarea_row = QHBoxLayout()
        btn_qarea = QPushButton("Export quant ion table (xlsx)…")
        btn_qarea.setToolTip(
            "Build an xlsx with one sheet per intensity metric:\n"
            "  Class | Compound | Adduct | Mode | RT (min) | m/z |\n"
            "  <pos sample columns…> | <neg sample columns…>\n"
            "Only spots with final_status = kept for the chosen quant\n"
            "ion mode (per class) are included.\n"
            "The 'area' sheet exists only when the alignment was run\n"
            "with Height (otherwise area is not stored).")
        btn_qarea.clicked.connect(self._utility_export_quant_areas)
        qarea_row.addWidget(btn_qarea)
        qarea_row.addStretch()
        gb_qarea_lay.addLayout(qarea_row)
        # アノテーションに通らなかった feature も含めた全 feature 表。
        nt_row = QHBoxLayout()
        self._btn_export_nontarget = QPushButton(
            "Export all features (non-target)…")
        self._btn_export_nontarget.setToolTip(
            "Write every aligned M+0 feature — including the ones that "
            "were never annotated — with a unique RT/mz id, the accepted "
            "annotation (kept only), every candidate with its status and, "
            "for rejected ones, which filter removed it. One sheet per "
            "polarity and intensity metric (pos_height, pos_area, …).")
        self._btn_export_nontarget.clicked.connect(
            self._export_nontarget_features)
        nt_row.addWidget(self._btn_export_nontarget)
        nt_row.addStretch()
        gb_qarea_lay.addLayout(nt_row)
        lay.addWidget(gb_qarea)

        # Parameters セクションは Session タブに統合済みのため削除。
        # Reset は Advanced ダイアログ側に残置。

        lay.addStretch()
        # Utility 本体は container の 1 サブタブになる。保持のみ。
        self._utility_reports_tab = utility_tab

        # ── Tasks タブ ──
        # MainWindow と並存。AP Tasks タブは MainWindow の操作の "second view"。
        # 全アクションは host (MainWindow) のメソッドに delegate し、表示は
        # _refresh_tasks_tab_from_host() で MainWindow の状態から再構築。
        tasks_tab = QWidget()
        tasks_lay = QVBoxLayout(tasks_tab)
        tasks_lay.setContentsMargins(8, 8, 8, 8)
        tasks_lay.setSpacing(6)

        # Output Directory (AP 側のミラー)
        ap_out_box = QGroupBox("Output Directory")
        ap_ob = QHBoxLayout(ap_out_box)
        self._ap_lbl_outdir = QLabel("(default: same as input file)")
        ap_ob.addWidget(self._ap_lbl_outdir, 1)
        ap_btn_outdir = QPushButton("Choose …")
        ap_btn_outdir.clicked.connect(self._ap_choose_outdir)
        ap_ob.addWidget(ap_btn_outdir)
        tasks_lay.addWidget(ap_out_box)

        # Tasks 表(AP 側のミラー、read-only display)
        ap_task_box = QGroupBox("Tasks")
        ap_tb = QVBoxLayout(ap_task_box)
        self._ap_task_table = QTableWidget(0, 5)
        self._ap_task_table.setHorizontalHeaderLabels(
            ["Class", "RT start", "RT end", "Ion Mode", "File"])
        self._ap_task_table.horizontalHeader().setStretchLastSection(True)
        # Allow user to edit cells; on cellChanged we propagate to host
        self._ap_task_table.cellChanged.connect(self._ap_on_task_cell_changed)
        ap_tb.addWidget(self._ap_task_table)

        ap_th = QHBoxLayout()
        for txt, slot in [
            ("Add Row",              self._ap_add_task_row),
            ("Delete Row",           self._ap_delete_task_row),
            ("Paste from Clipboard", self._ap_paste_tasks),
            ("Copy All Tasks",       self._ap_copy_all_tasks),
        ]:
            b = QPushButton(txt)
            b.clicked.connect(slot)
            ap_th.addWidget(b)
        ap_tb.addLayout(ap_th)
        tasks_lay.addWidget(ap_task_box)

        # Run Export (AP 側)
        self._ap_btn_run = QPushButton("▶  Run Export")
        self._ap_btn_run.setStyleSheet("font-size:14px;padding:8px;")
        self._ap_btn_run.clicked.connect(self._ap_run)
        tasks_lay.addWidget(self._ap_btn_run)

        self._ap_lbl_log = QLabel("")
        tasks_lay.addWidget(self._ap_lbl_log)
        tasks_lay.addStretch()
        # Tasks(LipidQuant relay TXT 出力)は Utility 内サブタブへ。
        # NNLS 定量が確定するまで外部検算の手段として機能は残す。
        self._relay_export_tab = tasks_tab

        # ── Quantification タブ(プレースホルダ)──────────────────
        quant_tab = self._build_quantification_tab()
        self._tabs.addTab(quant_tab, TAB_LABELS['quant'])

        # ── Utility container(サブタブ構成)─────────────────
        # Session / Reports・Export / LipidQuant relay export を 1 タブに束ねる。
        # 構築順の都合で Utility は最後に addTab するので、結果の並びは
        #   Annotation | Quantification | Utility
        # になる(fix19 で Peak Detection & Alignment が先頭に入る)。
        util_container = QWidget()
        _uc_lay = QVBoxLayout(util_container)
        _uc_lay.setContentsMargins(0, 0, 0, 0)
        _uc_lay.setSpacing(0)
        self._utility_subtabs = QTabWidget()
        try:
            self._utility_subtabs.addTab(
                self._session_tab_widget, SUBTAB_LABELS['session'])
        except Exception as _e:
            log.warning(f"session subtab failed: {_e}")
        try:
            self._utility_subtabs.addTab(
                self._build_external_import_tab(), SUBTAB_LABELS['import'])
        except Exception as _e:
            import traceback
            traceback.print_exc()
            log.warning(f"import subtab failed: {_e}")
        try:
            self._utility_subtabs.addTab(
                self._utility_reports_tab, SUBTAB_LABELS['reports'])
        except Exception as _e:
            log.warning(f"reports subtab failed: {_e}")
        try:
            self._utility_subtabs.addTab(
                self._relay_export_tab, SUBTAB_LABELS['relay'])
        except Exception as _e:
            log.warning(f"relay subtab failed: {_e}")
        _uc_lay.addWidget(self._utility_subtabs)
        self._tabs.addTab(util_container, TAB_LABELS['util'])

        # ── ① Detection タブを先頭に挿入 ──
        # 最終的な並び:
        #   Detection | Annotation | Quantification | Utility
        try:
            self._tabs.insertTab(0, self._build_alignment_tab(),
                                 TAB_LABELS['align'])
        except Exception as _e:
            import traceback
            traceback.print_exc()
            log.warning(f"alignment tab build failed: {_e}")

        # 起動時のデフォルトタブは ① Detection。
        # 解析は raw mzML のピーク検出から始めるのが通常の流れなので、
        # 最初に開くタブもそこに合わせる。
        # ① の構築に失敗した場合だけ Annotation に落ちる。
        try:
            _want = [TAB_LABELS['align'], TAB_LABELS['annot']]
            for _label in _want:
                _hit = next((_i for _i in range(self._tabs.count())
                             if self._tabs.tabText(_i) == _label), None)
                if _hit is not None:
                    self._tabs.setCurrentIndex(_hit)
                    break
        except Exception:
            pass

        # タブ切替時に relay export サブタブの内容を host から refresh。
        # Tasks はサブタブになったので、トップレベル(Utility 選択)と
        # サブタブ(relay 選択)の両方を監視する。
        try:
            self._tabs.currentChanged.connect(self._on_tabs_current_changed)
        except Exception:
            pass
        try:
            self._utility_subtabs.currentChanged.connect(
                self._on_utility_subtab_changed)
        except Exception:
            pass

    # ── Utility: Export Current View ───────────────────────────────
    def _utility_save_scatter(self, mode: str):
        """Save the main scatter as PNG. single-mode dialog,
        both Pos and Neg buttons save the same scatter."""
        path, _ = QFileDialog.getSaveFileName(
            self, f"Save {mode} scatter PNG",
            timestamped_filename(f"scatter_{mode}.png"), "PNG (*.png)")
        if not path:
            return
        try:
            self._fig.savefig(path, dpi=150, bbox_inches='tight')
            QMessageBox.information(self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))

    # ── ヘルパ ──────────────────────────────────
    def _utility_scatter_action(self, source: str, mode: str):
        """Scatter Plots セクションの統一ハンドラ。

        source: 'current' = 現在の描画(全マッチ),
                'curated' = ⑦ Select Quant Ion で選択されたクラス×モードのみ
        mode:   'pos' | 'neg' | 'merged'
        全モードを JointPlotDialog で統一処理。クラス色アノテーション
        を反映。Merged は両モードを同一軸に重ね合わせ。
          - Histogram ON  → Joint Plot (散布図 + 周辺ヒストグラム)
          - Histogram OFF → Scatter Preview (散布図のみ)
        source='curated' で _quant_ion_choices を使ってクラスを絞る。
        """
        is_curated = (source == 'curated')
        if is_curated:
            qi = getattr(self, '_quant_ion_choices', None) or {}
            if not qi or not any(v in ('pos', 'neg') for v in qi.values()):
                QMessageBox.warning(
                    self, "No Quant Ion Selection",
                    "Please use '⑦ Select Quant Ion…' first to choose\n"
                    "pos / neg mode per lipid class before opening\n"
                    "the Curated view.")
                return
        include_hist = (hasattr(self, '_chk_include_hist')
                        and self._chk_include_hist.isChecked())
        dlg = JointPlotDialog(
            host=self, parent=self,
            include_histograms=include_hist,
            initial_mode=mode,
            curated=is_curated)
        dlg.exec()

    # ── Utility: Joint Plots ──────────────────────────────────────
    # ── Utility: Adduct Pattern Figures ───────────────────────────
    def _utility_save_patterns_png(self):
        """IS adduct pattern heatmap をプレビューしてから保存。
        旧: 直接 File Save Dialog → 修正後: HeatmapPreviewDialog → ユーザーが保存。"""
        # 統合パターンを構築
        view: dict = {}
        runtime = self._runtime_is_fingerprints or {}
        source = dict(runtime)
        for cls, fp in CURATED_ADDUCT_FINGERPRINTS.items():
            source[cls] = fp
        if not source:
            QMessageBox.warning(self, "No data",
                                "No fingerprints available. Run IS Filter first.")
            return
        for cls, cls_fp in source.items():
            view[cls] = {}
            for adu in ADDUCT_SET_V54:
                spec = cls_fp.get(adu, {}) or {}
                view[cls][adu] = {
                    'expected': bool(spec.get('expected', False)),
                    'value': float(spec.get('weight', 0.0) or 0.0),
                }
        dlg = HeatmapPreviewDialog(
            host=self, view=view,
            title="IS Adduct Pattern Reference",
            default_filename="is_adduct_pattern_heatmap.png",
            parent=self)
        dlg.exec()

    def _utility_save_report_png(self):
        """Adduct observation report heatmap を
        プレビューしてから保存。
        AdductReportDialog と同じ強度グラデーション(行内 max 正規化、
        [0.25, 1.0] の範囲にマップ)を Utility 経由でも適用。
        両モードの _adduct_attribution を統合(active mode のみだと
        他モードの帰属スポットが抜け落ちる問題を修正)。"""
        # active + 他モードの attribution を統合
        active_attr = self._adduct_attribution or {}
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        other_attr = snap.get('_adduct_attribution') or {}
        if not active_attr and not other_attr:
            QMessageBox.warning(
                self, "No data",
                "No Adduct Filter results. Run Adduct Ion Filter first.")
            return

        # mode 情報を付与した統合リストを構築
        merged_items = []
        for coord, result in active_attr.items():
            merged_items.append((self._active_mode, coord, result))
        for coord, result in other_attr.items():
            merged_items.append((other_mode, coord, result))
        # mode, RT, m/z でソート(可読性のため)
        merged_items.sort(key=lambda t: (t[0], float(t[1][0]),
                                         float(t[1][1])))

        # AdductReportDialog の view 構築ロジックをミラー
        view: dict = {}
        for mode_v, coord, result in merged_items:
            rt, mz = coord
            obs = result.get('observations', {}) or {}
            obs_int = result.get('observation_intensities', {}) or {}
            winner = result.get('winner') or ('skip' if result.get('skipped')
                                              else 'undecided')
            row_label = f"[{mode_v}] {rt:.3f}/{mz:.4f} [{winner}]"
            # 行内 max で正規化(各 spot 単位でグラデーション)
            row_max = 0.0
            for adu in ADDUCT_SET_V54:
                if obs.get(adu, False):
                    v = float(obs_int.get(adu, 0.0) or 0.0)
                    if v > row_max:
                        row_max = v
            view[row_label] = {}
            for adu in ADDUCT_SET_V54:
                observed = bool(obs.get(adu, False))
                if observed:
                    if row_max > 0:
                        intensity = float(obs_int.get(adu, 0.0) or 0.0)
                        # 範囲 [0.25, 1.0] にマッピング(最弱でも 0.25 で見える)
                        value = 0.25 + 0.75 * (intensity / row_max)
                    else:
                        value = 1.0   # 強度情報無し → 一律塗りつぶし
                else:
                    value = 0.0
                view[row_label][adu] = {
                    'expected': observed,
                    'value': value,
                }
        dlg = HeatmapPreviewDialog(
            host=self, view=view,
            title="Adduct Observation Report",
            default_filename="adduct_observation_report.png",
            parent=self)
        dlg.exec()

    def _render_adduct_heatmap_png(self, view: dict, path_or_buf, title: str):
        """matplotlib で adduct heatmap を PNG 出力する。
        path_or_buf は str(ファイルパス)または file-like(BytesIO)を受け付ける。

        view: {row_label: {col_label: {'expected': bool, 'value': float}}}
        """
        rows = list(view.keys())
        cols = list(ADDUCT_SET_V54)
        n_rows = len(rows)
        n_cols = len(cols)
        if n_rows == 0:
            raise RuntimeError("No data rows to render")

        # matrix と expected mask を作成
        max_val = max(
            (float(view[r][c].get('value', 0.0)) for r in rows for c in cols),
            default=1.0)
        if max_val <= 0:
            max_val = 1.0
        mat = np.zeros((n_rows, n_cols))
        expected_mat = np.zeros((n_rows, n_cols), dtype=bool)
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                spec = view[r].get(c, {})
                mat[i, j] = float(spec.get('value', 0.0) or 0.0) / max_val
                expected_mat[i, j] = bool(spec.get('expected', False))

        fig_h = max(2.5, 0.35 * n_rows + 1.5)
        # colorbar の幅を見込んで横を少し広げる
        fig_w = max(7.5, 0.7 * n_cols + 2.8)
        fig = Figure(figsize=(fig_w, fig_h), constrained_layout=True)
        ax = fig.add_subplot(111)
        im = ax.imshow(mat, cmap='Blues', aspect='auto', vmin=0, vmax=1)
        # カラーバーを追加(右側)
        cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label('Relative intensity', fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        # 記号(○/×)
        for i in range(n_rows):
            for j in range(n_cols):
                text = "○" if expected_mat[i, j] else "×"
                color = ('white' if mat[i, j] > 0.55 and expected_mat[i, j]
                         else '#666')
                ax.text(j, i, text, ha='center', va='center',
                        color=color, fontsize=10)
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(cols, rotation=30, ha='right')
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(rows)
        ax.set_title(title)
        # file-like buffer 対応(format 明示)
        if isinstance(path_or_buf, str):
            fig.savefig(path_or_buf, dpi=150, bbox_inches='tight')
        else:
            fig.savefig(path_or_buf, format='png', dpi=150,
                        bbox_inches='tight')

    # ── Utility: Tables ───────────────────────────────────────────
    # ── Utility: Parameters ───────────────────────────────────────
    def _gather_current_params(self) -> dict:
        return {
            'version': '5.5.0',
            'ppm_tol': float(self._spin_ppm.value()),
            'void_rt': float(self._spin_void_rt.value()),
            'is_tol': float(self._spin_is_tol.value()),
            'adduct_ppm': float(self._spin_adduct_ppm.value()),
            'adduct_rt': float(self._spin_adduct_rt.value()),
            'adduct_int_floor': float(self._spin_adduct_int.value()),
            'adduct_vote_threshold': float(self._spin_adduct_vote.value()),
            'mad_k': float(self._spin_iqr.value()),
            'sigma_threshold': float(self._spin_sigma.value()),
            'loose_coherence': self._loose_coherence_enabled(),
            'two_series_classes':
                ', '.join(self._two_series_classes()),
            'sum_two_series': self._sum_two_series_enabled(),
        }

    def _apply_params(self, p: dict):
        try:
            if 'ppm_tol' in p:
                self._spin_ppm.setValue(float(p['ppm_tol']))
            if 'void_rt' in p:
                self._spin_void_rt.setValue(float(p['void_rt']))
            if 'is_tol' in p:
                self._spin_is_tol.setValue(float(p['is_tol']))
            if 'adduct_ppm' in p:
                self._spin_adduct_ppm.setValue(float(p['adduct_ppm']))
            if 'adduct_rt' in p:
                self._spin_adduct_rt.setValue(float(p['adduct_rt']))
            if 'adduct_int_floor' in p:
                self._spin_adduct_int.setValue(float(p['adduct_int_floor']))
            if 'adduct_vote_threshold' in p:
                self._spin_adduct_vote.setValue(
                    float(p['adduct_vote_threshold']))
            if 'mad_k' in p:
                self._spin_iqr.setValue(float(p['mad_k']))
            if 'sigma_threshold' in p:
                self._spin_sigma.setValue(float(p['sigma_threshold']))
            if 'loose_coherence' in p:
                self._cb_loose_coherence.setChecked(
                    bool(p['loose_coherence']))
            if 'two_series_classes' in p:
                _v = p['two_series_classes']
                if not isinstance(_v, str):
                    _v = ', '.join(str(x) for x in (_v or ()))
                self._ed_two_series.setText(_v)
            if 'sum_two_series' in p:
                self._cb_sum_two_series.setChecked(
                    bool(p['sum_two_series']))
        except Exception as e:
            log.warning(f"[Apply Params] failed: {e}")

    def _utility_reset_params(self):
        defaults = {
            'ppm_tol': 10.0, 'is_tol': 0.10,
            'adduct_ppm': 10.0, 'adduct_rt': 0.05, 'adduct_int_floor': 0.0,
            'mad_k': 3.0, 'sigma_threshold': 3.0,
            'loose_coherence': COH_LOOSE_ATTRIBUTION_DEFAULT,
            'two_series_classes':
                ', '.join(COH_TWO_SERIES_CLASSES_DEFAULT),
            'sum_two_series': QUANT_SUM_TWO_SERIES_DEFAULT,
        }
        self._apply_params(defaults)

    # ── 修正: 軸エイリアス再割当ヘルパー ─────────────────────
    def _refresh_axis_aliases(self):
        """active モードに従って _ax_sc / _ax_sc_other / _base_scatter を
        再割当する。Run All 内や _on_mode_radio_changed で使用。"""
        if not hasattr(self, '_ax_sc_pos') or not hasattr(self, '_ax_sc_neg'):
            return
        if self._active_mode == 'pos':
            self._ax_sc       = self._ax_sc_pos
            self._ax_sc_other = self._ax_sc_neg
            if hasattr(self, '_base_scatter_pos'):
                self._base_scatter = self._base_scatter_pos
        else:
            self._ax_sc       = self._ax_sc_neg
            self._ax_sc_other = self._ax_sc_pos
            if hasattr(self, '_base_scatter_neg'):
                self._base_scatter = self._base_scatter_neg

    # ── 修正: 各 Apply を両モード実行 ────────────────────────

        # active 軸に対応する annot を self._annot に
        if self._active_mode == 'pos':
            self._annot = getattr(self, '_annot_pos',
                                  getattr(self, '_annot', None))
        else:
            self._annot = getattr(self, '_annot_neg',
                                  getattr(self, '_annot', None))

    def _apply_to_all_modes(self, fn):
        """指定された step 関数を全ロード済みモードで実行する。

        現在 active なモードを保存し、他モードに切替えて実行、最後に元の
        モードに戻す。各 Apply ボタンが「特別指定しない限り両モードで実行」
        するために使用。
        """
        # ロード済みモードのリスト
        # pos → neg 固定順序(上プロット → 下プロット)。
        # ファイルの登録順に依存しない。
        loaded = set()
        for fe2 in (self._all_entries or []):
            m = getattr(fe2, 'ion_mode', None)
            if m:
                loaded.add(m)
        modes = [m for m in ('pos', 'neg') if m in loaded]
        if not modes:
            modes = [self._active_mode]

        original_mode = self._active_mode
        # iteration 中は redraw を抑制
        self._suppress_overlay_draw = True
        try:
            for run_mode in modes:
                if run_mode != self._active_mode:
                    # 現在の状態を snapshot 保存し、他モードに切替
                    self._snapshot_mode_state(self._active_mode)
                    if not self._load_mode_state(run_mode):
                        log.info(f"[apply all modes] skip {run_mode}: not loaded")
                        continue
                    self._active_mode = run_mode
                    self._refresh_axis_aliases()
                    # ラジオボタンも追従(signal を一時無効化)
                    try:
                        for rb in (self._mode_radio_pos,
                                    self._mode_radio_neg):
                            rb.blockSignals(True)
                        if run_mode == 'pos':
                            self._mode_radio_pos.setChecked(True)
                        else:
                            self._mode_radio_neg.setChecked(True)
                    except Exception:
                        pass
                    finally:
                        try:
                            for rb in (self._mode_radio_pos,
                                        self._mode_radio_neg):
                                rb.blockSignals(False)
                        except Exception:
                            pass
                # 実行(_draw_overlay は suppress により no-op)
                try:
                    fn()
                except Exception as exc:
                    log.warning(f"[apply all modes] {fn.__name__} on "
                          f"{run_mode} failed: {exc}")
        finally:
            # iteration 終了 → suppress 解除
            self._suppress_overlay_draw = False
        try:
            # 元のモードに戻す
            if self._active_mode != original_mode:
                self._snapshot_mode_state(self._active_mode)
                if self._load_mode_state(original_mode):
                    self._active_mode = original_mode
                    self._refresh_axis_aliases()
                    try:
                        for rb in (self._mode_radio_pos,
                                    self._mode_radio_neg):
                            rb.blockSignals(True)
                        if original_mode == 'pos':
                            self._mode_radio_pos.setChecked(True)
                        else:
                            self._mode_radio_neg.setChecked(True)
                    except Exception:
                        pass
                    finally:
                        try:
                            for rb in (self._mode_radio_pos,
                                        self._mode_radio_neg):
                                rb.blockSignals(False)
                        except Exception:
                            pass
            # 最後に 1 回だけ UI 更新
            classes = []
            if self._match_df is not None:
                try:
                    classes = self._match_df['lipid_class'].unique().tolist()
                except Exception:
                    pass
            try:
                self._rebuild_filters(classes)
                # 両モード(active + 他方)を描画。個別 step ボタン
                # 経由でも他方モードの neg プロットが空になる問題を修正。
                self._draw_both_overlays()
                if self._match_df is not None:
                    matched = self._match_df[
                        self._match_df['matched']]
                    self._rebuild_match_table(matched)
            except Exception as e:
                log.warning(f"[apply all modes] final UI update failed: {e}")
            # モード切替で active mode が変わった可能性があり、
            # active mode の _adduct_attribution に応じてボタン状態を再同期する。
            try:
                self._sync_adduct_button_state()
            except Exception:
                pass
        finally:
            self._suppress_overlay_draw = False

    # ── per-mode 状態管理 ──────────────────────────
    # 切替対象の状態変数名(クラス変数)。これに含まれる属性が mode 切替時に
    # snapshot/restore される。fingerprints などクラス横断のものは含めない。
    _PER_MODE_STATE_VARS = (
        '_rt', '_mz', '_sample_intensities', '_sample_col_names',
        '_match_df', '_matched_coords', '_reserved_coords',
        '_conflict_coords', '_conflict_map',
        '_outlier_coords', '_coherence_outliers',
        '_coherence_assignments', '_coherence_models',
        '_coherence_pairs', '_coherence_all_pairs',
        '_adduct_attribution', '_is_filter_choices', '_class_ref_rt',
        '_manual_winner_coords',
    )

    def _get_match_df_for_mode(self, mode: str):
        """指定 mode の match_df を返す。
        active mode なら self._match_df、それ以外は _mode_state[mode] から復元。
        見つからなければ None。"""
        try:
            if getattr(self, '_active_mode', None) == mode:
                return getattr(self, '_match_df', None)
            snap = (self._mode_state or {}).get(mode)
            if snap is None:
                return None
            return snap.get('_match_df')
        except Exception:
            return None

    def _utility_export_quant_areas(self):
        """⑦ Select Quant Ion で選択された (class, mode) の kept スポットの
        サンプル別エリア値を 1 シート xlsx で出力。"""
        qi = getattr(self, '_quant_ion_choices', None) or {}
        valid = {cls: m for cls, m in qi.items() if m in ('pos', 'neg')}
        if not valid:
            QMessageBox.warning(
                self, "No Quant Ion Selection",
                "Please use '⑦ Select Quant Ion…' first to choose\n"
                "pos / neg mode per lipid class before exporting\n"
                "the area table.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Quant Ion Table",
            timestamped_filename("quant_ion_table.xlsx"),
            "Excel (*.xlsx)")
        if not path:
            return
        if not path.lower().endswith('.xlsx'):
            path += '.xlsx'
        # mode -> FileEntry
        pos_fe = next(
            (fe for fe in (self._all_entries or [])
             if getattr(fe, 'ion_mode', None) == 'pos'), None)
        neg_fe = next(
            (fe for fe in (self._all_entries or [])
             if getattr(fe, 'ion_mode', None) == 'neg'), None)
        pos_sample_cols = list(pos_fe.sample_columns()) if pos_fe else []
        neg_sample_cols = list(neg_fe.sample_columns()) if neg_fe else []
        # 指標(height / area)ごとのサンプル値。1 極性で 1〜2 個。
        _tabs = {}
        for _m, _f in (('pos', pos_fe), ('neg', neg_fe)):
            if _f is None:
                continue
            try:
                _tabs[_m] = dict(fe_intensity_tables(_f))
            except Exception as _e:
                log.warning(f"intensity tables failed ({_m}): {_e}")
                _tabs[_m] = {}
        _metrics = []
        for _m in ('pos', 'neg'):
            for _k in (_tabs.get(_m) or {}):
                if _k not in _metrics:
                    _metrics.append(_k)
        # union (pos 先頭 → neg だけ追加)
        all_sample_cols = list(pos_sample_cols) + [
            c for c in neg_sample_cols if c not in pos_sample_cols]
        # mode -> match_df
        pos_mdf = self._get_match_df_for_mode('pos')
        neg_mdf = self._get_match_df_for_mode('neg')
        # 信頼度列のもとになる帰属（両モードぶん）を 1 回だけ作る
        _attr_q: dict = {}
        try:
            for (_m, _rt, _mz), _ent in (
                    self._build_merged_attribution() or {}).items():
                _attr_q.setdefault(str(_m), {})[
                    (round(float(_rt), 6), round(float(_mz), 6))] = _ent
        except Exception as _e:
            log.warning(f"attribution build failed (quant): {_e}")
        try:
            _sig_thr_q = float(self._spin_sigma.value())
        except Exception:
            _sig_thr_q = 3.0
        rows_by_metric = {k: [] for k in _metrics}
        skipped_cls: list[str] = []
        for cls, mode in valid.items():
            fe   = pos_fe   if mode == 'pos' else neg_fe
            mdf  = pos_mdf  if mode == 'pos' else neg_mdf
            scs  = pos_sample_cols if mode == 'pos' else neg_sample_cols
            if fe is None or mdf is None:
                skipped_cls.append(f"{cls} ({mode}: no data)")
                continue
            if not isinstance(mdf, pd.DataFrame) or mdf.empty:
                skipped_cls.append(f"{cls} ({mode}: empty match_df)")
                continue
            if ('lipid_class' not in mdf.columns
                    or 'final_status' not in mdf.columns):
                skipped_cls.append(f"{cls} ({mode}: invalid match_df)")
                continue
            sub = mdf[(mdf['lipid_class'] == cls)
                      & (mdf['final_status'] == STATUS_KEPT)]
            if sub.empty:
                skipped_cls.append(f"{cls} ({mode}: no kept rows)")
                continue
            for _, r in sub.iterrows():
                # その座標の帰属から信頼度を引く
                _e = None
                try:
                    _e = (_attr_q.get(mode) or {}).get(
                        (round(float(r['obs_rt']), 6),
                         round(float(r['obs_mz']), 6)))
                except Exception:
                    _e = None
                _lab, _sg, _by, _vs = annotation_confidence(_e, r, _sig_thr_q)
                _meta = {
                    'Class':     cls,
                    'Compound':  r.get('compound', ''),
                    'Adduct':    r.get('adduct', ''),
                    'Mode':      mode,
                    'RT (min)':  r.get('obs_rt', np.nan),
                    'm/z':       r.get('obs_mz', np.nan),
                    'confidence':  _lab,
                    'sigma':       _sg,
                    'decided_by':  _by,
                    'competitors': _vs,
                }
                # 指標ごとに同じ行を作り、値だけ差し替える
                for _k in _metrics:
                    _vals = (_tabs.get(mode) or {}).get(_k)
                    if _vals is None:
                        continue
                    row_dict = dict(_meta)
                    try:
                        midx = int(r['matched_idx'])
                        src = _vals.iloc[midx]
                        for col in scs:
                            if col in src.index:
                                row_dict[col] = src[col]
                    except Exception as _e:
                        log.warning(f"[export_quant_table] row build failed "
                              f"({cls} {mode} {_k}): {_e}")
                    rows_by_metric[_k].append(row_dict)
        if not any(rows_by_metric.values()):
            QMessageBox.warning(
                self, "Nothing to export",
                "No kept matches found for the selected quant ion classes.\n"
                "Run All first, then choose ⑦ Select Quant Ion, "
                "then export.")
            return
        # 二峰性クラス(DG)は 1 つの分子種が 2 スポットとして残る。
        # プリカーサー定量では 2 本の合計が求める量なので、ここで足す。
        # IS も同じ規則でまとまるので正規化の分母も 2 本の和になる。
        _sum_on = False
        try:
            _sum_on = bool(self._sum_two_series_enabled()
                           and self._two_series_classes())
        except Exception:
            _sum_on = False
        if _sum_on:
            _ts_cls = self._two_series_classes()
            _n_before = sum(len(v) for v in rows_by_metric.values())
            for _k in list(rows_by_metric.keys()):
                rows_by_metric[_k] = _quant_sum_two_series_rows(
                    rows_by_metric[_k], _ts_cls, all_sample_cols)
            _n_after = sum(len(v) for v in rows_by_metric.values())
            log.info(f"quant export summed two-series peaks: "
                  f"{_n_before} \u2192 {_n_after} rows "
                  f"(classes: {', '.join(_ts_cls)})")
        meta_cols = ['Class', 'Compound', 'Adduct', 'Mode', 'RT (min)',
                     'm/z']
        if _sum_on:
            meta_cols += ['peaks', 'RT (peaks)']
        meta_cols += list(ANNOTATION_CONFIDENCE_COLUMNS)
        cols = meta_cols + all_sample_cols
        _out = {}
        for _k in _metrics:
            if not rows_by_metric.get(_k):
                continue
            _df = pd.DataFrame(rows_by_metric[_k], columns=cols)
            # Class → Mode → RT 順でソート
            try:
                _df = _df.sort_values(
                    by=['Class', 'Mode', 'RT (min)']).reset_index(drop=True)
            except Exception:
                pass
            _out[_k] = _df
        df = _out[_metrics[0]] if _metrics and _metrics[0] in _out else \
            next(iter(_out.values()))
        try:
            with pd.ExcelWriter(path) as _xw:
                for _k, _df in _out.items():
                    _df.to_excel(_xw, sheet_name=_k, index=False)
        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))
            return
        msg = (f"Saved: {path}\n\nSheets: {', '.join(_out)}\n"
               f"Rows: {len(df)}   Classes: {len(valid)}")
        if len(_out) == 1:
            msg += ("\n\nOnly one intensity metric is available — the "
                    "alignment was run with that column, so the other "
                    "one was never stored.")
        if skipped_cls:
            msg += ("\n\nSkipped classes (no kept rows):\n  - "
                    + "\n  - ".join(skipped_cls))
        QMessageBox.information(self, "Saved", msg)

    def _snapshot_mode_state(self, mode: str):
        """現在の状態を mode キーに保存する。"""
        snap = {}
        for var in self._PER_MODE_STATE_VARS:
            snap[var] = getattr(self, var, None)
        self._mode_state[mode] = snap

    def _load_mode_state(self, mode: str):
        """mode キーから状態を復元する(ない場合は初期化)。"""
        snap = self._mode_state.get(mode)
        if snap is None:
            # 初期化(その mode の FileEntry を all_entries から探す)
            target_fe = None
            for fe2 in (self._all_entries or []):
                if getattr(fe2, 'ion_mode', None) == mode:
                    target_fe = fe2
                    break
            if target_fe is None:
                return False
            # ピーク表ロード
            try:
                rt = target_fe.rt_array()
                mz = target_fe.mz_array()
                if rt is None or mz is None:
                    return False
                import numpy as _np
                mask = ~(_np.isnan(rt) | _np.isnan(mz))
                self._rt = rt[mask]
                self._mz = mz[mask]
                try:
                    sample_cols = target_fe.sample_columns()
                    import pandas as _pd
                    self._sample_intensities = (
                        target_fe.df[sample_cols]
                        .apply(_pd.to_numeric, errors='coerce')
                        .values[mask].astype(float))
                    self._sample_col_names = list(sample_cols)
                except Exception:
                    self._sample_intensities = None
                    self._sample_col_names = []
            except Exception as e:
                log.warning(f"[mode switch] failed to load {mode}: {e}")
                return False
            # 解析関連はリセット
            self._match_df = None
            self._matched_coords = set()
            self._reserved_coords = set()
            self._conflict_coords = set()
            self._conflict_map = {}
            self._outlier_coords = {}
            self._coherence_outliers = {}
            self._coherence_assignments = {}
            self._coherence_models = {}
            self._coherence_pairs = []
            self._coherence_all_pairs = []
            self._adduct_attribution = {}
            self._manual_winner_coords = set()
            self._is_filter_choices = {}
            self._class_ref_rt = {}
            # FileEntry を切替えた扱いで self.fe を差替
            self.fe = target_fe
            # そのモードの保存済み設定があれば適用。
            # is_filter_choices / class_ref_rt / library_path /
            # pipeline_state などが復元され、auto-replay のモード切替時にも
            # IS Filter ダイアログをスキップできるようになる。
            saved = (self._all_initial_settings or {}).get(mode)
            if saved:
                try:
                    self._restore_settings(saved)
                except Exception as e:
                    log.warning(f"[_load_mode_state] restore {mode} failed: {e}")
            return True
        # snapshot から復元
        for var, val in snap.items():
            setattr(self, var, val)
        # FileEntry も切替
        for fe2 in (self._all_entries or []):
            if getattr(fe2, 'ion_mode', None) == mode:
                self.fe = fe2
                break
        return True

    def _on_mode_radio_changed(self, checked):
        """Pos/Neg selector で mode を切替える。

        現在のモードの状態を snapshot 保存し、他モードの状態を復元 or 初期化。
        切替後は filter rebuild + scatter redraw で UI を更新。
        """
        if not checked:
            return  # 排他的トグルの片方のみ反応
        target_mode = 'pos' if self._mode_radio_pos.isChecked() else 'neg'
        if target_mode == self._active_mode:
            return  # 既に target mode

        # 1) 現在の状態を snapshot
        self._snapshot_mode_state(self._active_mode)

        # 2) target mode の状態を復元 (or 初期化)
        ok = self._load_mode_state(target_mode)
        if not ok:
            QMessageBox.warning(
                self, "Mode switch failed",
                f"Could not load {target_mode.upper()} mode data.\n"
                f"Make sure the {target_mode.upper()} peak table is loaded "
                "in the main window.")
            # ラジオを元に戻す
            if self._active_mode == 'pos':
                self._mode_radio_pos.setChecked(True)
            else:
                self._mode_radio_neg.setChecked(True)
            return

        # 3) active mode 更新
        self._active_mode = target_mode

        # 4) UI 更新: filter / scatter / match table を再構築
        try:
            classes = []
            if self._match_df is not None:
                classes = self._match_df['lipid_class'].unique().tolist()
            self._rebuild_filters(classes)
        except Exception as e:
            log.warning(f"[mode switch] _rebuild_filters failed: {e}")

        # 軸エイリアスを先に再割当
        # _base_scatter_pos / _base_scatter_neg を物理位置別に作っており、
        # 両者は __init__ 時点で既に正しいモードのデータを保持している。
        # モード切替時に offsets を上書きする必要はない。
        self._refresh_axis_aliases()

        # Mode 切替時は scatter を再描画しない
        # 旧仕様では _draw_overlay を呼んでいたが、active 軸 vs secondary 軸の
        # 描画スタイル差(凡例、特殊マーカー)で feature spot 表示が swap して
        # 見える問題があった。Mode 切替は Match table の表示モード選択に
        # 役割を絞り、scatter は前回のレンダリングを維持する。
        try:
            if self._match_df is not None:
                matched = self._match_df[self._match_df['matched']]
                self._rebuild_match_table(matched)
            else:
                self._match_table.setRowCount(0)
        except Exception as e:
            log.warning(f"[mode switch] match table rebuild failed: {e}")

        # ⑤ ボタン状態も切替後の active mode の attribution に同期
        try:
            self._sync_adduct_button_state()
        except Exception:
            pass

        log.info(f"[mode switch] active mode: {self._active_mode}")

    # ── Parameters ダイアログ ─────────────────────────────
    def _open_parameters_dialog(self):
        """全ステップのパラメータを 1 個のダイアログでまとめて編集する。

        Spinbox 群は self._params_holder の子として保持されているため、
        ダイアログにはこれらを参照するラベル + 値を表示する。OK で各
        spinbox の値が確定する(変更は spinbox に直接反映される)。
        """
        dlg = ParametersDialog(
            host=self,
            parent=self,
        )
        dlg.exec()

    def _open_manage_files_dialog(self):
        """Open the Manage Files dialog.

        The dialog is modeless (can stay open alongside the main view).
        Forwards add/remove/inspect actions to MainWindow (host) via parent.
        """
        host = self._host  # MainWindow
        if host is None or not hasattr(host, '_file_entries'):
            QMessageBox.warning(
                self, "Cannot open dialog",
                "Host main window not available.")
            return
        # Reuse a single instance to keep state across opens
        dlg = getattr(self, '_manage_files_dlg', None)
        if dlg is None or not dlg.isVisible():
            dlg = ManageFilesDialog(host=host, parent=self)
            dlg.setAttribute(Qt.WA_DeleteOnClose, False)
            dlg.filesChanged.connect(self._refresh_data_summary)
            self._manage_files_dlg = dlg
        dlg.refresh()
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()

    def _refresh_data_summary(self):
        """Update the Data: summary label."""
        host = self._host
        if host is None or not hasattr(host, '_file_entries'):
            return
        entries = list(host._file_entries or [])
        if not entries:
            text = "(no data loaded)"
        elif len(entries) == 1:
            fe = entries[0]
            text = f"#{fe.tag.lstrip('#')} [{fe.ion_mode}] {fe.path.name}"
        else:
            n_pos = sum(1 for fe in entries if fe.ion_mode == 'pos')
            n_neg = sum(1 for fe in entries if fe.ion_mode == 'neg')
            text = (f"{len(entries)} files loaded "
                    f"({n_pos} pos, {n_neg} neg)")
        try:
            self._lbl_data_summary.setText(text)
        except Exception:
            pass
        try:
            self._refresh_import_summary()
        except Exception:
            pass
        # stub fe で起動していた場合、最初の実 fe に切り替える
        # stub 以外でも、ファイル追加で他方モードが新たに揃った場合
        # 他方モード scatter を更新する
        # pos が読み込まれていれば active mode を必ず pos に固定
        # (プロット上=pos、下=neg、処理順 pos→neg と整合させる)
        try:
            if entries:
                pos_fe = next(
                    (f for f in entries if getattr(f, 'ion_mode', None) == 'pos'),
                    None)
                cur_is_stub = (self.fe.__class__.__name__ == 'StubFileEntry')
                cur_mode = getattr(self.fe, 'ion_mode', None)
                if cur_is_stub:
                    # stub から実 fe へ:pos があれば pos、なければ最初の fe
                    self._load_fe(pos_fe if pos_fe is not None else entries[0])
                elif pos_fe is not None and cur_mode != 'pos':
                    # neg active 状態で pos が後から追加された → pos に切替
                    self._load_fe(pos_fe)
                else:
                    # 通常の差分更新(他方モードの scatter のみ)
                    self._refresh_other_mode_scatter()
        except Exception as _e:
            log.warning(f"[refresh -> _load_fe] failed: {_e}")

    def _refresh_other_mode_scatter(self):
        """他方モード(現在の active mode の反対側)の scatter plot を
        host._file_entries から最新データで描き直す。active mode 側は触らない。
        """
        if self.fe is None:
            return
        try:
            import numpy as _np
            host = self._host
            if host is not None and hasattr(host, '_file_entries'):
                self._all_entries = list(host._file_entries or [])
            other_rt, other_mz, _ = self._get_other_mode_data()
            target_scatter = (self._base_scatter_neg
                              if self.fe.ion_mode == 'pos'
                              else self._base_scatter_pos)
            target_ax = (self._ax_sc_neg if self.fe.ion_mode == 'pos'
                         else self._ax_sc_pos)
            if other_rt is not None and len(other_rt) > 0:
                target_scatter.set_offsets(
                    _np.column_stack([other_rt, other_mz]))
                o_rt_min = float(other_rt.min())
                o_rt_max = float(other_rt.max())
                o_mz_min = float(other_mz.min())
                o_mz_max = float(other_mz.max())
                o_rt_pad = max((o_rt_max - o_rt_min) * 0.02, 0.1)
                o_mz_pad = max((o_mz_max - o_mz_min) * 0.02, 1.0)
                o_xlim = (o_rt_min - o_rt_pad, o_rt_max + o_rt_pad)
                o_ylim = (o_mz_min - o_mz_pad, o_mz_max + o_mz_pad)
                target_ax.set_xlim(*o_xlim)
                target_ax.set_ylim(*o_ylim)
                # home リセット先も更新
                if self.fe.ion_mode == 'pos':
                    self._home_xlim_neg = o_xlim
                    self._home_ylim_neg = o_ylim
                else:
                    self._home_xlim_pos = o_xlim
                    self._home_ylim_pos = o_ylim
            else:
                target_scatter.set_offsets(_np.empty((0, 2)))
            if self._canvas is not None:
                self._canvas.draw_idle()
        except Exception as _e:
            log.warning(f"[_refresh_other_mode_scatter] failed: {_e}")

    def _load_fe(self, new_fe):
        """Stub fe を実 fe に差し替えて AP を再初期化する。

        StubFileEntry で起動した AP に最初の実 FileEntry を読み込んだとき、
        または Manage Files で files が変化したときに呼ぶ。scatter plot の
        データ更新、active mode 切替、window title 更新を行う。
        """
        if new_fe is None:
            return
        try:
            import numpy as _np
            old_mode = getattr(self.fe, 'ion_mode', None)
            self.fe = new_fe
            # Update title
            self.setWindowTitle(ap_window_title(new_fe.label))
            # Update rt / mz arrays
            rt = new_fe.rt_array()
            mz = new_fe.mz_array()
            if rt is None or mz is None:
                return
            mask = ~(_np.isnan(rt) | _np.isnan(mz))
            self._rt = rt[mask]
            self._mz = mz[mask]
            # Sample intensities
            try:
                sample_cols = new_fe.sample_columns()
                self._sample_intensities = (
                    new_fe.df[sample_cols].apply(
                        pd.to_numeric, errors="coerce")
                      .values[mask].astype(float))
                self._sample_col_names = list(sample_cols)
            except Exception:
                self._sample_intensities = None
                self._sample_col_names = []
            # Update _all_entries
            host = self._host
            if host is not None and hasattr(host, '_file_entries'):
                self._all_entries = list(host._file_entries or [])
            # Update scatter data on active mode side
            try:
                if new_fe.ion_mode == 'pos':
                    self._base_scatter_pos.set_offsets(
                        list(zip(self._rt, self._mz))
                        if len(self._rt) > 0 else _np.empty((0, 2)))
                else:
                    self._base_scatter_neg.set_offsets(
                        list(zip(self._rt, self._mz))
                        if len(self._rt) > 0 else _np.empty((0, 2)))
                # 他方モードの scatter も _get_other_mode_data から更新
                try:
                    other_rt, other_mz, _other_int = self._get_other_mode_data()
                except Exception:
                    other_rt, other_mz = None, None
                if new_fe.ion_mode == 'pos':
                    if other_rt is not None and len(other_rt) > 0:
                        self._base_scatter_neg.set_offsets(
                            _np.column_stack([other_rt, other_mz]))
                    else:
                        self._base_scatter_neg.set_offsets(_np.empty((0, 2)))
                else:
                    if other_rt is not None and len(other_rt) > 0:
                        self._base_scatter_pos.set_offsets(
                            _np.column_stack([other_rt, other_mz]))
                    else:
                        self._base_scatter_pos.set_offsets(_np.empty((0, 2)))
                # Rescale axes - active mode 側
                if len(self._rt) > 0:
                    rt_min, rt_max = float(self._rt.min()), float(self._rt.max())
                    mz_min, mz_max = float(self._mz.min()), float(self._mz.max())
                    rt_pad = max((rt_max - rt_min) * 0.02, 0.1)
                    mz_pad = max((mz_max - mz_min) * 0.02, 1.0)
                    new_xlim = (rt_min - rt_pad, rt_max + rt_pad)
                    new_ylim = (mz_min - mz_pad, mz_max + mz_pad)
                    if new_fe.ion_mode == 'pos':
                        self._ax_sc_pos.set_xlim(*new_xlim)
                        self._ax_sc_pos.set_ylim(*new_ylim)
                        # home リセット先も更新
                        self._home_xlim_pos = new_xlim
                        self._home_ylim_pos = new_ylim
                    else:
                        self._ax_sc_neg.set_xlim(*new_xlim)
                        self._ax_sc_neg.set_ylim(*new_ylim)
                        self._home_xlim_neg = new_xlim
                        self._home_ylim_neg = new_ylim
                    # active mode 用エイリアスも更新
                    self._home_xlim = new_xlim
                    self._home_ylim = new_ylim
                # Rescale axes - 他方モード側
                if other_rt is not None and len(other_rt) > 0:
                    o_rt_min, o_rt_max = float(other_rt.min()), float(other_rt.max())
                    o_mz_min, o_mz_max = float(other_mz.min()), float(other_mz.max())
                    o_rt_pad = max((o_rt_max - o_rt_min) * 0.02, 0.1)
                    o_mz_pad = max((o_mz_max - o_mz_min) * 0.02, 1.0)
                    o_xlim = (o_rt_min - o_rt_pad, o_rt_max + o_rt_pad)
                    o_ylim = (o_mz_min - o_mz_pad, o_mz_max + o_mz_pad)
                    if new_fe.ion_mode == 'pos':
                        self._ax_sc_neg.set_xlim(*o_xlim)
                        self._ax_sc_neg.set_ylim(*o_ylim)
                        # home リセット先も更新
                        self._home_xlim_neg = o_xlim
                        self._home_ylim_neg = o_ylim
                    else:
                        self._ax_sc_pos.set_xlim(*o_xlim)
                        self._ax_sc_pos.set_ylim(*o_ylim)
                        self._home_xlim_pos = o_xlim
                        self._home_ylim_pos = o_ylim
                # active mode ポインタ更新
                if new_fe.ion_mode == 'pos':
                    self._ax_sc       = self._ax_sc_pos
                    self._ax_sc_other = self._ax_sc_neg
                    self._base_scatter = self._base_scatter_pos
                else:
                    self._ax_sc       = self._ax_sc_neg
                    self._ax_sc_other = self._ax_sc_pos
                    self._base_scatter = self._base_scatter_neg
                self._active_mode = new_fe.ion_mode
                if self._canvas is not None:
                    self._canvas.draw_idle()
            except Exception as _e:
                log.warning(f"[_load_fe scatter update] {_e}")
        except Exception as e:
            log.warning(f"[_load_fe] failed: {e}")

    # ── AP relay export サブタブの delegate / refresh ────────────
    def _on_tabs_current_changed(self, index: int):
        """Utility タブに切り替わったとき、relay export サブタブが
        表示中なら host から refresh。(fix17: 旧 Tasks タブ)"""
        try:
            text = self._tabs.tabText(index)
        except Exception:
            text = ""
        # ① タブを開いたときに pyOpenMS を一度だけ確認する
        if text == TAB_LABELS['align']:
            try:
                self._al_ensure_pyopenms_checked()
            except Exception as e:
                log.warning(f"{e}")
            return
        if text != TAB_LABELS['util']:
            return
        try:
            sub = self._utility_subtabs
            if sub.tabText(sub.currentIndex()) == SUBTAB_LABELS['relay']:
                self._refresh_tasks_tab_from_host()
        except Exception:
            pass

    def _on_utility_subtab_changed(self, index: int):
        """Utility 内サブタブ切替。relay export に入ったら refresh。"""
        try:
            text = self._utility_subtabs.tabText(index)
        except Exception:
            return
        if text == SUBTAB_LABELS['relay']:
            self._refresh_tasks_tab_from_host()

    def _refresh_tasks_tab_from_host(self):
        """MainWindow の状態から AP Tasks タブを再構築。"""
        host = self._host
        if host is None:
            return
        # Output directory
        try:
            text = host._lbl_outdir.text() if hasattr(host, '_lbl_outdir') else ""
            self._ap_lbl_outdir.setText(text)
        except Exception:
            pass
        # Tasks 表(host の _task_table をコピー)
        try:
            src = getattr(host, '_task_table', None)
            if src is not None:
                # 編集シグナル抑止
                self._ap_task_table.blockSignals(True)
                try:
                    n_rows = src.rowCount()
                    self._ap_task_table.setRowCount(n_rows)
                    for r in range(n_rows):
                        for c in range(5):
                            # host の Ion Mode (col 3) / File (col 4)
                            # は QComboBox の cellWidget。それ以外は QTableWidgetItem。
                            widget = src.cellWidget(r, c)
                            if widget is not None and isinstance(widget, QComboBox):
                                txt = widget.currentText()
                            else:
                                it = src.item(r, c)
                                txt = it.text() if it is not None else ""
                            ap_it = QTableWidgetItem(txt)
                            self._ap_task_table.setItem(r, c, ap_it)
                finally:
                    self._ap_task_table.blockSignals(False)
        except Exception as _e:
            log.warning(f"[AP Tasks refresh] failed: {_e}")
        # Log
        try:
            log_text = host._lbl_log.text() if hasattr(host, '_lbl_log') else ""
            self._ap_lbl_log.setText(log_text)
        except Exception:
            pass

    def _ap_choose_outdir(self):
        host = self._host
        if host is None or not hasattr(host, '_choose_outdir'):
            return
        host._choose_outdir()
        self._refresh_tasks_tab_from_host()

    def _ap_add_task_row(self):
        host = self._host
        if host is None or not hasattr(host, '_add_task_row'):
            return
        host._add_task_row()
        self._refresh_tasks_tab_from_host()

    def _ap_delete_task_row(self):
        host = self._host
        if host is None or not hasattr(host, '_delete_task_row'):
            return
        # MainWindow の _delete_task_row は selected row を使うが、AP 側で選択
        # された行を host にも反映してから呼ぶ。
        try:
            ap_row = self._ap_task_table.currentRow()
            src = getattr(host, '_task_table', None)
            if src is not None and ap_row >= 0 and ap_row < src.rowCount():
                src.setCurrentCell(ap_row, 0)
        except Exception:
            pass
        host._delete_task_row()
        self._refresh_tasks_tab_from_host()

    def _ap_paste_tasks(self):
        host = self._host
        if host is None or not hasattr(host, '_paste_tasks'):
            return
        host._paste_tasks()
        self._refresh_tasks_tab_from_host()

    def _ap_copy_all_tasks(self):
        host = self._host
        if host is None or not hasattr(host, '_copy_all_tasks'):
            return
        host._copy_all_tasks()
        # クリップボード操作のみなので refresh は不要だが整合のため呼ぶ
        self._refresh_tasks_tab_from_host()

    def _ap_run(self):
        host = self._host
        if host is None or not hasattr(host, '_run'):
            return
        host._run()
        self._refresh_tasks_tab_from_host()

    def _ap_on_task_cell_changed(self, row: int, col: int):
        """AP Tasks 表でユーザーがセルを編集したら host 側にも反映。"""
        host = self._host
        if host is None:
            return
        src = getattr(host, '_task_table', None)
        if src is None:
            return
        try:
            ap_it = self._ap_task_table.item(row, col)
            if ap_it is None:
                return
            txt = ap_it.text()
            # host 側の表を同じサイズにしておく
            while src.rowCount() <= row:
                src.insertRow(src.rowCount())
            src.blockSignals(True)
            try:
                # col 3 (Ion Mode) / col 4 (File) は QComboBox なので
                # cellWidget の現在テキストを書き換える
                widget = src.cellWidget(row, col)
                if widget is not None and isinstance(widget, QComboBox):
                    idx = widget.findText(txt)
                    if idx >= 0:
                        widget.setCurrentIndex(idx)
                else:
                    src.setItem(row, col, QTableWidgetItem(txt))
            finally:
                src.blockSignals(False)
        except Exception as _e:
            log.warning(f"[AP task cell -> host] failed: {_e}")

    # ── Session menu delegate ─────────────────────────
    def _ap_save_session(self):
        """AP Session menu → host の _save_session を呼ぶ。"""
        host = self._host
        if host is None or not hasattr(host, '_save_session'):
            QMessageBox.warning(
                self, "Cannot save",
                "Host main window is not available.")
            return
        host._save_session()

    def _ap_load_session(self):
        """AP Session menu → host の _load_session を呼んだ後、AP 側を完全 refresh。

        既存 AP インスタンスを再利用する v6.x のフローでは、host が
        _preview_settings_by_mode に保存した設定(library_path / class_ref_rt /
        is_filter_choices / pipeline_state など)を AP に流し込まないと、
        Select Quant Ion 等の状態が反映されない。new PreviewRTDialog 作成時の
        __init__(initial_settings, all_initial_settings) と同等の処理を
        明示的に行う。
        """
        host = self._host
        if host is None or not hasattr(host, '_load_session'):
            QMessageBox.warning(
                self, "Cannot load",
                "Host main window is not available.")
            return
        host._load_session()
        # 新しい FileEntry に切替
        try:
            new_entries = list(getattr(host, '_file_entries', None) or [])
            if new_entries:
                new_pos_fe = next(
                    (fe for fe in new_entries
                     if getattr(fe, 'ion_mode', None) == 'pos'),
                    new_entries[0])
                self._load_fe(new_pos_fe)
        except Exception as _e:
            log.warning(f"[_ap_load_session force load_fe] failed: {_e}")
        # Data 行更新(_refresh_data_summary 内で stale 検出のロジックも走る)
        try:
            self._refresh_data_summary()
        except Exception as _e:
            log.warning(f"[_ap_load_session refresh data] failed: {_e}")
        # Tasks タブも refresh
        try:
            self._refresh_tasks_tab_from_host()
        except Exception as _e:
            log.warning(f"[_ap_load_session refresh tasks] failed: {_e}")
        # stale な mode_state を clear。Load Session 後の
        # Run All で各モードを fresh で初期化させる。
        # per-mode state も明示的に clear して、stale state が
        # 残留しないようにする(load 直後の neg 処理を確実にするため)。
        try:
            self._mode_state = {}
            self._match_df = None
            self._matched_coords = set()
            self._reserved_coords = set()
            self._conflict_coords = set()
            self._conflict_map = {}
            self._outlier_coords = {}
            self._coherence_outliers = {}
            self._coherence_assignments = {}
            self._coherence_models = {}
            self._coherence_pairs = []
            self._coherence_all_pairs = []
            self._adduct_attribution = {}
            self._manual_winner_coords = set()
        except Exception as _e:
            log.warning(f"[_ap_load_session state reset] failed: {_e}")
        # 全モード分の preview_settings を _all_initial_settings に流し込み
        try:
            host_all = getattr(host, '_preview_settings_by_mode', None) or {}
            self._all_initial_settings = dict(host_all)
        except Exception as _e:
            log.warning(f"[_ap_load_session all_initial_settings] failed: {_e}")
        # active mode の設定を _restore_settings で適用
        # (library_path / class_ref_rt / is_filter_choices / pipeline_state)
        try:
            cur_mode = getattr(self.fe, 'ion_mode', None) or self._active_mode
            saved = (self._all_initial_settings or {}).get(cur_mode)
            if saved:
                self._restore_settings(saved)
        except Exception as _e:
            log.warning(f"[_ap_load_session restore_settings] failed: {_e}")
        # 保存済み quant_ion / filter_exemptions / params をこのダイアログに適用
        try:
            fe_exempt = getattr(host, '_loaded_session_filter_exempt', None)
            qi = getattr(host, '_loaded_session_quant_ion', None)
            params = getattr(host, '_loaded_session_params', None)
            if params and hasattr(self, '_apply_params'):
                try:
                    self._apply_params(params)
                except Exception as _e:
                    log.warning(f"[_ap_load_session apply_params] failed: {_e}")
            if qi:
                try:
                    self._quant_ion_choices = dict(qi)
                except Exception:
                    pass
            if fe_exempt:
                try:
                    self._filter_class_exemptions = {
                        k: set(v) for k, v in fe_exempt.items()}
                except Exception:
                    pass
            # ⑦ Select Quant Ion など各 step ボタンの状態を再計算
            try:
                self._sync_step_button_states()
            except Exception as _e:
                log.warning(f"[_ap_load_session sync_btn] failed: {_e}")
            if hasattr(self, 'auto_replay_after_load'):
                self.auto_replay_after_load(
                    filter_exemptions=fe_exempt,
                    quant_ion_choices=qi)
        except Exception as e:
            log.warning(f"[_ap_load_session auto_replay] failed: {e}")

    def _ap_load_session_settings_only(self):
        """AP Session タブ: Settings-only ロード。
        Parameters + Filter exemptions のみ別セッションから適用し、
        現在のデータファイル / Tasks / ion 選択はそのまま保持。
        適用後に自動で Run All を実行する。"""
        host = self._host
        if host is None or not hasattr(host, '_load_session_settings_only'):
            QMessageBox.warning(
                self, "Cannot load",
                "Host main window is not available.")
            return
        # 前提条件の確認ポップアップ
        _confirm = QMessageBox(self)
        _confirm.setIcon(QMessageBox.Information)
        _confirm.setWindowTitle("Load Settings Only — Prerequisites")
        _confirm.setText(
            "Before applying settings, make sure the following are ready in "
            "the Analysis tab:")
        _confirm.setInformativeText(
            "  •  Data files are loaded (Load data…)\n"
            "  •  Lipid library is loaded (Load library…)\n\n"
            "Settings-only load will apply parameters and filter exemptions\n"
            "from the selected session file, then automatically execute Run All\n"
            "on the currently loaded data.\n\n"
            "Continue?")
        _confirm.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        _confirm.setDefaultButton(QMessageBox.Ok)
        if _confirm.exec() != QMessageBox.Ok:
            return
        ok = host._load_session_settings_only()
        if not ok:
            return
        # 自身に params / filter_exemptions を適用
        try:
            params = getattr(host, '_loaded_session_params', None)
            fe_exempt = getattr(host, '_loaded_session_filter_exempt', None)
            if params and hasattr(self, '_apply_params'):
                try:
                    self._apply_params(params)
                except Exception as _e:
                    log.warning(f"[_ap_load_session_settings_only apply_params] failed: {_e}")
            if fe_exempt:
                try:
                    self._filter_class_exemptions = {
                        k: set(v) for k, v in fe_exempt.items()}
                except Exception as _e:
                    log.warning(f"[_ap_load_session_settings_only fe_exempt] failed: {_e}")
            # step ボタン状態を再計算
            try:
                self._sync_step_button_states()
            except Exception:
                pass
        except Exception as e:
            log.warning(f"[_ap_load_session_settings_only apply] failed: {e}")
        # 自動で Run All を実行(library 未ロード時は _run_all が警告)
        try:
            self._run_all()
        except Exception as e:
            log.warning(f"[_ap_load_session_settings_only run_all] failed: {e}")

    def _open_match_table_dialog(self):
        """化合物テーブルダイアログを表示する。

        初回呼び出し時にインスタンスを生成し、self._match_table と
        self._lbl_match_count をダイアログに格納する。以降はダイアログを
        show + raise + activate する(modeless)。
        """
        if self._match_table_dialog is None:
            self._match_table_dialog = MatchTableDialog(
                self, self._match_table, self._lbl_match_count,
                title=f"Compounds Table — {self.fe.tag}")
        self._match_table_dialog.show()
        self._match_table_dialog.raise_()
        self._match_table_dialog.activateWindow()

    def _on_finished(self, _result):
        """ダイアログ終了時に設定/reserved_coords/coherence 帰属結果を親へ通知"""
        # 子の MatchTableDialog も明示的に閉じる
        if self._match_table_dialog is not None:
            try:
                self._match_table_dialog.hide()
                self._match_table_dialog.deleteLater()
            except Exception:
                pass
            self._match_table_dialog = None
        self.settings_saved.emit(self._current_preview_settings())
        self.reserved_updated.emit(self.fe.tag, set(self._reserved_coords))
        # match_df を MainWindow に伝播
        try:
            self.match_df_updated.emit(self.fe.tag, self._match_df)
        except Exception:
            pass
        # Coherence 帰属結果を loser クラス別の座標集合として渡す
        # エクスポート時、各クラスから「そのクラスが敗者だった座標」を除外
        loser_coords_by_class: dict[str, set] = {}
        for coord, a in self._coherence_assignments.items():
            loser_coords_by_class.setdefault(
                a['loser_class'], set()).add(coord)
        self.coherence_updated.emit(self.fe.tag, loser_coords_by_class)
        # AP 閉じ時はアプリ全体を終了
        # (Stage 5a の host.show() は MainWindow が用済みになったため廃止)
        try:
            QApplication.quit()
        except Exception as _e:
            log.warning(f"[AP close -> quit] failed: {_e}")

    # ── 設定の収集・復元 ──────────────────────────────────────────────

    def _current_preview_settings_for_mode(self, mode: str) -> dict:
        """指定モードの preview settings を返す。
        active モードなら live state、他モードなら _mode_state スナップショット
        からビルドする。両モード分を JSON に保存するために使う。"""
        if mode == self._active_mode:
            return self._current_preview_settings()
        snap = (self._mode_state or {}).get(mode, {}) or {}
        # pipeline_state を snap から導出
        try:
            pipeline_state = {
                "match_run": snap.get('_match_df') is not None,
                "is_filter_on": bool(snap.get('_is_filter_choices') or {}),
                "adduct_filter_on": bool(snap.get('_adduct_attribution') or {}),
                "rt_outlier_on": bool(
                    any(len(v) > 0
                        for v in (snap.get('_outlier_coords') or {}).values())),
                "coherence_on": bool(
                    snap.get('_coherence_models')
                    or snap.get('_coherence_assignments')),
            }
        except Exception:
            pipeline_state = {}
        # class_ref_rt は snap 内のものを使う
        class_ref_rt_snap = snap.get('_class_ref_rt') or {}
        return {
            # global UI 状態(両モード共通)
            "library_path": str(self._lib_path) if self._lib_path else None,
            "ion_mode":     mode,
            "ppm_tol":      self._spin_ppm.value(),
            "void_rt":      self._spin_void_rt.value(),
            "mad_k":        self._spin_iqr.value(),
            "is_tol":       self._spin_is_tol.value(),
            # per-mode 状態(snap 由来)
            "class_ref_rt": {
                cls: list(vals)
                for cls, vals in class_ref_rt_snap.items()
            },
            "is_filter_choices": snap.get('_is_filter_choices') or {},
            "pipeline_state": pipeline_state,
            # 他モードの手動帰属を snap から取得
            "manual_winner_coords": [
                [float(c[0]), float(c[1])]
                for c in (snap.get('_manual_winner_coords') or set())
            ],
            "manual_overrides": self._serialize_manual_overrides(
                snap.get('_match_df')),
        }

    def _current_preview_settings(self) -> dict:
        """現在のPreview設定をdictで返す（セッション保存用）"""
        # pipeline state flags
        try:
            pipeline_state = {
                "match_run": self._match_df is not None,
                "is_filter_on": bool(self._btn_is_filter.isChecked()
                                     if hasattr(self, '_btn_is_filter')
                                     else False),
                "adduct_filter_on": bool(
                    self._btn_adduct_filter.isChecked()
                    if hasattr(self, '_btn_adduct_filter') else False),
                "rt_outlier_on": bool(
                    any(len(v) > 0
                        for v in (self._outlier_coords or {}).values())),
                "coherence_on": bool(
                    self._coherence_models or self._coherence_assignments),
            }
        except Exception:
            pipeline_state = {}
        return {
            "library_path": str(self._lib_path) if self._lib_path else None,
            "ion_mode":     self.fe.ion_mode,
            "ppm_tol":      self._spin_ppm.value(),
            "void_rt":      self._spin_void_rt.value(),
            "mad_k":        self._spin_iqr.value(),
            "is_tol":       self._spin_is_tol.value(),
            "class_ref_rt": {
                cls: list(vals)   # (ref_rt, tol, source) → [ref_rt, tol, source]
                for cls, vals in self._class_ref_rt.items()
            },
            "is_filter_choices": self._is_filter_choices,
            "pipeline_state": pipeline_state,
            # 手動帰属(右クリック curation)を保存
            "manual_winner_coords": [
                [float(c[0]), float(c[1])]
                for c in (getattr(self, '_manual_winner_coords', None)
                          or set())
            ],
            "manual_overrides": self._serialize_manual_overrides(
                self._match_df),
        }

    def _serialize_manual_overrides(self, df) -> list:
        """match_df の manual_status != STATUS_NA 行を
        list-of-dict にシリアライズして JSON 保存可能にする。"""
        try:
            if df is None or df.empty:
                return []
            if 'manual_status' not in df.columns:
                return []
            sub = df[df['manual_status'] != STATUS_NA]
            if sub.empty:
                return []
            out = []
            for _, row in sub.iterrows():
                out.append({
                    'rt': round(float(row['obs_rt']), 6),
                    'mz': round(float(row['obs_mz']), 6),
                    'lipid_class': str(row['lipid_class']),
                    'compound': str(row.get('compound', '')),
                    'adduct': str(row.get('adduct', '')),
                    'status': str(row['manual_status']),
                })
            return out
        except Exception as e:
            log.warning(f"[_serialize_manual_overrides] failed: {e}")
            return []

    def auto_replay_after_load(self,
                                filter_exemptions: dict | None = None,
                                quant_ion_choices: dict | None = None) -> None:
        """Load Session 後に Preview RT を開いた時、
        保存時のパイプライン状態を自動再現する。
        filter_exemptions: 事前に適用するクラス別 Filter exempt 設定。
        quant_ion_choices: 事前に適用する ⑦ Select Quant Ion の選択。
        """
        # filter_exemptions を適用(クラス別 Filter exempt 設定)
        if filter_exemptions:
            try:
                self._filter_class_exemptions = {
                    k: set(v) for k, v in filter_exemptions.items()}
            except Exception:
                pass
        # Quant Ion 選択を事前に適用
        if quant_ion_choices:
            try:
                self._quant_ion_choices = dict(quant_ion_choices)
            except Exception:
                pass
        # 保存時に match_run=True で、ライブラリが存在するなら Run All を起動
        ps = getattr(self, '_restored_pipeline_state', None) or {}
        if not ps.get('match_run'):
            return
        if not self._lib_path or not Path(self._lib_path).exists():
            return
        try:
            # auto-replay フラグを立て、Run All 終了後に
            # 必ず下ろす。例外が起きてもクリアされるよう finally 経由。
            self._auto_replay_in_progress = True
            def _run_with_cleanup():
                try:
                    self._run_all()
                    # Run All 完了後、保存された手動帰属を再適用
                    try:
                        self._apply_restored_manual_attribution()
                    except Exception as e:
                        log.warning(f"[manual restore] failed: {e}")
                    # 確実な最終描画(両モード)。Run All 内の
                    # _draw_both_overlays の後に matplotlib の遅延処理が
                    # 残っていることがあるため、もう一度キックする。
                    try:
                        QTimer.singleShot(50, self._draw_both_overlays)
                    except Exception:
                        pass
                finally:
                    self._auto_replay_in_progress = False
            QTimer.singleShot(0, _run_with_cleanup)
        except Exception as e:
            self._auto_replay_in_progress = False
            log.warning(f"[auto_replay_after_load] failed: {e}")

    def _apply_restored_manual_attribution(self) -> None:
        """Run All 完了後、_all_initial_settings に
        保存された手動帰属(_manual_winner_coords + manual_overrides)を
        各モードの match_df に再適用する。"""
        all_settings = getattr(self, '_all_initial_settings', None) or {}
        if not all_settings:
            return
        any_change = False
        for mode, s in all_settings.items():
            if not s:
                continue
            mw_coords = s.get('manual_winner_coords') or []
            overrides = s.get('manual_overrides') or []
            if not mw_coords and not overrides:
                continue
            # 対象モードの match_df を取得
            if mode == getattr(self, '_active_mode', None):
                df = self._match_df
            else:
                snap = (self._mode_state or {}).get(mode, {}) or {}
                df = snap.get('_match_df')
            if df is None or df.empty:
                continue
            # manual_status 上書きを適用
            if overrides and 'manual_status' in df.columns:
                for ov in overrides:
                    try:
                        rt = round(float(ov.get('rt')), 6)
                        mz = round(float(ov.get('mz')), 6)
                        cls = ov.get('lipid_class')
                        status = ov.get('status')
                        if not status:
                            continue
                        mask = (
                            (df['obs_rt'].round(6) == rt) &
                            (df['obs_mz'].round(6) == mz) &
                            (df['lipid_class'] == cls) &
                            df['matched']
                        )
                        # compound / adduct も一致させるとより厳密
                        comp = ov.get('compound')
                        if comp and 'compound' in df.columns:
                            mask = mask & (df['compound'].astype(str) == str(comp))
                        add = ov.get('adduct')
                        if add and 'adduct' in df.columns:
                            mask = mask & (df['adduct'].astype(str) == str(add))
                        if mask.any():
                            df.loc[mask, 'manual_status'] = status
                            any_change = True
                    except Exception as e:
                        log.warning(f"[manual restore override] skip: {e}")
                df = _update_final_status_df(df)
            # _manual_winner_coords を復元
            coords_set = {
                (round(float(c[0]), 6), round(float(c[1]), 6))
                for c in mw_coords
            }
            if mode == getattr(self, '_active_mode', None):
                self._match_df = df
                self._manual_winner_coords = coords_set
            else:
                snap = (self._mode_state or {}).get(mode, {}) or {}
                snap['_match_df'] = df
                snap['_manual_winner_coords'] = coords_set
                self._mode_state[mode] = snap
        # 再描画
        if any_change:
            try:
                self._draw_overlay()
            except Exception:
                pass
            try:
                matched = self._match_df[self._match_df['matched']]
                self._rebuild_match_table(matched)
            except Exception:
                pass

    def _restore_settings(self, s: dict | None):
        """セッションから読み込んだ設定をウィジェットに反映する。
        マッチング自体は自動実行しない（ユーザーが Match & Overlay を押す）。
        """
        if not s:
            return
        # Ion mode は FileEntry から取得するので復元不要
        if s.get("ppm_tol") is not None:
            self._spin_ppm.setValue(float(s["ppm_tol"]))
        if "void_rt" in s:
            self._spin_void_rt.setValue(float(s["void_rt"]))
        if s.get("mad_k") is not None:
            self._spin_iqr.setValue(float(s["mad_k"]))
        if s.get("is_tol") is not None:
            self._spin_is_tol.setValue(float(s["is_tol"]))
        if s.get("class_ref_rt"):
            self._class_ref_rt = {
                cls: (float(vals[0]), float(vals[1]),
                      str(vals[2]), str(vals[3]))
                for cls, vals in s["class_ref_rt"].items()
                if len(vals) >= 4
            }
        # IS Filter プレビューの選択状態を復元
        if s.get("is_filter_choices"):
            self._is_filter_choices = s["is_filter_choices"]
        lib = s.get("library_path")
        if lib and Path(lib).exists():
            self._lib_path = Path(lib)
            self._lbl_lib.setText(self._lib_path.name)
            self._lbl_match_count.setText(
                "Library restored. Click 'Match & Overlay' to reload.")
        elif lib:
            # パスが存在しない場合はラベルに警告表示
            self._lbl_lib.setText(f"(not found) {Path(lib).name}")
        # pipeline state を保存(後で auto-replay で参照)
        self._restored_pipeline_state = s.get("pipeline_state") or {}

    # ── All / None 一括切替 ───────────────────────────────────────────

    def _on_curated_view_filter_toggled(self, _state=None):
        """Curated view filter checkbox の stateChanged
        ハンドラ。ON で ⑦ 選択を検証して curated_mode=True にし、クラス
        checkbox を同期。OFF で curated_mode=False。"""
        cb = self._class_checks.get(self._TAG_CURATED_VIEW)
        if cb is None:
            return
        checked = bool(cb.isChecked())
        if not checked:
            self._curated_mode = False
            try:
                self._draw_both_overlays()
            except Exception:
                pass
            return
        # ON: ⑦ 検証
        qi = getattr(self, '_quant_ion_choices', None) or {}
        if not qi:
            QMessageBox.warning(
                self, "No Quant Ion Selection",
                "Please use '⑦ Select Quant Ion…' first to choose\n"
                "pos / neg / skip per lipid class before enabling\n"
                "Curated view.")
            cb.blockSignals(True)
            cb.setChecked(False)
            cb.blockSignals(False)
            self._curated_mode = False
            return
        self._curated_mode = True
        # クラス checkbox を ⑦ 選択と同期
        try:
            non_class_tags = self._filter_modifier_tags()
            selected_classes = {
                str(c) for c, ch in qi.items() if ch in ('pos', 'neg')
            }
            for key, ccb in self._class_checks.items():
                if key in non_class_tags:
                    continue
                desired = str(key) in selected_classes
                ccb.blockSignals(True)
                ccb.setChecked(desired)
                self._user_filter_prefs[key] = desired
                ccb.blockSignals(False)
        except Exception as e:
            log.warning(f"[curated view filter] sync failed: {e}")
        try:
            self._draw_both_overlays()
        except Exception as e:
            log.warning(f"[curated view filter] redraw failed: {e}")

    def _select_all_classes(self):
        """クラス membership checkbox を全 ON にして再描画する。
        両モード(上下のスキャッタ)を同時に再描画する。
        "Show classes" ボタンなのでクラス系のみ操作し、
        Filters 修飾子(RT outliers / IS-rejected / Adduct-rejected /
        Coherence outliers / Low-confidence)の状態は維持する。
        Curated view checkbox は OFF にして curated
        mode を解除する(All を押すと全表示が期待される挙動)。"""
        non_class_tags = self._filter_modifier_tags()
        for key, cb in self._class_checks.items():
            if key in non_class_tags:
                continue
            cb.blockSignals(True)
            cb.setChecked(True)
            self._user_filter_prefs[key] = True
            cb.blockSignals(False)
        # Curated view checkbox を OFF + curated_mode フラグも解除
        try:
            self._curated_mode = False
            cv_cb = self._class_checks.get(self._TAG_CURATED_VIEW)
            if cv_cb is not None:
                cv_cb.blockSignals(True)
                cv_cb.setChecked(False)
                cv_cb.blockSignals(False)
        except Exception:
            pass
        self._draw_both_overlays()

    def _select_no_classes(self):
        """クラス membership checkbox を全 OFF にして再描画する。
        Filters 修飾子は変更しない(ラベルは "Show classes")。
        Curated view も OFF にする。"""
        non_class_tags = self._filter_modifier_tags()
        for key, cb in self._class_checks.items():
            if key in non_class_tags:
                continue
            cb.blockSignals(True)
            cb.setChecked(False)
            self._user_filter_prefs[key] = False
            cb.blockSignals(False)
        try:
            self._curated_mode = False
            cv_cb = self._class_checks.get(self._TAG_CURATED_VIEW)
            if cv_cb is not None:
                cv_cb.blockSignals(True)
                cv_cb.setChecked(False)
                cv_cb.blockSignals(False)
        except Exception:
            pass
        self._draw_both_overlays()

    def _apply_filter_class_exemptions(self) -> None:
        """各 Filter で exempt 指定されたクラスの
        rejected フラグを剥がし、final_status を再計算する。
        両モードの match_df に対して適用。"""
        exempt = getattr(self, '_filter_class_exemptions', None) or {}
        ex_is = exempt.get('is', set())
        ex_ad = exempt.get('adduct', set())
        ex_rt = exempt.get('rt_outlier', set())
        ex_co = exempt.get('coherence', set())
        if not (ex_is or ex_ad or ex_rt or ex_co):
            return

        def _fixup(df):
            if df is None or getattr(df, 'empty', True):
                return df
            if 'lipid_class' not in df.columns:
                return df
            cls = df['lipid_class']
            changed = False
            if ex_is and 'is_filter_status' in df.columns:
                m = cls.isin(ex_is) & (
                    df['is_filter_status'] == STATUS_REJ_WINDOW)
                if m.any():
                    df.loc[m, 'is_filter_status'] = STATUS_KEPT
                    changed = True
            if ex_ad and 'adduct_filter_status' in df.columns:
                m = cls.isin(ex_ad) & (
                    df['adduct_filter_status'] == STATUS_REJ_ADDUCT)
                if m.any():
                    df.loc[m, 'adduct_filter_status'] = STATUS_KEPT
                    changed = True
            if ex_rt and 'rt_outlier_status' in df.columns:
                m = cls.isin(ex_rt) & (
                    df['rt_outlier_status'] == STATUS_REJ_OUTLIER)
                if m.any():
                    df.loc[m, 'rt_outlier_status'] = STATUS_KEPT
                    changed = True
            if ex_co and 'coherence_status' in df.columns:
                m = cls.isin(ex_co) & df['coherence_status'].isin(
                    [STATUS_REJ_RESIDUAL, STATUS_REJ_LOSER])
                if m.any():
                    df.loc[m, 'coherence_status'] = STATUS_KEPT
                    changed = True
            if changed:
                df = _update_final_status_df(df)
            return df

        # active mode
        self._match_df = _fixup(self._match_df)
        # other mode snapshots
        mode_state = getattr(self, '_mode_state', None) or {}
        for mode_key, snap in mode_state.items():
            if not isinstance(snap, dict):
                continue
            mdf = snap.get('_match_df')
            if mdf is not None:
                snap['_match_df'] = _fixup(mdf)

        # _outlier_coords / _coherence_outliers セットからも exempt クラスを除く
        # (描画は status 列ベースだが、報告系の整合性のため)
        if ex_rt:
            for cls_name in list(self._outlier_coords.keys()):
                if cls_name in ex_rt:
                    self._outlier_coords[cls_name] = set()
        if ex_co:
            for cls_name in list(self._coherence_outliers.keys()):
                if cls_name in ex_co:
                    self._coherence_outliers[cls_name] = set()

    def _filter_modifier_tags(self) -> set:
        """クラス系でない(=Filters 修飾子)タグの集合。
        All/None で操作しない対象。Curated view も含める。"""
        return {
            self._TAG_RT_OUTLIERS,
            self._TAG_COHERENCE_OUTLIERS,
            self._TAG_LOWCONF,
            self._TAG_SHOW_IS_REJECTED,
            self._TAG_SHOW_ADDUCT_REJECTED,
            self._TAG_CURATED_VIEW,
            self._TAG_IS_ONLY,
        }

    def _is_only_enabled(self) -> bool:
        """IS only が ON か。UI が無い経路でも落ちないようにする。"""
        try:
            cb = self._class_checks.get(self._TAG_IS_ONLY)
            return bool(cb is not None and cb.isChecked())
        except Exception:
            return False

    @staticmethod
    def _restrict_to_is(df):
        """IS の行だけに絞る。is_IS 列が無ければ空にする。"""
        if df is None or len(df) == 0:
            return df
        if 'is_IS' not in df.columns:
            log.info("[IS only] match_df に is_IS 列が無いので表示できない")
            return df.iloc[0:0]
        try:
            return df[df['is_IS'].fillna(False).astype(bool)]
        except Exception as _e:
            log.warning(f"[IS only] 絞り込みに失敗: {_e}")
            return df.iloc[0:0]

    def _draw_both_overlays(self):
        """active 軸 + secondary 軸の両方を再描画する。
        フィルタ操作を両モードに同時反映するためのヘルパ。"""
        self._draw_overlay()
        try:
            self._draw_other_mode_overlay()
        except Exception as e:
            log.warning(f"[draw_both_overlays] secondary draw failed: {e}")
        try:
            self._canvas.draw_idle()
        except Exception:
            pass

    # ── ライブラリ操作 ────────────────────────────────────────────────
    def _load_library(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Lipid Library",
            "", "Excel files (*.xlsx *.xlsm);;All (*)")
        if path:
            self._lib_path = Path(path)
            self._lbl_lib.setText(self._lib_path.name)
            self._lbl_match_count.setText(
                "Library loaded. Click 'Match & Overlay'.")

    def _run_match(self):
        """Step 1: 通常マッチング（m/z誤差最小）+ 散布図描画。
        IS処理・RT外れ値・競合検出は行わない。"""
        if not self._lib_path:
            QMessageBox.warning(self, "No Library",
                                "Please load a library file first.")
            return
        try:
            lib_df = load_lipid_library(
                self._lib_path, self.fe.ion_mode)
        except Exception as e:
            QMessageBox.critical(self, "Library Error", str(e))
            return
        if lib_df.empty:
            QMessageBox.warning(self, "Empty",
                                "No entries found for the selected ion mode.")
            return

        # 現ライブラリに存在しないクラスのFix RT設定は破棄
        valid_classes = set(lib_df['lipid_class'].unique())
        self._class_ref_rt = {
            k: v for k, v in self._class_ref_rt.items() if k in valid_classes
        }

        # 通常マッチング（IS処理なし、手動Fix RTは適用）
        match_df, updated_ref_rt, _ = match_library_to_peaks(
            lib_df, self._mz, self._rt, self._spin_ppm.value(),
            void_rt=float(self._spin_void_rt.value()),
            class_ref_rt=self._class_ref_rt or None,
            is_tol=self._spin_is_tol.value(),
            apply_is_filter=False)
        self._match_df = match_df
        # manual Fix RT のみ保持（auto_is は Apply IS Filter で再計算）
        self._class_ref_rt = updated_ref_rt
        matched = match_df[match_df['matched']]

        # クラス色割り当て
        classes = match_df['lipid_class'].unique().tolist()
        colors  = _class_colors(len(classes))
        self._class_colors = dict(zip(classes, colors))

        # マッチ済み座標セット（Not matched 判定用）
        self._matched_coords = set(
            zip(matched['obs_rt'].round(6), matched['obs_mz'].round(6)))

        # 初期マッチング時点では外れ値・競合・reservedは未算出
        self._outlier_coords  = {}
        self._coherence_outliers    = {}
        self._coherence_models      = {}
        self._coherence_assignments = {}
        self._coherence_pairs       = []
        self._coherence_all_pairs   = []
        self._conflict_coords = set()
        self._conflict_map    = {}
        self._reserved_coords = set()
        # 上流再実行で Adduct Ion Filter 状態は無効化
        self._adduct_attribution = {}

        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        self._update_match_count_label(matched, lib_df)

        # 次ステップのボタンを有効化
        self._btn_is_filter.setEnabled(True)
        self._btn_conflict.setEnabled(True)
        self._btn_adduct_filter.setEnabled(True)
        self._btn_adduct_patterns.setEnabled(True)
        self._btn_adduct_report.setEnabled(True)
        self._btn_rt_outlier.setEnabled(True)
        self._btn_coherence.setEnabled(True)
        # トグルボタンを uncheck 状態に同期
        self._sync_adduct_button_state()

    def _toggle_is_filter(self):
        """③ IS Filter トグルハンドラ(両モード一括)。

        ボタンの新しい checked 状態に応じて:
          - ON: dialog を開いて IS Filter を適用(全モード対象)
          - OFF: 全モードで IS Filter を取り消し、非 IS マッチに戻す
        """
        if self._match_df is None or not self._lib_path:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            self._set_btn_checked_silent(self._btn_is_filter, False)
            return
        target_state = bool(self._btn_is_filter.isChecked())
        if target_state:
            self._apply_to_all_modes(self._apply_is_filter)
            # Filter Classes exempt 反映
            self._apply_filter_class_exemptions()
            # apply_is_filter で dialog を cancel された場合も同期し直す
            self._sync_step_button_states()
            self._draw_overlay()
        else:
            self._apply_to_all_modes(self._undo_is_filter)

    def _undo_is_filter(self):
        """③ IS Filter の取り消し(active mode で実行)。

        - lib_df を再ロードし、apply_is_filter=False で match_library_to_peaks を
          再実行(非 IS 通常マッチに戻す)
        - _reserved_coords / _runtime_is_fingerprints / _is_filter_choices をクリア
        - _class_ref_rt から auto_is エントリを削除(manual は尊重)
        - 後段 ④⑤⑥ の state もカスケードクリア
        - 散布図 / テーブル / フィルタを再構築
        """
        if not self._lib_path or self._match_df is None:
            return
        try:
            lib_df = load_lipid_library(
                self._lib_path, self.fe.ion_mode)
        except Exception as e:
            log.warning(f"[IS Filter undo] library reload failed: {e}")
            return

        # auto_is 由来のエントリだけ削除(manual はユーザー明示なので保持)
        self._class_ref_rt = {
            cls: vals for cls, vals in (self._class_ref_rt or {}).items()
            if not (len(vals) >= 3 and vals[2] == "auto_is")
        }

        # 非 IS マッチを再実行
        match_df, updated_ref_rt, _ = match_library_to_peaks(
            lib_df, self._mz, self._rt, self._spin_ppm.value(),
            void_rt=float(self._spin_void_rt.value()),
            class_ref_rt=self._class_ref_rt or None,
            apply_is_filter=False)
        self._match_df = match_df
        self._class_ref_rt = updated_ref_rt
        matched = match_df[match_df['matched']]

        # ③ 関連 state クリア
        self._reserved_coords = set()
        self._runtime_is_fingerprints = {}
        self._runtime_is_fingerprints_by_mode = {}
        self._is_filter_choices = {}

        # match coords 再計算
        self._matched_coords = set(
            zip(matched['obs_rt'].round(6), matched['obs_mz'].round(6)))

        # ② Conflict / ④ Outlier / ⑤ Adduct / ⑥ Coherence をカスケードクリア
        self._conflict_coords = set()
        self._conflict_map = {}
        self._adduct_attribution = {}
        self._outlier_coords = {}
        self._coherence_outliers = {}
        self._coherence_assignments = {}
        self._coherence_models = {}

        # UI 更新
        classes = match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        self._update_match_count_label(matched, lib_df)
        self._sync_step_button_states()
        log.info(f"[IS Filter] reverted (mode={self.fe.ion_mode})")

    def _apply_is_filter(self):
        """Step 2: ISベース2パスマッチングに切り替え、reserved_coords確立

        2 パスマッチングまで到達したかどうかを
        self._is_filter_applied_ok に残す。Run All はこれを見て、③ が
        適用されなかった場合に ④〜⑥ へ進まず中止する。① Match Overlay が
        _class_ref_rt を消したあとに ③ が早期 return すると、④⑤⑥ が基準 RT
        なしで走って黙って壊れるため。
        """
        self._is_filter_applied_ok = False
        if not self._lib_path or self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            return
        try:
            lib_df = load_lipid_library(
                self._lib_path, self.fe.ion_mode)
        except Exception as e:
            QMessageBox.critical(self, "Library Error", str(e))
            return

        # ── IS プレビュー: IS エントリのみ仮マッチング ───────────────
        is_lib_df = lib_df[lib_df['is_IS']].reset_index(drop=True)
        if is_lib_df.empty:
            QMessageBox.warning(
                self, "No IS Entries",
                "No IS entries found in the library. "
                "IS filter cannot be applied.")
            return

        # IS のみを通常マッチングして検出状況を取得
        # Fix RT の manual 設定があればそれを尊重（手動RT窓内でIS検出）
        # IS Filter Preview は 1 IS エントリ = 1 行 + Multiple
        # candidates ドロップダウンの仕様。multi_label=False で抑制。
        # 強度加味の自動選択(既定 OFF)。手動選択が最優先なのは不変。
        _use_int = False
        try:
            _use_int = bool(self._cb_is_auto_intensity.isChecked())
        except Exception:
            pass
        _peak_int = None
        if self._sample_intensities is not None:
            try:
                _si = np.asarray(self._sample_intensities, dtype=float)
                _si = np.where(_si > 0, _si, np.nan)
                with np.errstate(all='ignore'):
                    _peak_int = np.nanmean(_si, axis=1)
                _peak_int = np.where(np.isfinite(_peak_int), _peak_int, 0.0)
            except Exception:
                _peak_int = None

        is_match_df, _, _ = match_library_to_peaks(
            is_lib_df, self._mz, self._rt, self._spin_ppm.value(),
            void_rt=float(self._spin_void_rt.value()),
            class_ref_rt=self._class_ref_rt or None,
            apply_is_filter=False,
            multi_label=False,
            is_peak_intensities=_peak_int,
            is_use_intensity=_use_int)

        # 各 IS の ppm 内全候補ピークを列挙して is_match_df に付加
        # (混雑領域での誤マッチを防ぐ Candidate Picker 機能のデータ源)
        ppm_tol = self._spin_ppm.value()
        candidates_per_row = []
        for _, r in is_match_df.iterrows():
            cands = _collect_match_candidates(
                float(r['theoretical_mz']),
                self._mz, self._rt, ppm_tol,
                sample_intensities=self._sample_intensities,
                void_rt=float(self._spin_void_rt.value()),
            )
            candidates_per_row.append(cands)
        is_match_df = is_match_df.copy()
        is_match_df['candidates'] = candidates_per_row

        # auto-replay 中は ISPreviewDialog をスキップして
        # 復元済み _is_filter_choices を直接使う(ユーザー操作不要)
        # Run All 中も同じ扱いにする。Run All を 2 回目に押したときに
        #   ここでダイアログが開くと、① が消した _class_ref_rt を ③ が
        #   埋め直す前にキャンセルされ得る。保存済みの選択があるなら黙って
        #   使い、Run All が何回押されても同じ結果になるようにする。
        #   保存済み選択が無い初回はダイアログを出す(IS のピン留めは人が
        #   決めるべきなので、そこは自動化しない)。
        if ((getattr(self, '_auto_replay_in_progress', False)
             or getattr(self, '_run_all_in_progress', False))
                and self._is_filter_choices):
            log.info("[IS Filter] Run All / auto-replay: skip preview dialog, "
                  "use saved choices")
            # auto-replay ではダイアログも soft block 警告も出ない。
            #   session を読んで Run All するとき、複数候補を持つ IS が
            #   ピン留めされていなければ自動選択が黙って決めてしまう
            #   (S6 で pos/HexCer が実際にこれで誤選択されていた)。
            #   ダイアログは出さずに警告だけ出す。
            try:
                _multi = collect_multi_candidate_is(
                    self._is_filter_choices, is_lib_df, self._mz,
                    self._spin_ppm.value())
            except Exception as _e:
                _multi = []
                log.warning(f"[IS Filter] multi-candidate check failed: {_e}")
            # pos / neg を続けて処理するので、モードごとに持つ。
            # 単一の属性にすると後続モードの結果で上書きされてしまう。
            _mode = getattr(self.fe, 'ion_mode', None) or self._active_mode
            if not isinstance(getattr(self, '_pending_multi_candidate_is', None), dict):
                self._pending_multi_candidate_is = {}
            if not isinstance(getattr(self, '_multi_cand_warned', None), dict):
                self._multi_cand_warned = {}
            self._pending_multi_candidate_is[_mode] = list(_multi)
            if _multi:
                log.info(f"[IS Filter] {_mode}: multi-candidate IS not pinned: "
                      + ", ".join(_multi))
                if not self._multi_cand_warned.get(_mode):
                    self._multi_cand_warned[_mode] = True
                    QMessageBox.warning(
                        self, f"Multi-candidate IS not pinned [{_mode}]",
                        f"[{_mode}] {len(_multi)} internal standard(s) have several\n"
                        "candidate peaks inside the ppm window, but none was\n"
                        "chosen manually:\n\n"
                        + "\n".join(f"  - {s}" for s in _multi)
                        + "\n\nThese are decided by auto-pick. In a crowded region it\n"
                          "can grab the wrong peak and shift the RT window of the\n"
                          "class.\n\n"
                          "Run ③ IS Filter manually, pick the right peak in the\n"
                          "Candidate Picker, and save the session again.")

        else:
            # プレビューダイアログを表示（前回の選択を復元）
            prev_choices = self._is_filter_choices
            # ion_mode をタイトルに表示
            # Candidate Picker で EIC を出すために生データのパスを渡す。
            # 取り込み経路によっては存在しないので、その場合は空リスト。
            try:
                _raw = resolve_entry_raw_paths(self.fe)
            except Exception as _e:
                log.warning(f"raw path resolve failed: {_e}")
                _raw = []
            dlg = ISPreviewDialog(is_match_df, prev_choices,
                                  sample_col_names=self._sample_col_names,
                                  ion_mode=self.fe.ion_mode,
                                  raw_paths=_raw,
                                  parent=self)
            if not dlg.exec():
                return  # Cancel → 何もしない

            # 複数候補があるが未確認の IS があれば soft block 警告
            warn, unconfirmed = dlg.has_unconfirmed_multi_candidates()
            if warn:
                msg = (
                    f"{len(unconfirmed)} internal standard(s) have multiple "
                    f"candidate peaks within ppm tolerance:\n\n"
                    + "\n".join(f"  - {s}" for s in unconfirmed)
                    + "\n\nWithout explicit selection, the m/z-closest "
                    "candidate will be used,\nwhich may pick an incorrect "
                    "peak in crowded regions.\n\n"
                    "Strongly recommend reviewing each candidate before "
                    "applying the IS filter."
                )
                mb = QMessageBox(self)
                mb.setIcon(QMessageBox.Warning)
                mb.setWindowTitle("Multi-candidate IS not reviewed")
                mb.setText(msg)
                btn_review = mb.addButton("Cancel and Review",
                                          QMessageBox.RejectRole)
                btn_proceed = mb.addButton("Proceed Anyway",
                                           QMessageBox.AcceptRole)
                mb.setDefaultButton(btn_review)
                mb.exec()
                if mb.clickedButton() == btn_review:
                    # ISPreviewDialog を再度開いてレビューしてもらう
                    # (再帰的に _apply_is_filter を呼び出す)
                    return self._apply_is_filter()
                # Proceed Anyway → そのまま続行

            self._is_filter_choices = dlg.result_choices()

        # 承認された IS 情報で2パスマッチング実行
        match_df, updated_ref_rt, is_missing = match_library_to_peaks(
            lib_df, self._mz, self._rt, self._spin_ppm.value(),
            void_rt=float(self._spin_void_rt.value()),
            class_ref_rt=self._class_ref_rt or None,
            is_tol=self._spin_is_tol.value(),
            apply_is_filter=True,
            is_choices=self._is_filter_choices,
            is_peak_intensities=_peak_int,
            is_use_intensity=_use_int)
        self._match_df = match_df
        self._class_ref_rt = updated_ref_rt
        self._is_filter_applied_ok = True   # ここまで来たら適用済み
        matched = match_df[match_df['matched']]

        # クラス色は維持（既に割り当て済み）、変更があった場合のみ補完
        for cls in match_df['lipid_class'].unique():
            if cls not in self._class_colors:
                self._class_colors[cls] = '#888888'

        self._matched_coords = set(
            zip(matched['obs_rt'].round(6), matched['obs_mz'].round(6)))

        # reserved_coords: IS優先クラス（自動または手動）で帰属された座標
        is_bearing = {
            cls for cls, vals in self._class_ref_rt.items()
            if vals[2] in ("auto_is", "manual")
        }
        reserved_df = matched[matched['lipid_class'].isin(is_bearing)]
        self._reserved_coords = set(
            zip(reserved_df['obs_rt'].round(6),
                reserved_df['obs_mz'].round(6)))

        # Review で再計算されるまで、RT outlier / Coherence outlier はリセット
        # 新パイプライン順序 (② Conflict → ③ IS) では
        # Conflict Search が IS Filter より先に実行されているため、
        # _conflict_coords をクリアせず、IS Filter 後の新しい match_df で
        # 衝突を再計算する。
        self._outlier_coords  = {}
        self._coherence_outliers = {}
        try:
            self._conflict_coords, self._conflict_map = calc_conflicts(
                self._match_df)
        except Exception as exc:
            log.warning(f"[IS Filter] re-calc_conflicts failed: {exc}")
            self._conflict_coords = set()
            self._conflict_map = {}

        # ② IS Filter 完了後、IS 由来アダクトパターンを生成
        # 内蔵 CURATED_ADDUCT_FINGERPRINTS が空(オプション 2)の場合に
        # ⑤ Adduct Ion Filter が参照する。
        try:
            self._compute_is_fingerprints(lib_df)
        except Exception as exc:
            log.warning(f"[Adduct Ion Filter] _compute_is_fingerprints failed: {exc}")
            self._runtime_is_fingerprints = {}
            self._runtime_is_fingerprints_by_mode = {}

        # IS Filter 再実行で Adduct Ion Filter 結果は無効化
        self._adduct_attribution = {}
        self._sync_adduct_button_state()

        classes = match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        self._update_match_count_label(matched, lib_df)

        if is_missing:
            QMessageBox.warning(
                self, "IS Not Detected",
                "The following IS-bearing classes fell back to "
                "standard m/z-closest matching because their IS "
                "entries did not match any peak:\n\n"
                + ", ".join(sorted(is_missing)))

    # ── ランタイム IS 由来アダクトパターン生成 ──────────────────
    def _get_other_mode_data(self):
        """他モード(現在の dialog が pos なら neg)のピーク表データを取得。

        まず self._all_entries(__init__ で受け取った FileEntry リスト)から
        探し、それでも見つからなければ親 MainWindow の self._file_entries を
        探す。rt / mz / mean intensity 配列を返す。見つからない or
        読み込めない場合は (None, None, None) を返す。
        """
        other_mode = 'neg' if self.fe.ion_mode == 'pos' else 'pos'
        # 1) self._all_entries(__init__ で渡された全 FileEntry)
        file_entries = list(getattr(self, '_all_entries', []) or [])
        if not file_entries:
            # 2) 親 MainWindow の _file_entries (_host 経由)
            parent = self._host
            file_entries = list(getattr(parent, '_file_entries', None) or [])
        if not file_entries:
            return None, None, None
        other_fe = next(
            (f for f in file_entries
             if getattr(f, 'ion_mode', None) == other_mode),
            None
        )
        if other_fe is None:
            return None, None, None
        try:
            rt_arr = other_fe.rt_array()
            mz_arr = other_fe.mz_array()
            if rt_arr is None or mz_arr is None:
                return None, None, None
            mask = ~(np.isnan(rt_arr) | np.isnan(mz_arr))
            rt_v = rt_arr[mask]
            mz_v = mz_arr[mask]
            sample_cols = other_fe.sample_columns()
            if not sample_cols:
                return rt_v, mz_v, None
            ints = (other_fe.df[sample_cols]
                    .apply(pd.to_numeric, errors="coerce")
                    .values[mask].astype(float))
            mean_int = np.nanmean(ints, axis=1)
            return rt_v, mz_v, mean_int
        except Exception:
            return None, None, None

    def _rebuild_is_fingerprints_now(self):
        """現モードの IS アダクト指紋を、今の state で作り直す。

        _compute_is_fingerprints はクロスモード検索に pos↔neg RT シフトを
        使うが、それは相手モードの _class_ref_rt が埋まっていて初めて
        求まる。③ の中で呼ばれるため、Run All の Phase 1 で先に走る
        モードの指紋はシフト既定値(±0.05 min)で作られてしまい、
        あとから走るモードと条件が揃わない。

        Run All の Phase 2 冒頭でこれを呼ぶことで、両モードの ③ が
        終わった状態、つまり両モードの class_ref_rt が揃った状態で
        指紋を作り直す。
        """
        if not self._lib_path or self._match_df is None:
            return
        try:
            lib_df = load_lipid_library(self._lib_path, self.fe.ion_mode)
        except Exception as exc:
            log.warning("[Run All] fingerprint rebuild: library reload failed: "
                  f"{exc}")
            return
        try:
            self._compute_is_fingerprints(lib_df)
        except Exception as exc:
            log.warning(f"[Run All] fingerprint rebuild failed: {exc}")

    def _compute_is_fingerprints(self, lib_df):
        """② IS Filter 後に、IS 由来のアダクト fingerprint を生成。

        対象クラス: self._class_ref_rt のうち source == 'auto_is' のもの
                    (= IS が観測されたクラス)。

        各クラスについて、IS の neutral mass を起点に 9 アダクトの理論 m/z
        を計算し、自モード / 他モードのピーク表で ppm + RT 窓内に観測がある
        かを確認。観測ありなら expected=True、weight=最強アダクトに対する
        強度比、観測なしなら expected=False、weight=0.0。

        結果は self._runtime_is_fingerprints_by_mode に
        (クラス, 極性) キーで保存される:
          {(class_name, mode): {adduct: {'expected': bool, 'weight': float,
                                         'ratio_target': None,
                                         'source_mode': str}}}
        表示用に、クラス名だけをキーにしたビュー
        self._runtime_is_fingerprints も作り直す。

        fix51 の変更点:
          - キーに極性を含めた。クラス名だけだと、両極性に IS があるクラス
            (検証データでは HexCer / Cer) が後に走った極性で上書きされる。
          - クロスモード検索に ⑤ と同じ pos↔neg RT 補正を入れた。
          - expected / weight の判定を ⑤ が採点で使う基準で行う(自己検査)。
            中心 m/z は IS の観測ピークから逆算した neutral、RT は IS の
            実測 RT、許容誤差は ⑤ のもの。緩い条件でしか見つからない
            エントリは載せない。
            検証データでは IS 自身の同極性副アダクト 21 件のうち 7 件が
            「指紋にあるのに同条件で再検索すると見つからない」状態だった。
        """
        # 既存の fingerprints を保持(両モード処理で上書きされないように)
        if not hasattr(self, '_runtime_is_fingerprints') or self._runtime_is_fingerprints is None:
            self._runtime_is_fingerprints = {}
        if (not hasattr(self, '_runtime_is_fingerprints_by_mode')
                or self._runtime_is_fingerprints_by_mode is None):
            self._runtime_is_fingerprints_by_mode = {}
        if not self._class_ref_rt:
            return
        ppm_tol = float(self._spin_ppm.value())
        rt_win = float(self._spin_is_tol.value())
        # expected の判定は ⑤ が採点で使う許容誤差で行う。
        # 指紋は「⑤ が探して見つけられるもの」だけを載せる約束にする。
        try:
            ppm_chk = float(self._spin_adduct_ppm.value())
        except Exception:
            ppm_chk = ppm_tol
        try:
            rt_chk = float(self._spin_adduct_rt.value())
        except Exception:
            rt_chk = rt_win
        # クロスモード検索に ⑤ と同じ pos↔neg RT 補正を入れる。
        # fix50 までは補正なしで探していたため、⑤ 側と基準が違っていた。
        try:
            _shift, _sstd, _nis = self._compute_pos_neg_rt_shift()
        except Exception:
            _shift, _sstd, _nis = 0.0, 0.0, 0
        cross_win = max(0.01, 3.0 * _sstd) if _nis > 0 else rt_win
        cross_chk = max(0.01, 3.0 * _sstd) if _nis > 0 else rt_chk
        _dropped: list[str] = []

        # IS のみ抽出(同じクラスに複数 IS がある場合は最初のもの)
        is_entries = lib_df[lib_df['is_IS']]
        if is_entries.empty:
            return

        # 自モードの mean intensity 配列
        my_mean_int: np.ndarray | None
        if self._sample_intensities is not None:
            my_mean_int = np.nanmean(self._sample_intensities, axis=1)
        else:
            my_mean_int = None

        # 他モードのデータを 1 度だけ取得してキャッシュ
        other_rt, other_mz, other_mean_int = self._get_other_mode_data()

        # 極性は表示中のモードで決める。self.fe は _load_mode_state で
        # 差し替わる約束だが、差し替えに失敗した場合に備えて _active_mode を
        # 優先する。
        my_mode = (getattr(self, '_active_mode', None)
                   or getattr(self.fe, 'ion_mode', None) or 'pos')
        pos_set = set(ADDUCT_SET_POS_V54)

        for cls, vals in self._class_ref_rt.items():
            if vals[2] != 'auto_is':
                continue  # IS が検出されたクラスのみ
            ref_rt = float(vals[0])

            # 当該クラスの IS エントリを取得(複数あれば最初)
            cls_is = is_entries[is_entries['lipid_class'] == cls]
            if cls_is.empty:
                continue
            is_row = cls_is.iloc[0]
            try:
                exact_mass = float(is_row['exact_mass'])
            except (KeyError, TypeError, ValueError):
                continue

            # 自己検査の基準は ⑤ が採点で使うものに揃える。
            #   ⑤ は「観測ピークの m/z から逆算した neutral + offset」を
            #   「衝突ピークの実測 RT」で探す。理論質量とクラス基準 RT で
            #   検査すると、⑤ が見つけられないエントリを通してしまう。
            neutral_chk, ref_obs = exact_mass, ref_rt
            try:
                _mdf = self._match_df
                if _mdf is not None and 'is_IS' in _mdf.columns:
                    _g = _mdf[(_mdf['matched'] == True)
                              & (_mdf['is_IS'] == True)
                              & (_mdf['lipid_class'].astype(str) == str(cls))]
                    if not _g.empty:
                        _r = _g.iloc[0]
                        _a0 = ADDUCT_ALIAS.get(str(_r['adduct']),
                                               str(_r['adduct']))
                        _o0 = ADDUCT_OFFSETS.get(_a0)
                        if _o0 is not None:
                            neutral_chk = float(_r['obs_mz']) - _o0
                            ref_obs = float(_r['obs_rt'])
            except Exception:
                pass

            # 9 アダクトそれぞれを観測判定
            adduct_results: dict[str, dict[str, object]] = {}
            intensities: dict[str, float] = {}
            _n_loose = 0          # 緩い基準で見つかった本数（ログ用）
            for adduct in ADDUCT_SET_V54:
                offset = ADDUCT_OFFSETS.get(adduct)
                if offset is None:
                    continue
                mz_theo = exact_mass + offset
                adduct_mode = 'pos' if adduct in pos_set else 'neg'

                if adduct_mode == my_mode:
                    rt_arr = self._rt
                    mz_arr = self._mz
                    int_arr = my_mean_int
                    target_rt = ref_rt
                    tol_win, tol_chk = rt_win, rt_chk
                else:
                    rt_arr = other_rt
                    mz_arr = other_mz
                    int_arr = other_mean_int
                    # ⑤ と同じ向きに補正する
                    target_rt = (ref_rt + _shift if my_mode == 'pos'
                                 else ref_rt - _shift)
                    tol_win, tol_chk = cross_win, cross_chk

                if rt_arr is None or mz_arr is None or len(rt_arr) == 0:
                    adduct_results[adduct] = {
                        'expected': False, 'weight': 0.0,
                        'ratio_target': None,
                    }
                    continue

                ppm_err = np.abs(mz_arr - mz_theo) / mz_theo * 1e6
                rt_err = np.abs(rt_arr - target_rt)
                hit_mask = (ppm_err <= ppm_tol) & (rt_err <= tol_win)

                # 自己検査。⑤ が採点で使う基準（観測由来の neutral と
                # IS の実測 RT、⑤ の許容誤差）でも見つかるか。
                mz_chk = neutral_chk + offset
                if adduct_mode == my_mode:
                    rt_chk_target = ref_obs
                else:
                    rt_chk_target = (ref_obs + _shift if my_mode == 'pos'
                                     else ref_obs - _shift)
                chk_mask = (
                    (np.abs(mz_arr - mz_chk) / mz_chk * 1e6 <= ppm_chk)
                    & (np.abs(rt_arr - rt_chk_target) <= tol_chk))

                if np.any(hit_mask):
                    _n_loose += 1
                if not np.any(chk_mask):
                    if np.any(hit_mask):
                        # 緩い条件では見つかるが ⑤ の条件では見つからない
                        _dropped.append(f'{cls}[{my_mode}]/{adduct}')
                    adduct_results[adduct] = {
                        'expected': False, 'weight': 0.0,
                        'ratio_target': None,
                    }
                    continue

                # 候補のうち最強強度を採用(intensity 不明なら ppm 最近接)
                cand_idx = np.where(chk_mask)[0]
                if int_arr is not None:
                    intensity = float(np.nanmax(int_arr[cand_idx]))
                else:
                    intensity = 1.0  # 強度情報なし、placeholder

                adduct_results[adduct] = {
                    'expected': True, 'weight': intensity,
                    'ratio_target': None,
                    'source_mode': my_mode,
                }
                intensities[adduct] = intensity

            # weight を「最強アダクト = 1.0」基準に正規化
            if intensities:
                max_int = max(intensities.values())
                if max_int > 0:
                    for adu, spec in adduct_results.items():
                        if spec['expected']:
                            spec['weight'] = (
                                intensities.get(adu, 0.0) / max_int)

            # (クラス, 極性) で保存する。表示用のマージ済みビューは
            # 下でまとめて作り直す。
            self._runtime_is_fingerprints_by_mode[(cls, my_mode)] = \
                adduct_results
            _n_exp = sum(1 for _s in adduct_results.values()
                         if _s.get('expected'))
            if _n_exp != _n_loose:
                log.info(f"[IS Filter] fingerprint {cls}[{my_mode}]: "
                      f"expected {_n_exp} 本 (緩い基準では {_n_loose} 本)")

        # 表示用ビューを作り直す(同極性を優先、無ければ他極性を流用)
        self._rebuild_fingerprint_view()
        if _dropped:
            log.info("[IS Filter] fingerprint entries dropped by self-check "
                  f"({len(_dropped)}): " + ", ".join(_dropped))

    def _rebuild_fingerprint_view(self):
        """(クラス, 極性) のストアから、表示用の {クラス: spec} を作る。

        Patterns ダイアログや Utility の heatmap は従来どおりクラス名だけを
        キーにした辞書を期待するため、表示用のビューを別に持つ。判定に
        使うのは `_fingerprints_for_mode()` の方。
        """
        by = getattr(self, '_runtime_is_fingerprints_by_mode', None) or {}
        view: dict = {}
        cur = getattr(self, '_active_mode', None) or 'pos'
        for (cls, m), spec in by.items():
            if m == cur:
                view[cls] = spec
        for (cls, m), spec in by.items():
            view.setdefault(cls, spec)
        self._runtime_is_fingerprints = view

    def _fingerprints_for_mode(self, mode: str) -> dict:
        """指定極性で ⑤ が使う指紋を返す。

        その極性で IS が検出されたクラスの指紋を優先し、片極性しか
        無いクラスは他極性のものを流用する(どのアダクトを出すかは分子の
        性質で、run の違いは基準 RT にしか効かないため)。流用したものは
        spec の 'source_mode' で区別できる。
        """
        by = getattr(self, '_runtime_is_fingerprints_by_mode', None) or {}
        out: dict = {}
        for (cls, m), spec in by.items():
            if m == mode:
                out[cls] = spec
        for (cls, m), spec in by.items():
            out.setdefault(cls, spec)
        if not out:
            # fix51 以前の状態(ビューだけある)に備えたフォールバック
            out = dict(getattr(self, '_runtime_is_fingerprints', None) or {})
        return out

    def _toggle_conflict_search(self):
        """② Conflict Search トグルハンドラ(両モード一括)。

        ボタンの新しい checked 状態に応じて apply / undo を全モードで実行する。
        """
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            self._set_btn_checked_silent(self._btn_conflict, False)
            return
        target_state = bool(self._btn_conflict.isChecked())
        if target_state:
            self._apply_to_all_modes(self._conflict_search)
        else:
            self._apply_to_all_modes(self._undo_conflict_search)

    def _undo_conflict_search(self):
        """② Conflict Search の取り消し(active mode で実行)。

        - self._conflict_coords / _conflict_map をクリア
        - ⑤ adduct attribution / adduct_filter_status は ② に依存しているため
          カスケードクリア
        - 散布図 / テーブルを再描画
        """
        self._conflict_coords = set()
        self._conflict_map = {}
        # ⑤ Adduct Ion Filter は conflicts に依存 → カスケードクリア
        self._adduct_attribution = {}
        if self._match_df is not None:
            match_df = self._match_df.copy()
            if 'adduct_filter_status' in match_df.columns:
                matched_mask = match_df['matched']
                match_df.loc[matched_mask, 'adduct_filter_status'] = STATUS_KEPT
                match_df.loc[~matched_mask, 'adduct_filter_status'] = STATUS_NA
            match_df = _update_final_status_df(match_df)
            self._match_df = match_df
            matched = match_df[match_df['matched']]
            classes = match_df['lipid_class'].unique().tolist()
            self._rebuild_filters(classes)
            self._draw_overlay()
            self._rebuild_match_table(matched)
        self._sync_step_button_states()
        log.info(f"[Conflict Search] reverted (mode={self.fe.ion_mode})")

    def _conflict_search(self):
        """Step ② 競合マッチ検出。冪等(同じ入力で同じ結果)。

        旧来のトグル動作を一旦廃止していたが、
        ボタンを setCheckable に戻し、_toggle_conflict_search 経由で apply
        専用の役割に。undo は _undo_conflict_search に分離されている。
        """
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            return

        # 常に再計算(冪等)
        self._conflict_coords, self._conflict_map = calc_conflicts(
            self._match_df)

        # 衝突集合が変わるため Adduct Ion Filter 結果は無効化
        self._adduct_attribution = {}
        self._sync_step_button_states()

        # 右パネルのチェックボックスを含めて再構築し、散布図とテーブルを更新
        matched = self._match_df[self._match_df['matched']]
        classes = self._match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)

    # ── Step 4 Adduct Ion Filter ───────────────────────────────
    def _compute_pos_neg_rt_shift(self) -> tuple[float, float, int]:
        """pos と neg の _class_ref_rt(IS から推定された
        クラス別 reference RT)を比較して、両モードに共通する IS class の
        RT 差から系統的 pos→neg shift を推定する。

        戻り値:
          (median_shift_min, sigma_shift_min, n_used_is)
          median_shift > 0 は neg が pos より遅く溶出していることを意味する
          (= neg_rt = pos_rt + shift)。
          n_used_is = 0 のとき shift=0, sigma=0.05 (fallback)。

        σ を np.std から **中央値 + MAD の頑健推定** に替え、
        中央値から 3σ を超える組を外れ値として除いてから再推定する。
        fix51 までは外れ値 1 つで σ が膨らみ、クロスモード窓が
        ±0.29 min まで広がって RT をほぼ見ない状態になっていた
        (検証データの 13 組のうち +0.018 / +0.164 / +0.330 が外れ値)。
        σ には POS_NEG_SHIFT_SIGMA_FLOOR を下限として置く。
        """
        # 自モードの _class_ref_rt は self._class_ref_rt、他モードは _mode_state
        active_mode = self._active_mode
        other_mode = 'neg' if active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}

        active_ref = self._class_ref_rt or {}
        other_ref = other_state.get('_class_ref_rt') or {}

        if active_mode == 'pos':
            pos_ref, neg_ref = active_ref, other_ref
        else:
            pos_ref, neg_ref = other_ref, active_ref

        if not pos_ref or not neg_ref:
            # n=0 になった理由が分かるようにする
            log.info("[Adduct Ion Filter] pos→neg shift: class_ref_rt が片側に "
                  f"無い (pos {len(pos_ref)} / neg {len(neg_ref)} クラス) "
                  "→ shift=0、クロスモード窓は既定値")
            return 0.0, 0.05, 0

        shifts = []
        for cls in pos_ref:
            if cls not in neg_ref:
                continue
            pos_entry = pos_ref[cls]
            neg_entry = neg_ref[cls]
            # auto_is 由来のものだけを使う(manual は人為設定なので除外)
            if (len(pos_entry) >= 3
                    and len(neg_entry) >= 3
                    and pos_entry[2] == 'auto_is'
                    and neg_entry[2] == 'auto_is'):
                shifts.append(float(neg_entry[0]) - float(pos_entry[0]))

        if not shifts:
            # 共通クラスはあるが auto_is 由来の組が無い場合
            _common = sorted(set(pos_ref) & set(neg_ref))
            log.info("[Adduct Ion Filter] pos→neg shift: 共通クラス "
                  f"{len(_common)} 件だが auto_is 由来の組が無い"
                  + (f" ({', '.join(_common[:8])})" if _common else "")
                  + " → shift=0、クロスモード窓は既定値")
            return 0.0, 0.05, 0

        # 中央値 + MAD の頑健推定。3σ を超える組は外れ値として除く。
        _arr = np.asarray(shifts, dtype=float)
        _med0 = float(np.median(_arr))
        _mad0 = float(np.median(np.abs(_arr - _med0))) * 1.4826
        if not np.isfinite(_mad0) or _mad0 <= 0.0:
            _mad0 = POS_NEG_SHIFT_SIGMA_FLOOR
        _thr = 3.0 * max(_mad0, POS_NEG_SHIFT_SIGMA_FLOOR)
        _keep = np.abs(_arr - _med0) <= _thr
        _out = _arr[~_keep]
        _kept = _arr[_keep] if bool(_keep.any()) else _arr
        median = float(np.median(_kept))
        _mad = float(np.median(np.abs(_kept - median))) * 1.4826
        sigma = max(_mad if np.isfinite(_mad) else 0.0,
                    POS_NEG_SHIFT_SIGMA_FLOOR)

        detail = ", ".join(f"{s:+.4f}" for s in sorted(_kept))
        log.info(f"[Adduct Ion Filter] IS class shifts (neg − pos, "
              f"n={len(_kept)}/{len(_arr)}): [{detail}]")
        if len(_out):
            log.info("[Adduct Ion Filter] shift outliers excluded ("
                  f"|Δ−median| > {_thr:.4f} min): "
                  + ", ".join(f"{s:+.4f}" for s in sorted(_out)))
        return median, sigma, int(len(_kept))

    def _apply_adduct_filter(self):
        """Step 4: Adduct Ion Filter 本体。

        衝突候補に IS あり クラスを 1 つ以上含むスポットで、観測アダクト
        パターンから 1 クラスを採用する。判定不能(または衝突候補に IS なし
        クラスが混在)の場合は何もせず、後続の Coherence Filter に流す。

        前提: ① Match Overlay、② IS Filter、③ Conflict Search が既に実行
              されていること。② IS Filter 後に self._runtime_is_fingerprints
              が _compute_is_fingerprints() で生成されている必要がある。
        """
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            return
        if not self._conflict_coords:
            # 衝突なし or Conflict Search 未実行 → 何もしない
            self._adduct_attribution = {}
            return

        # 期待パターンを統合(内蔵 CURATED_ADDUCT_FINGERPRINTS が空の場合は
        # ランタイム IS 派生のみ。両方ある場合は内蔵を優先)
        # その極性で作られた指紋を使う(クラス名だけのビューは
        # 両極性に IS があるクラスで上書きが起きるため使わない)。
        fingerprints: dict = self._fingerprints_for_mode(
            getattr(self, '_active_mode', None)
            or getattr(self.fe, 'ion_mode', None) or 'pos')
        for cls, fp in CURATED_ADDUCT_FINGERPRINTS.items():
            fingerprints[cls] = fp  # 内蔵で上書き

        if not fingerprints:
            QMessageBox.information(
                self, "No Fingerprints",
                "No adduct fingerprints available.\n"
                "Please run 'Apply IS Filter' first to derive runtime "
                "fingerprints from IS observations, or populate "
                "CURATED_ADDUCT_FINGERPRINTS in LipidZoner.py.")
            return

        # ピーク表データ(自モード + 他モード)を準備
        # 強度配列も取得して heatmap グラデーション用に流す
        my_int = None
        if self._sample_intensities is not None:
            try:
                my_int = np.nanmean(self._sample_intensities, axis=1)
            except Exception:
                my_int = None
        if self.fe.ion_mode == 'pos':
            pos_rt, pos_mz = self._rt, self._mz
            pos_int = my_int
            other_rt, other_mz, other_int = self._get_other_mode_data()
            neg_rt, neg_mz, neg_int = other_rt, other_mz, other_int
        else:
            neg_rt, neg_mz = self._rt, self._mz
            neg_int = my_int
            other_rt, other_mz, other_int = self._get_other_mode_data()
            pos_rt, pos_mz, pos_int = other_rt, other_mz, other_int

        # Adduct Ion Filter 専用パラメータを使用
        try:
            ppm_tol = float(self._spin_adduct_ppm.value())
        except Exception:
            ppm_tol = float(self._spin_ppm.value())
        try:
            rt_win = float(self._spin_adduct_rt.value())
        except Exception:
            rt_win = float(self._spin_is_tol.value())
        try:
            int_floor = float(self._spin_adduct_int.value())
        except Exception:
            int_floor = 0.0
        # vote_threshold を spinbox から読む(default=1.0)
        try:
            vote_threshold = float(self._spin_adduct_vote.value())
        except Exception:
            vote_threshold = 1.0

        # 診断ログ: fingerprint 状態を出力
        fp_classes = sorted(fingerprints.keys())
        log.info(f"[Adduct Ion Filter] mode={self.fe.ion_mode}, "
              f"params: ppm={ppm_tol}, rt={rt_win}, int={int_floor}, "
              f"vote={vote_threshold}")
        log.info(f"[Adduct Ion Filter] fingerprints classes "
              f"({len(fp_classes)}): {fp_classes}")
        # PC/PE が両方含まれているか確認
        for must in ('PC', 'PE'):
            if must in fingerprints:
                fp = fingerprints[must]
                expected_adducts = [
                    a for a, s in fp.items() if s.get('expected')]
                log.info(f"  [{must}] expected adducts: {expected_adducts}")
            else:
                log.info(f"  [{must}] FINGERPRINT MISSING - "
                      "may not have been generated by the IS Filter")

        # IS から pos↔neg RT shift を推定
        shift, shift_std, n_is = self._compute_pos_neg_rt_shift()
        cross_rt_win = max(0.01, 3.0 * shift_std) if n_is > 0 else rt_win
        if n_is > 0:
            log.info(f"[Adduct Ion Filter] IS-derived pos→neg RT shift: "
                  f"{shift:+.4f} ± {shift_std:.4f} min  (n={n_is})  "
                  f"→ cross-mode tolerance = ±{cross_rt_win:.4f} min")
        else:
            log.info(f"[Adduct Ion Filter] no common IS in pos/neg, "
                  f"cross-mode tolerance fallback = ±{cross_rt_win:.4f} min")

        # スコアリング実行
        attribution = calc_adduct_filter(
            self._conflict_coords, self._conflict_map, fingerprints,
            pos_rt, pos_mz, neg_rt, neg_mz,
            ppm_tol, rt_win, int_floor, vote_threshold,
            conflict_mode=self.fe.ion_mode,
            cross_mode_rt_shift=shift,
            cross_mode_rt_win=cross_rt_win,
            pos_int=pos_int, neg_int=neg_int,
        )
        self._adduct_attribution = attribution

        # undecided / decided 詳細ログ
        n_decided = sum(
            1 for r in attribution.values()
            if not r.get('skipped') and r.get('winner') is not None)
        n_undecided = sum(
            1 for r in attribution.values()
            if not r.get('skipped') and r.get('winner') is None)
        n_skipped = sum(
            1 for r in attribution.values() if r.get('skipped'))
        log.info(f"[Adduct Ion Filter] result: decided={n_decided}, "
              f"undecided={n_undecided}, skipped={n_skipped}")

        undecided_samples = []
        for coord, result in attribution.items():
            if result.get('skipped'):
                continue
            if result.get('winner') is None:
                undecided_samples.append((coord, result))
            if len(undecided_samples) >= 8:
                break
        if undecided_samples:
            log.info(f"[Adduct Ion Filter] undecided detail "
                  f"(showing up to 8):")
            for coord, result in undecided_samples:
                rt, mz = coord
                scores = result.get('scores', {})
                disc_map = result.get('discriminator_map', {})
                per_cls = result.get('per_class_observations', {})
                cls_to_adu = result.get('class_to_adduct', {})
                log.info(
                    f"  ───  RT={rt:.3f}  mz={mz:.4f}  ───")
                log.info(
                    f"    candidates: "
                    + ", ".join(f"{c}({cls_to_adu.get(c, '?')})"
                                for c in scores))
                if not disc_map:
                    log.info(
                        "    no discriminating adduct: every candidate expects the "
                        "same adduct → chemically indistinguishable, left to the "
                        "Coherence Filter")
                else:
                    log.info(
                        "    discriminating adducts vs observation:")
                    for adu, expects_cls in disc_map.items():
                        # この adduct について、各候補がどう観測したか
                        obs_str = ", ".join(
                            f"{c}={'O' if per_cls.get(c, {}).get(adu) else 'X'}"
                            for c in scores)
                        log.info(
                            f"      {adu} (expected: {'/'.join(expects_cls)}) "
                            f"→ {obs_str}")
                log.info(f"    scores: {scores}")
                # 同点 or 0 vs 0 の理由を明記
                top = max(scores.values()) if scores else 0
                tied = [c for c, s in scores.items() if s == top]
                if len(tied) > 1:
                    log.info(
                        f"    → undecided: {len(tied)} candidates tied at the top "
                        f"score {top} ({', '.join(tied)})")
                else:
                    second = sorted(scores.values(), reverse=True)
                    diff = (second[0] - second[1]) if len(second) > 1 else 0
                    log.info(
                        f"    → undecided: top - second = {diff} < "
                        f"vote_threshold={vote_threshold}")

        # match_df の adduct_filter_status 列を更新
        match_df = self._match_df.copy()
        if 'adduct_filter_status' not in match_df.columns:
            match_df['adduct_filter_status'] = STATUS_NA
        # まず matched 行を kept、unmatched は n/a に揃える
        matched_mask = match_df['matched']
        match_df.loc[matched_mask, 'adduct_filter_status'] = STATUS_KEPT
        match_df.loc[~matched_mask, 'adduct_filter_status'] = STATUS_NA

        # 衝突解決の loser を rejected_adduct に
        n_decided = 0
        n_undecided = 0
        n_skipped = 0
        for coord, result in attribution.items():
            if result.get('skipped'):
                n_skipped += 1
                continue
            winner = result.get('winner')
            if winner is None:
                n_undecided += 1
                continue
            n_decided += 1
            rt, mz = coord
            coord_mask = (
                (match_df['obs_rt'].round(6) == rt)
                & (match_df['obs_mz'].round(6) == mz)
                & matched_mask
            )
            for cls in result.get('losers', set()):
                cls_mask = coord_mask & (match_df['lipid_class'] == cls)
                match_df.loc[cls_mask, 'adduct_filter_status'] = (
                    STATUS_REJ_ADDUCT)

        # final_status を再計算
        match_df = _update_final_status_df(match_df)
        self._match_df = match_df

        # 新順序では RT Outlier(④)は Adduct(⑤)より上流。よって Adduct
        # 適用時に RT Outlier 結果(_outlier_coords / rt_outlier_status)を
        # クリアしてはならない(以前の版が適用側の消去を消し忘れていたバグ。
        # calc_rt_outliers と calc_adduct_filter は互いに独立)。クリアするのは
        # 下流の Coherence のみ。rt_outlier_status 列は match_df 上で保持される。
        self._coherence_outliers = {}

        # 散布図とテーブル更新
        matched = match_df[match_df['matched']]
        # 右パネルの再構築(adduct_filter_status の有無を反映)
        classes = match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)

        # ログ
        log.info(f"[Adduct Ion Filter] decided: {n_decided}, "
              f"undecided: {n_undecided}, "
              f"skipped (no fingerprint): {n_skipped}, "
              f"total conflicts: {len(self._conflict_coords)}")

        # 適用後はトグルボタンを checked に同期
        self._sync_adduct_button_state()

    # ── トグル動作のための apply / undo / router ──
    def _set_btn_checked_silent(self, btn, value: bool):
        """signal を blocking しつつ btn.setChecked(value) を実行する小ヘルパ。"""
        if btn is None:
            return
        was_blocked = btn.blockSignals(True)
        try:
            btn.setChecked(bool(value))
        finally:
            btn.blockSignals(was_blocked)

    def _sync_adduct_button_state(self):
        """互換: ⑤ ボタンのみを同期(従来 API)。内部で _sync_step_button_states
        を呼ぶことで ②③④⑤⑥ 全てを一度に同期する。"""
        self._sync_step_button_states()

    def _on_match_clicked(self):
        """① Match Overlay ボタン押下時のハンドラ。

        - 常に _run_match を両モードで実行(checked/unchecked に関わらず)
        - checked 状態は _sync_step_button_states で _match_df 有無から導出
        - クリックによる Qt の自動 toggle で一時的に状態がずれる場合があるため
          実行後に必ず同期し直す
        """
        self._apply_to_all_modes(self._run_match)
        try:
            self._sync_step_button_states()
        except Exception:
            pass

    def _sync_step_button_states(self):
        """①②③④⑤⑥ 各ステップボタンの checked 状態を、
        対応する内部 state の有無に同期する。signal は blockSignals で抑制。

        判定基準:
          ① Match Overlay     : self._match_df is not None and not empty
          ② Conflict Search   : bool(self._conflict_coords)
          ③ IS Filter         : bool(self._reserved_coords)
          ⑤ Adduct Ion Filter : bool(self._adduct_attribution)
          ④ RT Outlier Filter : 任意の lipid_class で outliers が存在
          ⑥ Coherence Filter  : models / assignments / outliers のいずれかが存在
        """
        # ① Match Overlay
        mdf = getattr(self, '_match_df', None)
        match_active = (mdf is not None and not mdf.empty)
        self._set_btn_checked_silent(
            getattr(self, '_btn_match', None), match_active)
        # ② Conflict Search
        self._set_btn_checked_silent(
            getattr(self, '_btn_conflict', None),
            bool(getattr(self, '_conflict_coords', None)))
        # ③ IS Filter
        self._set_btn_checked_silent(
            getattr(self, '_btn_is_filter', None),
            bool(getattr(self, '_reserved_coords', None)))
        # ⑤ Adduct Ion Filter
        self._set_btn_checked_silent(
            getattr(self, '_btn_adduct_filter', None),
            bool(getattr(self, '_adduct_attribution', None)))
        # ④ RT Outlier Filter
        outliers = getattr(self, '_outlier_coords', {}) or {}
        self._set_btn_checked_silent(
            getattr(self, '_btn_rt_outlier', None),
            any(len(v) > 0 for v in outliers.values()))
        # ⑥ Coherence Filter
        coh_outliers = getattr(self, '_coherence_outliers', {}) or {}
        coh_active = (
            bool(getattr(self, '_coherence_models', None))
            or bool(getattr(self, '_coherence_assignments', None))
            or any(len(v) > 0 for v in coh_outliers.values())
        )
        self._set_btn_checked_silent(
            getattr(self, '_btn_coherence', None),
            coh_active)
        # ⑦ Select Quant Ion 完了状態を同期
        self._set_btn_checked_silent(
            getattr(self, '_btn_quant_ion', None),
            bool(getattr(self, '_quant_ion_choices', None)))

    def _toggle_adduct_filter(self):
        """⑤ Adduct Ion Filter ボタンのトグルハンドラ。

        ボタンが checked になった瞬間に apply、uncheck になった瞬間に undo を
        両モード一括で実行する。state 不整合を避けるため、ボタンの新しい
        checked 状態に応じて常に同じ方向の操作を全モードで行う。
        """
        # 前提チェック: match_df / IS Filter / Conflict Search が走っているか
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            self._btn_adduct_filter.blockSignals(True)
            self._btn_adduct_filter.setChecked(False)
            self._btn_adduct_filter.blockSignals(False)
            return

        target_state = bool(self._btn_adduct_filter.isChecked())
        if target_state:
            # 適用: 既存ロジックを両モードで
            self._apply_to_all_modes(self._apply_adduct_filter)
            # Filter Classes exempt 反映
            self._apply_filter_class_exemptions()
            self._draw_overlay()
        else:
            # 取り消し: 全モードで _undo_adduct_filter
            self._apply_to_all_modes(self._undo_adduct_filter)

    def _undo_adduct_filter(self):
        """⑤ Adduct Ion Filter の効果を取り消す。

        - self._adduct_attribution をクリア
        - match_df の adduct_filter_status を STATUS_NA(matched 行は KEPT
          のままだが、最終的に列リセットで揃える)に戻す
        - final_status を再計算(_update_final_status_df)
        - 後続の outlier / coherence 結果はリセット
        - フィルタ・散布図・テーブルを再描画

        IS Filter で adduct fingerprints をリセット済みの状態でも安全に
        呼べる(state がなければ no-op に近い)。
        """
        # state クリア
        self._adduct_attribution = {}

        if self._match_df is None:
            return

        match_df = self._match_df.copy()
        if 'adduct_filter_status' in match_df.columns:
            # matched 行は KEPT、unmatched は NA に揃える(初期状態と同等)
            matched_mask = match_df['matched']
            match_df.loc[matched_mask, 'adduct_filter_status'] = STATUS_KEPT
            match_df.loc[~matched_mask, 'adduct_filter_status'] = STATUS_NA
        # final_status 再計算
        match_df = _update_final_status_df(match_df)
        self._match_df = match_df

        # 後続フィルタ結果(Coherence)はリセット。
        # RT Outlier は新順序では上流(④)のため保持する。
        self._coherence_outliers = {}

        # UI 更新
        matched = match_df[match_df['matched']]
        classes = match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        log.info(f"[Adduct Ion Filter] reverted (mode={self.fe.ion_mode})")

        # 取り消し後はトグルボタンを uncheck に同期
        self._sync_adduct_button_state()

    def _open_adduct_patterns(self):
        """[Patterns…] ダイアログを開く。

        Adduct Ion Filter が参照する期待アダクトパターンを heatmap で
        表示する。データソースは内蔵 CURATED_ADDUCT_FINGERPRINTS および
        ② IS Filter で生成された self._runtime_is_fingerprints。
        """
        dlg = AdductPatternsDialog(
            curated=CURATED_ADDUCT_FINGERPRINTS,
            runtime=self._runtime_is_fingerprints,
            adduct_columns=ADDUCT_SET_V54,
            parent=self,
        )
        dlg.exec()

    def _collect_rejected_rows(self, status_col: str, status_val: str) -> list:
        """両モードの match_df から指定ステータスの行を抽出。
        戻り値: [{mode, lipid_class, compound, adduct, obs_rt, obs_mz, ppm}, ...]
        """
        out = []
        # active mode + 他モード snapshot
        contexts = [(self._active_mode, self._match_df)]
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        contexts.append((other_mode, snap.get('_match_df')))
        for mode, mdf in contexts:
            if mdf is None or mdf.empty:
                continue
            if status_col not in mdf.columns:
                continue
            sub = mdf[(mdf['matched']) & (mdf[status_col] == status_val)]
            for _, row in sub.iterrows():
                out.append({
                    'mode': mode,
                    'lipid_class': row.get('lipid_class', ''),
                    'compound': row.get('compound', ''),
                    'adduct': row.get('adduct', ''),
                    'obs_rt': float(row.get('obs_rt', 0.0)),
                    'obs_mz': float(row.get('obs_mz', 0.0)),
                    'ppm': float(row.get('delta_ppm', 0.0)),
                })
        out.sort(key=lambda r: (r['mode'], r['lipid_class'],
                                r['obs_rt'], r['obs_mz']))
        return out

    def _open_is_filter_report(self):
        """③ IS Filter で REJ_WINDOW された行を表示。
        AdvancedDialog から host 側に移動(背面散布図で
        ハイライトを見やすくするため Advanced を閉じてから呼び出す)。"""
        rows = self._collect_rejected_rows(
            'is_filter_status', STATUS_REJ_WINDOW)
        if not rows:
            QMessageBox.information(
                self, "No Rejected Spots",
                "No spots were rejected by ③ IS Filter.\n"
                "Either step ③ has not been run yet, or all candidates "
                "fell within the IS RT window.")
            return
        dlg = RejectedSpotsDialog(
            title="③ IS Filter — Rejected Spots",
            rows=rows,
            why="Outside IS RT window (beyond reference RT ± IS tolerance).",
            parent=self)
        dlg.spotSelected.connect(self._highlight_coherence_spot)
        # ダイアログ close 時にハイライトを消去
        dlg.finished.connect(self._clear_coherence_highlight)
        dlg.exec()

    def _open_rt_outlier_report(self):
        """④ RT Outlier Filter で外れ値判定された行を表示。"""
        rows = self._collect_rejected_rows(
            'rt_outlier_status', STATUS_REJ_OUTLIER)
        if not rows:
            QMessageBox.information(
                self, "No Outlier Spots",
                "No spots were flagged as outliers by ④ RT Outlier Filter.\n"
                "Step ④ may not have been executed yet.")
            return
        dlg = RejectedSpotsDialog(
            title="④ RT Outlier Filter — Rejected Spots",
            rows=rows,
            why="Deviates from class median RT by more than MAD × k.",
            parent=self)
        dlg.spotSelected.connect(self._highlight_coherence_spot)
        # ダイアログ close 時にハイライトを消去
        dlg.finished.connect(self._clear_coherence_highlight)
        dlg.exec()

    def _collect_all_events(self) -> list:
        """全 filtering ステップの結果を 1 リストに統合。

        各 event は以下のキーを持つ:
          mode, step, status, lipid_class, compound, adduct,
          obs_rt, obs_mz, detail, final

        ⑤ / ⑥ の行を **match_df の status 列から** 組み立てる。
        fix48 まではここだけが `_adduct_attribution` /
        `_coherence_assignments` という別の辞書を読んでいたため、
        Spot Status Inspector(status 列を読む)と食い違っていた
        (検証データで 292 件中 71 件)。辞書は「なぜそう判定したか」の
        詳細(予測 RT・残差・σ)を引くためだけに使う。

        あわせて ⑥ の判定 B(同一化合物が複数ピークに当たったときの
        残差最小選択 = rejected_residual)も出すようにした。fix48 まで
        Summary に一切現れず、「Summary に出ていないのに消えている」
        状態だった(pos 209 行中 134 行、neg 93 行中 75 行)。
        """
        events: list = []
        # 詳細引き当て用。key = (mode, rt_r6, mz_r6)
        _attr = self._build_merged_attribution(
            include_methods=('adduct', 'coherence', 'undecided'))

        def _attr_at(_mode, _rt, _mz):
            return _attr.get((_mode, round(float(_rt), 6),
                              round(float(_mz), 6)))

        # ── ③/④ Reject events from match_df (both modes) ──
        mode_contexts = [(self._active_mode, self._match_df)]
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        mode_contexts.append((other_mode, snap.get('_match_df')))
        for mode, mdf in mode_contexts:
            if mdf is None or mdf.empty:
                continue
            sub = mdf[mdf['matched']]
            # ③ IS Filter
            if 'is_filter_status' in sub.columns:
                rej = sub[sub['is_filter_status'] == STATUS_REJ_WINDOW]
                for _, row in rej.iterrows():
                    events.append({
                        'mode': mode,
                        'step': '③ IS Filter',
                        'status': 'Rejected',
                        'lipid_class': str(row.get('lipid_class', '')),
                        'compound': str(row.get('compound', '')),
                        'adduct': str(row.get('adduct', '')),
                        'obs_rt': float(row.get('obs_rt', 0.0)),
                        'obs_mz': float(row.get('obs_mz', 0.0)),
                        'detail': 'Outside IS RT window',
                        'final': str(row.get('final_status', '')),
                    })
            # ④ RT Outlier
            if 'rt_outlier_status' in sub.columns:
                rej = sub[sub['rt_outlier_status'] == STATUS_REJ_OUTLIER]
                for _, row in rej.iterrows():
                    events.append({
                        'mode': mode,
                        'step': '④ RT Outlier',
                        'status': 'Rejected',
                        'lipid_class': str(row.get('lipid_class', '')),
                        'compound': str(row.get('compound', '')),
                        'adduct': str(row.get('adduct', '')),
                        'obs_rt': float(row.get('obs_rt', 0.0)),
                        'obs_mz': float(row.get('obs_mz', 0.0)),
                        'detail': 'Deviates from class median RT by MAD×k',
                        'final': str(row.get('final_status', '')),
                    })

        # ── ⑤ Adduct + ⑥ Coherence (match_df の status 列から) ──
        # ここが辞書ではなく status 列を読むようになった。
        for mode, mdf in mode_contexts:
            if mdf is None or mdf.empty:
                continue
            sub = mdf[mdf['matched']]
            for col, rej, step, _mk in (
                    ('adduct_filter_status', STATUS_REJ_ADDUCT,
                     '⑤ Adduct Filter', 'adduct'),
                    ('coherence_status', STATUS_REJ_LOSER,
                     '⑥ Coherence Filter', 'coherence'),
                    ('coherence_status', STATUS_REJ_RESIDUAL,
                     '⑥ Coherence (residual)', 'residual')):
                if col not in sub.columns:
                    continue
                for _, row in sub[sub[col] == rej].iterrows():
                    _rt = float(row.get('obs_rt', 0.0))
                    _mz = float(row.get('obs_mz', 0.0))
                    a = _attr_at(mode, _rt, _mz)
                    if _mk == 'residual':
                        detail = ('another peak of the same compound fits '
                                  'the coherence model better')
                    elif a is not None:
                        detail = (f"lost to {a.get('winner_class', '')} "
                                  f"{a.get('winner_compound', '')}")
                        if _mk == 'coherence':
                            try:
                                detail += (f"  (σ={float(a.get('confidence', 0.0)):.2f}"
                                           f", {int(a.get('n_candidates', 0))} candidates)")
                            except Exception:
                                pass
                    else:
                        detail = 'rejected at a conflict spot'
                    events.append({
                        'mode': mode, 'step': step, 'status': 'Rejected',
                        'lipid_class': str(row.get('lipid_class', '')),
                        'compound': str(row.get('compound', '')),
                        'adduct': str(row.get('adduct', '')),
                        'obs_rt': _rt, 'obs_mz': _mz,
                        'detail': detail,
                        'final': str(row.get('final_status', '')),
                    })

        # ── 勝者 / △ Undecided は辞書から(座標の最終判定を併記) ──
        merged = _attr
        for key, a in merged.items():
            if isinstance(key, tuple) and len(key) == 3:
                mode = str(key[0]); rt = float(key[1]); mz = float(key[2])
            else:
                mode = str(a.get('mode', ''))
                rt = float(key[0]); mz = float(key[1])
            method = a.get('method', '')
            if method == 'adduct':
                step = '⑤ Adduct Filter'
                detail = (f"winner over "
                          f"{int(a.get('n_candidates', 2)) - 1} other "
                          f"candidate(s)")
            elif method == 'coherence':
                step = '⑥ Coherence Filter'
                _ls = a.get('losers') or []
                detail = (f"σ={a.get('confidence', 0.0):.2f}, winner over "
                          + ', '.join(str(l.get('compound'))
                                      for l in _ls[:4])
                          + (f" +{len(_ls) - 4}" if len(_ls) > 4 else ''))
            elif method == 'undecided':
                step = '△ Undecided'
                events.append({
                    'mode': mode, 'step': step, 'status': 'Undecided',
                    'lipid_class': str(a.get('winner_class', '')),
                    'compound': str(a.get('winner_compound', '')),
                    'adduct': str(a.get('winner_adduct', '')),
                    'obs_rt': rt, 'obs_mz': mz,
                    'detail': 'No method could decide the conflict',
                    'final': self._final_status_at(
                        mode, rt, mz, a.get('winner_compound')),
                })
                continue
            else:
                continue
            _fin = self._final_status_at(
                mode, rt, mz, a.get('winner_compound'))
            events.append({
                'mode': mode,
                'step': step,
                'status': ('Attributed' if _fin == STATUS_KEPT
                           else 'Attributed → rejected later'),
                'lipid_class': str(a.get('winner_class', '')),
                'compound': str(a.get('winner_compound', '')),
                'adduct': str(a.get('winner_adduct', '')),
                'obs_rt': rt,
                'obs_mz': mz,
                'detail': detail,
                'final': _fin,
            })

        # ── ✋ Manual attribution events (両モード) ──
        # 手動帰属したスポットも Summary に表示する。
        manual_contexts = [(
            self._active_mode, self._match_df,
            getattr(self, '_manual_winner_coords', None) or set(),
        )]
        snap2 = self._mode_state.get(other_mode, {}) or {}
        manual_contexts.append((
            other_mode, snap2.get('_match_df'),
            snap2.get('_manual_winner_coords') or set(),
        ))
        for mode, mdf, mcs in manual_contexts:
            if mdf is None or mdf.empty or not mcs:
                continue
            if 'manual_status' not in mdf.columns:
                continue
            for coord in mcs:
                rt = float(coord[0])
                mz = float(coord[1])
                same_coord = (
                    (mdf['obs_rt'].round(6) == round(rt, 6))
                    & (mdf['obs_mz'].round(6) == round(mz, 6))
                    & mdf['matched']
                )
                winner_rows = mdf[
                    same_coord
                    & (mdf['manual_status'] == STATUS_KEPT)
                ]
                if winner_rows.empty:
                    continue
                wrow = winner_rows.iloc[0]
                loser_rows = mdf[
                    same_coord
                    & (mdf['manual_status'] == STATUS_REJ_MANUAL)
                ]
                if not loser_rows.empty:
                    loser_cls = ' / '.join(
                        sorted(loser_rows['lipid_class'].unique().tolist()))
                    detail = f"manually attributed (vs {loser_cls})"
                else:
                    detail = "manually attributed"
                events.append({
                    'mode': mode,
                    'step': '✋ Manual',
                    'status': 'Attributed',
                    'lipid_class': str(wrow.get('lipid_class', '')),
                    'compound': str(wrow.get('compound', '')),
                    'adduct': str(wrow.get('adduct', '')),
                    'obs_rt': rt,
                    'obs_mz': mz,
                    'detail': detail,
                    'final': str(wrow.get('final_status', '')),
                })

        events.sort(key=lambda e: (e['mode'], e['step'],
                                   e['obs_rt'], e['obs_mz']))
        return events

    def _final_status_at(self, mode: str, rt: float, mz: float,
                         compound=None) -> str:
        """指定座標(と化合物)の final_status を match_df から引く。

        Summary の行に「この帰属は最終的にどうなったか」を併記するため。
        winner として記録されていても、後段の ⑥ 判定 B(residual)や
        手動操作で rejected になっていることがある。
        """
        if mode == self._active_mode:
            mdf = self._match_df
        else:
            mdf = (self._mode_state.get(mode, {}) or {}).get('_match_df')
        if mdf is None or getattr(mdf, 'empty', True):
            return ''
        try:
            m = (mdf['matched']
                 & (mdf['obs_rt'].round(6) == round(float(rt), 6))
                 & (mdf['obs_mz'].round(6) == round(float(mz), 6)))
            if compound is not None:
                m = m & (mdf['compound'].astype(str) == str(compound))
            s = mdf[m]
            if s.empty:
                return ''
            return str(s.iloc[0].get('final_status', ''))
        except Exception:
            return ''

    def _open_full_report(self):
        """全 filter ステップの結果を統合表示。"""
        events = self._collect_all_events()
        if not events:
            QMessageBox.information(
                self, "No Results",
                "No filtering results to display.\n"
                "Run ③ IS Filter / ④ RT Outlier / ⑤ Adduct / ⑥ Coherence "
                "first.")
            return
        dlg = FullReportDialog(events=events, parent=self)
        dlg.spotSelected.connect(self._highlight_coherence_spot)
        # ダイアログ close 時にハイライトを消去
        dlg.finished.connect(self._clear_coherence_highlight)
        dlg.exec()

    def _open_adduct_report(self):
        """[Report…] ダイアログを開く。

        ⑤ Adduct Ion Filter で計算された衝突スポットの観測アダクトパターンを
        heatmap で表示し、PNG / xlsx エクスポートを提供する。

        ⑤ 決着分(+ 未決 △)の attribution table も
        表示するよう拡張。
        """
        if not self._adduct_attribution:
            QMessageBox.information(
                self, "No Adduct Filter results",
                "Please run 'Apply' for ⑤ Adduct Ion Filter first.\n"
                "The report shows observation results of that step.")
            return
        # ⑤ 決着分 + undecided の attribution dict を構築
        adduct_attribution_table = self._build_merged_attribution(
            include_methods=('adduct', 'undecided'))
        # heatmap 用に両モードの adduct_attribution を統合
        # (キーを (mode, rt, mz) の 3-tuple に正規化)
        merged_heatmap_attr: dict = {}
        for coord, r in (self._adduct_attribution or {}).items():
            merged_heatmap_attr[(self._active_mode, coord[0], coord[1])] = r
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        other_attr = snap.get('_adduct_attribution') or {}
        for coord, r in other_attr.items():
            merged_heatmap_attr[(other_mode, coord[0], coord[1])] = r

        # heatmap と attribution table に共通 ID を付与
        # ソート優先順位: (skipped, undecided, mode, rt, mz)
        # → 決着済みを先頭、未決を末尾に配置(heatmap の natural display 順と一致)
        all_keys: set = set(merged_heatmap_attr.keys())
        all_keys |= set(adduct_attribution_table.keys())

        def _key_sort(k):
            res_h = merged_heatmap_attr.get(k) or {}
            res_t = adduct_attribution_table.get(k) or {}
            # heatmap 由来の skipped / winner を優先(⑤ の判定状態)
            skipped = bool(res_h.get('skipped', False))
            # winner が None かつ method='undecided' なら未決
            undecided = (res_t.get('method') == 'undecided'
                         or res_h.get('winner') is None)
            if isinstance(k, tuple) and len(k) == 3:
                mode_v = str(k[0]); rt_v = float(k[1]); mz_v = float(k[2])
            else:
                mode_v = ''; rt_v = float(k[0]); mz_v = float(k[1])
            return (skipped, undecided, mode_v, rt_v, mz_v)
        sorted_keys = sorted(all_keys, key=_key_sort)
        key_to_id = {k: (i + 1) for i, k in enumerate(sorted_keys)}

        dlg = AdductReportDialog(
            attribution=merged_heatmap_attr,
            adduct_columns=ADDUCT_SET_V54,
            attribution_table=adduct_attribution_table,
            sigma_threshold=self._spin_sigma.value(),
            key_to_id=key_to_id,
            parent=self,
        )
        # 行クリックで散布図ハイライト
        try:
            dlg.spotSelected.connect(self._highlight_coherence_spot)
        except Exception:
            pass
        # ダイアログ close 時にハイライトを消去
        try:
            dlg.finished.connect(self._clear_coherence_highlight)
        except Exception:
            pass
        dlg.exec()

    # ── Run All ───────────────────────────────────────────────
    def _run_all(self):
        """① 〜 ⑥ を順に走らせる(やり直し型、毎回最初から)。

        両モードがロードされている場合、各 step を全モードで
        実行する(現モード → 他モード → 元モードに戻す)。ユーザーは「Run All」
        1 回で pos/neg 両方の解析を完了できる。

        二段階実行にした。
          Phase 1 … 両モードで ①〜③（Match / Conflict / IS Filter）
          Phase 2 … 両モードで ④〜⑥（RT Outlier / Adduct / Coherence）
        従来は「モードごとに ①〜⑥」だったので、先に走るモードの ⑤ は
        相手モードの _class_ref_rt がまだ空で、pos↔neg RT シフトを推定
        できずクロスモード窓が既定値 ±0.05 min に落ちていた。Phase 1 で
        両モードの基準 RT を揃えてから ⑤ に入れば、どちらのモードでも
        IS 由来の窓（検証データでは ±0.015 min）が使える。

        Run All は何回押しても同じ結果になる:
          - ① Match Overlay が ②〜⑥ の state をカスケードクリアするので、
            トグル式の ④ ⑥ も常に「適用」側に入る
          - ③ は保存済み IS 選択があればプレビューダイアログを出さない
          - ③ が適用されなかった場合は ④〜⑥ に進まず中止する
        """
        if not self._lib_path:
            QMessageBox.warning(
                self, "No Library",
                "Please load a library first (Load library… button).")
            return

        steps = [
            # Conflict Search を IS Filter の前に移動
            ("① Match Overlay",      self._run_match),
            ("② Conflict Search",    self._conflict_search),
            ("③ IS Filter",          self._apply_is_filter),
            ("④ RT Outlier Filter",  self._rt_outlier_search),
            ("⑤ Adduct Ion Filter",  self._apply_adduct_filter),
            ("⑥ Coherence Filter",   self._coherence_search),
        ]
        # 実行対象モード: all_entries に含まれる全 ion_mode
        loaded_modes = set()
        for fe2 in (self._all_entries or []):
            m = getattr(fe2, 'ion_mode', None)
            if m:
                loaded_modes.add(m)
        # pos → neg 固定順序(上プロット → 下プロット)。
        # ファイルの登録順に依存しない。
        modes_to_run = [m for m in ('pos', 'neg') if m in loaded_modes]
        if not modes_to_run:
            modes_to_run = [self._active_mode]
        # アクティブモードは最後に戻ってくる
        original_mode = self._active_mode

        # Run All 中も中間描画を抑制
        self._suppress_overlay_draw = True
        self._run_all_in_progress = True
        try:
            self._btn_run_all.setEnabled(False)
            QApplication.processEvents()

            # 二段階実行。Phase 1 で両モードの _class_ref_rt を
            # 揃えてから Phase 2 に入る。こうしないと先に走るモードの ⑤ が
            # pos↔neg シフトを推定できない。
            # Phase 2 の先頭で指紋を作り直す。Phase 1 で先に走った
            # モードの指紋は相手モードの class_ref_rt がまだ無く、
            # pos↔neg シフトが既定値で作られているため。
            _fp_step = ("③' IS 指紋の再構築", self._rebuild_is_fingerprints_now)
            phases = (("Phase 1 (①〜③)", tuple(steps[:3])),
                      ("Phase 2 (④〜⑥)", (_fp_step,) + tuple(steps[3:])))
            _abort = False
            for phase_label, phase_steps in phases:
                if _abort:
                    break
                log.info(f"[Run All] ===== {phase_label} =====")
                for run_mode in modes_to_run:
                    if run_mode != self._active_mode:
                        log.info(f"[Run All] switching to {run_mode}")
                        self._snapshot_mode_state(self._active_mode)
                        if self._load_mode_state(run_mode):
                            self._active_mode = run_mode
                            self._refresh_axis_aliases()  # 修正
                            # ラジオも同期
                            try:
                                for rb in (self._mode_radio_pos,
                                           self._mode_radio_neg):
                                    rb.blockSignals(True)
                                if run_mode == 'pos':
                                    self._mode_radio_pos.setChecked(True)
                                else:
                                    self._mode_radio_neg.setChecked(True)
                            except Exception:
                                pass
                            finally:
                                try:
                                    for rb in (self._mode_radio_pos,
                                               self._mode_radio_neg):
                                        rb.blockSignals(False)
                                except Exception:
                                    pass
                        else:
                            log.info(f"[Run All] skip mode {run_mode}: "
                                  "data not loaded")
                            continue

                    log.info(f"[Run All] === Mode {run_mode.upper()} "
                          f"/ {phase_label} ===")
                    for label, fn in phase_steps:
                        log.info(f"[Run All] {label} …")
                        try:
                            fn()
                        except Exception as exc:
                            log.warning(f"[Run All] {label} raised: {exc}")
                            QMessageBox.warning(
                                self, f"{label} failed",
                                f"Step '{label}' failed in {run_mode} mode:"
                                f"\n{exc}\n"
                                "Run All will continue with the next step.")
                        QApplication.processEvents()

                        # ③ が適用されなかったら ④〜⑥ に進まない。
                        #   ① が _class_ref_rt を消したあとで ③ が早期
                        #   return すると、④⑤⑥ が基準 RT なしで走って
                        #   黙って壊れる(検証データでは ④ の外れ値が
                        #   29 → 236、⑥ の帰属が 38 → 115 になる)。
                        # バインドメソッドは毎回新しいオブジェクトに
                        # なるので `is` では比較できない。ラベルで判定する。
                        if (label.startswith("③")
                                and not getattr(self,
                                                '_is_filter_applied_ok',
                                                False)):
                            log.info("[Run All] ③ IS Filter が適用されなかった"
                                  f"ため中止 (mode={run_mode})")
                            QMessageBox.warning(
                                self, "Run All aborted",
                                "\u2462 IS Filter was not applied in "
                                f"{run_mode} mode.\n\n"
                                "Steps \u2463-\u2465 rely on the per-class "
                                "reference RT that \u2462 establishes. "
                                "Running them without it silently produces "
                                "wrong filtering, so Run All stopped here.\n\n"
                                "Run \u2462 IS Filter manually, check the "
                                "Candidate Picker, then press Run All again.")
                            _abort = True
                            break
                    if _abort:
                        break

            # 元のモードに戻す
            if self._active_mode != original_mode:
                log.info(f"[Run All] returning to original mode {original_mode}")
                self._snapshot_mode_state(self._active_mode)
                if self._load_mode_state(original_mode):
                    self._active_mode = original_mode
                    self._refresh_axis_aliases()  # 修正
                    # blockSignals 必須
                    try:
                        for rb in (self._mode_radio_pos,
                                   self._mode_radio_neg):
                            rb.blockSignals(True)
                        if original_mode == 'pos':
                            self._mode_radio_pos.setChecked(True)
                        else:
                            self._mode_radio_neg.setChecked(True)
                    except Exception:
                        pass
                    finally:
                        try:
                            for rb in (self._mode_radio_pos,
                                       self._mode_radio_neg):
                                rb.blockSignals(False)
                        except Exception:
                            pass
                    classes = []
                    if self._match_df is not None:
                        classes = self._match_df['lipid_class'].unique().tolist()
                    self._rebuild_filters(classes)
                    self._draw_overlay()
                    if self._match_df is not None:
                        matched = self._match_df[self._match_df['matched']]
                        self._rebuild_match_table(matched)

            log.info("[Run All] All steps completed for all modes.")
        finally:
            self._run_all_in_progress = False
            self._btn_run_all.setEnabled(True)
            # 最終 UI 更新(suppress 解除して 1 回だけ描画)
            self._suppress_overlay_draw = False
            try:
                classes = []
                if self._match_df is not None:
                    classes = self._match_df['lipid_class'].unique().tolist()
                self._rebuild_filters(classes)
                # 両モード(active + 他方)を描画。
                # _draw_overlay だけだと _mode_state に保存された他方モード
                # の match_df が _draw_other_mode_overlay に渡らず、neg
                # プロットが空になる問題を修正。
                self._draw_both_overlays()
                if self._match_df is not None:
                    matched = self._match_df[self._match_df['matched']]
                    self._rebuild_match_table(matched)
            except Exception as e:
                log.warning(f"[Run All] final UI update failed: {e}")

    def _loose_coherence_enabled(self) -> bool:
        """緩い帰属が有効か。UI が無い経路でも落ちないようにする。"""
        try:
            return bool(self._cb_loose_coherence.isChecked())
        except Exception:
            return bool(COH_LOOSE_ATTRIBUTION_DEFAULT)



    def _two_series_classes(self) -> tuple:
        """二峰性として扱うクラス。UI が無い経路でも落ちない。"""
        try:
            _s = str(self._ed_two_series.text() or '')
        except Exception:
            return tuple(COH_TWO_SERIES_CLASSES_DEFAULT)
        _out = tuple(_t.strip() for _t in _s.replace(';', ',').split(',')
                     if _t.strip())
        return _out

    def _sum_two_series_enabled(self) -> bool:
        """定量エクスポートで 2 本を合算するか。"""
        try:
            return bool(self._cb_sum_two_series.isChecked())
        except Exception:
            return bool(QUANT_SUM_TWO_SERIES_DEFAULT)

    def _rt_outlier_search(self):
        """Step 4: RT外れ値検出のトグル（ON/OFF切替）。×マーカー更新

        外れ値判定結果は self._outlier_coords(描画用)と
        self._match_df['rt_outlier_status'] 列(エクスポート判定用)の
        両方に反映する。final_status は AND で再計算される。
        """
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            return

        # トグル: 既に検出済みならクリア、未検出なら実行
        has_outliers = any(len(v) > 0 for v in self._outlier_coords.values())
        if has_outliers:
            self._outlier_coords = {}
            # ④ undo 時に⑥ Coherence の outliers もクリア
            # (Coherence outliers は RT outliers を学習データから除く前提)
            self._coherence_outliers = {}
        else:
            # Advanced Setting → Filter Classes → ③ で
            # チェックを外したクラスは「③ が効いていない」ものとして
            # 扱い、IS があっても ④ の対象にする。
            _ex_is = (getattr(self, '_filter_class_exemptions', None)
                      or {}).get('is', set())
            self._outlier_coords = calc_rt_outliers(
                self._match_df, self._spin_iqr.value(),
                is_exempt_classes=_ex_is,
                class_ref_rt=dict(getattr(self, '_class_ref_rt', None) or {}))

        # rt_outlier_status 列を更新
        self._match_df = _apply_rt_outlier_status(
            self._match_df, self._outlier_coords)

        # Filter Classes exempt 反映
        self._apply_filter_class_exemptions()

        # 右パネルのチェックボックスを含めて再構築し、散布図とテーブルを更新
        matched = self._match_df[self._match_df['matched']]
        classes = self._match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        # ④ ボタン状態を同期(⑥ も連動)
        self._sync_step_button_states()

    def _coherence_search(self):
        """Step ⑥ Coherence Engine のトグル(ON/OFF 切替)。

        トグル時の処理:
          - ON 済み(モデル or 帰属 or 外れ値のいずれかが存在)→ 全 Coherence
            状態をクリア。衝突表示は △(未解決)に戻る。
          - OFF → 以下の 3 段階を実行:
              1. link × クラス別の pooled 回帰 RT=α·C+β·U+γ_class を fit
                 (学習データ: 非衝突 + 非 RT 外れ値のマッチ点、IS は anchor)
              2. 有効な衝突ペア(ペア選択ダイアログの ON 分)に対して
                 残差最小のクラスへ自動帰属。信頼度 σ を併記。
              3. 非衝突 + 帰属後の元衝突点について、自クラストレンドから
                 MAD×k を超える残差の点を Coherence 外れ値としてフラグ。
        """
        if self._match_df is None:
            QMessageBox.warning(self, "Run Match First",
                                "Please run 'Match & Overlay' first.")
            return

        # トグル判定: いずれかの Coherence 状態が存在すれば ON とみなす
        has_any = (bool(self._coherence_models)
                   or bool(self._coherence_assignments)
                   or any(len(v) > 0 for v in self._coherence_outliers.values()))

        if has_any:
            # クリア
            self._coherence_models = {}
            self._coherence_assignments = {}
            # status 列もクリア
            if self._match_df is not None:
                self._match_df['coherence_status'] = np.where(
                    self._match_df.get('matched', False) == True,
                    STATUS_KEPT, STATUS_NA)
                self._match_df = _update_final_status_df(self._match_df)
            self._coherence_outliers = {}
            self._btn_coherence_pairs.setEnabled(False)
            self._btn_coherence_report.setEnabled(False)
        else:
            # ライブラリから衝突ペアを自動検出(未検出の場合)
            if not self._coherence_all_pairs:
                try:
                    lib_df = load_lipid_library(
                        self._lib_path, self.fe.ion_mode)
                except Exception as e:
                    QMessageBox.warning(
                        self, "Library Error",
                        f"Failed to reload library for pair detection:\n{e}")
                    return
                self._coherence_all_pairs = _detect_collision_pairs(
                    lib_df, self._spin_ppm.value())
                # 初回はデフォルト全 ON
                self._coherence_pairs = list(self._coherence_all_pairs)

            enabled = {frozenset(p) for p in self._coherence_pairs}

            # ⑤ Adduct Ion Filter で勝者確定済みの衝突座標は
            # 学習データから除外しない(=非衝突点扱い)。⑤ winner 行(final_status=
            # 'kept')が pooled regression に含まれることで、PC/PE 等の γ 推定が
            # ⑤ の確定情報を活用してより精密になる。⑤ で undecided/skipped の
            # 衝突のみを ⑥ の学習除外対象とする。
            resolved_by_adduct = {
                coord for coord, result in (self._adduct_attribution or {}).items()
                if not result.get('skipped')
                and result.get('winner') is not None
            }
            undecided_conflict_coords = (
                self._conflict_coords - resolved_by_adduct)
            n_resolved = len(resolved_by_adduct)
            n_undecided = len(undecided_conflict_coords)
            log.info(f"[Coherence Filter] training: ⑤ winner spots "
                  f"included = {n_resolved}, undecided conflicts "
                  f"excluded = {n_undecided}")

            # 1. pooled fit(非衝突・非 RT 外れ値で学習。⑤ winner は含む)
            self._coherence_models = _fit_coherence_model(
                self._match_df,
                conflict_coords=undecided_conflict_coords,
                rt_outlier_coords=self._outlier_coords,
                min_points=5,
                mad_k=self._spin_iqr.value(),
                two_series_classes=self._two_series_classes(),
            )

            # 2. 衝突帰属: ⑤ で未決の conflict のみ coherence で再帰属
            #    ⑤ で決着済みの spot は _adduct_attribution に既に勝者が
            #    記録されているため、⑥ で重複処理しない。
            self._coherence_assignments = resolve_conflicts_by_coherence(
                undecided_conflict_coords,
                self._conflict_map,
                self._coherence_models,
                enabled,
                sigma_threshold=self._spin_sigma.value(),
                loose=self._loose_coherence_enabled(),
            )

            # σ の内訳を要約して出す。外挿している側が勝った /
            # 負けたスポットがどれだけあるかが分かるようにする。
            try:
                _as = list(self._coherence_assignments.values())
                if _as:
                    _low = sum(1 for a in _as if a.get('low_confidence'))
                    _ext = sum(
                        1 for a in _as
                        if max(float(a.get('leverage_winner') or 0.0),
                               float(a.get('leverage_loser') or 0.0)) > 3.0)
                    log.info(f"[Coherence Filter] attributions {len(_as)}: "
                          f"low-confidence {_low}, "
                          f"extrapolated candidate (h>3) {_ext}")
            except Exception as _e:
                log.warning(f"[Coherence Filter] sigma summary failed: {_e}")

            # 3. Coherence 外れ値(帰属適用後の match_df で判定)
            resolved_df = self._build_resolved_match_df()
            self._coherence_outliers = calc_coherence_outliers(
                resolved_df, self._coherence_models,
                mad_k=self._spin_iqr.value())

            # status 列に Coherence 結果を反映
            self._match_df = _apply_coherence_status(
                self._match_df,
                self._coherence_models,
                self._coherence_assignments)

            # Filter Classes exempt 反映
            self._apply_filter_class_exemptions()

            self._btn_coherence_pairs.setEnabled(True)
            self._btn_coherence_report.setEnabled(True)

        # 表示を更新
        matched = self._match_df[self._match_df['matched']]
        classes = self._match_df['lipid_class'].unique().tolist()
        self._rebuild_filters(classes)
        self._draw_overlay()
        self._rebuild_match_table(matched)
        # ⑥ ボタン状態を同期
        self._sync_step_button_states()



    def _build_resolved_match_df(self) -> pd.DataFrame:
        """帰属適用後の match_df を生成する(物理的にはコピー)。

        帰属確定した衝突点については、敗者クラスの行を除外した
        新しい DataFrame を返す。元の match_df は変更しない。
        """
        if self._match_df is None:
            return pd.DataFrame()
        df = self._match_df.copy()
        if not self._coherence_assignments:
            return df
        # 敗者クラス行を特定して除外
        drop_mask = np.zeros(len(df), dtype=bool)
        coords_arr = list(zip(df['obs_rt'].round(6), df['obs_mz'].round(6)))
        for i, coord in enumerate(coords_arr):
            a = self._coherence_assignments.get(coord)
            if a is None:
                continue
            if df.iloc[i]['lipid_class'] == a['loser_class']:
                drop_mask[i] = True
        return df[~drop_mask].reset_index(drop=True)

    def _open_coherence_pairs(self):
        """Coherence Pairs 選択ダイアログを開く。
        ペアの ON/OFF を変更したら Coherence を再実行して反映。"""
        if not self._coherence_all_pairs:
            QMessageBox.information(
                self, "No Pairs",
                "No colliding class pairs were detected in the library.")
            return
        dlg = CoherencePairsDialog(
            self._coherence_all_pairs, self._coherence_pairs, parent=self)
        if dlg.exec():
            self._coherence_pairs = dlg.selected_pairs()
            # 現在 Coherence が ON なら再実行(クリア→再実行で反映)
            has_any = (bool(self._coherence_models)
                       or bool(self._coherence_assignments)
                       or any(len(v) > 0 for v in self._coherence_outliers.values()))
            if has_any:
                # 一旦クリアして再実行(トグル 2 回)
                self._coherence_search()
                self._coherence_search()

    def _open_coherence_report(self):
        """Coherence Attribution Report ダイアログを開く。

        ⑥ Coherence で帰属したスポットのみを表示する。
        ⑤ で決着したスポットは ⑤ Report に分離。
        """
        # ⑥ のみ(両モード統合)
        merged = self._build_merged_attribution(include_methods=('coherence',))
        if not merged:
            QMessageBox.information(
                self, "No Attributions",
                "No coherence attributions yet. "
                "Run ⑥ Coherence Filter first.")
            return
        dlg = CoherenceReportDialog(
            merged,
            self._spin_sigma.value(),
            parent=self)
        # 選択されたスポットへジャンプするシグナル
        dlg.spotSelected.connect(self._highlight_coherence_spot)
        # ダイアログ close 時にハイライトを消去
        dlg.finished.connect(self._clear_coherence_highlight)
        dlg.exec()

    def _build_merged_attribution(self,
                                  include_methods: tuple = (
                                      'coherence', 'adduct', 'undecided')
                                  ) -> dict:
        """帰属辞書を構築。
        include_methods で表示対象の method を絞れる。
          - ⑥ Coherence Report: ('coherence',)
          - ⑤ Adduct Report: ('adduct', 'undecided')

        各エントリは 'method' キーで由来を示す:
          'adduct'    ─ ⑤ Adduct Ion Filter で決着
          'coherence' ─ ⑥ Coherence Filter で帰属
          'undecided' ─ どちらでも決着できなかった △ 状態

        返り値の各エントリは以下のキーを持つ:
          'method', 'mode', 'confidence', 'winner_class', 'winner_compound',
          'winner_adduct', 'loser_class', 'loser_compound', 'loser_adduct',
          'pred_winner', 'pred_loser', 'res_winner', 'res_loser',
          'low_confidence', and detail fields (winner_C, winner_alpha, ...)

        Key は (mode, rt, mz) の 3-tuple で、両モードのエントリを共存させる。
        """
        merged: dict = {}

        # 処理対象のモード一覧を構築
        mode_contexts: list[dict] = []
        # active モード
        mode_contexts.append({
            'mode': self._active_mode,
            'adduct_attribution': self._adduct_attribution or {},
            'coherence_assignments': self._coherence_assignments or {},
            'conflict_coords': self._conflict_coords or set(),
            'conflict_map': self._conflict_map or {},
            'coherence_models': self._coherence_models or {},
        })
        # 他モード(snapshot から)
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        if snap:
            mode_contexts.append({
                'mode': other_mode,
                'adduct_attribution': snap.get('_adduct_attribution') or {},
                'coherence_assignments': snap.get('_coherence_assignments') or {},
                'conflict_coords': snap.get('_conflict_coords') or set(),
                'conflict_map': snap.get('_conflict_map') or {},
                'coherence_models': snap.get('_coherence_models') or {},
            })

        for ctx in mode_contexts:
            self._populate_merged_for_mode(merged, ctx, include_methods)
        return merged

    def _populate_merged_for_mode(
            self, merged: dict, ctx: dict,
            include_methods: tuple = (
                'coherence', 'adduct', 'undecided')) -> None:
        """1 モード分の attribution を merged dict に追加する
        (key は (mode, rt, mz) の 3-tuple)。"""
        mode = ctx['mode']
        adduct_attribution = ctx['adduct_attribution']
        coherence_assignments = ctx['coherence_assignments']
        conflict_coords = ctx['conflict_coords']
        conflict_map = ctx['conflict_map']
        coh_models = ctx['coherence_models']

        # ⑥ Coherence の帰属はそのまま採用(method='coherence')
        if 'coherence' in include_methods:
            for coord, a in (coherence_assignments or {}).items():
                key = (mode, coord[0], coord[1])
                entry = dict(a)
                entry['method'] = 'coherence'
                entry['mode'] = mode
                merged[key] = entry

        # ⑤ Adduct Ion Filter の勝者 → Coherence-format に変換
        if 'adduct' not in include_methods:
            adduct_attribution = {}
        for coord, result in (adduct_attribution or {}).items():
            key = (mode, coord[0], coord[1])
            if key in merged:
                continue  # ⑥ が再帰属したものは ⑥ を優先
            if result.get('skipped'):
                continue
            winner_cls = result.get('winner')
            if winner_cls is None:
                continue  # undecided はスキップ
            losers = result.get('losers') or set()
            scores = result.get('scores') or {}
            # 勝者・敗者の代表 compound / adduct を conflict_map から取得
            entries = conflict_map.get(coord, [])
            winner_compound = ""
            winner_adduct = ""
            loser_compound = ""
            loser_adduct = ""
            loser_cls = next(iter(losers), "") if losers else ""
            for e in entries:
                if e.get('lipid_class') == winner_cls and not winner_compound:
                    winner_compound = e.get('compound', '')
                    winner_adduct = e.get('adduct', '')
                if loser_cls and e.get('lipid_class') == loser_cls \
                        and not loser_compound:
                    loser_compound = e.get('compound', '')
                    loser_adduct = e.get('adduct', '')
            # confidence: スコア差を proxy として採用
            top = scores.get(winner_cls, 0.0)
            others = [s for c, s in scores.items() if c != winner_cls]
            second = max(others) if others else 0.0
            margin = top - second
            # ⑤ 決着分にも coherence model 由来の pred/res
            # モデル詳細(α,β,γ,rmse,C,U)も保存
            entry = {
                'method':          'adduct',
                'adduct_method':   str(result.get('method', '') or ''),
                'mode':            mode,
                # エクスポートの competitors 用に負けた候補を全部持つ
                'losers': [
                    {'class': e.get('lipid_class', ''),
                     'compound': e.get('compound', ''),
                     'adduct': e.get('adduct', '')}
                    for e in entries
                    if e.get('lipid_class') != winner_cls
                ],
                'confidence':      float(margin),
                'winner_class':    winner_cls,
                'winner_compound': winner_compound,
                'winner_adduct':   winner_adduct,
                'loser_class':     loser_cls,
                'loser_compound':  loser_compound,
                'loser_adduct':    loser_adduct,
                'pred_winner':     float('nan'),
                'pred_loser':      float('nan'),
                'res_winner':      float('nan'),
                'res_loser':       float('nan'),
                'low_confidence':  False,
                'adduct_scores':   dict(scores),
            }
            try:
                Cw, Uw = _parse_c_and_u(winner_compound)
                if Cw is not None:
                    link_w = _link_kind(winner_compound)
                    Mw = coh_models.get((winner_cls, link_w))
                    if Mw is not None:
                        # 近い系列の予測を使う
                        pred_w = float(_coherence_best_pred(
                            Mw, Cw, Uw, coord[0])[0])
                        entry.update({
                            'pred_winner': pred_w,
                            'res_winner':  float(coord[0]) - pred_w,
                            'winner_C':    int(Cw), 'winner_U':    int(Uw),
                            'winner_link': link_w,
                            'winner_alpha':float(Mw['alpha']),
                            'winner_beta': float(Mw['beta']),
                            'winner_gamma':float(Mw['gamma']),
                            'winner_rmse': float(Mw.get('rmse', 0.0)),
                            'winner_fit_mode': Mw.get('fit_mode', 'pooled'),
                        })
                if loser_cls and loser_compound:
                    Cl, Ul = _parse_c_and_u(loser_compound)
                    if Cl is not None:
                        link_l = _link_kind(loser_compound)
                        Ml = coh_models.get((loser_cls, link_l))
                        if Ml is not None:
                            # 近い系列の予測を使う
                            pred_l = float(_coherence_best_pred(
                                Ml, Cl, Ul, coord[0])[0])
                            entry.update({
                                'pred_loser': pred_l,
                                'res_loser':  float(coord[0]) - pred_l,
                                'loser_C':    int(Cl), 'loser_U':    int(Ul),
                                'loser_link': link_l,
                                'loser_alpha':float(Ml['alpha']),
                                'loser_beta': float(Ml['beta']),
                                'loser_gamma':float(Ml['gamma']),
                                'loser_rmse': float(Ml.get('rmse', 0.0)),
                                'loser_fit_mode': Ml.get('fit_mode', 'pooled'),
                            })
            except Exception:
                pass
            merged[key] = entry

        # ⑤ も ⑥ も決着できなかった衝突(△ のまま)
        if 'undecided' not in include_methods:
            conflict_coords = set()
        for coord in (conflict_coords or set()):
            key = (mode, coord[0], coord[1])
            if key in merged:
                continue
            candidates = conflict_map.get(coord, [])
            if not candidates:
                continue
            # 各候補について coherence model 予測を計算
            scored = []
            for cand in candidates:
                cls = cand.get('lipid_class')
                compound = cand.get('compound', '')
                C, U = _parse_c_and_u(compound)
                if C is None:
                    continue
                link = _link_kind(compound)
                model = coh_models.get((cls, link))
                if model is None:
                    continue
                # 近い系列の予測を使う
                pred = float(_coherence_best_pred(
                    model, C, U, coord[0])[0])
                res = float(coord[0]) - pred
                scored.append({
                    'class': cls, 'compound': compound,
                    'adduct': cand.get('adduct', ''),
                    'pred': pred, 'res': res, 'abs_res': abs(res),
                    'rmse': float(model.get('rmse', 0.0)),
                })
            scored.sort(key=lambda x: x['abs_res'])
            if len(scored) >= 2:
                w = scored[0]; l = scored[1]
                rmse_ref = (w['rmse'] + l['rmse']) / 2.0
                margin = l['abs_res'] - w['abs_res']
                conf = (margin / rmse_ref) if rmse_ref > 0 else 0.0
                # 詳細値も埋める
                undec_entry = {
                    'method':          'undecided',
                    'mode':            mode,
                    'confidence':      float(conf),
                    'winner_class':    w['class'],
                    'winner_compound': w['compound'],
                    'winner_adduct':   w['adduct'],
                    'loser_class':     l['class'],
                    'loser_compound':  l['compound'],
                    'loser_adduct':    l['adduct'],
                    'pred_winner':     float(w['pred']),
                    'pred_loser':      float(l['pred']),
                    'res_winner':      float(w['res']),
                    'res_loser':       float(l['res']),
                    'low_confidence':  True,
                    'rmse_ref':        float(rmse_ref),
                }
                for side, cand in [('winner', w), ('loser', l)]:
                    Cx, Ux = _parse_c_and_u(cand['compound'])
                    if Cx is None: continue
                    link_x = _link_kind(cand['compound'])
                    Mx = coh_models.get((cand['class'], link_x))
                    if Mx is None: continue
                    undec_entry.update({
                        f'{side}_C':    int(Cx), f'{side}_U':    int(Ux),
                        f'{side}_link': link_x,
                        f'{side}_alpha':float(Mx['alpha']),
                        f'{side}_beta': float(Mx['beta']),
                        f'{side}_gamma':float(Mx['gamma']),
                        f'{side}_rmse': float(Mx.get('rmse', 0.0)),
                        f'{side}_fit_mode': Mx.get('fit_mode', 'pooled'),
                    })
                merged[key] = undec_entry
            else:
                # model 不足 — 候補クラスのみ列挙
                cls_set = sorted(set(c.get('lipid_class', '') for c in candidates))
                cmp_set = sorted(set(c.get('compound', '') for c in candidates))
                merged[key] = {
                    'method':          'undecided',
                    'mode':            mode,
                    'confidence':      0.0,
                    'winner_class':    cls_set[0] if cls_set else '?',
                    'winner_compound': ' / '.join(cmp_set[:3]),
                    'winner_adduct':   '',
                    'loser_class':     ' / '.join(cls_set[1:]) if len(cls_set) > 1 else '',
                    'loser_compound':  '',
                    'loser_adduct':    '',
                    'pred_winner':     float('nan'),
                    'pred_loser':      float('nan'),
                    'res_winner':      float('nan'),
                    'res_loser':       float('nan'),
                    'low_confidence':  True,
                }

    def _clear_coherence_highlight(self, *args, **kwargs):
        """Report 系ダイアログ close 時に呼ぶ。
        _is_coherence_highlight タグ付き artist を全削除して再描画。
        QDialog.finished シグナル(int 引数)からも呼べるよう *args 受容。"""
        removed = False
        for artist in list(self._overlay_artists):
            if getattr(artist, '_is_coherence_highlight', False):
                try:
                    artist.remove()
                    removed = True
                except Exception:
                    pass
                try:
                    self._overlay_artists.remove(artist)
                except Exception:
                    pass
        if removed and self._canvas is not None:
            self._canvas.draw_idle()

    def _highlight_coherence_spot(
            self, rt_r6: float, mz_r6: float, mode: str = ""):
        """Report 行クリック時の散布図ハイライト。

        - mode に応じて pos / neg どちらの軸でハイライトするか決定
        - 大きな赤エッジの中空円(○)を一時マーカーとして描画
        - 前回のハイライトは自動クリア
        - 該当軸を spot 周辺にズーム
        """
        # 既存のハイライトを削除(_is_coherence_highlight タグ付き artist)
        for artist in list(self._overlay_artists):
            if getattr(artist, '_is_coherence_highlight', False):
                try: artist.remove()
                except Exception: pass
                try: self._overlay_artists.remove(artist)
                except Exception: pass

        # 描画対象の軸を決定
        target_ax = None
        m = (mode or "").lower()
        if m == 'pos':
            target_ax = getattr(self, '_ax_sc_pos', None)
        elif m == 'neg':
            target_ax = getattr(self, '_ax_sc_neg', None)
        if target_ax is None:
            target_ax = self._ax_sc
        if target_ax is None:
            return

        # 大きな赤エッジ中空円 + 内側に十字風の薄いリング
        sc1 = target_ax.scatter(
            [rt_r6], [mz_r6],
            s=260, marker='o', facecolors='none',
            edgecolors='#E0241E', linewidths=2.4, zorder=10)
        sc1._is_coherence_highlight = True
        self._overlay_artists.append(sc1)
        sc2 = target_ax.scatter(
            [rt_r6], [mz_r6],
            s=80, marker='o', facecolors='none',
            edgecolors='#E0241E', linewidths=1.2,
            alpha=0.7, zorder=10)
        sc2._is_coherence_highlight = True
        self._overlay_artists.append(sc2)

        # 該当軸を spot 周辺にズーム(±0.3 min, ±30 m/z)
        try:
            target_ax.set_xlim(rt_r6 - 0.3, rt_r6 + 0.3)
            target_ax.set_ylim(mz_r6 - 30.0, mz_r6 + 30.0)
        except Exception:
            pass
        if self._canvas is not None:
            self._canvas.draw_idle()


    def _update_match_count_label(self, matched: pd.DataFrame, lib_df: pd.DataFrame):
        """マッチ件数ラベルを更新（IS-derived RT 窓があれば併記）"""
        base = (f"{matched['compound'].nunique()}/{len(lib_df)} compounds matched"
                f"  ({self._spin_ppm.value()} ppm)")
        if self._class_ref_rt:
            parts = []
            for c, vals in sorted(self._class_ref_rt.items()):
                ref, tol = vals[0], vals[1]
                start = ref - tol
                end   = ref + tol
                parts.append(f"{c} ({start:.2f}–{end:.2f})")
            base += f"   |   IS RT: {', '.join(parts)}"
        self._lbl_match_count.setText(base)

    # _open_fix_rt 撤廃
    # (IS Filter の手動 RT 入力で代替可能なため不要)

    def _clear_overlay(self):
        """両モードを完全にクリアする(以前は active mode のみ)。

        修正点:
          - active mode と他モードの両方で全 state をクリア
          - _runtime_is_fingerprints / _adduct_attribution もクリア
          - ⑤ Adduct Ion Filter / Patterns / Report サブボタンも disable
          - secondary 軸の overlay artists を削除し base_scatter を再表示
          - 散布図タイトルを "raw peaks only" に戻す
          - _sync_step_button_states でボタン checked 状態を全同期
        """
        # ── active mode の state をクリア ──────────────────────
        self._match_df        = None
        self._matched_coords  = set()
        self._outlier_coords  = {}
        self._coherence_outliers    = {}
        self._coherence_models      = {}
        self._coherence_assignments = {}
        self._coherence_pairs       = []
        self._coherence_all_pairs   = []
        self._conflict_coords = set()
        self._conflict_map    = {}
        self._reserved_coords = set()
        self._class_ref_rt    = {}
        self._is_filter_choices = {}
        self._adduct_attribution = {}
        # 両モード共通の fingerprints もクリア
        self._runtime_is_fingerprints = {}
        self._runtime_is_fingerprints_by_mode = {}

        # ── 他モードの snapshot をクリア(_mode_state) ──────────
        for mode_key in list(self._mode_state.keys()):
            snap = self._mode_state.get(mode_key, {}) or {}
            for var in self._PER_MODE_STATE_VARS:
                if var in ('_rt', '_mz',
                           '_sample_intensities', '_sample_col_names'):
                    continue  # ピーク表データ自体は保持
                # state を「未マッチ」相当の初期値で上書き
                if var in ('_match_df',):
                    snap[var] = None
                elif var in ('_matched_coords', '_reserved_coords',
                             '_conflict_coords'):
                    snap[var] = set()
                elif var in ('_class_ref_rt', '_conflict_map',
                             '_outlier_coords', '_coherence_outliers',
                             '_coherence_assignments', '_coherence_models',
                             '_adduct_attribution', '_is_filter_choices'):
                    snap[var] = {} if var != '_class_ref_rt' else {}
                elif var in ('_coherence_pairs', '_coherence_all_pairs'):
                    snap[var] = []
            self._mode_state[mode_key] = snap

        # ── テーブル + フィルタ UI クリア ──────────────────────
        self._match_table.setRowCount(0)
        self._lbl_match_count.setText("Overlay cleared.")
        for cb in self._class_checks.values():
            cb.deleteLater()
        self._class_checks.clear()
        # Clear はチェックボックス永続辞書も全リセット
        # (これがないと「Clear → Match」後にクラス checkbox が前回の
        # uncheck 状態のまま復元され、散布図が灰色になるバグ)
        self._user_filter_prefs = {}
        # Clear で Curated view も解除する。これがないと
        # _quant_ion_choices だけが空リセットされ、_curated_mode が True の
        # まま残るため、Clear→Match 後に active_classes が空になり
        # 散布図が全点グレーになる。
        self._curated_mode = False
        while self._filter_grid.count():
            item = self._filter_grid.takeAt(0)
            if item and item.widget():
                item.widget().deleteLater()
        self._filter_container.setFixedHeight(22)  # 空の最小高さ
        self._filter_container.updateGeometry()

        # ── ステップボタンを全て disable ─────────────────────
        # ⑤ Adduct Ion Filter とサブボタンも追加で無効化
        for btn in (
            getattr(self, '_btn_is_filter', None),
            getattr(self, '_btn_conflict', None),
            getattr(self, '_btn_rt_outlier', None),
            getattr(self, '_btn_coherence', None),
            getattr(self, '_btn_adduct_filter', None),
            getattr(self, '_btn_adduct_patterns', None),
            getattr(self, '_btn_adduct_report', None),
            getattr(self, '_btn_coherence_pairs', None),
            getattr(self, '_btn_coherence_report', None),
        ):
            if btn is not None:
                btn.setEnabled(False)
        # トグル checked 状態を全部 false に同期
        try:
            self._sync_step_button_states()
        except Exception:
            pass

        # ── active 軸の overlay artists 削除 ────────────────
        for artist in self._overlay_artists:
            try: artist.remove()
            except Exception: pass
        self._overlay_artists.clear()
        # ── secondary 軸の overlay artists も削除 ──────────────
        if hasattr(self, '_other_overlay_artists'):
            for artist in self._other_overlay_artists:
                try: artist.remove()
                except Exception: pass
            self._other_overlay_artists = []

        # ── 両軸の base_scatter を再表示 + タイトル戻し ─────────
        try:
            self._base_scatter.set_visible(True)
        except Exception:
            pass
        for ax_attr, mode_label in (('_ax_sc_pos', 'pos'),
                                     ('_ax_sc_neg', 'neg')):
            ax = getattr(self, ax_attr, None)
            if ax is not None:
                try:
                    ax.set_title(f"RT vs m/z  ({mode_label}) — raw peaks only")
                except Exception:
                    pass
        # secondary base_scatter も再表示
        for base_attr in ('_base_scatter_pos', '_base_scatter_neg'):
            base = getattr(self, base_attr, None)
            if base is not None:
                try:
                    base.set_visible(True)
                except Exception:
                    pass
        self._canvas.draw_idle()

    # ── チェックボックス構築 ──────────────────────────────────────────
    def _rebuild_filters(self, classes: list[str]):
        """2グループ構成のチェックボックスを再構築する:

        グループ1 (所属ラベル, OR条件):
          クラス名 (PC, PE, ...) / Not matched / Conflicts
        グループ2 (修飾ラベル, フィルタ):
          Include RT outliers

        Mode セレクタ撤廃に伴い、両モードの class を union
        して表示する。引数 classes(active mode 由来)に他モードの classes を
        合算する。
        _suppress_overlay_draw 中は rebuild をスキップ
        (両モード iteration 後の最終 rebuild で済むため)。
        """
        # multi-mode iteration 中はスキップ
        if getattr(self, '_suppress_overlay_draw', False):
            return
        # 他モードのクラスを統合
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}
        other_match_df = other_state.get('_match_df')
        other_classes: list[str] = []
        if other_match_df is not None and not other_match_df.empty:
            try:
                other_classes = (
                    other_match_df['lipid_class'].dropna().unique().tolist())
            except Exception:
                other_classes = []
        # union を保持(active mode の順序を優先)
        seen = set(classes)
        for c in other_classes:
            if c not in seen:
                classes = list(classes) + [c]
                seen.add(c)
        # 再構築前のチェック状態をスナップショットして復元可能にする
        prev_checked: dict[str, bool] = {
            key: cb.isChecked()
            for key, cb in self._class_checks.items()
        }
        # deleteLater 前に setVisible(False) で即座に隠す。
        # deleteLater は次の event loop で実行されるため、その間に
        # 古い widget と新しい widget が同時に描画される可能性がある
        # (特に modal dialog 表示中)。先に hide することで visual overlap を防ぐ。
        for cb in self._class_checks.values():
            try:
                cb.setVisible(False)
            except Exception:
                pass
            cb.deleteLater()
        self._class_checks.clear()
        while self._filter_grid.count():
            item = self._filter_grid.takeAt(0)
            if item and item.widget():
                w = item.widget()
                try:
                    w.setVisible(False)
                except Exception:
                    pass
                w.deleteLater()

        COLS = 3
        # ─── グループ1: 所属ラベル (OR条件) ───────────────────────
        membership_items = list(classes) + [self._TAG_NOT_MATCHED]
        if self._conflict_coords:
            membership_items.append(self._TAG_CONFLICTS)

        # ─── グループ2: 修飾ラベル (フィルタ) ─────────────────────
        has_outliers  = any(len(v) > 0 for v in self._outlier_coords.values())
        has_coh_out   = any(len(v) > 0 for v in self._coherence_outliers.values())
        has_lowconf   = any(a['low_confidence']
                            for a in self._coherence_assignments.values())
        filter_items = []
        if has_outliers:
            filter_items.append(self._TAG_RT_OUTLIERS)
        if has_coh_out:
            filter_items.append(self._TAG_COHERENCE_OUTLIERS)
        if self._coherence_assignments:
            # 低信頼度の有無に関わらず帰属があれば Highlight チェックを出す
            filter_items.append(self._TAG_LOWCONF)
        # IS Filter で reject された候補の表示切替(あれば)
        # mode 切替直後など match_df が None の場合に備えガード
        if (self._match_df is not None
                and 'is_filter_status' in self._match_df.columns):
            has_is_rej = (self._match_df['is_filter_status']
                          == STATUS_REJ_WINDOW).any()
            if has_is_rej:
                filter_items.append(self._TAG_SHOW_IS_REJECTED)
        # Adduct Filter で reject された候補の表示切替(あれば)
        if (self._match_df is not None
                and 'adduct_filter_status' in self._match_df.columns):
            has_adduct_rej = (self._match_df['adduct_filter_status']
                              == STATUS_REJ_ADDUCT).any()
            if has_adduct_rej:
                filter_items.append(self._TAG_SHOW_ADDUCT_REJECTED)

        # Curated view (Restrict display)
        # ⑦ 未実行でも常時表示(クリック時に警告)
        restrict_items = [self._TAG_CURATED_VIEW,
                          self._TAG_IS_ONLY]

        row = 0

        # グループ1ヘッダ
        hdr1 = QLabel("Include (any of):")
        hdr1.setStyleSheet(
            "font-size: 10px; color: #666; font-weight: bold; "
            "padding-top: 2px;")
        self._filter_grid.addWidget(hdr1, row, 0, 1, COLS)
        row += 1

        # グループ1チェックボックス
        group1_start_row = row
        for i, key in enumerate(membership_items):
            self._filter_grid.addWidget(
                self._make_filter_checkbox(key, prev_checked),
                group1_start_row + i // COLS, i % COLS)
        group1_rows = (len(membership_items) + COLS - 1) // COLS
        row = group1_start_row + group1_rows

        # グループ2: View filters (Expand display + Restrict display)
        # "Filters:" → "View filters:" にリネーム、
        # Expand / Restrict のサブヘッダで内訳を視覚的に区分。
        if filter_items or restrict_items:
            hdr2 = QLabel("View filters:")
            hdr2.setStyleSheet(
                "font-size: 10px; color: #666; font-weight: bold; "
                "padding-top: 6px;")
            self._filter_grid.addWidget(hdr2, row, 0, 1, COLS)
            row += 1

            # Expand display サブヘッダ + 既存フィルタ
            if filter_items:
                sub_exp = QLabel("Expand display:")
                sub_exp.setStyleSheet(
                    "font-size: 9px; color: #888; font-style: italic; "
                    "padding-left: 8px;")
                self._filter_grid.addWidget(sub_exp, row, 0, 1, COLS)
                row += 1
                group2_start_row = row
                for i, key in enumerate(filter_items):
                    self._filter_grid.addWidget(
                        self._make_filter_checkbox(key, prev_checked),
                        group2_start_row + i // COLS, i % COLS)
                group2_rows = (len(filter_items) + COLS - 1) // COLS
                row = group2_start_row + group2_rows

            # Restrict display サブヘッダ + Curated view
            if restrict_items:
                sub_res = QLabel("Restrict display:")
                sub_res.setStyleSheet(
                    "font-size: 9px; color: #888; font-style: italic; "
                    "padding-left: 8px; padding-top: 4px;")
                self._filter_grid.addWidget(sub_res, row, 0, 1, COLS)
                row += 1
                group3_start_row = row
                for i, key in enumerate(restrict_items):
                    self._filter_grid.addWidget(
                        self._make_filter_checkbox(key, prev_checked),
                        group3_start_row + i // COLS, i % COLS)
                group3_rows = (len(restrict_items) + COLS - 1) // COLS
                row = group3_start_row + group3_rows

        # 高さ調整: ヘッダ + チェックボックス行数
        total_rows = row
        self._filter_container.setFixedHeight(total_rows * 24 + 12)
        self._filter_container.updateGeometry()

    def _make_filter_checkbox(
        self, key: str, prev_checked: dict[str, bool]
    ) -> QCheckBox:
        """チェックボックス1つを生成する（ラベル・色・ツールチップ設定済み）"""
        tip = ""
        if key == self._TAG_NOT_MATCHED:
            label, color = "Not matched", "#B4B2A9"
        elif key == self._TAG_RT_OUTLIERS:
            label, color = "Include RT outliers", "#F0997B"
            tip = ("If checked, RT outliers (× markers) are included in the\n"
                   "display. Uncheck to hide outlier spots.")
        elif key == self._TAG_COHERENCE_OUTLIERS:
            label, color = "Include Coherence outliers", "#D25C1F"
            tip = ("If checked, Coherence outliers (+ markers) are included\n"
                   "in the display. Uncheck to hide trend-based outlier spots.\n"
                   "(Coherence outliers are spots whose residual from the\n"
                   "pooled RT = α·C + β·U + γ_class regression within their\n"
                   "class×link group exceeds MAD × k.)")
        elif key == self._TAG_LOWCONF:
            label, color = "Highlight low-confidence", "#B22222"
            tip = ("If checked, attributed conflict spots with confidence\n"
                   "below the σ threshold are highlighted with a dark red\n"
                   "edge. Uncheck to render them the same as high-confidence\n"
                   "attributions (no visual distinction).")
        elif key == self._TAG_CONFLICTS:
            label, color = "Conflicts", "#E8B84B"
            tip = ("If checked, conflict spots (△ markers) are shown\n"
                   "when any of the conflicting classes is also checked.")
        elif key == self._TAG_SHOW_IS_REJECTED:
            label, color = "Show IS-rejected", "#A0A0A0"
            tip = ("If checked, candidates outside the IS RT window are\n"
                   "displayed as light gray spots. Uncheck to hide them\n"
                   "(default).")
        elif key == self._TAG_SHOW_ADDUCT_REJECTED:
            label, color = "Show Adduct-rejected", "#A088A8"
            tip = ("If checked, candidates rejected by Adduct Ion Filter\n"
                   "(loser of conflict resolution) are displayed as light\n"
                   "purple spots. Uncheck to hide them (default).")
        elif key == self._TAG_IS_ONLY:
            # IS のスポットだけを出す
            label, color = "IS only", "#7B4FA3"
            tip = ("If checked, show only internal standard spots\n"
                   "(rows with is_IS = True). Everything else is hidden,\n"
                   "including the gray 'Not matched' background, so the IS\n"
                   "positions are easy to read off the RT axis.\n\n"
                   "Class checkboxes still apply, so you can look at the IS\n"
                   "of one class at a time. Combine with Curated view to\n"
                   "restrict further.\n\n"
                   "The match table is not filtered — only the plots.")
        elif key == self._TAG_CURATED_VIEW:
            # Curated view を View filters > Restrict
            # display 内のチェックボックスとして提供
            label, color = "Curated view", "#028090"
            tip = ("If checked, restrict the display to lipid classes\n"
                   "selected as 'pos' or 'neg' in ⑦ Select Quant Ion.\n"
                   "Class checkboxes auto-sync to that selection.\n"
                   "Skip / unselected classes are hidden.")
        else:
            color = self._class_colors.get(key, '#888888')
            # Matched features は両モード共通表示なので、
            # active モードだけでなく snapshot の他モード _class_ref_rt も
            # 合わせて鍵マーク判定する。いずれかのモードで IS RT window が
            # 立っていれば鍵を表示。
            ref_active = (
                self._class_ref_rt.get(key)
                if self._class_ref_rt else None)
            other_mode = (
                'neg' if getattr(self, '_active_mode', 'pos') == 'pos'
                else 'pos')
            other_snap = (self._mode_state or {}).get(other_mode, {}) or {}
            other_ref_rt = other_snap.get('_class_ref_rt') or {}
            ref_other = other_ref_rt.get(key)
            if ref_active or ref_other:
                label = f"{key} 🔒"
                tip_lines = ["IS-derived RT window:"]
                if ref_active:
                    ref, tol = ref_active[0], ref_active[1]
                    tip_lines.append(
                        f"  {self._active_mode}: "
                        f"{ref - tol:.3f} – {ref + tol:.3f} min  "
                        f"(ref={ref:.3f}, ±{tol:.3f})")
                if ref_other:
                    ref, tol = ref_other[0], ref_other[1]
                    tip_lines.append(
                        f"  {other_mode}: "
                        f"{ref - tol:.3f} – {ref + tol:.3f} min  "
                        f"(ref={ref:.3f}, ±{tol:.3f})")
                tip = "\n".join(tip_lines)
            else:
                label = key

        cb = QCheckBox(label)
        # 前回の状態を維持。新規キー(初回登場や新たに追加されたタグ)の
        # デフォルト値はキー種別で分ける:
        #   - 外れ値フィルタ系(RT outlier / Coherence outlier): デフォルト OFF
        #     → ステップ実行で表示が段階的にクリーンになる UX
        #   - Not matched: デフォルト OFF
        #     帰属プロセスの集中視認性向上のため
        #   - Conflicts: デフォルト ON
        #     衝突 △ マーカーは帰属で重要なので revert to ON
        #   - それ以外(クラス, Highlight low-confidence): ON
        if key in (self._TAG_RT_OUTLIERS, self._TAG_COHERENCE_OUTLIERS,
                   self._TAG_SHOW_IS_REJECTED,
                   self._TAG_SHOW_ADDUCT_REJECTED,
                   self._TAG_LOWCONF,                # B-3: 既定 OFF へ
                   self._TAG_NOT_MATCHED,
                   self._TAG_IS_ONLY,                # 既定 OFF
                   self._TAG_CURATED_VIEW):          #
            # Curated view は ⑦ Select Quant Ion 実行後に明示 ON する運用
            default_state = False
        else:
            default_state = True
        # 永続辞書(_user_filter_prefs)を最優先で参照。
        # モード切替で消えたクラスの状態も復元できる。
        if key in self._user_filter_prefs:
            cb.setChecked(self._user_filter_prefs[key])
        elif key in prev_checked:
            cb.setChecked(prev_checked[key])
        else:
            cb.setChecked(default_state)
        # チェックボックス個別の最大幅を制限（3列×95pxで右パネル内に収める）
        cb.setMaximumWidth(115)  # "Include RT outliers" が入るよう少し広げる
        cb.setMinimumWidth(0)
        if tip:
            cb.setToolTip(tip)
        cb.setStyleSheet(
            f"QCheckBox {{ font-size: 11px; }}"
            f"QCheckBox::indicator:checked {{"
            f"  background-color: {color};"
            f"  border: 1px solid #555; border-radius: 2px; }}")
        # フィルタ操作を両モードに同時反映する
        # 状態変更時に永続辞書を更新
        def _on_state_change(_state, _key=key, _cb=cb):
            self._user_filter_prefs[_key] = _cb.isChecked()
            self._draw_both_overlays()
        cb.stateChanged.connect(_on_state_change)
        # Curated view 専用の追加ハンドラ
        if key == self._TAG_CURATED_VIEW:
            cb.stateChanged.connect(self._on_curated_view_filter_toggled)
        self._class_checks[key] = cb
        return cb

    # ── 散布図描画 ────────────────────────────────────────────────────
    def _draw_other_mode_overlay(self):
        """secondary モード(_ax_sc_other)に Not matched 灰色 +
        library overlay クラス色 + 特殊マーカー(△ conflict、× RT outlier、
        P coherence outlier)を描画する。

        active 軸の _draw_overlay と視覚的にほぼ同等のレンダリングになる。
        legend、低信頼度ハイライト、attr 帰属は省略(secondary は参照表示)。
        """
        if not hasattr(self, '_ax_sc_other') or self._ax_sc_other is None:
            return

        # 前回の overlay artist を削除
        if not hasattr(self, '_other_overlay_artists'):
            self._other_overlay_artists = []
        for art in self._other_overlay_artists:
            try:
                art.remove()
            except Exception:
                pass
        self._other_overlay_artists = []

        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}

        # 他モードのピーク表データ取得(_mode_state 優先)
        other_rt = other_state.get('_rt')
        other_mz = other_state.get('_mz')
        if other_rt is None or other_mz is None or len(other_rt) == 0:
            try:
                other_rt, other_mz, _ = self._get_other_mode_data()
            except Exception:
                other_rt, other_mz = None, None

        if other_rt is None or other_mz is None or len(other_rt) == 0:
            self._ax_sc_other.set_title(
                f"RT vs m/z  ({other_mode})  — no data")
            return

        # 他モードの状態を取得
        other_match_df = other_state.get('_match_df')
        # match_df が無い(Match Overlay 未実行)場合は base_scatter
        # を表示したまま raw peaks のみを見せる。match_df があるときだけ hide
        # して overlay artists で塗り直す(従来挙動)。
        other_base = (self._base_scatter_neg
                      if other_mode == 'neg'
                      else getattr(self, '_base_scatter_pos', None))
        if other_base is not None:
            try:
                if other_match_df is not None and not other_match_df.empty:
                    other_base.set_visible(False)
                else:
                    other_base.set_visible(True)
            except Exception:
                pass
        other_conflict_coords = other_state.get('_conflict_coords') or set()
        other_outlier_coords = other_state.get('_outlier_coords') or {}
        other_coh_outliers = other_state.get('_coherence_outliers') or {}
        other_class_colors = other_state.get('_class_colors') or {}
        active_class_colors = getattr(self, '_class_colors', {}) or {}

        # secondary 軸でも ⑤ Adduct Ion Filter / ⑥ Coherence
        # で帰属済みの座標は △ として描画しない(active 軸と同じ挙動)。
        other_adduct_attribution = (
            other_state.get('_adduct_attribution') or {})
        other_coh_assignments = (
            other_state.get('_coherence_assignments') or {})
        other_adduct_resolved = {
            coord for coord, result in other_adduct_attribution.items()
            if result.get('winner') is not None
            and not result.get('skipped', False)
        }
        # 他モードの手動 winner も「解決済み」扱い
        other_manual_winner_coords = (
            other_state.get('_manual_winner_coords') or set())
        # 未解決の衝突 = 全衝突 − ⑤ 帰属済み − ⑥ 帰属済み − 手動 winner
        other_unresolved_conflict_coords = (
            other_conflict_coords
            - other_adduct_resolved
            - set(other_coh_assignments.keys())
            - other_manual_winner_coords
        )

        # matched_coords を計算
        # n_matched は表示中クラスでフィルタした件数も併記
        matched_coords: set = set()
        n_matched = 0
        n_shown_other = 0
        if other_match_df is not None and not other_match_df.empty:
            m = other_match_df[other_match_df['matched']]
            if 'final_status' in m.columns:
                m = m[m['final_status'] == STATUS_KEPT]
            if not m.empty:
                matched_coords = set(zip(
                    m['obs_rt'].round(6).values,
                    m['obs_mz'].round(6).values))
                if self._is_only_enabled():
                    m = self._restrict_to_is(m)
                n_matched = len(m)
                # チェック中クラスでフィルタ
                non_class_tags = self._filter_modifier_tags() | {
                    self._TAG_NOT_MATCHED, self._TAG_CONFLICTS}
                active_classes_set = {
                    cls for cls, cb in self._class_checks.items()
                    if cls not in non_class_tags and cb.isChecked()
                }
                if active_classes_set:
                    n_shown_other = len(
                        m[m['lipid_class'].isin(active_classes_set)])

        # Not matched(灰色)を描画
        # Not matched チェックボックスも secondary に反映
        nm_cb = self._class_checks.get(self._TAG_NOT_MATCHED)
        show_nm = nm_cb is None or nm_cb.isChecked()
        _is_only_other = self._is_only_enabled()
        if _is_only_other:
            show_nm = False
        if show_nm:
            if matched_coords:
                nm_mask = np.array([
                    (round(float(r), 6), round(float(mz), 6))
                    not in matched_coords
                    for r, mz in zip(other_rt, other_mz)])
                if nm_mask.any():
                    sc_nm = self._ax_sc_other.scatter(
                        other_rt[nm_mask], other_mz[nm_mask],
                        s=4, color='#B4B2A9', alpha=0.35, zorder=2,
                        edgecolors='none')
                    self._other_overlay_artists.append(sc_nm)
            else:
                sc_nm = self._ax_sc_other.scatter(
                    other_rt, other_mz,
                    s=4, color='#B4B2A9', alpha=0.35, zorder=2,
                    edgecolors='none')
                self._other_overlay_artists.append(sc_nm)

        # secondary 軸でも Include RT outliers /
        # Include Coherence toggle に応じて rejected 行を all_matched に含める。
        # toggle 状態は active panel の class_checks から読む(active 軸と共有)。
        _show_out_cb_other = self._class_checks.get(self._TAG_RT_OUTLIERS)
        _show_coh_cb_other = self._class_checks.get(self._TAG_COHERENCE_OUTLIERS)
        _show_out_other = (_show_out_cb_other is None
                           or _show_out_cb_other.isChecked())
        _show_coh_other = (_show_coh_cb_other is None
                           or _show_coh_cb_other.isChecked())
        _show_is_rej_cb_other = self._class_checks.get(
            self._TAG_SHOW_IS_REJECTED)
        _show_is_rej_other = (_show_is_rej_cb_other is not None
                              and _show_is_rej_cb_other.isChecked())
        _show_add_rej_cb_other = self._class_checks.get(
            self._TAG_SHOW_ADDUCT_REJECTED)
        _show_add_rej_other = (_show_add_rej_cb_other is not None
                               and _show_add_rej_cb_other.isChecked())

        # secondary 軸でも rejected 行は active_classes に限定。
        _NON_CLASS_TAGS_OTHER = (
            self._TAG_NOT_MATCHED, self._TAG_RT_OUTLIERS,
            self._TAG_COHERENCE_OUTLIERS, self._TAG_CONFLICTS,
            self._TAG_LOWCONF,
            self._TAG_SHOW_IS_REJECTED, self._TAG_SHOW_ADDUCT_REJECTED,
            self._TAG_IS_ONLY,
        )
        _active_classes_other = {
            cls for cls, cb in self._class_checks.items()
            if cls not in _NON_CLASS_TAGS_OTHER and cb.isChecked()
        }

        # クラス別に library overlay + 特殊マーカーを描画
        if other_match_df is not None and not other_match_df.empty:
            matched = other_match_df[other_match_df['matched']]
            # / Curated mode 中は
            # other_mode で _quant_ion_choices[cls] == other_mode の
            # クラスのみ残す
            if getattr(self, '_curated_mode', False):
                qi = getattr(self, '_quant_ion_choices', None) or {}
                curated_other = {
                    str(c) for c, ch in qi.items() if ch == other_mode
                }
                if curated_other:
                    matched = matched[
                        matched['lipid_class'].astype(str).isin(curated_other)
                    ]
                else:
                    matched = matched.iloc[0:0]
            # IS only — IS のスポットだけに絞る
            if _is_only_other:
                matched = self._restrict_to_is(matched)
            if 'final_status' in matched.columns:
                _kept_m = matched['final_status'] == STATUS_KEPT
                _act_cls_m_other = matched['lipid_class'].isin(
                    _active_classes_other)
                _extra_m = pd.Series(False, index=matched.index)
                if (_show_is_rej_other
                        and 'is_filter_status' in matched.columns):
                    _extra_m = _extra_m | (
                        _act_cls_m_other
                        & (matched['is_filter_status'] == STATUS_REJ_WINDOW))
                if (_show_add_rej_other
                        and 'adduct_filter_status' in matched.columns):
                    _extra_m = _extra_m | (
                        _act_cls_m_other
                        & (matched['adduct_filter_status']
                           == STATUS_REJ_ADDUCT))
                if (_show_out_other
                        and 'rt_outlier_status' in matched.columns):
                    _extra_m = _extra_m | (
                        _act_cls_m_other
                        & (matched['rt_outlier_status']
                           == STATUS_REJ_OUTLIER))
                if (_show_coh_other
                        and 'coherence_status' in matched.columns):
                    _extra_m = _extra_m | (
                        _act_cls_m_other
                        & (matched['coherence_status']
                           == STATUS_REJ_RESIDUAL))
                matched = matched[_kept_m | _extra_m]
            if not matched.empty:
                # 行ごとに座標と特殊状態をマーク
                rt_r6 = matched['obs_rt'].round(6).values
                mz_r6 = matched['obs_mz'].round(6).values
                cls_arr = matched['lipid_class'].values
                n_rows = len(matched)

                # 衝突座標フラグ
                is_conf_mask = np.array([
                    (float(rt_r6[i]), float(mz_r6[i]))
                    in other_unresolved_conflict_coords
                    for i in range(n_rows)
                ]) if other_unresolved_conflict_coords else np.zeros(
                    n_rows, dtype=bool)
                # status 列を「単一の真実のソース」とする。
                # 旧 _outlier_coords / _coherence_outliers セットは status 列と
                # desync する場合があるため、描画では他モードの match_df の
                # rt_outlier_status / coherence_status 列を直接見る。
                if 'rt_outlier_status' in matched.columns:
                    _rt_out_st = matched['rt_outlier_status'].values
                    is_out_mask = (_rt_out_st == STATUS_REJ_OUTLIER)
                else:
                    is_out_mask = np.zeros(n_rows, dtype=bool)
                # Coherence outlier フラグ(残差外れ値のみ; loser は
                # final_status != kept で既に除外されている)
                if 'coherence_status' in matched.columns:
                    _coh_st = matched['coherence_status'].values
                    is_coh_mask = (_coh_st == STATUS_REJ_RESIDUAL)
                else:
                    is_coh_mask = np.zeros(n_rows, dtype=bool)

                # active panel の class_checks フィルタを
                # secondary 軸にも適用する。クラスのチェックボックスが OFF
                # なら描画スキップ。チェックボックスが存在しない他モード固有
                # クラスはデフォルト ON 扱い。
                # / Curated mode 中は
                # other mode 側も _quant_ion_choices[cls] == other_mode の
                # クラスだけ True
                _curated_on = getattr(self, '_curated_mode', False)
                cls_show_map = {}
                if _curated_on:
                    qi_other = getattr(
                        self, '_quant_ion_choices', None) or {}
                    curated_other = {
                        str(c) for c, ch in qi_other.items()
                        if ch == other_mode
                    }
                    for cls in set(matched['lipid_class'].unique()):
                        cls_show_map[cls] = str(cls) in curated_other
                else:
                    for cls in set(matched['lipid_class'].unique()):
                        cb = self._class_checks.get(cls)
                        cls_show_map[cls] = cb is None or cb.isChecked()
                # tag チェックボックスの状態取得
                show_conf_cb = self._class_checks.get(self._TAG_CONFLICTS)
                show_out_cb = self._class_checks.get(self._TAG_RT_OUTLIERS)
                show_coh_cb = self._class_checks.get(
                    self._TAG_COHERENCE_OUTLIERS)
                show_conf_g = (show_conf_cb is None
                               or show_conf_cb.isChecked())
                show_out_g = (show_out_cb is None
                              or show_out_cb.isChecked())
                show_coh_g = (show_coh_cb is None
                              or show_coh_cb.isChecked())

                # 手動 winner mask
                is_manual_mask = np.array([
                    (float(rt_r6[i]), float(mz_r6[i]))
                    in other_manual_winner_coords
                    for i in range(n_rows)
                ]) if other_manual_winner_coords else np.zeros(
                    n_rows, dtype=bool)

                # クラスごとに描画
                for cls, grp in matched.groupby('lipid_class'):
                    if not cls_show_map.get(cls, True):
                        continue  # フィルタで OFF
                    idx_pos = [matched.index.get_loc(i) for i in grp.index]
                    g_conf = is_conf_mask[idx_pos]
                    g_out = is_out_mask[idx_pos]
                    g_coh = is_coh_mask[idx_pos]
                    g_manual = is_manual_mask[idx_pos]
                    color = (
                        other_class_colors.get(cls)
                        or active_class_colors.get(cls)
                        or '#888888')

                    # 排他的に分類: RT outlier > Coherence outlier > conflict > normal
                    # 手動 winner も通常 (normal) として描画
                    normal_g = grp[~g_out & ~g_coh & ~g_conf]
                    conflict_g = grp[~g_out & ~g_coh & g_conf]
                    coh_out_g = grp[~g_out & g_coh]
                    rt_out_g = grp[g_out]

                    if len(normal_g):
                        sc = self._ax_sc_other.scatter(
                            normal_g['obs_rt'].values,
                            normal_g['obs_mz'].values,
                            s=22, color=color, alpha=0.55, zorder=5,
                            edgecolors='none')
                        self._other_overlay_artists.append(sc)
                    if show_conf_g and len(conflict_g):
                        sc_conf = self._ax_sc_other.scatter(
                            conflict_g['obs_rt'].values,
                            conflict_g['obs_mz'].values,
                            s=60, marker='^', facecolors='none',
                            edgecolors=color, linewidths=1.5,
                            zorder=5, alpha=0.9)
                        self._other_overlay_artists.append(sc_conf)
                    # 手動 winner は normal と統合(別マーカー廃止)
                    if show_coh_g and len(coh_out_g):
                        sc_coh = self._ax_sc_other.scatter(
                            coh_out_g['obs_rt'].values,
                            coh_out_g['obs_mz'].values,
                            s=55, marker='P', color=color, alpha=0.7,
                            zorder=6, edgecolors='#D25C1F',
                            linewidths=0.8)
                        self._other_overlay_artists.append(sc_coh)
                    if show_out_g and len(rt_out_g):
                        sc_ro = self._ax_sc_other.scatter(
                            rt_out_g['obs_rt'].values,
                            rt_out_g['obs_mz'].values,
                            s=45, marker='x', color=color, alpha=0.5,
                            zorder=6, linewidths=1.2)
                        self._other_overlay_artists.append(sc_ro)

        # タイトル更新
        # Curated mode 中は "Curated" prefix
        _curated_prefix2 = (
            "Curated " if getattr(self, '_curated_mode', False) else "")
        if self._is_only_enabled():
            _curated_prefix2 += "IS only "
        if n_matched > 0:
            if n_shown_other == n_matched:
                title = (f"{_curated_prefix2}RT vs m/z  ({other_mode}) — "
                         f"{n_shown_other} matched (library overlay)")
            else:
                title = (f"{_curated_prefix2}RT vs m/z  ({other_mode}) — "
                         f"{n_shown_other} / {n_matched} matched "
                         "(library overlay)")
        else:
            title = f"RT vs m/z  ({other_mode}) — raw peaks only"
        self._ax_sc_other.set_title(title)


    def _draw_overlay(self):
        # 両モード処理中の中間描画を抑制
        if getattr(self, '_suppress_overlay_draw', False):
            return
        for artist in self._overlay_artists:
            try: artist.remove()
            except Exception: pass
        self._overlay_artists.clear()

        if self._match_df is None:
            self._base_scatter.set_visible(True)
            self._canvas.draw_idle()
            return

        show_not_matched  = self._class_checks.get(self._TAG_NOT_MATCHED, None)
        show_outliers_cb  = self._class_checks.get(self._TAG_RT_OUTLIERS, None)
        show_coh_cb       = self._class_checks.get(self._TAG_COHERENCE_OUTLIERS, None)
        show_conflicts_cb = self._class_checks.get(self._TAG_CONFLICTS, None)
        show_lowconf_cb   = self._class_checks.get(self._TAG_LOWCONF, None)
        show_nm     = show_not_matched  is None or show_not_matched.isChecked()
        # IS only の間は Not matched の灰色を隠す。IS だけ見たいのに
        # 背景の灰色が残っていては意味がない(チェック状態は変えない)。
        _is_only = self._is_only_enabled()
        if _is_only:
            show_nm = False
        show_out    = show_outliers_cb  is None or show_outliers_cb.isChecked()
        show_coh    = show_coh_cb       is None or show_coh_cb.isChecked()
        show_conf   = show_conflicts_cb is None or show_conflicts_cb.isChecked()
        show_lowconf = show_lowconf_cb  is None or show_lowconf_cb.isChecked()

        # アクティブな所属クラス(OR 条件の判定対象)
        # Show *-rejected 系タグもメンバーシップ集合から除外
        _NON_CLASS_TAGS = (
            self._TAG_NOT_MATCHED, self._TAG_RT_OUTLIERS,
            self._TAG_COHERENCE_OUTLIERS, self._TAG_CONFLICTS,
            self._TAG_LOWCONF,
            self._TAG_SHOW_IS_REJECTED, self._TAG_SHOW_ADDUCT_REJECTED,
            self._TAG_IS_ONLY,
        )
        # / Curated mode 中は
        # Show classes の個別チェックボックスを無視し、_quant_ion_choices に
        # 基づき active_classes を上書きする。
        if getattr(self, '_curated_mode', False):
            qi = getattr(self, '_quant_ion_choices', None) or {}
            active_classes = {
                str(c) for c, ch in qi.items()
                if ch == self._active_mode
            }
        else:
            active_classes = {
                cls for cls, cb in self._class_checks.items()
                if cls not in _NON_CLASS_TAGS and cb.isChecked()
            }

        # final_status='kept' のみデフォルト表示
        # Show IS-rejected / Show Adduct-rejected がチェックされたら、
        # それぞれの rejected ステータスも表示に含める
        show_is_rej_cb = self._class_checks.get(self._TAG_SHOW_IS_REJECTED)
        show_is_rej = show_is_rej_cb is not None and show_is_rej_cb.isChecked()
        show_adduct_rej_cb = self._class_checks.get(
            self._TAG_SHOW_ADDUCT_REJECTED)
        show_adduct_rej = (show_adduct_rej_cb is not None
                           and show_adduct_rej_cb.isChecked())
        m = self._match_df['matched']
        if 'final_status' in self._match_df.columns:
            kept_mask = m & (self._match_df['final_status'] == STATUS_KEPT)
            # Show/Include 系 toggle で all_matched に
            # 含める rejected 行は active_classes(チェック中クラス)に
            # 限定する。SE のみチェックしている時に Include RT outliers を
            # ON しても LPI の rejected 行が拾われない、というのが期待挙動。
            _active_cls_mask = self._match_df['lipid_class'].isin(active_classes)
            extra_mask = pd.Series(False, index=self._match_df.index)
            if show_is_rej:
                extra_mask = extra_mask | (
                    m & _active_cls_mask
                    & (self._match_df.get('is_filter_status', STATUS_NA)
                         == STATUS_REJ_WINDOW))
            if show_adduct_rej and 'adduct_filter_status' in self._match_df.columns:
                extra_mask = extra_mask | (
                    m & _active_cls_mask
                    & (self._match_df['adduct_filter_status']
                         == STATUS_REJ_ADDUCT))
            # Include RT outliers / Include Coherence
            # は "表示する" toggle。final_status='rejected' の RT-outlier /
            # Coherence-residual 行を all_matched に含めるためここで extra_mask
            # に追加する(従来は display_mask での「隠す」処理のみで、
            # all_matched から既に除外されているため可視化不可だった)
            if show_out and 'rt_outlier_status' in self._match_df.columns:
                extra_mask = extra_mask | (
                    m & _active_cls_mask
                    & (self._match_df['rt_outlier_status']
                         == STATUS_REJ_OUTLIER))
            if show_coh and 'coherence_status' in self._match_df.columns:
                extra_mask = extra_mask | (
                    m & _active_cls_mask
                    & (self._match_df['coherence_status']
                         == STATUS_REJ_RESIDUAL))
            all_matched = self._match_df[kept_mask | extra_mask]
        else:
            # 後方互換: status 列がない古い match_df
            all_matched = self._match_df[m]

        # / Curated mode 中は
        # active モードで _quant_ion_choices[cls] == self._active_mode の
        # クラスのみ残す
        if getattr(self, '_curated_mode', False):
            qi = getattr(self, '_quant_ion_choices', None) or {}
            curated_set = {
                str(c) for c, ch in qi.items()
                if ch == self._active_mode
            }
            if curated_set:
                all_matched = all_matched[
                    all_matched['lipid_class'].astype(str).isin(curated_set)
                ]
            else:
                all_matched = all_matched.iloc[0:0]

        # IS only — IS のスポットだけに絞る
        if _is_only:
            all_matched = self._restrict_to_is(all_matched)

        if not all_matched.empty:
            # 行単位のマスク計算
            rt_r6 = all_matched['obs_rt'].round(6).values
            mz_r6 = all_matched['obs_mz'].round(6).values
            cls_arr = all_matched['lipid_class'].values
            n_rows = len(all_matched)

            # Coherence 帰属適用: 敗者クラスの行は表示から除外
            # (匿属後の元衝突点は勝者クラスに帰属したので、敗者側は現れない)
            is_loser_mask = np.zeros(n_rows, dtype=bool)
            if self._coherence_assignments:
                for i in range(n_rows):
                    a = self._coherence_assignments.get(
                        (float(rt_r6[i]), float(mz_r6[i])))
                    if a is not None and cls_arr[i] == a['loser_class']:
                        is_loser_mask[i] = True

            # 衝突判定(帰属未解決のもののみを △ として扱う)
            # Adduct Ion Filter で帰属できた座標も「解決済み」扱い
            # 手動 winner も「解決済み」扱い(△ 表示しない)
            adduct_resolved_coords = {
                coord for coord, result in self._adduct_attribution.items()
                if result.get('winner') is not None
                and not result.get('skipped', False)
            }
            manual_winner_coords = (
                getattr(self, '_manual_winner_coords', None) or set())
            unresolved_conf_coords = (
                self._conflict_coords
                - set(self._coherence_assignments.keys())
                - adduct_resolved_coords
                - manual_winner_coords
                if self._conflict_coords else set())
            is_conf_mask = np.array([
                (float(rt_r6[i]), float(mz_r6[i])) in unresolved_conf_coords
                for i in range(n_rows)
            ]) if unresolved_conf_coords else np.zeros(n_rows, dtype=bool)
            # 手動 winner マスク
            is_manual_winner_mask = np.array([
                (float(rt_r6[i]), float(mz_r6[i])) in manual_winner_coords
                for i in range(n_rows)
            ]) if manual_winner_coords else np.zeros(n_rows, dtype=bool)

            # RT 外れ値判定
            # status 列を「単一の真実のソース」とする。
            # _outlier_coords / _coherence_outliers セットは Coherence の再計算や
            # 手動変更で status 列と desync することがあるため、描画では
            # rt_outlier_status / coherence_status 列を直接見る。
            if 'rt_outlier_status' in all_matched.columns:
                _rt_out_st = all_matched['rt_outlier_status'].values
                is_out_mask = (_rt_out_st == STATUS_REJ_OUTLIER)
            else:
                is_out_mask = np.zeros(n_rows, dtype=bool)
            # Coherence 外れ値判定(帰属後の自クラストレンドから外れた点)
            # 注: STATUS_REJ_LOSER 行は final_status != kept なので
            # そもそも all_matched に含まれない。Coherence 外れ値視覚は
            # 残差外れ値 (STATUS_REJ_RESIDUAL) のみを対象とする
            # (元の _coherence_outliers セットと同じ意味)。
            if 'coherence_status' in all_matched.columns:
                _coh_st = all_matched['coherence_status'].values
                is_coh_mask = (_coh_st == STATUS_REJ_RESIDUAL)
            else:
                is_coh_mask = np.zeros(n_rows, dtype=bool)

            # 所属判定 + OR 条件
            own_class_mask = np.array(
                [c in active_classes for c in cls_arr])
            display_mask = own_class_mask | (is_conf_mask & show_conf)
            # フィルタ適用
            if not show_out:
                display_mask = display_mask & ~is_out_mask
            if not show_coh:
                display_mask = display_mask & ~is_coh_mask
            # 敗者行は常に非表示(帰属された分子種は勝者クラスとして描画される)
            display_mask = display_mask & ~is_loser_mask

            draw_df = all_matched[display_mask].copy()
            draw_is_conf = is_conf_mask[display_mask]
            draw_is_out  = is_out_mask[display_mask]
            draw_is_coh  = is_coh_mask[display_mask]

            # 帰属された(勝者側の)スポット判定 + 低信頼度フラグ
            is_attr_winner = np.zeros(len(draw_df), dtype=bool)
            is_lowconf     = np.zeros(len(draw_df), dtype=bool)
            if self._coherence_assignments:
                d_rt = draw_df['obs_rt'].round(6).values
                d_mz = draw_df['obs_mz'].round(6).values
                d_cls = draw_df['lipid_class'].values
                for i in range(len(draw_df)):
                    a = self._coherence_assignments.get(
                        (float(d_rt[i]), float(d_mz[i])))
                    if a is not None and d_cls[i] == a['winner_class']:
                        is_attr_winner[i] = True
                        if a['low_confidence']:
                            is_lowconf[i] = True

            # クラスごとに描画(色はクラスで決定)
            # 手動 winner mask を draw 用に絞る
            draw_is_manual = is_manual_winner_mask[display_mask]

            # 優先順位: RT 外れ値(×) > Coherence 外れ値(P) > 未解決衝突(△) > 通常(●)
            for cls, grp in draw_df.groupby('lipid_class'):
                idx = grp.index
                grp_pos = [draw_df.index.get_loc(i) for i in idx]
                g_conf  = draw_is_conf[grp_pos]
                g_out   = draw_is_out[grp_pos]
                g_coh   = draw_is_coh[grp_pos]
                g_attr  = is_attr_winner[grp_pos]
                g_low   = is_lowconf[grp_pos]
                g_manual = draw_is_manual[grp_pos]

                color = self._class_colors.get(cls, '#888888')
                # 排他的に分類(優先順位に従って)
                # 手動 winner は通常 (normal) と同じ描画。
                # △ から除外されるだけで、見た目は class 色の通常スポット。
                normal      = grp[~g_out & ~g_coh & ~g_conf & ~(g_attr & g_low)]
                attr_low    = grp[~g_out & ~g_coh & ~g_conf &  g_attr & g_low]
                conflict    = grp[~g_out & ~g_coh &  g_conf]
                coh_out     = grp[~g_out &  g_coh]
                rt_out      = grp[ g_out]

                if len(normal):
                    # 複数候補が重なって見えるよう alpha を下げる
                    # IS-rejected / Adduct-rejected は色を抑えて参考表示
                    if 'final_status' in normal.columns:
                        kept_only = normal[normal['final_status'] == STATUS_KEPT]
                        rej_is = normal[normal.get('is_filter_status', STATUS_NA)
                                        == STATUS_REJ_WINDOW]
                        if 'adduct_filter_status' in normal.columns:
                            rej_adduct = normal[
                                normal['adduct_filter_status']
                                == STATUS_REJ_ADDUCT]
                        else:
                            rej_adduct = normal.iloc[0:0]
                    else:
                        kept_only = normal
                        rej_is = normal.iloc[0:0]
                        rej_adduct = normal.iloc[0:0]
                    if len(kept_only):
                        sc = self._ax_sc.scatter(
                            kept_only['obs_rt'].values, kept_only['obs_mz'].values,
                            s=30, color=color, alpha=0.5, zorder=5,
                            edgecolors='none', label=cls)
                        self._overlay_artists.append(sc)
                    if len(rej_is):
                        # IS Filter で reject されたものは薄いグレーで参考表示
                        sc_rej = self._ax_sc.scatter(
                            rej_is['obs_rt'].values, rej_is['obs_mz'].values,
                            s=20, color='#C8C8C8', alpha=0.4, zorder=4,
                            edgecolors='none')
                        self._overlay_artists.append(sc_rej)
                    if len(rej_adduct):
                        # Adduct Filter で reject されたものは薄い紫で参考表示
                        sc_rej_adu = self._ax_sc.scatter(
                            rej_adduct['obs_rt'].values,
                            rej_adduct['obs_mz'].values,
                            s=20, color='#A088A8', alpha=0.4, zorder=4,
                            edgecolors='none')
                        self._overlay_artists.append(sc_rej_adu)

                if len(attr_low):
                    # 低信頼度帰属: 濃赤縁取り
                    edge = '#B22222' if show_lowconf else 'none'
                    lw   = 1.5      if show_lowconf else 0.0
                    # 複数候補対応で alpha を下げる
                    sc_al = self._ax_sc.scatter(
                        attr_low['obs_rt'].values, attr_low['obs_mz'].values,
                        s=36, color=color, alpha=0.5, zorder=6,
                        edgecolors=edge, linewidths=lw)
                    self._overlay_artists.append(sc_al)

                if len(conflict):
                    sc_conf = self._ax_sc.scatter(
                        conflict['obs_rt'].values, conflict['obs_mz'].values,
                        s=60, marker='^', facecolors='none', edgecolors=color,
                        linewidths=1.5, zorder=5, alpha=0.9)
                    self._overlay_artists.append(sc_conf)

                # 手動 winner は normal と統合(別マーカー廃止)

                if len(coh_out):
                    sc_coh = self._ax_sc.scatter(
                        coh_out['obs_rt'].values, coh_out['obs_mz'].values,
                        s=55, marker='P', color=color, alpha=0.7, zorder=6,
                        edgecolors='#D25C1F', linewidths=0.8)
                    self._overlay_artists.append(sc_coh)

                if len(rt_out):
                    sc_ro = self._ax_sc.scatter(
                        rt_out['obs_rt'].values, rt_out['obs_mz'].values,
                        s=45, marker='x', color=color, alpha=0.5, zorder=6,
                        linewidths=1.2)
                    self._overlay_artists.append(sc_ro)

        # Not matched スポット（ベース散布図の代替として描画）
        self._base_scatter.set_visible(False)
        if show_nm and self._matched_coords:
            nm_mask = np.array([
                (round(float(r), 6), round(float(m), 6))
                not in self._matched_coords
                for r, m in zip(self._rt, self._mz)])
            sc_nm = self._ax_sc.scatter(
                self._rt[nm_mask], self._mz[nm_mask],
                s=4, color='#B4B2A9', alpha=0.35, zorder=2,
                edgecolors='none')
            self._overlay_artists.append(sc_nm)

        # 凡例

        # active 軸のタイトルを「表示中クラスに絞った matched 件数」で更新
        try:
            n_total = 0   # matched & kept 全件
            n_shown = 0   # チェック中クラスでフィルタした表示件数
            if self._match_df is not None and not self._match_df.empty:
                m = self._match_df[self._match_df['matched']]
                if 'final_status' in m.columns:
                    m = m[m['final_status'] == STATUS_KEPT]
                if self._is_only_enabled():
                    m = self._restrict_to_is(m)
                n_total = len(m)
                # 表示中クラスでフィルタ(空集合は全て非表示扱い)
                non_class_tags = self._filter_modifier_tags() | {
                    self._TAG_NOT_MATCHED, self._TAG_CONFLICTS}
                active_classes_set = {
                    cls for cls, cb in self._class_checks.items()
                    if cls not in non_class_tags and cb.isChecked()
                }
                if active_classes_set:
                    m_shown = m[m['lipid_class'].isin(active_classes_set)]
                    n_shown = len(m_shown)
                else:
                    n_shown = 0
            active_axis = (self._ax_sc_pos
                           if self._active_mode == 'pos'
                           else self._ax_sc_neg)
            # Curated mode 中は "Curated" prefix
            _curated_prefix = (
                "Curated " if getattr(self, '_curated_mode', False) else "")
            if self._is_only_enabled():
                _curated_prefix += "IS only "
            if n_total > 0:
                # 表示中 / 全件 を併記
                if n_shown == n_total:
                    t = (f"{_curated_prefix}RT vs m/z  ({self._active_mode}) — "
                         f"{n_shown} matched (library overlay)")
                else:
                    t = (f"{_curated_prefix}RT vs m/z  ({self._active_mode}) — "
                         f"{n_shown} / {n_total} matched (library overlay)")
            else:
                t = f"{_curated_prefix}RT vs m/z  ({self._active_mode})"
            active_axis.set_title(t)
        except Exception as e:
            log.warning(f"[draw_overlay title] failed: {e}")
        # secondary モードの library overlay も描画
        try:
            self._draw_other_mode_overlay()
        except Exception as e:
            log.warning(f"[draw_other_mode_overlay] failed: {e}")
        self._canvas.draw_idle()

    def _combined_matched_rows(self) -> pd.DataFrame:
        """両モードの matched 行を 'mode' 列付きで統合。

        active mode は self._match_df、他モードは _mode_state スナップショット
        から取得する。片モードしかロードされていない場合はそれだけ返す。
        """
        parts: list[pd.DataFrame] = []
        # active mode
        if self._match_df is not None and not self._match_df.empty:
            m = self._match_df[self._match_df['matched']].copy()
            if not m.empty:
                m['mode'] = self._active_mode
                parts.append(m)
        # 他モード(snapshot から)
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}
        other_match_df = other_state.get('_match_df')
        if other_match_df is not None and not other_match_df.empty:
            m = other_match_df[other_match_df['matched']].copy()
            if not m.empty:
                m['mode'] = other_mode
                parts.append(m)
        if not parts:
            return pd.DataFrame()
        return pd.concat(parts, ignore_index=True)

    def _combined_conflict_state(self) -> tuple[set, dict]:
        """両モードの conflict_coords / conflict_map を統合。"""
        coords = set(self._conflict_coords or set())
        cmap = dict(self._conflict_map or {})
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}
        other_coords = other_state.get('_conflict_coords') or set()
        other_cmap = other_state.get('_conflict_map') or {}
        coords |= set(other_coords)
        for k, v in other_cmap.items():
            cmap.setdefault(k, v)
        return coords, cmap

    def _rebuild_match_table(self, matched: pd.DataFrame = None):
        """両モードの matched 行を統合表示する。

        引数 matched は後方互換のため残しているが無視する。常に
        _combined_matched_rows() の結果を表示する。
        """
        self._match_table.setRowCount(0)
        combined = self._combined_matched_rows()
        if combined.empty:
            return
        all_conflict_coords, all_conflict_map = self._combined_conflict_state()
        # 両モードの class colors を統合(active 優先、なければ snapshot)
        active_colors = getattr(self, '_class_colors', {}) or {}
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        other_state = self._mode_state.get(other_mode, {}) or {}
        other_colors = other_state.get('_class_colors') or {}
        merged_colors = {**other_colors, **active_colors}

        for _, row in combined.sort_values(
                ['mode', 'lipid_class', 'delta_ppm']).iterrows():
            r = self._match_table.rowCount()
            self._match_table.insertRow(r)
            q_color = QColor(merged_colors.get(row['lipid_class'], '#FFFFFF'))
            q_color.setAlpha(80)
            coord = (round(float(row['obs_rt']), 6),
                     round(float(row['obs_mz']), 6))
            row_mode = row.get('mode', self._active_mode)
            is_conflict = coord in all_conflict_coords
            data = {k: row[k] for k in
                    ['obs_rt', 'obs_mz', 'compound',
                     'lipid_class', 'adduct', 'delta_ppm']}
            data['mode'] = row_mode  # 右クリック curation 用
            cls_label = f"⚠ {row['lipid_class']}" if is_conflict else row['lipid_class']

            # 競合相手の文字列（自分以外のエントリを列挙）
            conflicts_with = ""
            if is_conflict and coord in all_conflict_map:
                others = [
                    f"{e['lipid_class']} {e['compound']} {e['adduct']}"
                    for e in all_conflict_map[coord]
                    if not (e['lipid_class'] == row['lipid_class']
                            and e['compound'] == row['compound'])
                ]
                conflicts_with = "; ".join(others)

            for c, val in enumerate([
                row_mode, cls_label, row['compound'], row['adduct'],
                f"{row['delta_ppm']:.2f}", f"{row['obs_rt']:.3f}",
                conflicts_with,
            ]):
                item = QTableWidgetItem(str(val))
                item.setBackground(q_color)
                item.setData(Qt.UserRole, data)
                self._match_table.setItem(r, c, item)

    # ── 手動キュレーション ─────────────────────────────────
    def _on_match_table_context_menu(self, pos):
        """マッチテーブルの右クリックメニュー。

        操作:
          - Reject this candidate    → manual_status = STATUS_REJ_MANUAL
          - Keep (override)          → manual_status = STATUS_KEPT(自動 reject を上書き)
          - Restore default          → manual_status = STATUS_NA(自動判定に戻す)
        """
        if self._match_df is None or 'manual_status' not in self._match_df.columns:
            return
        idx_at_pos = self._match_table.indexAt(pos)
        if not idx_at_pos.isValid():
            return
        # 選択行(複数選択対応、なければクリック行のみ)
        sel_rows = sorted({i.row() for i in
                           self._match_table.selectionModel().selectedRows()})
        if not sel_rows:
            sel_rows = [idx_at_pos.row()]

        menu = QMenu(self._match_table)
        act_reject  = menu.addAction("Reject this candidate")
        act_keep    = menu.addAction("Keep (override)")
        act_restore = menu.addAction("Restore default")
        menu.addSeparator()
        info = menu.addAction(f"Selected rows: {len(sel_rows)}")
        info.setEnabled(False)

        chosen = menu.exec(self._match_table.viewport().mapToGlobal(pos))
        if chosen is None:
            return
        if chosen == act_reject:
            new_status = STATUS_REJ_MANUAL
        elif chosen == act_keep:
            new_status = STATUS_KEPT
        elif chosen == act_restore:
            new_status = STATUS_NA
        else:
            return

        # 選択行の (compound, obs_rt, obs_mz, mode) で match_df のキー特定
        # _rebuild_match_table がソートして行を作るので、UserRole から拾う
        # 行の mode に応じて active / snapshot を更新する
        affected_modes: set = set()
        for r in sel_rows:
            item = self._match_table.item(r, 0)
            if item is None:
                continue
            data = item.data(Qt.UserRole)
            if not data:
                continue
            row_mode = data.get('mode', self._active_mode)
            if row_mode == self._active_mode:
                target_df = self._match_df
            else:
                snap = self._mode_state.get(row_mode, {}) or {}
                target_df = snap.get('_match_df')
            if target_df is None:
                continue
            mask = (
                (target_df['compound']    == data['compound']) &
                (target_df['lipid_class'] == data['lipid_class']) &
                (target_df['adduct']      == data['adduct']) &
                (target_df['obs_rt'].round(6) ==
                    round(float(data['obs_rt']), 6)) &
                (target_df['obs_mz'].round(6) ==
                    round(float(data['obs_mz']), 6))
            )
            target_df.loc[mask, 'manual_status'] = new_status
            affected_modes.add(row_mode)

        # final_status を再計算 + 表示更新(影響を受けたモードごとに)
        for m in affected_modes:
            if m == self._active_mode:
                self._match_df = _update_final_status_df(self._match_df)
            else:
                snap = self._mode_state.get(m, {}) or {}
                if snap.get('_match_df') is not None:
                    snap['_match_df'] = _update_final_status_df(
                        snap['_match_df'])
                    self._mode_state[m] = snap
        # 通知 to MainWindow(active mode のみ)
        try:
            self.match_df_updated.emit(self.fe.tag, self._match_df)
        except Exception:
            pass
        self._draw_overlay()
        self._rebuild_match_table()

    # ── イベントハンドラ ──────────────────────────────────────────────
    def _on_scroll(self, event):
        # 2 段散布図対応 - active 軸 / 他軸どちらでもズーム可能
        # m/z (y) 軸も両軸同時に同じ factor で zoom することで
        #                  上下散布図の m/z スケールが乖離しないようにする。
        valid_axes = []
        if hasattr(self, '_ax_sc_pos') and self._ax_sc_pos is not None:
            valid_axes.append(self._ax_sc_pos)
        if hasattr(self, '_ax_sc_neg') and self._ax_sc_neg is not None:
            valid_axes.append(self._ax_sc_neg)
        if not valid_axes:
            valid_axes = [self._ax_sc] if self._ax_sc is not None else []
        if event.inaxes not in valid_axes:
            return
        if self._toolbar:
            self._toolbar.push_current()
        factor = 0.85 if event.button == 'up' else 1.0 / 0.85
        ax_active = event.inaxes
        xd, yd = event.xdata, event.ydata

        # X 軸は sharex で連動。active 軸に set_xlim すれば他軸も追従。
        ax_active.set_xlim(
            [xd + (x - xd) * factor for x in ax_active.get_xlim()])

        # Y 軸(m/z)は独立。両軸とも同じ factor で zoom する。
        # active 軸: カーソル位置をアンカー
        # 他軸: 現範囲の中点をアンカー(自分の m/z 範囲を中心に保つ)
        for ax in valid_axes:
            ylim = ax.get_ylim()
            if ax is ax_active:
                anchor = yd
            else:
                anchor = (ylim[0] + ylim[1]) / 2.0
            ax.set_ylim([anchor + (y - anchor) * factor for y in ylim])

        self._canvas.draw_idle()

    def _on_mouse_press(self, event):
        # ダブルクリックで home リセット
        # 両軸を同時に home reset(m/z 軸も同期)
        # 右クリックで手動帰属メニュー
        valid_axes = []
        if hasattr(self, '_ax_sc_pos') and self._ax_sc_pos is not None:
            valid_axes.append(self._ax_sc_pos)
        if hasattr(self, '_ax_sc_neg') and self._ax_sc_neg is not None:
            valid_axes.append(self._ax_sc_neg)
        if not valid_axes:
            valid_axes = [self._ax_sc] if self._ax_sc is not None else []
        if event.inaxes not in valid_axes:
            return
        # 右クリック(button=3)で手動帰属メニュー
        if event.button == 3 and not event.dblclick:
            self._open_plot_curation_menu(event)
            return
        # ダブルクリックで home リセット
        if not event.dblclick:
            return
        # axis ごとに自分の home に戻す
        for ax in valid_axes:
            if ax is getattr(self, '_ax_sc_pos', None):
                xlim = getattr(self, '_home_xlim_pos', self._home_xlim)
                ylim = getattr(self, '_home_ylim_pos', self._home_ylim)
            elif ax is getattr(self, '_ax_sc_neg', None):
                xlim = getattr(self, '_home_xlim_neg', self._home_xlim)
                ylim = getattr(self, '_home_ylim_neg', self._home_ylim)
            else:
                xlim = self._home_xlim
                ylim = self._home_ylim
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        if self._toolbar:
            self._toolbar.push_current()
        self._canvas.draw_idle()

    def _raw_paths_for_mode(self, mode: str) -> list:
        """指定モードのエントリに対応する mzML パスを返す。

        散布図は pos / neg の 2 軸あるので、クリックされた軸に
        対応するエントリから引く必要がある。見つからなければ空リスト。
        """
        for fe2 in (self._all_entries or []):
            if getattr(fe2, 'ion_mode', None) != mode:
                continue
            try:
                paths = resolve_entry_raw_paths(fe2)
            except Exception as e:
                log.warning(f"raw path resolve failed ({mode}): {e}")
                paths = []
            if paths:
                return paths
        return []

    def _spot_sample_intensities(self, mode: str, coord) -> list:
        """スポット座標に対応する行の、サンプル別強度を返す。

        EIC を出すサンプルを「最も強く出ている 1 本」に決めるために使う。
        取れなければ空リスト(その場合は先頭サンプルになる)。
        """
        for fe2 in (self._all_entries or []):
            if getattr(fe2, 'ion_mode', None) != mode:
                continue
            try:
                df = fe2.df
                rt = pd.to_numeric(df[fe2.rt_col], errors='coerce').values
                mz = pd.to_numeric(df[fe2.mz_col], errors='coerce').values
                d = np.abs(rt - coord[0]) + np.abs(mz - coord[1])
                i = int(np.nanargmin(d))
                cols = fe2.sample_columns()
                return [float(v) for v in pd.to_numeric(
                    df.iloc[i][cols], errors='coerce').fillna(0.0).values]
            except Exception as e:
                log.warning(f"spot intensities failed ({mode}): {e}")
                return []
        return []

    def _open_spot_eic(self, mode: str, coord, same_coord):
        """スポットの EIC ダイアログを開く。"""
        raw = self._raw_paths_for_mode(mode)
        if not raw:
            QMessageBox.information(
                self, "No raw data",
                "This data has no raw data (mzML) attached, so no EIC can\n"
                "be shown. It happens when an external table was imported,\n"
                "or when the mzML files were moved.")
            return
        # ラベル: そのスポットに乗っている化合物名(複数クラスなら並べる)
        try:
            names = []
            for _i, r in same_coord.iterrows():
                nm = f"{r['lipid_class']} {r['compound']}"
                if nm not in names:
                    names.append(nm)
            label = " / ".join(names[:3]) + ("  …" if len(names) > 3 else "")
        except Exception:
            label = ""
        names_col = []
        for fe2 in (self._all_entries or []):
            if getattr(fe2, 'ion_mode', None) == mode:
                try:
                    names_col = list(fe2.sample_columns())
                except Exception:
                    names_col = []
                break
        dlg = SpotEICDialog(
            mz=coord[1], rt=coord[0], label=label,
            raw_paths=raw, sample_col_names=names_col,
            sample_intensities=self._spot_sample_intensities(mode, coord),
            ion_mode=mode, parent=self)
        dlg.exec()

    def _open_plot_curation_menu(self, event):
        """scatter plot 上の右クリックで手動帰属メニューを表示。

        クリック位置の近傍のスポット(matched 行)を取得し、状況に応じたメニュー:
          - 衝突 (△ 未解決): 候補クラス間で winner を選ぶ
          - matched/winner (●): 「reject this annotation」or 「restore default」
          - 外れ値 (×/P): 「include anyway」(manual_status=KEPT で復活)

        manual_status を更新して final_status を再計算する。
        """
        # クリックされた軸とそのモードを判定
        if event.inaxes is self._ax_sc_pos:
            target_mode = 'pos'
        elif event.inaxes is self._ax_sc_neg:
            target_mode = 'neg'
        else:
            return

        # 対象モードの match_df を取得
        if target_mode == self._active_mode:
            target_df = self._match_df
            target_conflict_map = self._conflict_map
        else:
            snap = self._mode_state.get(target_mode, {}) or {}
            target_df = snap.get('_match_df')
            target_conflict_map = snap.get('_conflict_map') or {}
        if target_df is None:
            return

        matched = target_df[target_df['matched']]
        if matched.empty:
            return

        # クリック近傍の最近接スポットを探す(画面座標相対距離 < 0.03)
        xlim = event.inaxes.get_xlim()
        ylim = event.inaxes.get_ylim()
        xr = max(xlim[1] - xlim[0], 1e-6)
        yr = max(ylim[1] - ylim[0], 1e-6)
        dx = (matched['obs_rt'].values - event.xdata) / xr
        dy = (matched['obs_mz'].values - event.ydata) / yr
        dist = np.sqrt(dx**2 + dy**2)
        best = int(np.argmin(dist))
        if dist[best] >= 0.03:
            return  # 近接スポット無し
        clicked_row = matched.iloc[best]
        coord = (round(float(clicked_row['obs_rt']), 6),
                 round(float(clicked_row['obs_mz']), 6))

        # 同 coord の全 matched 行を取得(複数候補クラスがあるかも)
        same_coord_mask = (
            (target_df['obs_rt'].round(6) == coord[0]) &
            (target_df['obs_mz'].round(6) == coord[1]) &
            target_df['matched']
        )
        same_coord = target_df[same_coord_mask]

        # メニュー構築
        menu = QMenu(self)
        title_action = menu.addAction(
            f"Spot [{target_mode}] RT={coord[0]:.3f} m/z={coord[1]:.4f}")
        title_action.setEnabled(False)
        menu.addSeparator()

        # このスポットの EIC を出す。生データが無い経路
        # (外部テーブル取り込み)では理由を出して無効化する。
        _raw = self._raw_paths_for_mode(target_mode)
        act_eic = menu.addAction("Show EIC…")
        if _raw:
            act_eic.setToolTip(
                "Draw the EIC at this spot's m/z.\n"
                "The single most intense sample is chosen automatically.\n"
                "Takes 15-30 sec per mzML.")
        else:
            act_eic.setEnabled(False)
            act_eic.setText("Show EIC…  (no raw mzML)")
        menu.addSeparator()

        # 衝突 (複数クラス) の場合: winner 選択
        unique_classes = same_coord['lipid_class'].unique().tolist()
        winner_actions: dict = {}
        if len(unique_classes) >= 2:
            sub = menu.addMenu("Set as winner (resolve conflict)")
            for cls in sorted(unique_classes):
                act = sub.addAction(cls)
                winner_actions[id(act)] = cls

        # 単独 or 全候補: 個別 row 操作(matched=True 各行)
        sub_each = menu.addMenu("Per-row action")
        row_actions: dict = {}
        for idx, r in same_coord.iterrows():
            label = f"{r['lipid_class']} {r['compound']} {r['adduct']}"
            sub_r = sub_each.addMenu(label[:60])
            for action_label, status in (
                ("Keep this (override)", STATUS_KEPT),
                ("Reject this", STATUS_REJ_MANUAL),
                ("Restore default", STATUS_NA),
            ):
                act = sub_r.addAction(action_label)
                row_actions[id(act)] = (idx, status)

        menu.addSeparator()
        # 全候補が誤りのとき、スポットごと一括除外/復元する
        reject_spot_act = menu.addAction(
            "Reject entire spot (none is correct)")
        restore_spot_act = menu.addAction(
            "Restore entire spot (clear manual)")
        menu.addSeparator()
        # スポットの全 status を確認するメニュー項目
        inspect_act = menu.addAction("Inspect spot status…")
        menu.addSeparator()
        cancel_act = menu.addAction("Cancel")

        # 表示位置(マウスのスクリーン座標)
        try:
            global_pos = self._canvas.mapToGlobal(
                self._canvas.mapFromGlobal(QCursor.pos()))
        except Exception:
            global_pos = QCursor.pos()
        chosen = menu.exec(global_pos)
        # モーダルメニューがマウス release を奪うと matplotlib の
        # pan/zoom ドラッグが終了せず、以後カーソル移動で図が動いてしまう。
        # 進行中の pan/zoom を明示終了してドラッグ状態を解除する。
        self._end_stuck_toolbar_drag()
        if chosen is None or chosen is cancel_act:
            return

        # このスポットの EIC を出す
        if chosen is act_eic:
            self._open_spot_eic(target_mode, coord, same_coord)
            return

        if id(chosen) in winner_actions:
            # winner 選択
            winner_cls = winner_actions[id(chosen)]
            self._apply_manual_winner(target_mode, coord, winner_cls)
        elif id(chosen) in row_actions:
            row_idx, new_status = row_actions[id(chosen)]
            self._apply_manual_row_status(target_mode, row_idx, new_status)
        elif chosen is reject_spot_act:
            # 同 coord の全候補を手動 reject(全部誤りのとき)
            self._apply_spot_status_all(
                target_mode, coord, STATUS_REJ_MANUAL)
        elif chosen is restore_spot_act:
            # 同 coord の手動指定を既定に戻す
            self._apply_spot_status_all(target_mode, coord, STATUS_NA)
        elif chosen is inspect_act:
            # spot status inspector を開く
            self._open_spot_inspector(target_mode, coord, same_coord)

    def _end_stuck_toolbar_drag(self):
        """右クリックのモーダルメニューで release イベントを奪われ、
        matplotlib のツールバー pan/zoom ドラッグが解除されずに残る問題を
        修正する。進行中の pan/zoom があれば motion への drag 接続を切り、
        状態を破棄する(ラバーバンドも消す)。pan/zoom 中でなければ何もしない。"""
        tb = getattr(self, '_toolbar', None)
        if tb is None:
            return
        try:
            active = (getattr(tb, '_pan_info', None) is not None
                      or getattr(tb, '_zoom_info', None) is not None)
            if not active:
                return
            idd = getattr(tb, '_id_drag', None)
            if idd is not None:
                try:
                    tb.canvas.mpl_disconnect(idd)
                except Exception:
                    pass
            if hasattr(tb, '_pan_info'):
                tb._pan_info = None
            if hasattr(tb, '_zoom_info'):
                tb._zoom_info = None
            try:
                tb.remove_rubberband()
            except Exception:
                pass
            try:
                self._canvas.draw_idle()
            except Exception:
                pass
        except Exception as e:
            log.warning(f"[_end_stuck_toolbar_drag] {e}")

    def _open_spot_inspector(self, mode: str, coord: tuple, rows_df):
        """指定 spot の全 status を表示するダイアログを開く。
        rows_df は同 coord の全 matched 行(複数候補クラス含む)。
        kept / rejected の両方が含まれるので、なぜ rejected なのに表示されている(
        または逆)かを診断できる。"""
        dlg = SpotInspectorDialog(
            mode=mode, coord=coord, rows_df=rows_df, parent=self)
        dlg.exec()

    def _apply_manual_winner(self, mode: str, coord: tuple, winner_cls: str):
        """衝突座標で winner クラスを手動指定。

        - winner クラスの行: manual_status = STATUS_KEPT
        - 他クラスの同 coord 行: manual_status = STATUS_REJ_MANUAL
        - final_status 再計算 → エクスポートで accepted として扱われる
        """
        if mode == self._active_mode:
            df = self._match_df
        else:
            snap = self._mode_state.get(mode, {}) or {}
            df = snap.get('_match_df')
            if df is None:
                return
        same = (
            (df['obs_rt'].round(6) == coord[0]) &
            (df['obs_mz'].round(6) == coord[1]) &
            df['matched']
        )
        winner = same & (df['lipid_class'] == winner_cls)
        loser = same & (df['lipid_class'] != winner_cls)
        df.loc[winner, 'manual_status'] = STATUS_KEPT
        df.loc[loser, 'manual_status'] = STATUS_REJ_MANUAL
        df = _update_final_status_df(df)
        if mode == self._active_mode:
            self._match_df = df
            # 手動帰属座標を記録
            if not hasattr(self, '_manual_winner_coords') \
                    or self._manual_winner_coords is None:
                self._manual_winner_coords = set()
            self._manual_winner_coords.add(coord)
        else:
            snap['_match_df'] = df
            mw = snap.get('_manual_winner_coords') or set()
            mw.add(coord)
            snap['_manual_winner_coords'] = mw
            self._mode_state[mode] = snap
        # 通知 + UI 更新
        try:
            self.match_df_updated.emit(self.fe.tag, self._match_df)
        except Exception:
            pass
        self._draw_overlay()
        self._rebuild_match_table()
        log.info(f"[Manual] Set {winner_cls} as winner at {mode} "
              f"(RT={coord[0]:.3f}, m/z={coord[1]:.4f})")

    def _apply_manual_row_status(self, mode: str, row_idx, new_status: str):
        """個別行の manual_status を更新する。"""
        if mode == self._active_mode:
            df = self._match_df
        else:
            snap = self._mode_state.get(mode, {}) or {}
            df = snap.get('_match_df')
            if df is None:
                return
        df.loc[row_idx, 'manual_status'] = new_status
        df = _update_final_status_df(df)
        if mode == self._active_mode:
            self._match_df = df
        else:
            snap['_match_df'] = df
            self._mode_state[mode] = snap
        try:
            self.match_df_updated.emit(self.fe.tag, self._match_df)
        except Exception:
            pass
        self._draw_overlay()
        self._rebuild_match_table()
        log.info(f"[Manual] Row {row_idx} → manual_status={new_status} ({mode})")

    def _apply_spot_status_all(self, mode: str, coord: tuple,
                               new_status: str):
        """スポット(同 coord)の全候補行に manual_status を一括適用。
        new_status=STATUS_REJ_MANUAL で全候補を除外(全部誤りのとき)、
        STATUS_NA で手動指定を既定に戻す。manual_status はセッションに
        自動保存される(_serialize_manual_overrides)。"""
        if mode == self._active_mode:
            df = self._match_df
        else:
            snap = self._mode_state.get(mode, {}) or {}
            df = snap.get('_match_df')
            if df is None:
                return
        same = (
            (df['obs_rt'].round(6) == coord[0]) &
            (df['obs_mz'].round(6) == coord[1]) &
            df['matched']
        )
        if not bool(same.any()):
            return
        df.loc[same, 'manual_status'] = new_status
        df = _update_final_status_df(df)
        if mode == self._active_mode:
            self._match_df = df
        else:
            snap['_match_df'] = df
            self._mode_state[mode] = snap
        try:
            self.match_df_updated.emit(self.fe.tag, self._match_df)
        except Exception:
            pass
        self._draw_overlay()
        self._rebuild_match_table()
        _verb = ('rejected' if new_status == STATUS_REJ_MANUAL
                 else 'restored')
        log.info(f"[Manual] Spot {_verb} at {mode} "
              f"RT={coord[0]:.3f} m/z={coord[1]:.4f} "
              f"({int(same.sum())} candidate rows)")

    def _on_hover(self, event):
        # event.inaxes に対応する annot を使い、座標を
        # その軸の data coords で正しく解釈させる(ポップアップ位置ズレ修正)
        # 加えて、active 軸でホバーした場合だけ self._match_df を参照する。
        # 他軸では他モードの match_df を参照(snapshot/_mode_state 経由)。
        valid_axes = []
        if hasattr(self, '_ax_sc_pos') and self._ax_sc_pos is not None:
            valid_axes.append(self._ax_sc_pos)
        if hasattr(self, '_ax_sc_neg') and self._ax_sc_neg is not None:
            valid_axes.append(self._ax_sc_neg)
        if not valid_axes:
            valid_axes = [self._ax_sc] if self._ax_sc is not None else []

        # 先に両 annot を hide(切替時の残像防止)
        for a in (getattr(self, '_annot_pos', None),
                  getattr(self, '_annot_neg', None)):
            if a is not None:
                a.set_visible(False)

        if event.inaxes not in valid_axes:
            self._canvas.draw_idle()
            return

        # event.inaxes に対応する match_df と annot を取得
        if event.inaxes is self._ax_sc_pos:
            target_annot = self._annot_pos
            target_mode = 'pos'
        elif event.inaxes is self._ax_sc_neg:
            target_annot = self._annot_neg
            target_mode = 'neg'
        else:
            self._canvas.draw_idle()
            return

        # 自モードの match_df は self._match_df、他モードは _mode_state から
        if target_mode == self._active_mode:
            target_match_df = self._match_df
            target_conflict_map = self._conflict_map
        else:
            other_state = self._mode_state.get(target_mode, {}) or {}
            target_match_df = other_state.get('_match_df')
            target_conflict_map = other_state.get('_conflict_map') or {}

        if target_match_df is None:
            self._canvas.draw_idle()
            return
        # 描画されているスポット(=kept)のみを対象に最近接探索。
        # rejected 行も含めると、画面に出ていない rejected の情報が popup される
        # ため、ユーザーが見ている spot と内容がずれる。
        matched = target_match_df[target_match_df['matched']]
        if 'final_status' in matched.columns:
            matched = matched[matched['final_status'] == STATUS_KEPT]
        if matched.empty:
            self._canvas.draw_idle()
            return

        xlim = event.inaxes.get_xlim()
        ylim = event.inaxes.get_ylim()
        xr = max(xlim[1] - xlim[0], 1e-6)
        yr = max(ylim[1] - ylim[0], 1e-6)
        dx = (matched['obs_rt'].values - event.xdata) / xr
        dy = (matched['obs_mz'].values - event.ydata) / yr
        best = np.argmin(dx**2 + dy**2)
        if (dx[best]**2 + dy[best]**2) ** 0.5 < 0.03:
            r = matched.iloc[best]
            coord = (round(float(r['obs_rt']), 6),
                     round(float(r['obs_mz']), 6))
            target_annot.xy = (r['obs_rt'], r['obs_mz'])

            text = (f"{r['lipid_class']}  {r['compound']}\n"
                    f"adduct: {r['adduct']}\n"
                    f"obs m/z: {r['obs_mz']:.5f}\n"
                    f"Δppm: {r['delta_ppm']:.2f}\n"
                    f"RT: {r['obs_rt']:.3f} min")

            if coord in target_conflict_map:
                others = [
                    e for e in target_conflict_map[coord]
                    if not (e['lipid_class'] == r['lipid_class']
                            and e['compound'] == r['compound'])
                ]
                if others:
                    text += "\n⚠ Conflicts with:"
                    for e in others:
                        text += (f"\n  • {e['lipid_class']} "
                                 f"{e['compound']} {e['adduct']}")

            target_annot.set_text(text)
            target_annot.set_visible(True)
        self._canvas.draw_idle()

    def _on_table_select(self):
        # 既存ハイライトを削除
        for artist in list(self._overlay_artists):
            if getattr(artist, '_is_highlight', False):
                try: artist.remove()
                except Exception: pass
                self._overlay_artists.remove(artist)

        selected = self._match_table.selectedItems()
        if not selected:
            self._canvas.draw_idle()
            return
        data = selected[0].data(Qt.UserRole)
        if not data:
            return

        # 行の mode に応じて pos / neg どちらの軸で
        # ハイライトするか決定する。
        row_mode = data.get('mode', self._active_mode)
        target_ax = (
            getattr(self, '_ax_sc_pos', None) if row_mode == 'pos'
            else getattr(self, '_ax_sc_neg', None))
        if target_ax is None:
            target_ax = self._ax_sc
        sc = target_ax.scatter(
            [data['obs_rt']], [data['obs_mz']],
            s=120, color='none', edgecolors='#E24B4A',
            linewidths=2.0, zorder=8)
        sc._is_highlight = True
        self._overlay_artists.append(sc)

        self._annot.xy = (data['obs_rt'], data['obs_mz'])
        coord = (round(float(data['obs_rt']), 6),
                 round(float(data['obs_mz']), 6))
        text = (f"[{row_mode}] {data['lipid_class']}  {data['compound']}\n"
                f"adduct: {data['adduct']}\n"
                f"obs m/z: {data['obs_mz']:.5f}\n"
                f"Δppm: {data['delta_ppm']:.2f}\n"
                f"RT: {data['obs_rt']:.3f} min")
        # 両モードの conflict_map を統合チェック
        _, all_cmap = self._combined_conflict_state()
        if coord in all_cmap:
            others = [
                e for e in all_cmap[coord]
                if not (e['lipid_class'] == data['lipid_class']
                        and e['compound'] == data['compound'])
            ]
            if others:
                text += "\n⚠ Conflicts with:"
                for e in others:
                    text += f"\n  • {e['lipid_class']} {e['compound']} {e['adduct']}"
        self._annot.set_text(text)
        self._annot.set_visible(True)
        self._canvas.draw_idle()

    def _on_select(self, xmin, xmax, source_mode: str = ""):
        """pos / neg どちらの軸からの選択も受け付ける。
        他軸の選択ハイライトはクリアする(同時に2 つ赤帯が出るのを避ける)。
        SpanSelector の選択直後に消える矩形ではなく
        永続的な axvspan を描画し、ユーザーが選択範囲を視認できるようにする。
        """
        # 極小幅クリック(誤クリック)は無視
        if abs(xmax - xmin) < 1e-6:
            return
        self._selected_range = (round(xmin, 3), round(xmax, 3))
        self._lbl_range.setText(
            f"Selected RT: {self._selected_range[0]} – "
            f"{self._selected_range[1]}")
        # 他軸の SpanSelector visual + 永続矩形をクリア
        other_mode = 'neg' if source_mode == 'pos' else 'pos'
        other_span_attr = '_span_neg' if source_mode == 'pos' else '_span_pos'
        other_span = getattr(self, other_span_attr, None)
        if other_span is not None:
            try:
                other_span.set_visible(False)
            except Exception:
                pass
        # 既存の永続矩形を全モード分削除
        for m, art in list(self._persistent_range_artists.items()):
            try:
                art.remove()
            except Exception:
                pass
        self._persistent_range_artists.clear()
        # 選択された軸に永続的な axvspan を描画
        src_ax = (self._ax_sc_pos if source_mode == 'pos'
                  else self._ax_sc_neg)
        if src_ax is not None:
            try:
                art = src_ax.axvspan(
                    self._selected_range[0],
                    self._selected_range[1],
                    alpha=0.25, facecolor='red', edgecolor='none',
                    zorder=1)
                self._persistent_range_artists[source_mode] = art
            except Exception as e:
                log.warning(f"[_on_select] axvspan failed: {e}")
        try:
            if self._canvas is not None:
                self._canvas.draw_idle()
        except Exception:
            pass

    def _clear_range_selection(self):
        """scatter plot の RT 範囲選択を消去する。
        誤って範囲を引いた場合などに、_selected_range と両軸の span 可視を
        リセットする。
        永続的な axvspan ハイライトも削除する。
        """
        self._selected_range = None
        try:
            self._lbl_range.setText("Select a range on the scatter plot")
        except Exception:
            pass
        # 両軸の span visual をクリア
        for attr in ('_span_pos', '_span_neg'):
            sp = getattr(self, attr, None)
            if sp is None:
                continue
            try:
                sp.set_visible(False)
            except Exception:
                pass
            # extents もリセット(interactive=False では set_visible だけで
            # 描画は消えるが、内部 state も初期化しておく)
            try:
                sp.extents = (0.0, 0.0)
            except Exception:
                pass
        # 永続的な axvspan を削除
        for m, art in list(self._persistent_range_artists.items()):
            try:
                art.remove()
            except Exception:
                pass
        self._persistent_range_artists.clear()
        try:
            if self._canvas is not None:
                self._canvas.draw_idle()
        except Exception:
            pass

    def _open_add_task(self):
        if self._selected_range is None:
            QMessageBox.warning(self, "No Range",
                                "Please select a range first.")
            return
        dlg = AddTaskPopup(
            self._selected_range[0], self._selected_range[1],
            self._all_entries, parent=self)
        dlg.task_ready.connect(self.task_ready)
        dlg.exec()

    def _open_quant_ion_selector(self):
        """⑦ Select Quant Ion トグル。

        ボタンが checkable のため click 時に Qt が自動で checked を反転する。
        - 新しい state が True(灰→青、つまり「ON する」操作): ダイアログを開く
        - 新しい state が False(青→灰、つまり「OFF する」操作): 確認ダイアログを
          経て _quant_ion_choices をクリア。ユーザーが Cancel なら state を
          青に戻す。
        """
        # Qt の auto-toggle 後の新 state
        try:
            new_state = bool(self._btn_quant_ion.isChecked())
        except Exception:
            new_state = True

        # 青 → 灰: fix16 — 既に選択がある場合はクリアせず、
        #   「選択を保持して再編集 / Clear all / Cancel」の3択にする。
        #   既定は Re-edit(保持)。微調整のたびに全部選び直す手間をなくす。
        if not new_state:
            if not (self._quant_ion_choices or {}):
                # 既に空ならそのまま OFF を確定
                return
            box = QMessageBox(self)
            box.setWindowTitle("Select Quant Ion")
            box.setText("A Quant Ion selection already exists.")
            box.setInformativeText(
                "Re-edit keeps your current choices (tweak only a few).\n"
                "Clear all removes the selection (disables Curated view / "
                "Auto-add Tasks filters).")
            b_reedit = box.addButton(
                "Re-edit (keep selection)", QMessageBox.AcceptRole)
            b_clear = box.addButton("Clear all", QMessageBox.DestructiveRole)
            box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(b_reedit)
            box.exec()
            clicked = box.clickedButton()
            if clicked is b_clear:
                self._quant_ion_choices = {}
                try:
                    self._draw_overlay()
                except Exception:
                    pass
                return  # button は unchecked のまま(OFF 確定)
            elif clicked is b_reedit:
                # 選択を保持したままダイアログを開き直す(下の ON 経路へ)
                try:
                    self._set_btn_checked_silent(self._btn_quant_ion, True)
                except Exception:
                    pass
                new_state = True  # fall through to open dialog
            else:  # Cancel: 何もしない(checked に戻す)
                try:
                    self._set_btn_checked_silent(self._btn_quant_ion, True)
                except Exception:
                    pass
                return

        # 灰 → 青(ON する)経路: 従来通りダイアログを開く
        # 両モードの matched & kept をクラスごとに集計
        per_class: dict = {}   # cls -> {'pos': count, 'neg': count}
        contexts = [(self._active_mode, self._match_df)]
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        contexts.append((other_mode, snap.get('_match_df')))
        for mode, mdf in contexts:
            if mdf is None or mdf.empty:
                continue
            sub = mdf[mdf['matched']]
            if 'final_status' in sub.columns:
                sub = sub[sub['final_status'] == STATUS_KEPT]
            if sub.empty:
                continue
            for cls, grp in sub.groupby('lipid_class'):
                if not str(cls):
                    continue
                per_class.setdefault(str(cls), {'pos': 0, 'neg': 0})
                per_class[str(cls)][mode] = len(grp)
        if not per_class:
            QMessageBox.information(
                self, "No data",
                "No matched & kept lipid classes available.\n"
                "Run the annotation pipeline first.")
            # データなし時は button を OFF に戻す
            try:
                self._set_btn_checked_silent(self._btn_quant_ion, False)
            except Exception:
                pass
            return
        dlg = QuantIonSelectorDialog(
            per_class=per_class,
            initial=self._quant_ion_choices,
            parent=self)
        # クラス名クリック時に散布図をフォーカス
        dlg.classFocusRequested.connect(self._focus_on_class)
        if dlg.exec():
            self._quant_ion_choices = dlg.result_choices()
        # ボタンの checked 状態を _quant_ion_choices と同期
        # (Cancel 時も含めて auto-toggle の挙動を上書き)
        try:
            self._set_btn_checked_silent(
                self._btn_quant_ion,
                bool(self._quant_ion_choices))
        except Exception:
            pass
        # Curated view 表示中なら再描画して反映
        try:
            self._draw_overlay()
        except Exception:
            pass

    def _focus_on_class(self, cls: str):
        """指定クラスのみ散布図に表示する。
        cls が空文字なら全クラス表示に戻す。
        Filters 修飾子(RT outliers / IS-rejected 等)の checked 状態は維持。
        Not matched は focus 時は明示的に uncheck。
        Conflicts は維持(△ で focus クラスの衝突状況を見たいため)。"""
        if not cls:
            # 全クラスを再表示
            self._select_all_classes()
            return
        filter_modifiers = self._filter_modifier_tags()
        # Not matched は focus 時に uncheck する(Conflicts は維持)
        force_off_tags = {self._TAG_NOT_MATCHED}
        for key, cb in self._class_checks.items():
            if key in filter_modifiers:
                continue  # Filters 修飾子は触らない
            if key in force_off_tags:
                target = False
            else:
                target = (key == cls)
            cb.blockSignals(True)
            cb.setChecked(target)
            self._user_filter_prefs[key] = target
            cb.blockSignals(False)
        self._draw_both_overlays()

    def _auto_add_tasks_from_pipeline(self):
        """アノテーションパイプライン後に確定した
        各脂質クラスから Task を自動生成する。

        - 両モード(active + 他モード snapshot)の match_df を走査
        - matched & final_status='kept' の行を class ごとに集約
        - obs_rt の (min, max) ± padding を RT 範囲として task_ready emit
        - 対応する FileEntry が _all_entries にない mode はスキップ
        - ⑦ Select Quant Ion の選択でフィルタ
        """
        # mode → FileEntry.tag のマップを構築
        fe_by_mode: dict = {}
        for fe in (self._all_entries or []):
            m = getattr(fe, 'ion_mode', None)
            if m and m not in fe_by_mode:
                fe_by_mode[m] = fe.tag

        # mode → match_df のマップを構築
        mode_dfs: dict = {}
        if self._match_df is not None and not self._match_df.empty:
            mode_dfs[self._active_mode] = self._match_df
        other_mode = 'neg' if self._active_mode == 'pos' else 'pos'
        snap = self._mode_state.get(other_mode, {}) or {}
        other_df = snap.get('_match_df')
        if other_df is not None and not other_df.empty:
            mode_dfs[other_mode] = other_df

        if not mode_dfs:
            QMessageBox.warning(
                self, "No data",
                "No match_df available. Run Match & Overlay first.")
            return

        PAD = 0.05  # min, 0.05min ≒ 3sec padding
        candidates: list = []   # [(mode, file_tag, cls, rt_start, rt_end, n)]
        for mode, df in mode_dfs.items():
            file_tag = fe_by_mode.get(mode)
            if file_tag is None:
                continue
            sub = df[df['matched']]
            if 'final_status' in sub.columns:
                sub = sub[sub['final_status'] == STATUS_KEPT]
            if sub.empty:
                continue
            for cls, grp in sub.groupby('lipid_class'):
                if not str(cls):
                    continue
                rt_vals = grp['obs_rt'].astype(float).values
                if len(rt_vals) == 0:
                    continue
                rt_min = float(rt_vals.min()) - PAD
                rt_max = float(rt_vals.max()) + PAD
                # 単一スポットなら最低 0.1 min 幅を確保
                if rt_max - rt_min < 0.1:
                    mid = (rt_max + rt_min) / 2.0
                    rt_min = mid - 0.05
                    rt_max = mid + 0.05
                candidates.append(
                    (mode, file_tag, str(cls),
                     round(rt_min, 3), round(rt_max, 3), len(rt_vals)))

        if not candidates:
            QMessageBox.information(
                self, "No tasks to add",
                "No confirmed lipid classes found.\n"
                "Run the annotation pipeline (Match → Conflict → IS Filter → "
                "Adduct Ion Filter → RT Outlier → Coherence Filter) "
                "and ensure at least one row has final_status='kept'.")
            return

        # Quant Ion 選択でフィルタ。
        # _quant_ion_choices[cls] == mode のクラス×モードのみ採用。
        # 'skip' / 未選択 / 異なるモードは除外する。
        quant_choices = getattr(self, '_quant_ion_choices', None) or {}
        if not quant_choices or not any(
                v in ('pos', 'neg') for v in quant_choices.values()):
            QMessageBox.warning(
                self, "No Quant Ion Selection",
                "Please use '⑦ Select Quant Ion…' first to choose\n"
                "pos / neg mode per lipid class before using\n"
                "Auto-add Tasks.")
            return
        filtered = [
            c for c in candidates
            if quant_choices.get(c[2]) == c[0]
        ]
        n_dropped = len(candidates) - len(filtered)
        candidates = filtered
        if not candidates:
            QMessageBox.information(
                self, "No tasks to add",
                "No tasks match the current Quant Ion selection.\n"
                "Open '⑦ Select Quant Ion…' to review your choices.")
            return

        candidates.sort(key=lambda t: (t[0], t[2]))

        # 確認ダイアログ
        n = len(candidates)
        n_pos = sum(1 for c in candidates if c[0] == 'pos')
        n_neg = sum(1 for c in candidates if c[0] == 'neg')
        # 上位プレビュー
        preview_lines = []
        for c in candidates[:12]:
            preview_lines.append(
                f"  [{c[0]}] {c[2]:<8s}  RT {c[3]:.3f} – {c[4]:.3f} "
                f"({c[5]} spots)")
        if len(candidates) > 12:
            preview_lines.append(f"  … and {len(candidates) - 12} more")
        # フィルタで除外された候補数も併記
        msg = (f"About to add {n} task(s) "
               f"(Pos: {n_pos}, Neg: {n_neg}).\n"
               f"(Filtered by ⑦ Select Quant Ion; "
               f"{n_dropped} class/mode combination(s) excluded.)\n\n"
               + "\n".join(preview_lines)
               + "\n\nProceed?")
        ok = QMessageBox.question(
            self, "Auto-add Tasks", msg,
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes)
        if ok != QMessageBox.Yes:
            return

        # task_ready emit
        for mode, file_tag, cls, rt_start, rt_end, _n in candidates:
            self.task_ready.emit(dict(
                cls=cls, rt_start=rt_start, rt_end=rt_end,
                ion=mode, file_tag=file_tag,
            ))

    # _open_set_fix_rt 撤廃
    # (IS Filter の手動 RT 入力で代替可能なため不要)


# ════════════════════════════════════════════════════════════════════
#  AddTaskPopup
# ════════════════════════════════════════════════════════════════════
# ── 極性ごとウィジェットのプロパティを生成する ────────────────
#   self._al_sp_noise → self._al_w[self._al_cur_pol]['sp_noise']
#   こうすることで、_al_* を参照している既存の 392 箇所(44 ウィジェット)
#   に手を入れずに 2 列化できる。対象の列は _al_cur_pol が決める。
def _al_make_widget_property(_name: str):
    def _get(self):
        try:
            return self._al_w[getattr(self, '_al_cur_pol', 'pos')][_name]
        except (AttributeError, KeyError):
            raise AttributeError(f"_al_{_name} is not built yet")

    def _set(self, value):
        if not hasattr(self, '_al_w'):
            self._al_w = {'pos': {}, 'neg': {}}
        self._al_w[getattr(self, '_al_cur_pol', 'pos')][_name] = value

    return property(_get, _set, doc=f"per-polarity widget: {_name}")


for _al_wname in PreviewRTDialog._AL_PER_POL_WIDGETS:
    setattr(PreviewRTDialog, '_al_' + _al_wname,
            _al_make_widget_property(_al_wname))
del _al_wname



class AddTaskPopup(QDialog):
    task_ready = Signal(dict)

    def __init__(self, rt_start, rt_end, file_entries, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Add Task")
        lay = QFormLayout(self)

        self._le_class = QLineEdit()
        self._le_class.setPlaceholderText("e.g. TG, PC, SM …")
        lay.addRow("Lipid Class:", self._le_class)

        self._le_rts = QLineEdit(str(rt_start))
        self._le_rte = QLineEdit(str(rt_end))
        lay.addRow("RT start:", self._le_rts)
        lay.addRow("RT end:",   self._le_rte)

        self._cmb_ion = QComboBox()
        self._cmb_ion.addItems(["pos", "neg"])
        lay.addRow("Ion Mode:", self._cmb_ion)

        self._cmb_file = QComboBox()
        for fe in file_entries:
            self._cmb_file.addItem(fe.label, fe.tag)
        lay.addRow("File:", self._cmb_file)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self._accept)
        bb.rejected.connect(self.reject)
        lay.addRow(bb)

    def _accept(self):
        cls = self._le_class.text().strip()
        if not cls:
            QMessageBox.warning(self, "Input Error",
                                "Lipid class is required.")
            return
        self.task_ready.emit(dict(
            cls=cls,
            rt_start=float(self._le_rts.text()),
            rt_end=float(self._le_rte.text()),
            ion=self._cmb_ion.currentText(),
            file_tag=self._cmb_file.currentData(),
        ))
        self.accept()


# ════════════════════════════════════════════════════════════════════
#  ISPreviewDialog
# ════════════════════════════════════════════════════════════════════

# ════════════════════════════════════════════════════════════════════
#  ManageFilesDialog
# ════════════════════════════════════════════════════════════════════
class ManageFilesDialog(QDialog):
    """Modeless dialog to add / remove / inspect loaded MS-DIAL data files.

    Provides the same operations that previously lived in MainWindow's
    "Input Files" panel: add (multi-file), remove selected, inspect contents.
    The host (MainWindow) supplies the actual list of FileEntry objects via
    callbacks; this dialog only acts as a UI surface.

    Signals (the host wires them up):
      filesChanged  — emitted when the underlying file list changed
    """

    filesChanged = Signal()

    def __init__(self, host, parent=None):
        super().__init__(parent or host)
        self.setWindowTitle("Manage Data Files")
        self.resize(680, 360)
        self._host = host

        lay = QVBoxLayout(self)
        intro = QLabel(
            "Manage MS-DIAL data files. Add new files, inspect existing "
            "ones, or remove unused entries.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color:#555; padding:2px;")
        lay.addWidget(intro)

        # Input mode radio (v1.1.0-alpha.6)
        # Alignment result = 既存 (1 file with cross-sample columns)
        # Per-sample peak list = 新規 (N files, simple alignment in LipidZoner)
        mode_box = QGroupBox("Input Mode")
        mode_lay = QHBoxLayout(mode_box)
        self._rb_alignment_mode = QRadioButton("Alignment result file")
        self._rb_alignment_mode.setToolTip(
            "Existing mode: MS-DIAL Alignment result export (1 TXT per ion mode)")
        self._rb_per_sample_mode = QRadioButton("Per-sample peak list folder")
        self._rb_per_sample_mode.setToolTip(
            "New mode (v1.1.0): MS-DIAL Peak list result export "
            "(N TXT per ion mode, isotopes preserved)")
        # 既定を Per-sample peak list folder にする
        self._rb_per_sample_mode.setChecked(True)
        self._mode_group = QButtonGroup(self)
        self._mode_group.addButton(self._rb_alignment_mode)
        self._mode_group.addButton(self._rb_per_sample_mode)
        # Per-sample を先頭に表示
        mode_lay.addWidget(self._rb_per_sample_mode)
        mode_lay.addWidget(self._rb_alignment_mode)
        mode_lay.addStretch()
        lay.addWidget(mode_box)

        # File list table
        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(["Tag", "Ion mode", "Path"])
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectRows)
        self._table.setSelectionMode(QTableWidget.SingleSelection)
        lay.addWidget(self._table, stretch=1)

        # Action buttons
        btn_row = QHBoxLayout()
        btn_add = QPushButton("Add file…")
        btn_add.clicked.connect(self._on_add)
        btn_remove = QPushButton("Remove selected")
        btn_remove.clicked.connect(self._on_remove)
        btn_inspect = QPushButton("Inspect…")
        btn_inspect.clicked.connect(self._on_inspect)
        btn_row.addWidget(btn_add)
        btn_row.addWidget(btn_remove)
        btn_row.addWidget(btn_inspect)
        btn_row.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btn_row.addWidget(btn_close)
        lay.addLayout(btn_row)

        self.refresh()

    def refresh(self):
        """Reload the table from the host's _file_entries list."""
        entries = list(getattr(self._host, '_file_entries', []) or [])
        self._table.setRowCount(len(entries))
        for r, fe in enumerate(entries):
            self._table.setItem(r, 0, QTableWidgetItem(str(fe.tag)))
            self._table.setItem(r, 1, QTableWidgetItem(str(fe.ion_mode)))
            self._table.setItem(r, 2, QTableWidgetItem(str(fe.path)))
        self._table.resizeColumnsToContents()

    def _selected_tag(self):
        rows = self._table.selectionModel().selectedRows()
        if not rows:
            return None
        item = self._table.item(rows[0].row(), 0)
        return item.text() if item else None

    def _on_add(self):
        """Delegate to host's add method based on input mode selection.

        - Alignment result mode → host._add_file (file picker)
        - Per-sample peak list mode → host._add_per_sample_folder (folder picker)
        """
        host = self._host
        if self._rb_per_sample_mode.isChecked():
            # Per-sample peak list mode (v1.1.0-alpha.6+)
            if hasattr(host, '_add_per_sample_folder'):
                try:
                    host._add_per_sample_folder()
                except Exception as e:
                    QMessageBox.warning(self, "Add failed", str(e))
            else:
                QMessageBox.warning(self, "Not implemented",
                    "Per-sample peak list loading is not available in this build.")
        else:
            # Alignment result mode (existing)
            if hasattr(host, '_add_file'):
                try:
                    host._add_file()
                except Exception as e:
                    QMessageBox.warning(self, "Add failed", str(e))
        self.refresh()
        self.filesChanged.emit()

    def _on_remove(self):
        """Delegate to host's _remove_file using currently-selected tag."""
        tag = self._selected_tag()
        if not tag:
            QMessageBox.information(
                self, "No selection",
                "Select a file row first, then click Remove selected.")
            return
        host = self._host
        # Set MainWindow's combobox to the matching tag, then call _remove_file
        try:
            cb = getattr(host, '_cmb_files', None)
            if cb is not None:
                idx = cb.findData(tag)
                if idx >= 0:
                    cb.setCurrentIndex(idx)
            if hasattr(host, '_remove_file'):
                host._remove_file()
        except Exception as e:
            QMessageBox.warning(self, "Remove failed", str(e))
        self.refresh()
        self.filesChanged.emit()

    def _on_inspect(self):
        """Delegate to host's _inspect_file (opens InspectDialog)."""
        tag = self._selected_tag()
        if not tag:
            QMessageBox.information(
                self, "No selection",
                "Select a file row first, then click Inspect.")
            return
        host = self._host
        try:
            cb = getattr(host, '_cmb_files', None)
            if cb is not None:
                idx = cb.findData(tag)
                if idx >= 0:
                    cb.setCurrentIndex(idx)
            # host のメソッド名は _inspect_file ではなく _inspect
            if hasattr(host, '_inspect'):
                host._inspect()
        except Exception as e:
            QMessageBox.warning(self, "Inspect failed", str(e))


class ISPreviewDialog(QDialog):
    """
    IS フィルタ適用前の確認ダイアログ。

    ライブラリに登録されている全 IS エントリについて:
      - Use チェックボックス（採否）
      - Manual RT（空欄なら自動、数値入力で上書き）
      - Status（OK / Warning / Not detected）
    を表示し、ユーザーが手動で選択・編集できる。

    判定:
      - Not detected: ppm tol内にマッチなし → Use 自動OFF
      - Warning: 同クラス内の複数ISで obs_rt の標準偏差 > 0.1 min
      - OK: 上記以外
    """

    WARNING_RT_STD = 0.10  # クラス内IS間のRT標準偏差の警告閾値 [min]

    def __init__(
        self,
        is_match_df: pd.DataFrame,
        prev_choices: dict | None = None,
        sample_col_names: list[str] | None = None,
        ion_mode: str | None = None,
        raw_paths: "list[str] | None" = None,
        parent=None,
    ):
        """
        is_match_df: ISエントリのみを仮マッチングした結果
          必須列: lipid_class, compound, adduct, theoretical_mz,
                  obs_mz, obs_rt, delta_ppm, matched
        prev_choices: 前回のセッションでの選択状態（任意）
          形式: {cls: [{"compound", "adduct", "use", "manual_rt",
                        "selected_peak_idx"}, ...]}
        sample_col_names: サンプル列名(CandidatePickerDialog 用)
        ion_mode: 'pos' / 'neg' — タイトルバーに表示
        """
        super().__init__(parent)
        # pos/neg を一目で区別できるようにタイトルに含める
        title_suffix = (
            f"  -  {str(ion_mode).upper()} mode"
            if ion_mode else "")
        self.setWindowTitle(f"IS Filter Preview{title_suffix}")
        self.resize(1020, 520)
        self._result: dict = {}
        self._sample_col_names = list(sample_col_names) if sample_col_names else []
        # Candidate Picker の EIC 用に生データのパスを持ち回る
        self._raw_paths = [str(p) for p in (raw_paths or [])]
        # 行ごとの候補リストを保持(Candidates セルクリック時に参照)
        # row_idx → (theoretical_mz, candidates_list, lipid_class, compound, adduct)
        self._row_candidates: dict[int, dict] = {}

        prev_choices = prev_choices or {}

        lay = QVBoxLayout(self)
        info = QLabel(
            "Review the detected internal standards (IS) before applying "
            "the IS filter.\n"
            "Uncheck rows to exclude an IS from filtering. "
            "Enter a Manual RT to override the detected value.\n"
            "Classes with at least one checked IS will use IS-based "
            "two-pass matching."
        )
        info.setStyleSheet("color:#555; font-size:11px;")
        lay.addWidget(info)

        # ─── テーブル構築 ─────────────────────────────────────────
        # Use / Class / Compound / Adduct / Theor m/z / Obs m/z /
        # Δppm / Obs RT / Manual RT / RT tol / Status / Candidates
        # 「RT tol」を Manual RT の右に追加。クラスごとに
        # ③ の RT 窓幅(既定 ±0.1 min)を変えられるようにする。
        # HexCer のように 1 クラス内で O2 系と O3 系が 0.26 min 離れて
        # 溶出する場合、既定の窓では片方が丸ごと落ちるため。
        cols = ["Use", "Class", "IS Compound", "Adduct", "Theor. m/z",
                "Obs. m/z", "Δppm", "Obs. RT", "Manual RT", "RT tol",
                "Status", "Candidates"]
        self._table = QTableWidget(0, len(cols))
        self._table.setHorizontalHeaderLabels(cols)
        self._table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        hdr = self._table.horizontalHeader()
        hdr.setSectionResizeMode(2, QHeaderView.Stretch)  # Compound列を伸ばす
        # Candidates セルクリックで CandidatePickerDialog を開く
        self._table.cellClicked.connect(self._on_cell_clicked)
        lay.addWidget(self._table)

        # クラス内 obs_rt の標準偏差を計算（Warning判定用）
        class_rt_std: dict[str, float] = {}
        matched_only = is_match_df[is_match_df['matched']]
        for cls, grp in matched_only.groupby('lipid_class'):
            if len(grp) >= 2:
                class_rt_std[cls] = float(grp['obs_rt'].std())
            else:
                class_rt_std[cls] = 0.0

        # ─── 行を追加 ─────────────────────────────────────────────
        # prev_choices を (cls, compound, adduct) でルックアップできる辞書に
        prev_lookup: dict[tuple, dict] = {}
        for cls, rows in prev_choices.items():
            for row in rows:
                key = (cls, row.get("compound", ""), row.get("adduct", ""))
                prev_lookup[key] = row

        # IS全件（matched/unmatchedとも）を表示
        for _, r in is_match_df.sort_values(
                ['lipid_class', 'compound']).iterrows():
            row_idx = self._table.rowCount()
            self._table.insertRow(row_idx)

            matched = bool(r['matched'])
            cls = r['lipid_class']
            cmp = r['compound']
            add = r['adduct']

            # Status 判定
            if not matched:
                status = "Not detected"
                status_color = QColor("#E26666")
                default_use = False
            elif class_rt_std.get(cls, 0.0) > self.WARNING_RT_STD:
                status = "Warning: RT scatter"
                status_color = QColor("#D8A847")
                default_use = True
            else:
                status = "OK"
                status_color = QColor("#4AAE5A")
                default_use = True

            # 前回の選択を復元（あれば）
            prev = prev_lookup.get((cls, cmp, add))
            use = prev["use"] if prev else default_use
            manual_rt = prev.get("manual_rt") if prev else None
            # クラス別 RT tol(前回値)。クラス単位なので、同じ
            # クラスの他の行に入っていればそれを引き継ぐ。
            rt_tol_prev = prev.get("rt_tol") if prev else None
            if not rt_tol_prev:
                for _r in (prev_choices.get(cls) or []):
                    if _r.get("rt_tol"):
                        rt_tol_prev = _r["rt_tol"]
                        break
            # 候補リスト
            cand_list = r.get('candidates', []) if hasattr(r, 'get') else []
            if not isinstance(cand_list, list):
                cand_list = list(cand_list) if cand_list is not None else []
            n_candidates = len(cand_list)
            # 複数候補 (>=2) ならステータスを "Multiple candidates" に上書き
            if matched and n_candidates >= 2:
                status = "Multiple candidates"
                status_color = QColor("#E68A00")  # オレンジ系
                # default_use はそのまま ON(原則0、ユーザー判断尊重)

            # Use チェックボックス
            use_cb = QCheckBox()
            use_cb.setChecked(bool(use))
            # Not detected の場合はチェック不可
            if not matched:
                use_cb.setEnabled(False)
                use_cb.setChecked(False)
            use_widget = QWidget()
            use_lay = QHBoxLayout(use_widget)
            use_lay.setContentsMargins(0, 0, 0, 0)
            use_lay.addWidget(use_cb)
            use_lay.setAlignment(Qt.AlignCenter)
            self._table.setCellWidget(row_idx, 0, use_widget)

            # その他のセル
            for c, val in enumerate([
                cls, cmp, add,
                f"{r['theoretical_mz']:.5f}",
                f"{r['obs_mz']:.5f}" if matched else "",
                f"{r['delta_ppm']:.2f}" if matched else "",
                f"{r['obs_rt']:.3f}" if matched else "",
            ], start=1):
                item = QTableWidgetItem(str(val))
                self._table.setItem(row_idx, c, item)

            # Manual RT スピンボックス（編集可）
            manual_spin = QDoubleSpinBox()
            manual_spin.setRange(0.0, 999.0)
            manual_spin.setDecimals(3)
            manual_spin.setSingleStep(0.05)
            manual_spin.setSuffix(" min")
            manual_spin.setSpecialValueText("(auto)")  # 値0で"(auto)"表示
            manual_spin.setValue(
                float(manual_rt) if manual_rt is not None else 0.0)
            self._table.setCellWidget(row_idx, 8, manual_spin)

            # RT tol スピンボックス(クラス単位)
            tol_spin = QDoubleSpinBox()
            tol_spin.setRange(0.0, 5.0)
            tol_spin.setDecimals(3)
            tol_spin.setSingleStep(0.01)
            tol_spin.setSuffix(" min")
            tol_spin.setSpecialValueText("(default)")   # 0 なら既定値を使う
            tol_spin.setValue(float(rt_tol_prev) if rt_tol_prev else 0.0)
            tol_spin.setToolTip(
                "Half width of the RT window ③ uses for this class.\n"
                "(default) uses the IS RT tolerance from Advanced.\n"
                "Every row of the same class shares one value.\n\n"
                "Widen it when one class elutes in separate places (e.g.\n"
                "the O2 and O3 series of HexCer are 0.26 min apart).")
            tol_spin.setProperty("lipid_class", cls)
            tol_spin.valueChanged.connect(self._on_rt_tol_changed)
            self._table.setCellWidget(row_idx, 9, tol_spin)

            # Status セル（色付き）
            status_item = QTableWidgetItem(status)
            status_item.setForeground(status_color)
            self._table.setItem(row_idx, 10, status_item)

            # Candidates セル: 候補数を表示、>=2 ならクリック可能
            if n_candidates == 0:
                cand_text = "-"
                cand_bg = QColor("#FCE6E6")
            elif n_candidates == 1:
                cand_text = "1"
                cand_bg = None
            else:
                cand_text = f"{n_candidates} ▼"
                cand_bg = QColor("#FFE0B2")
            cand_item = QTableWidgetItem(cand_text)
            cand_item.setTextAlignment(Qt.AlignCenter)
            if cand_bg:
                cand_item.setBackground(QBrush(cand_bg))
            if n_candidates >= 2:
                cand_item.setToolTip(
                    f"{n_candidates} candidate peaks within ppm tolerance.\n"
                    "Click to open Candidate Picker.")
            self._table.setItem(row_idx, 11, cand_item)

            # 行データを保持(後で取得 + Candidate Picker 起動用)
            use_cb.setProperty("row_data", {
                "lipid_class": cls, "compound": cmp, "adduct": add})
            # 候補情報を行 index で保持
            self._row_candidates[row_idx] = {
                "theoretical_mz": float(r['theoretical_mz']),
                "candidates":     cand_list,
                "lipid_class":    cls,
                "compound":       cmp,
                "adduct":         add,
                "n_candidates":   n_candidates,
                "selected_peak_idx": (prev.get("selected_peak_idx")
                                      if prev else None),
            }

        # 列幅調整
        self._table.resizeColumnsToContents()
        self._table.setColumnWidth(0, 40)    # Use
        self._table.setColumnWidth(8, 100)   # Manual RT
        self._table.setColumnWidth(9, 100)   # RT tol
        self._table.setColumnWidth(11, 100)  # Candidates

        # ─── ボタン ──────────────────────────────────────────────
        bb = QHBoxLayout()
        bb.addStretch(1)
        btn_cancel = QPushButton("Cancel")
        btn_cancel.clicked.connect(self.reject)
        bb.addWidget(btn_cancel)

        btn_apply = QPushButton("Apply")
        btn_apply.setDefault(True)
        btn_apply.clicked.connect(self._on_apply)
        bb.addWidget(btn_apply)
        lay.addLayout(bb)

    def _on_rt_tol_changed(self, value: float):
        """RT tol はクラス単位なので、同じクラスの行を全部そろえる。

        表は IS 化合物ごとの行だが、③ の RT 窓はクラスに 1 つしかない。
        行ごとにバラバラの値を入れられると、どれが効いているのか
        画面から分からなくなるため、編集した瞬間に同期する。
        """
        if getattr(self, '_syncing_rt_tol', False):
            return
        src = self.sender()
        if src is None:
            return
        cls = src.property("lipid_class")
        if not cls:
            return
        self._syncing_rt_tol = True
        try:
            for row in range(self._table.rowCount()):
                wdg = self._table.cellWidget(row, 9)
                if wdg is None or wdg is src:
                    continue
                if wdg.property("lipid_class") == cls and wdg.value() != value:
                    wdg.setValue(value)
        finally:
            self._syncing_rt_tol = False

    def _on_cell_clicked(self, row: int, col: int):
        """Candidates 列(col=11)がクリックされたら CandidatePickerDialog を開く。

        RT tol 列を挿したので 10 → 11 にずれた。
        """
        if col != 11:
            return
        info = self._row_candidates.get(row)
        if info is None or info["n_candidates"] < 2:
            return  # 候補が 1 以下ならクリック無効

        title = (f"{info['lipid_class']} {info['compound']} "
                 f"candidates ({info['n_candidates']})")
        # 既存の Manual RT 値があれば、対応する候補を初期選択
        manual_spin = self._table.cellWidget(row, 8)
        cur_rt = manual_spin.value() if manual_spin else 0.0
        initial_pidx = info.get("selected_peak_idx")
        if initial_pidx is None and cur_rt > 0:
            # Manual RT に既に値があれば、それに最も近い候補を初期選択
            best = None; best_d = float('inf')
            for c in info["candidates"]:
                d = abs(c['obs_rt'] - cur_rt)
                if d < best_d:
                    best, best_d = c, d
            if best is not None:
                initial_pidx = best['peak_idx']

        dlg = CandidatePickerDialog(
            title=title,
            theoretical_mz=info["theoretical_mz"],
            candidates=info["candidates"],
            sample_col_names=self._sample_col_names,
            initial_peak_idx=initial_pidx,
            raw_paths=self._raw_paths,
            parent=self,
        )
        if dlg.exec():
            sel = dlg.selected_candidate()
            if sel is not None:
                # 選択候補の RT を Manual RT に書き込む
                manual_spin.setValue(float(sel['obs_rt']))
                # 選択 peak_idx を保持(Apply 時に結果に含める)
                info["selected_peak_idx"] = sel['peak_idx']

    def _on_apply(self):
        """Use/Manual RT の選択結果を収集して accept"""
        self._result = {}
        for row in range(self._table.rowCount()):
            use_widget = self._table.cellWidget(row, 0)
            use_cb = use_widget.findChild(QCheckBox)
            row_data = use_cb.property("row_data")
            cls = row_data["lipid_class"]
            cmp = row_data["compound"]
            add = row_data["adduct"]

            manual_spin = self._table.cellWidget(row, 8)
            manual_val = manual_spin.value()
            manual_rt = manual_val if manual_val > 0 else None

            # クラス別 RT tol。0 は「既定値を使う」の意味なので None。
            tol_spin = self._table.cellWidget(row, 9)
            tol_val = tol_spin.value() if tol_spin is not None else 0.0
            rt_tol = float(tol_val) if tol_val > 0 else None

            info = self._row_candidates.get(row, {})
            selected_peak_idx = info.get("selected_peak_idx")

            self._result.setdefault(cls, []).append({
                "compound":          cmp,
                "adduct":            add,
                "use":               use_cb.isChecked(),
                "manual_rt":         manual_rt,
                "rt_tol":            rt_tol,
                "selected_peak_idx": selected_peak_idx,
            })
        self.accept()

    def has_unconfirmed_multi_candidates(self) -> tuple[bool, list[str]]:
        """複数候補がある IS のうち、ユーザーが Manual RT で
        明示的に選択していないものをリストアップする。
        Apply 押下時の soft block 警告に使う。

        戻り値: (warning_needed, [説明文字列, ...])
        """
        unconfirmed = []
        for row, info in self._row_candidates.items():
            if info.get("n_candidates", 0) < 2:
                continue
            # selected_peak_idx で確認済みなら OK
            if info.get("selected_peak_idx") is not None:
                continue
            # Manual RT に値が入っているなら確認済みとみなす
            manual_spin = self._table.cellWidget(row, 8)
            manual_val = manual_spin.value() if manual_spin else 0.0
            if manual_val > 0:
                continue
            unconfirmed.append(
                f"{info['lipid_class']} {info['compound']} "
                f"({info['n_candidates']} candidates)")
        return (len(unconfirmed) > 0, unconfirmed)

    def result_choices(self) -> dict:
        """{cls: [{compound, adduct, use, manual_rt}, ...]} を返す"""
        return self._result


# ════════════════════════════════════════════════════════════════════
#  Coherence Engine 関連ダイアログ
# ════════════════════════════════════════════════════════════════════


class CoherencePairsDialog(QDialog):
    """衝突クラスペアの ON/OFF を選択するダイアログ。

    ライブラリから自動検出された全衝突ペアをチェックボックスで一覧表示する。
    デフォルトは全 ON。ユーザーは特定ペアの帰属を無効化できる(無効ペアの
    衝突は △ 未解決のまま残る)。
    """

    def __init__(
        self,
        all_pairs: list[tuple[str, str]],
        enabled_pairs: list[tuple[str, str]],
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Coherence Pairs")
        self.resize(360, 400)
        self._all_pairs = list(all_pairs)
        enabled_set = {tuple(sorted(p)) for p in enabled_pairs}

        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(
            "Enable/disable per-pair conflict attribution.\n"
            "Disabled pairs remain as unresolved conflicts (△)."))

        self._list = QListWidget()
        self._list.setSelectionMode(QListWidget.NoSelection)
        for pair in self._all_pairs:
            key = tuple(sorted(pair))
            text = f"{key[0]} ↔ {key[1]}"
            item = QListWidgetItem(text)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(
                Qt.Checked if key in enabled_set else Qt.Unchecked)
            item.setData(Qt.UserRole, key)
            self._list.addItem(item)
        lay.addWidget(self._list, stretch=1)

        # 一括 ON/OFF
        tool_row = QHBoxLayout()
        btn_all_on = QPushButton("Select All")
        btn_all_off = QPushButton("Clear All")
        btn_all_on.clicked.connect(lambda: self._set_all(Qt.Checked))
        btn_all_off.clicked.connect(lambda: self._set_all(Qt.Unchecked))
        tool_row.addWidget(btn_all_on)
        tool_row.addWidget(btn_all_off)
        tool_row.addStretch()
        lay.addLayout(tool_row)

        # OK / Cancel
        btns = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    def _set_all(self, state):
        for i in range(self._list.count()):
            self._list.item(i).setCheckState(state)

    def selected_pairs(self) -> list[tuple[str, str]]:
        """OK 押下後に呼び出して ON 状態のペアリストを取得。"""
        result = []
        for i in range(self._list.count()):
            item = self._list.item(i)
            if item.checkState() == Qt.Checked:
                result.append(item.data(Qt.UserRole))
        return result


class CoherenceReportDialog(QDialog):
    """Coherence 帰属の詳細レポートダイアログ。

    全帰属を信頼度降順(低信頼度が上)でテーブル表示する。
    行クリックで parent の散布図へジャンプするための spotSelected シグナルを発火。
    """

    spotSelected = Signal(float, float, str)  # (obs_rt, obs_mz, mode)

    def __init__(
        self,
        assignments: dict,
        sigma_threshold: float,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Coherence Attribution Report")
        # 詳細パネルを下部に配置するため縦サイズ拡大
        self.resize(1100, 780)
        self._assignments = assignments
        self._sigma_threshold = sigma_threshold

        lay = QVBoxLayout(self)

        # 概要ヘッダ
        n = len(assignments)
        n_low = sum(1 for a in assignments.values() if a.get('low_confidence'))
        # ⑤/⑥/△ 由来別カウント
        n_adduct = sum(
            1 for a in assignments.values() if a.get('method') == 'adduct')
        n_coherence = sum(
            1 for a in assignments.values() if a.get('method') == 'coherence')
        n_undecided = sum(
            1 for a in assignments.values() if a.get('method') == 'undecided')
        method_summary = ""
        if n_adduct or n_coherence or n_undecided:
            method_summary = (
                f"    Adduct (⑤): {n_adduct}    "
                f"Coherence (⑥): {n_coherence}    "
                f"Undecided (△): {n_undecided}")
        summary = (
            f"Total attributions: {n}    "
            f"Low-confidence (σ<{sigma_threshold:.1f}): {n_low}    "
            f"High-confidence: {n - n_low}"
            + method_summary)
        lbl = QLabel(summary)
        lbl.setStyleSheet("font-weight: bold; padding: 6px;")
        lay.addWidget(lbl)

        # テーブル
        headers = [
            "Mode", "Method", "σ", "obs RT", "obs m/z",
            "Winner class", "Winner compound",
            "Loser class", "Loser compound",
            "pred winner", "pred loser",
            "res winner", "res loser",
        ]
        self._tbl = QTableWidget()
        self._tbl.setColumnCount(len(headers))
        self._tbl.setHorizontalHeaderLabels(headers)
        self._tbl.setEditTriggers(QTableWidget.NoEditTriggers)
        self._tbl.setSelectionBehavior(QTableWidget.SelectRows)
        self._tbl.setAlternatingRowColors(True)

        # 信頼度昇順(低信頼度が上)
        rows = []
        for coord, a in assignments.items():
            rows.append((coord, a))
        rows.sort(key=lambda x: x[1].get('confidence', 0.0))

        # NaN を「-」として安全に表示するヘルパ
        def _fmt_float(v, fmt):
            try:
                fv = float(v)
                if not np.isfinite(fv):
                    return "-"
                return f"{fv:{fmt}}"
            except (TypeError, ValueError):
                return "-"

        self._tbl.setRowCount(len(rows))
        self._row_keys = []   # 行 index -> key (mode, rt, mz) or (rt, mz)
        for r, (key, a) in enumerate(rows):
            self._row_keys.append(key)
            # key が 3-tuple (mode, rt, mz) の場合の対応
            if isinstance(key, tuple) and len(key) == 3:
                mode_str = str(key[0])
                rt_v, mz_v = float(key[1]), float(key[2])
            else:
                mode_str = a.get('mode', '')
                rt_v, mz_v = float(key[0]), float(key[1])
            method = a.get('method', 'coherence')
            # 'undecided' を Method 列に明示表示
            if method == 'adduct':
                method_label = '⑤ Adduct'
            elif method == 'coherence':
                method_label = '⑥ Coherence'
            elif method == 'undecided':
                method_label = '△ Undecided'
            else:
                method_label = method
            cells = [
                mode_str,
                method_label,
                _fmt_float(a.get('confidence', 0.0), '.2f'),
                f"{rt_v:.3f}",
                f"{mz_v:.4f}",
                a.get('winner_class', ''),
                a.get('winner_compound', ''),
                a.get('loser_class', ''),
                a.get('loser_compound', ''),
                _fmt_float(a.get('pred_winner'), '.3f'),
                _fmt_float(a.get('pred_loser'), '.3f'),
                _fmt_float(a.get('res_winner'), '+.3f'),
                _fmt_float(a.get('res_loser'), '+.3f'),
            ]
            for c, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                if a.get('low_confidence'):
                    # 低信頼度行を濃赤でハイライト
                    item.setForeground(QBrush(QColor('#B22222')))
                self._tbl.setItem(r, c, item)
        self._tbl.resizeColumnsToContents()
        self._tbl.cellClicked.connect(self._on_cell_clicked)
        lay.addWidget(self._tbl, stretch=1)

        # 行選択時に計算詳細を表示する panel
        self._detail_label = QLabel("Click a row to display calculation details below.")
        self._detail_label.setStyleSheet(
            "font-size: 10px; color: #555; padding: 4px;")
        lay.addWidget(self._detail_label)
        self._detail_text = QTextEdit()
        self._detail_text.setReadOnly(True)
        self._detail_text.setFixedHeight(220)
        self._detail_text.setStyleSheet(
            "QTextEdit { background:#FAFAFA; color:#202020; "
            "font-family:'Consolas','Courier New',monospace; font-size:11px; "
            "border:1px solid #CCC; padding:4px; }")
        lay.addWidget(self._detail_text)

        # 閉じるボタン
        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        btns.button(QDialogButtonBox.Close).clicked.connect(self.accept)
        lay.addWidget(btns)

    def _on_cell_clicked(self, row: int, col: int):
        if 0 <= row < len(self._row_keys):
            key = self._row_keys[row]
            # key は (mode, rt, mz) または (rt, mz)
            if isinstance(key, tuple) and len(key) == 3:
                mode_v = str(key[0])
                rt_v, mz_v = float(key[1]), float(key[2])
            else:
                rt_v, mz_v = float(key[0]), float(key[1])
                # entry から mode を取得(なければ active 想定で空)
                a = self._assignments.get(key, {}) or {}
                mode_v = str(a.get('mode', ''))
            # mode 情報も emit してハイライト軸を判別可能に
            self.spotSelected.emit(rt_v, mz_v, mode_v)
            # 詳細パネルを更新
            self._update_detail_panel(key)

    def _update_detail_panel(self, key):
        """選択 spot の計算詳細を整形表示する。
        key は (mode, rt, mz) または (rt, mz) を受け付ける。"""
        a = self._assignments.get(key)
        if not a:
            self._detail_text.setPlainText("(no detail)")
            return
        # 表示用の coord (rt, mz) を抽出
        if isinstance(key, tuple) and len(key) == 3:
            coord = (float(key[1]), float(key[2]))
        else:
            coord = (float(key[0]), float(key[1]))
        method = a.get('method', 'coherence')
        method_label = {
            'adduct': '⑤ Adduct Ion Filter',
            'coherence': '⑥ Coherence Filter',
            'undecided': '△ Undecided',
        }.get(method, method)

        lines = []
        lines.append(f"━━ Spot: RT = {coord[0]:.3f} min,  m/z = {coord[1]:.4f} ━━")
        lines.append(f"  Method:  {method_label}")
        lines.append(f"  σ (confidence): {a.get('confidence', 0.0):.3f}")
        lines.append("")

        # Winner block
        wc = a.get('winner_class', '')
        wcomp = a.get('winner_compound', '')
        wadu = a.get('winner_adduct', '')
        lines.append(f"[Winner]  {wc}  {wcomp}  {wadu}")
        if 'winner_C' in a:
            wC = a['winner_C']; wU = a['winner_U']
            wA = a['winner_alpha']; wB = a['winner_beta']
            wG = a['winner_gamma']; wRMSE = a['winner_rmse']
            wMode = a.get('winner_fit_mode', '?')
            lines.append(f"  link={a.get('winner_link','?')}  C={wC}  U={wU}  "
                         f"(fit_mode: {wMode})")
            lines.append(f"  α = {wA:+.4f}  β = {wB:+.4f}  "
                         f"γ = {wG:.4f}  rmse = {wRMSE:.4f}")
            pw = a.get('pred_winner')
            rw = a.get('res_winner')
            if pw is not None and isinstance(pw, (int, float)) and np.isfinite(pw):
                lines.append(f"  pred = {wA:+.4f} × {wC} + {wB:+.4f} × {wU} "
                             f"+ {wG:.4f} = {pw:.4f} min")
                if rw is not None and np.isfinite(rw):
                    lines.append(f"  res  = obs RT − pred = "
                                 f"{coord[0]:.3f} − {pw:.4f} = {rw:+.4f} min")
        lines.append("")

        # Loser block
        lc = a.get('loser_class', '')
        lcomp = a.get('loser_compound', '')
        ladu = a.get('loser_adduct', '')
        lines.append(f"[Loser]   {lc}  {lcomp}  {ladu}")
        if 'loser_C' in a:
            lC = a['loser_C']; lU = a['loser_U']
            lA = a['loser_alpha']; lB = a['loser_beta']
            lG = a['loser_gamma']; lRMSE = a['loser_rmse']
            lMode = a.get('loser_fit_mode', '?')
            lines.append(f"  link={a.get('loser_link','?')}  C={lC}  U={lU}  "
                         f"(fit_mode: {lMode})")
            lines.append(f"  α = {lA:+.4f}  β = {lB:+.4f}  "
                         f"γ = {lG:.4f}  rmse = {lRMSE:.4f}")
            pl = a.get('pred_loser')
            rl = a.get('res_loser')
            if pl is not None and isinstance(pl, (int, float)) and np.isfinite(pl):
                lines.append(f"  pred = {lA:+.4f} × {lC} + {lB:+.4f} × {lU} "
                             f"+ {lG:.4f} = {pl:.4f} min")
                if rl is not None and np.isfinite(rl):
                    lines.append(f"  res  = obs RT − pred = "
                                 f"{coord[0]:.3f} − {pl:.4f} = {rl+0:+.4f} min")

        # Confidence calc
        rmse_ref = a.get('rmse_ref')
        if rmse_ref is not None:
            rw = a.get('res_winner', float('nan'))
            rl = a.get('res_loser', float('nan'))
            try:
                if np.isfinite(rw) and np.isfinite(rl):
                    margin = abs(rl) - abs(rw)
                    lines.append("")
                    lines.append(f"[σ calculation]")
                    lines.append(f"  margin = |res_loser| − |res_winner| = "
                                 f"|{rl:+.4f}| − |{rw:+.4f}| = {margin:+.4f}")
                    lines.append(f"  rmse_ref = (rmse_w + rmse_l) / 2 "
                                 f"= {rmse_ref:.4f}")
                    sigma = a.get('confidence', 0.0)
                    lines.append(f"  σ = margin / rmse_ref = "
                                 f"{margin:+.4f} / {rmse_ref:.4f} = "
                                 f"{sigma:.3f}")
            except Exception:
                pass

        # Adduct scores (⑤ method only)
        if method == 'adduct':
            scores = a.get('adduct_scores', {})
            if scores:
                lines.append("")
                lines.append(f"[Adduct evidence (⑤)]")
                for cls, sc in sorted(scores.items(), key=lambda x: -x[1]):
                    lines.append(f"  {cls}: score = {sc:.2f}")

        self._detail_text.setPlainText("\n".join(lines))


# ════════════════════════════════════════════════════════════════════
#  AdductHeatmapWidget
# ════════════════════════════════════════════════════════════════════

class AdductHeatmapWidget(QTableWidget):
    """クラス × アダクトのヒートマップ表示 PySide6 ウィジェット。

    旧 adduct_heatmap_widget.py(外部モジュール)を
    LipidZoner.py に inline 化。.exe 配布時の依存が減り自己完結化される。

    Args:
        adduct_columns: アダクト名のリスト/タプル(列の順序)。
                        指定なしなら DEFAULT_ADDUCT_COLUMNS を使用。
        editable     : True なら expected フラグをセルクリックで toggle 可。
                        既定 False(読み取り専用)。
        parent       : 親ウィジェット(オプション)。

    Signals:
        cellEdited(cls: str, adduct: str): 編集モードでセルが変更されたとき発火。
    """

    cellEdited = Signal(str, str)

    def __init__(self, adduct_columns=None, editable: bool = False,
                 parent=None):
        super().__init__(0, 0, parent)
        self._adduct_columns: tuple = tuple(
            adduct_columns or DEFAULT_ADDUCT_COLUMNS)
        self._editable: bool = bool(editable)
        self._classes: list = []
        self._data: dict = {}

        # テーブル基本設定
        self.setColumnCount(len(self._adduct_columns))
        self.setHorizontalHeaderLabels(list(self._adduct_columns))
        self.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self.verticalHeader().setVisible(True)
        self.verticalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)

        if not self._editable:
            self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        else:
            self.cellChanged.connect(self._on_cell_changed)

        self.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.setSelectionMode(QAbstractItemView.SingleSelection)

        # セルクリックでの toggle(編集モード時のみ)
        if self._editable:
            self.cellDoubleClicked.connect(self._on_cell_double_clicked)

    # ── データ操作 API ──────────────────────────────────────────

    def set_data(self, data: dict):
        """ヒートマップを再描画する。

        Args:
            data: {class_name: {adduct: {'value': float, 'expected': bool,
                                          'tooltip': str (optional)}}}
        """
        # 編集中シグナルが発火しないように一旦切る
        self.blockSignals(True)
        try:
            self._data = data or {}
            # 旧実装は sorted() で行ラベルを辞書順に並べていた
            # ため、"#1" < "#10" < "#2" のような文字列ソートで ID 順が崩れていた。
            # 呼び出し側が渡した dict の挿入順を尊重するよう list() に変更。
            self._classes = list(self._data.keys())
            self.setRowCount(len(self._classes))
            self.setVerticalHeaderLabels(self._classes)
            self.setColumnCount(len(self._adduct_columns))
            self.setHorizontalHeaderLabels(list(self._adduct_columns))

            # 値の最大値を取得(色スケール用)
            max_val = 0.0
            for cls_dict in self._data.values():
                for spec in cls_dict.values():
                    try:
                        v = float(spec.get('value', 0.0) or 0.0)
                    except (TypeError, ValueError):
                        v = 0.0
                    if v > max_val:
                        max_val = v
            if max_val <= 0:
                max_val = 1.0

            # セルを埋める
            for r, cls in enumerate(self._classes):
                row_data = self._data.get(cls, {})
                for c, adu in enumerate(self._adduct_columns):
                    spec = row_data.get(adu, {})
                    self._render_cell(r, c, cls, adu, spec, max_val)

            self.resizeColumnsToContents()
            self.resizeRowsToContents()
        finally:
            self.blockSignals(False)

    def get_classes(self) -> list:
        return list(self._classes)

    def get_columns(self) -> tuple:
        return self._adduct_columns

    def get_cell(self, cls: str, adduct: str) -> dict:
        """指定セルのデータを返す。"""
        return self._data.get(cls, {}).get(adduct, {})

    def get_data(self) -> dict:
        """現在のデータ全体を返す(deep copy ではない)。"""
        return self._data

    # ── 内部実装 ────────────────────────────────────────────────

    def _render_cell(self, row: int, col: int, cls: str, adu: str,
                     spec: dict, max_val: float):
        """1 セルの表示を更新する。"""
        try:
            value = float(spec.get('value', 0.0) or 0.0)
        except (TypeError, ValueError):
            value = 0.0
        expected = bool(spec.get('expected', False))
        tooltip_extra = str(spec.get('tooltip', ''))

        # 表示文字: ○ if expected, × otherwise
        text = "○" if expected else "×"
        item = QTableWidgetItem(text)
        item.setTextAlignment(Qt.AlignCenter)

        # フォント: ○/× をやや大きめに
        font = QFont()
        font.setPointSize(11)
        item.setFont(font)

        # 背景色: expected セルは値の濃淡で青系、非 expected セルは薄灰
        if expected and max_val > 0:
            intensity = min(max(value / max_val, 0.0), 1.0)
            # 白(255,255,255) から濃青(60, 120, 220) へグラデーション
            r_c = int(255 - intensity * (255 - 60))
            g_c = int(255 - intensity * (255 - 120))
            b_c = int(255 - intensity * (255 - 220))
            item.setBackground(QBrush(QColor(r_c, g_c, b_c)))
            # コントラスト維持: 濃い色のとき文字を白に
            if intensity > 0.5:
                item.setForeground(QBrush(QColor(255, 255, 255)))
        else:
            item.setBackground(QBrush(QColor(245, 245, 245)))
            item.setForeground(QBrush(QColor(160, 160, 160)))

        # tooltip
        tooltip = (f"{cls} × {adu}\n"
                   f"expected: {expected}\n"
                   f"value: {value:.4g}")
        if tooltip_extra:
            tooltip += f"\n{tooltip_extra}"
        item.setToolTip(tooltip)

        # 編集可否
        if not self._editable:
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
        else:
            # 編集モードでは UserRole に元の dict を持たせる(toggle 用)
            item.setData(Qt.UserRole, dict(spec))

        self.setItem(row, col, item)

    def _on_cell_double_clicked(self, row: int, col: int):
        """編集モード時、ダブルクリックで expected フラグを toggle。"""
        if row < 0 or col < 0:
            return
        if row >= len(self._classes) or col >= len(self._adduct_columns):
            return
        cls = self._classes[row]
        adu = self._adduct_columns[col]
        if cls not in self._data:
            self._data[cls] = {}
        if adu not in self._data[cls]:
            self._data[cls][adu] = {'expected': False, 'value': 0.0}
        cur = self._data[cls][adu]
        cur['expected'] = not bool(cur.get('expected', False))
        # 再描画
        max_val = 1.0
        for cls_dict in self._data.values():
            for spec in cls_dict.values():
                v = float(spec.get('value', 0.0) or 0.0)
                if v > max_val:
                    max_val = v
        self._render_cell(row, col, cls, adu, cur, max_val)
        self.cellEdited.emit(cls, adu)

    def _on_cell_changed(self, row: int, col: int):
        """編集モード時、セル値変更を検知(現在は未使用、将来拡張用)。"""
        pass


# ════════════════════════════════════════════════════════════════════
#  AdductPatternsDialog / AdductReportDialog
#  Adduct Ion Filter の参照パターンと観測結果の可視化ダイアログ
# ════════════════════════════════════════════════════════════════════


class AdductPatternsDialog(QDialog):
    """[Patterns…] ダイアログ。

    Adduct Ion Filter が参照する期待アダクトパターンを heatmap で表示する。
    データソース:
      - 内蔵 CURATED_ADDUCT_FINGERPRINTS(初期は空、必要に応じてキュレーション)
      - ランタイム IS 派生(現セッションの ② IS Filter で生成された
        self._runtime_is_fingerprints)

    Source 切替で「内蔵のみ / ランタイムのみ / 統合(ランタイムを内蔵で上書き)」を
    表示できる。

    現バージョンは表示専用(編集は将来の v5.4.x 拡張で対応)。
    """

    def __init__(
        self,
        curated: dict,
        runtime: dict,
        adduct_columns=None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Adduct Ion Pattern Reference")
        self.resize(900, 460)

        self._curated = curated or {}
        self._runtime = runtime or {}
        self._adduct_columns = tuple(
            adduct_columns or DEFAULT_ADDUCT_COLUMNS)

        lay = QVBoxLayout(self)

        # ── ヘッダ説明 ──
        info = QLabel(
            "Active adduct fingerprint patterns referenced by Adduct Ion "
            "Filter.\n"
            "Cell symbol: ○ = expected, × = not expected.  "
            "Cell color: weight (darker = larger).\n"
            "Source switch shows: Internal CURATED only / Runtime IS-derived "
            "only / Combined (Internal overrides Runtime)."
        )
        info.setWordWrap(True)
        lay.addWidget(info)

        # ── ソース切替 ──
        src_row = QHBoxLayout()
        src_row.addWidget(QLabel("Source:"))
        self._radio_internal = QRadioButton("Internal (CURATED)")
        self._radio_runtime = QRadioButton("Runtime (IS-derived)")
        self._radio_combined = QRadioButton("Combined")
        self._radio_combined.setChecked(True)
        for rb in (self._radio_internal, self._radio_runtime,
                   self._radio_combined):
            rb.toggled.connect(self._refresh)
            src_row.addWidget(rb)
        src_row.addStretch()
        # 件数ラベル
        self._lbl_count = QLabel("")
        src_row.addWidget(self._lbl_count)
        lay.addLayout(src_row)

        # ── ヒートマップ ──
        if AdductHeatmapWidget is None:
            lay.addWidget(QLabel(
                "Error: adduct_heatmap_widget module not loaded."))
            self._heatmap = None
        else:
            self._heatmap = AdductHeatmapWidget(
                adduct_columns=self._adduct_columns,
                editable=False,
                parent=self,
            )
            lay.addWidget(self._heatmap, stretch=1)

        # ── 閉じるボタン ──
        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        lay.addWidget(btns)

        self._refresh()

    def _build_view_data(self) -> dict:
        """ソース切替に応じた表示用データを構築。

        セルの 'value' には内蔵では weight、ランタイムでは weight を使用。
        合体時は内蔵 weight が優先(無ければランタイム weight)。
        """
        if self._radio_internal.isChecked():
            source = self._curated
        elif self._radio_runtime.isChecked():
            source = self._runtime
        else:
            # combined: ランタイム → 内蔵 で上書き
            source = dict(self._runtime)
            for cls, fp in self._curated.items():
                source[cls] = fp

        # AdductHeatmapWidget の入力形式に変換
        view: dict = {}
        for cls, cls_fp in source.items():
            view[cls] = {}
            for adu in self._adduct_columns:
                spec = cls_fp.get(adu, {}) or {}
                view[cls][adu] = {
                    'expected': bool(spec.get('expected', False)),
                    'value': float(spec.get('weight', 0.0) or 0.0),
                }
        return view

    def _refresh(self):
        if self._heatmap is None:
            return
        view = self._build_view_data()
        self._heatmap.set_data(view)
        self._lbl_count.setText(
            f"{len(view)} class(es) with patterns")


class AdductReportDialog(QDialog):
    """[Report…] ダイアログ。

    Adduct Ion Filter が判定した衝突スポットごとに、9 アダクトの観測有無を
    heatmap で表示する。各行は 1 衝突スポット (RT, m/z)、列は 9 アダクト、
    セルは ○(観測あり)/ ×(観測なし)。

    Attribution Table セクションを追加(⑤ 決着 + undecided)。
    行クリックで散布図ハイライト + 計算詳細表示。

    入力: PreviewRTDialog._adduct_attribution
      {coord: {'winner', 'losers', 'scores', 'observations', 'skipped'}}

    エクスポート: PNG (heatmap 画像) / xlsx (集計表)
    """

    spotSelected = Signal(float, float, str)

    def __init__(
        self,
        attribution: dict,
        adduct_columns=None,
        attribution_table: dict = None,
        sigma_threshold: float = 3.0,
        key_to_id: dict = None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Adduct Observation Report")
        # attribution table を加えるため縦サイズ拡大
        self.resize(1150, 780)

        self._attribution = attribution or {}
        self._adduct_columns = tuple(
            adduct_columns or DEFAULT_ADDUCT_COLUMNS)
        self._attribution_table = attribution_table or {}
        self._sigma_threshold = sigma_threshold
        # 共通 ID マップ + heatmap 行 → key の逆引き
        self._key_to_id = key_to_id or {}
        self._heatmap_row_to_key: list = []

        lay = QVBoxLayout(self)

        # ── 説明 ──
        info = QLabel(
            "Observed adduct patterns at each conflict spot.  "
            "Rows: conflict spots.  Columns: 9 adducts.\n"
            "Cell: ○ = adduct observed at this RT and theoretical m/z, "
            "× = not observed."
        )
        info.setWordWrap(True)
        lay.addWidget(info)

        # ── サマリ統計 ──
        n_total = len(self._attribution)
        n_decided = sum(
            1 for r in self._attribution.values()
            if r.get('winner') is not None and not r.get('skipped'))
        n_undecided = sum(
            1 for r in self._attribution.values()
            if r.get('winner') is None and not r.get('skipped'))
        n_skipped = sum(
            1 for r in self._attribution.values()
            if r.get('skipped'))
        summary = QLabel(
            f"Total conflicts: {n_total}    "
            f"Decided: {n_decided}    "
            f"Undecided: {n_undecided}    "
            f"Skipped (no fingerprint): {n_skipped}"
        )
        summary.setStyleSheet("color: #444; font-weight: bold;")
        lay.addWidget(summary)

        # ── ヒートマップ ──
        if AdductHeatmapWidget is None or not self._attribution:
            if not self._attribution:
                lay.addWidget(QLabel(
                    "No conflict spots to display. "
                    "Run Adduct Ion Filter first."))
            else:
                lay.addWidget(QLabel(
                    "Error: adduct_heatmap_widget module not loaded."))
            self._heatmap = None
        else:
            self._heatmap = AdductHeatmapWidget(
                adduct_columns=self._adduct_columns,
                editable=False,
                parent=self,
            )
            view = self._build_view_data()
            self._heatmap.set_data(view)
            lay.addWidget(self._heatmap, stretch=1)
            # heatmap 行クリック → 散布図ハイライト
            try:
                self._heatmap.cellClicked.connect(self._on_heatmap_row_click)
                self._heatmap.verticalHeader().sectionClicked.connect(
                    self._on_heatmap_row_click)
            except Exception:
                pass

        # Attribution Results table + 計算詳細パネル
        if self._attribution_table:
            sep = QLabel("Attribution Results (⑤ Adduct + Undecided)")
            sep.setStyleSheet(
                "font-weight:bold; color:#333; padding:8px 0 2px 0; "
                "border-top:1px solid #CCC;")
            lay.addWidget(sep)

            # ID 列を先頭に追加(heatmap と共通)
            headers = ["#", "Mode", "Method", "σ", "obs RT", "obs m/z",
                       "Winner class", "Winner compound",
                       "Loser class", "Loser compound",
                       "pred winner", "pred loser",
                       "res winner", "res loser"]
            self._attr_tbl = QTableWidget()
            self._attr_tbl.setColumnCount(len(headers))
            self._attr_tbl.setHorizontalHeaderLabels(headers)
            self._attr_tbl.setEditTriggers(QTableWidget.NoEditTriggers)
            self._attr_tbl.setSelectionBehavior(QTableWidget.SelectRows)
            self._attr_tbl.setAlternatingRowColors(True)

            def _fmt_float(v, fmt):
                try:
                    fv = float(v)
                    if not np.isfinite(fv):
                        return "-"
                    return f"{fv:{fmt}}"
                except (TypeError, ValueError):
                    return "-"
            self._attr_fmt_float = _fmt_float

            # heatmap と同じ union を反復して
            # 行数と ID を完全一致させる
            all_keys = set(
                (self._attribution_table or {}).keys()) | set(
                (self._attribution or {}).keys())
            rows = []
            for k in sorted(
                    all_keys,
                    key=lambda kk: self._key_to_id.get(kk, 1e9)):
                t_entry = self._attribution_table.get(k)
                if t_entry is None:
                    # heatmap には存在するが table 用には作られなかった key
                    # (例: ⑥ Coherence モデルが無く scored が 2 未満の coord)
                    # → 最低限の placeholder で表示
                    h_entry = self._attribution.get(k) or {}
                    winner = h_entry.get('winner')
                    skipped = bool(h_entry.get('skipped'))
                    mode_v = ''
                    if isinstance(k, tuple) and len(k) == 3:
                        mode_v = str(k[0])
                    t_entry = {
                        'method': ('adduct' if winner and not skipped
                                   else 'undecided' if not skipped
                                   else 'skipped'),
                        'mode': mode_v,
                        'confidence': 0.0,
                        'winner_class': winner or '',
                        'winner_compound': '',
                        'loser_class': '',
                        'loser_compound': '',
                        'pred_winner': float('nan'),
                        'pred_loser': float('nan'),
                        'res_winner': float('nan'),
                        'res_loser': float('nan'),
                        'low_confidence': False,
                    }
                rows.append((k, t_entry))

            self._attr_tbl.setRowCount(len(rows))
            self._attr_row_keys = []
            for r, (key, a) in enumerate(rows):
                self._attr_row_keys.append(key)
                if isinstance(key, tuple) and len(key) == 3:
                    mode_str = str(key[0])
                    rt_v, mz_v = float(key[1]), float(key[2])
                else:
                    mode_str = a.get('mode', '')
                    rt_v, mz_v = float(key[0]), float(key[1])
                method = a.get('method', 'adduct')
                method_label = {
                    'adduct': '⑤ Adduct',
                    'coherence': '⑥ Coherence',
                    'undecided': '△ Undecided',
                }.get(method, method)
                # 共通 ID
                idn = self._key_to_id.get(key, r + 1)
                cells = [
                    f"#{idn}",
                    mode_str, method_label,
                    _fmt_float(a.get('confidence', 0.0), '.2f'),
                    f"{rt_v:.3f}", f"{mz_v:.4f}",
                    a.get('winner_class', ''),
                    a.get('winner_compound', ''),
                    a.get('loser_class', ''),
                    a.get('loser_compound', ''),
                    _fmt_float(a.get('pred_winner'), '.3f'),
                    _fmt_float(a.get('pred_loser'), '.3f'),
                    _fmt_float(a.get('res_winner'), '+.3f'),
                    _fmt_float(a.get('res_loser'), '+.3f'),
                ]
                for c, txt in enumerate(cells):
                    item = QTableWidgetItem(txt)
                    if a.get('low_confidence'):
                        item.setForeground(QBrush(QColor('#B22222')))
                    self._attr_tbl.setItem(r, c, item)
            self._attr_tbl.resizeColumnsToContents()
            self._attr_tbl.cellClicked.connect(self._on_attr_cell_clicked)
            lay.addWidget(self._attr_tbl, stretch=1)

            # 計算詳細パネル
            self._attr_detail_text = QTextEdit()
            self._attr_detail_text.setReadOnly(True)
            self._attr_detail_text.setFixedHeight(180)
            self._attr_detail_text.setStyleSheet(
                "QTextEdit { background:#FAFAFA; color:#202020; "
                "font-family:'Consolas','Courier New',monospace; "
                "font-size:11px; border:1px solid #CCC; padding:4px; }")
            self._attr_detail_text.setPlainText(
                "Click a row to display calculation details here.")
            lay.addWidget(self._attr_detail_text)

        # ── エクスポート + 閉じるボタン ──
        btn_row = QHBoxLayout()
        btn_png = QPushButton("Save PNG…")
        btn_png.clicked.connect(self._save_png)
        btn_xlsx = QPushButton("Save xlsx…")
        btn_xlsx.clicked.connect(self._save_xlsx)
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.reject)
        btn_row.addWidget(btn_png)
        btn_row.addWidget(btn_xlsx)
        btn_row.addStretch()
        btn_row.addWidget(btn_close)
        lay.addLayout(btn_row)

    def _on_heatmap_row_click(self, row, col=None):
        """heatmap の行クリック → spotSelected 発火 +
        attribution table の同じ ID 行を選択 + 詳細パネル更新。
        verticalHeader.sectionClicked のときは引数 1 つ (row のみ)。
        """
        try:
            row_idx = int(row)
        except Exception:
            return
        if not (0 <= row_idx < len(self._heatmap_row_to_key)):
            return
        key = self._heatmap_row_to_key[row_idx]
        if isinstance(key, tuple) and len(key) == 3:
            mode_v = str(key[0]); rt_v, mz_v = float(key[1]), float(key[2])
        else:
            rt_v, mz_v = float(key[0]), float(key[1]); mode_v = ''
        self.spotSelected.emit(rt_v, mz_v, mode_v)
        # attribution table 側で対応する row を選択 + 詳細更新
        if hasattr(self, '_attr_row_keys') and key in self._attr_row_keys:
            try:
                attr_row = self._attr_row_keys.index(key)
                self._attr_tbl.selectRow(attr_row)
                self._attr_tbl.scrollToItem(
                    self._attr_tbl.item(attr_row, 0))
                a = self._attribution_table.get(key, {}) or {}
                self._update_attr_detail(a, (rt_v, mz_v))
            except Exception:
                pass

    def _on_attr_cell_clicked(self, row: int, col: int):
        """Attribution Table 行クリック時の処理"""
        if not (0 <= row < len(self._attr_row_keys)):
            return
        key = self._attr_row_keys[row]
        if isinstance(key, tuple) and len(key) == 3:
            mode_v = str(key[0])
            rt_v, mz_v = float(key[1]), float(key[2])
        else:
            a = self._attribution_table.get(key, {}) or {}
            mode_v = str(a.get('mode', ''))
            rt_v, mz_v = float(key[0]), float(key[1])
        self.spotSelected.emit(rt_v, mz_v, mode_v)
        # 詳細パネル更新
        a = self._attribution_table.get(key, {}) or {}
        self._update_attr_detail(a, (rt_v, mz_v))

    def _update_attr_detail(self, a: dict, coord: tuple):
        """⑤ attribution table 行の計算詳細表示"""
        if not a:
            self._attr_detail_text.setPlainText("(no detail)")
            return
        method = a.get('method', 'adduct')
        method_label = {
            'adduct': '⑤ Adduct Ion Filter',
            'undecided': '△ Undecided',
        }.get(method, method)
        lines = []
        lines.append(f"━━ Spot: RT = {coord[0]:.3f} min,  m/z = {coord[1]:.4f} ━━")
        lines.append(f"  Method:  {method_label}")
        lines.append(f"  σ (margin/confidence): {a.get('confidence', 0.0):.3f}")
        lines.append("")
        lines.append(f"[Winner]  {a.get('winner_class','')}  "
                     f"{a.get('winner_compound','')}  "
                     f"{a.get('winner_adduct','')}")
        if 'winner_C' in a:
            lines.append(f"  link={a.get('winner_link','?')}  "
                         f"C={a['winner_C']}  U={a['winner_U']}  "
                         f"(fit_mode: {a.get('winner_fit_mode','?')})")
            lines.append(f"  α={a['winner_alpha']:+.4f}  "
                         f"β={a['winner_beta']:+.4f}  "
                         f"γ={a['winner_gamma']:.4f}  "
                         f"rmse={a['winner_rmse']:.4f}")
            pw = a.get('pred_winner')
            rw = a.get('res_winner')
            if pw is not None and isinstance(pw, (int, float)) \
                    and np.isfinite(pw):
                lines.append(f"  pred = {pw:.4f}  res = {rw:+.4f}")
        lines.append("")
        lines.append(f"[Loser]   {a.get('loser_class','')}  "
                     f"{a.get('loser_compound','')}  "
                     f"{a.get('loser_adduct','')}")
        if 'loser_C' in a:
            lines.append(f"  link={a.get('loser_link','?')}  "
                         f"C={a['loser_C']}  U={a['loser_U']}  "
                         f"(fit_mode: {a.get('loser_fit_mode','?')})")
            lines.append(f"  α={a['loser_alpha']:+.4f}  "
                         f"β={a['loser_beta']:+.4f}  "
                         f"γ={a['loser_gamma']:.4f}  "
                         f"rmse={a['loser_rmse']:.4f}")
            pl = a.get('pred_loser')
            rl = a.get('res_loser')
            if pl is not None and isinstance(pl, (int, float)) \
                    and np.isfinite(pl):
                lines.append(f"  pred = {pl:.4f}  res = {rl:+.4f}")
        # adduct scores
        if method == 'adduct':
            scores = a.get('adduct_scores', {})
            if scores:
                lines.append("")
                lines.append("[Adduct evidence (⑤)]")
                for cls, sc in sorted(scores.items(), key=lambda x: -x[1]):
                    lines.append(f"  {cls}: score = {sc:.2f}")
        self._attr_detail_text.setPlainText("\n".join(lines))

    def _row_label(self, key, result: dict) -> str:
        """heatmap の行ラベルを生成。key は (mode, rt, mz) または (rt, mz)。"""
        # 3-tuple key (mode, rt, mz) 対応
        if isinstance(key, tuple) and len(key) == 3:
            mode, rt, mz = key
            mode_prefix = f"[{mode}] "
        else:
            rt, mz = key
            mode_prefix = ""
        # 共通 ID prefix(下段 attribution table と一致)
        idn = self._key_to_id.get(key)
        id_prefix = f"#{idn}  " if idn is not None else ""
        winner = result.get('winner')
        skipped = result.get('skipped', False)
        marker = ""
        if skipped:
            marker = " [skipped]"
        elif winner is None:
            marker = " [undecided]"
        else:
            marker = f" [{winner}]"
        return f"{id_prefix}{mode_prefix}{rt:.3f} / {mz:.4f}{marker}"

    def _build_view_data(self) -> dict:
        """attribution → AdductHeatmapWidget の入力形式に変換。

        旧実装では self._attribution.items() のみ反復していた
        ため、下段 attribution_table と key 集合が異なり ID が飛び飛びになる
        問題があった。両者の union を反復し、observations は self._attribution
        から(あれば)取得、無ければ空 obs として描画する。これで heatmap の
        行と attribution table の行が 1 対 1 で対応する。
        """
        view: dict = {}
        # union of heatmap + table keys, sorted by common ID
        all_keys = set(self._attribution.keys()) | set(
            (self._attribution_table or {}).keys())
        sorted_keys = sorted(
            all_keys, key=lambda k: self._key_to_id.get(k, 1e9))
        # heatmap 行 index → key の逆引きを構築
        self._heatmap_row_to_key = []
        for key in sorted_keys:
            result = self._attribution.get(key) or {}
            # table 由来のキー(_attribution に無い)は table エントリから
            # 表示用の補完情報を組み立てる
            if not result and self._attribution_table:
                t_entry = self._attribution_table.get(key) or {}
                if t_entry:
                    result = {
                        'winner': t_entry.get('winner_class'),
                        'observations': {},
                        'observation_intensities': {},
                        'skipped': False,
                    }
            row_key = self._row_label(key, result)
            obs = result.get('observations', {}) or {}
            # 観測ピーク強度(adduct → max intensity)
            obs_int = result.get('observation_intensities', {}) or {}
            # 行内 max で正規化(各 spot 単位でグラデーション)
            row_max = 0.0
            for adu in self._adduct_columns:
                if obs.get(adu, False):
                    v = float(obs_int.get(adu, 0.0) or 0.0)
                    if v > row_max:
                        row_max = v
            row_data: dict = {}
            for adu in self._adduct_columns:
                observed = bool(obs.get(adu, False))
                if observed:
                    if row_max > 0:
                        # log 圧縮で広 dynamic range の差を見やすく
                        intensity = float(obs_int.get(adu, 0.0) or 0.0)
                        # 範囲 [0.25, 1.0] にマッピング(最弱でも 0.25 で見える)
                        value = 0.25 + 0.75 * (intensity / row_max)
                    else:
                        value = 1.0   # 強度情報無し → 従来通り
                else:
                    value = 0.0
                row_data[adu] = {
                    'expected': observed,
                    'value': value,
                }
            view[row_key] = row_data
            self._heatmap_row_to_key.append(key)
        return view

    def _save_png(self):
        """heatmap を PNG にエクスポート。"""
        if self._heatmap is None:
            QMessageBox.warning(self, "No data", "No heatmap to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save heatmap PNG",
            timestamped_filename("adduct_observation_report.png"),
            "PNG (*.png)")
        if not path:
            return
        # QTableWidget 全体を QPixmap にレンダリング
        pixmap = self._heatmap.grab()
        if pixmap.save(path, "PNG"):
            QMessageBox.information(
                self, "Saved", f"Saved: {path}")
        else:
            QMessageBox.warning(
                self, "Save failed",
                f"Could not save PNG to:\n{path}")

    def _save_xlsx(self):
        """attribution を xlsx にエクスポート。"""
        if not self._attribution:
            QMessageBox.warning(self, "No data", "No data to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save observation report xlsx",
            timestamped_filename("adduct_observation_report.xlsx"),
            "Excel (*.xlsx)")
        if not path:
            return

        # 表データを構築
        rows = []
        for key, result in self._attribution.items():
            if isinstance(key, tuple) and len(key) == 3:
                mode, rt, mz = key
            else:
                rt, mz = key
                mode = ''
            obs = result.get('observations', {}) or {}
            scores = result.get('scores', {}) or {}
            row: dict = {
                'mode': mode,
                'RT': rt,
                'm/z': mz,
                'winner': result.get('winner') or '',
                'losers': '/'.join(sorted(result.get('losers', set()))),
                'skipped': result.get('skipped', False),
            }
            for adu in self._adduct_columns:
                row[f"obs:{adu}"] = bool(obs.get(adu, False))
            for cls, sc in scores.items():
                row[f"score:{cls}"] = float(sc)
            rows.append(row)
        df = pd.DataFrame(rows)
        try:
            df.to_excel(path, index=False)
            QMessageBox.information(
                self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(
                self, "Save failed",
                f"Could not save xlsx:\n{e}")


class RejectedSpotsDialog(QDialog):
    """③ IS Filter / ④ RT Outlier Filter で除外された
    スポットを一覧表示するダイアログ。

    機能:
      - rejected 行のテーブル表示(Mode / Class / Compound / Adduct / RT / m/z / ppm)
      - 行クリックで散布図ハイライト(spotSelected シグナルを emit)
      - xlsx エクスポートボタン
    """

    spotSelected = Signal(float, float, str)  # (obs_rt, obs_mz, mode)

    def __init__(
        self,
        title: str,
        rows: list,
        why: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(900, 560)
        self._rows = list(rows or [])
        self._why = why or ""

        lay = QVBoxLayout(self)

        # 概要ヘッダ
        n = len(self._rows)
        n_pos = sum(1 for r in self._rows if r.get('mode') == 'pos')
        n_neg = sum(1 for r in self._rows if r.get('mode') == 'neg')
        summary = (
            f"Rejected spots: {n}    "
            f"Pos: {n_pos}    Neg: {n_neg}\n"
            f"Reason: {self._why}")
        lbl = QLabel(summary)
        lbl.setStyleSheet(
            "font-weight: bold; padding: 6px; "
            "background:#FDF2F2; color:#7A2222; "
            "border:1px solid #E0241E;")
        lbl.setWordWrap(True)
        lay.addWidget(lbl)

        # テーブル
        headers = [
            "Mode", "Class", "Compound", "Adduct",
            "obs RT", "obs m/z", "Δppm",
        ]
        self._tbl = QTableWidget()
        self._tbl.setColumnCount(len(headers))
        self._tbl.setHorizontalHeaderLabels(headers)
        self._tbl.setEditTriggers(QTableWidget.NoEditTriggers)
        self._tbl.setSelectionBehavior(QTableWidget.SelectRows)
        self._tbl.setAlternatingRowColors(True)

        self._tbl.setRowCount(n)
        for r, row in enumerate(self._rows):
            cells = [
                str(row.get('mode', '')),
                str(row.get('lipid_class', '')),
                str(row.get('compound', '')),
                str(row.get('adduct', '')),
                f"{float(row.get('obs_rt', 0.0)):.3f}",
                f"{float(row.get('obs_mz', 0.0)):.4f}",
                f"{float(row.get('ppm', 0.0)):+.2f}",
            ]
            for c, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                # rejected を赤系で示す
                item.setForeground(QBrush(QColor('#7A2222')))
                self._tbl.setItem(r, c, item)
        self._tbl.resizeColumnsToContents()
        self._tbl.cellClicked.connect(self._on_cell_clicked)
        lay.addWidget(self._tbl, stretch=1)

        # ボタン群
        btns = QHBoxLayout()
        btn_export = QPushButton("Export xlsx…")
        btn_export.clicked.connect(self._on_export)
        btns.addWidget(btn_export)
        btns.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btns.addWidget(btn_close)
        lay.addLayout(btns)

    def _on_cell_clicked(self, row: int, col: int):
        if 0 <= row < len(self._rows):
            r = self._rows[row]
            try:
                rt = float(r.get('obs_rt', 0.0))
                mz = float(r.get('obs_mz', 0.0))
                mode = str(r.get('mode', ''))
                self.spotSelected.emit(rt, mz, mode)
            except (TypeError, ValueError):
                pass

    def _on_export(self):
        if not self._rows:
            QMessageBox.information(
                self, "No data", "Nothing to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save rejected spots",
            timestamped_filename("rejected_spots.xlsx"),
            "Excel files (*.xlsx)")
        if not path:
            return
        try:
            df = pd.DataFrame(self._rows)
            df.to_excel(path, index=False)
            QMessageBox.information(
                self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(
                self, "Save failed",
                f"Could not save xlsx:\n{e}")


class FullReportDialog(QDialog):
    """全 filter ステップの結果を 1 つの表に統合表示。

    含まれる event:
      ③ IS Filter — Rejected (REJ_WINDOW)
      ④ RT Outlier — Rejected (REJ_OUTLIER)
      ⑤ Adduct Filter — Attributed (winner) / Rejected (loser)
      ⑥ Coherence Filter — Attributed (winner) / Rejected (loser)
      ⑥ Coherence (residual) — Rejected (同一化合物の別ピークが採用された)
      △ Undecided — どの方法でも決着できなかった衝突

    ⑤ / ⑥ の行は match_df の status 列から作られる。FINAL 列に
    その行の final_status を併記するので、Spot Status Inspector と
    突き合わせられる。

    機能:
      - 概要ヘッダ(全件数 / ステップ別カウント)
      - ステップ別 filter checkbox(表示行を絞れる)
      - テーブル: Mode / Step / Status / Class / Compound / Adduct / RT / m/z / Detail
      - 行クリックで散布図ハイライト(spotSelected シグナル)
      - xlsx エクスポート
    """

    spotSelected = Signal(float, float, str)  # (obs_rt, obs_mz, mode)

    def __init__(self, events: list, parent=None):
        super().__init__(parent)
        self.setWindowTitle("All Filtering Results")
        self.resize(1100, 700)
        self._all_events = list(events or [])
        self._visible_events: list = list(self._all_events)

        lay = QVBoxLayout(self)

        # 概要ヘッダ
        n_total = len(self._all_events)
        from collections import Counter
        step_counts = Counter(e['step'] for e in self._all_events)
        mode_counts = Counter(e['mode'] for e in self._all_events)
        summary = (
            f"Total events: {n_total}    "
            f"Pos: {mode_counts.get('pos', 0)}    "
            f"Neg: {mode_counts.get('neg', 0)}\n"
            + "    ".join(
                f"{step}: {count}"
                for step, count in sorted(step_counts.items())))
        lbl = QLabel(summary)
        lbl.setStyleSheet(
            "font-weight: bold; padding: 6px; "
            "background:#F2F4F8; color:#2A3850; "
            "border:1px solid #5B7AB0;")
        lbl.setWordWrap(True)
        lay.addWidget(lbl)

        # ステップ別フィルタ checkbox
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Show steps:"))
        self._step_checks: dict = {}
        for step in sorted(step_counts.keys()):
            cb = QCheckBox(step)
            cb.setChecked(True)
            cb.stateChanged.connect(self._refresh_table)
            self._step_checks[step] = cb
            filter_row.addWidget(cb)
        filter_row.addStretch()
        lay.addLayout(filter_row)

        # テーブル
        # FINAL 列。Summary の判定と Spot Status Inspector の
        # 表示が食い違って見えないよう、その行の最終状態を併記する。
        self._headers = [
            "Mode", "Step", "Status",
            "Class", "Compound", "Adduct",
            "obs RT", "obs m/z", "FINAL", "Detail",
        ]
        self._tbl = QTableWidget()
        self._tbl.setColumnCount(len(self._headers))
        self._tbl.setHorizontalHeaderLabels(self._headers)
        self._tbl.setEditTriggers(QTableWidget.NoEditTriggers)
        self._tbl.setSelectionBehavior(QTableWidget.SelectRows)
        self._tbl.setAlternatingRowColors(True)
        self._tbl.cellClicked.connect(self._on_cell_clicked)
        lay.addWidget(self._tbl, stretch=1)

        # 行選択時の detail パネル
        self._detail_label = QLabel(
            "Click a row to display details here.")
        self._detail_label.setStyleSheet(
            "font-size: 10px; color: #555; padding: 4px;")
        lay.addWidget(self._detail_label)
        self._detail_text = QTextEdit()
        self._detail_text.setReadOnly(True)
        self._detail_text.setFixedHeight(120)
        self._detail_text.setStyleSheet(
            "QTextEdit { background:#FAFAFA; color:#202020; "
            "font-family:'Consolas','Courier New',monospace; "
            "font-size:11px; border:1px solid #CCC; padding:4px; }")
        lay.addWidget(self._detail_text)

        # ボタン群
        btns = QHBoxLayout()
        btn_export = QPushButton("Export xlsx…")
        btn_export.clicked.connect(self._on_export)
        btns.addWidget(btn_export)
        btns.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btns.addWidget(btn_close)
        lay.addLayout(btns)

        # 初期描画
        self._refresh_table()

    def _refresh_table(self):
        """ステップフィルタを反映してテーブルを再描画。"""
        active_steps = {
            step for step, cb in self._step_checks.items() if cb.isChecked()
        }
        self._visible_events = [
            e for e in self._all_events if e['step'] in active_steps
        ]
        self._tbl.setRowCount(len(self._visible_events))
        for r, e in enumerate(self._visible_events):
            cells = [
                str(e.get('mode', '')),
                str(e.get('step', '')),
                str(e.get('status', '')),
                str(e.get('lipid_class', '')),
                str(e.get('compound', '')),
                str(e.get('adduct', '')),
                f"{float(e.get('obs_rt', 0.0)):.3f}",
                f"{float(e.get('obs_mz', 0.0)):.4f}",
                str(e.get('final', '')),
                str(e.get('detail', '')),
            ]
            for c, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                # status により色付け
                status = e.get('status', '')
                if status == 'Rejected':
                    item.setForeground(QBrush(QColor('#7A2222')))
                elif status == 'Undecided':
                    item.setForeground(QBrush(QColor('#A06820')))
                elif status.endswith('rejected later'):
                    item.setForeground(QBrush(QColor('#A06820')))
                # FINAL 列は kept / rejected で塗り分ける
                if self._headers[c] == 'FINAL':
                    if txt == 'kept':
                        item.setBackground(QColor('#DCF3D6'))
                    elif txt == 'rejected':
                        item.setBackground(QColor('#F8DCDC'))
                self._tbl.setItem(r, c, item)
        self._tbl.resizeColumnsToContents()

    def _on_cell_clicked(self, row: int, col: int):
        if 0 <= row < len(self._visible_events):
            e = self._visible_events[row]
            try:
                rt = float(e.get('obs_rt', 0.0))
                mz = float(e.get('obs_mz', 0.0))
                mode = str(e.get('mode', ''))
                self.spotSelected.emit(rt, mz, mode)
            except (TypeError, ValueError):
                pass
            # detail panel
            lines = [
                f"━━ Spot: RT = {e.get('obs_rt', 0.0):.3f} min, "
                f"m/z = {e.get('obs_mz', 0.0):.4f} ━━",
                f"  Mode:     {e.get('mode', '')}",
                f"  Step:     {e.get('step', '')}",
                f"  Status:   {e.get('status', '')}",
                f"  Class:    {e.get('lipid_class', '')}",
                f"  Compound: {e.get('compound', '')}",
                f"  Adduct:   {e.get('adduct', '')}",
                f"  Detail:   {e.get('detail', '')}",
            ]
            self._detail_text.setPlainText("\n".join(lines))

    def _on_export(self):
        if not self._visible_events:
            QMessageBox.information(
                self, "No data", "Nothing to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save full filtering report",
            timestamped_filename("all_filtering_results.xlsx"),
            "Excel files (*.xlsx)")
        if not path:
            return
        try:
            df = pd.DataFrame(self._visible_events)
            df.to_excel(path, index=False)
            QMessageBox.information(
                self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(
                self, "Save failed",
                f"Could not save xlsx:\n{e}")


# JointPlotDialog で使用する図のフォーマット設定
DEFAULT_FORMAT_SETTINGS: dict = {
    'title': '',           # 空文字 = 自動生成
    'show_title': True,
    'title_size': 11,
    'label_size': 10,
    'tick_size': 9,
    'legend_size': 8,
    'marker_size': 24,
    'marker_alpha': 0.7,
    'save_dpi': 150,
    'fig_width_in': 8.0,
    'fig_height_in': 6.0,
}


class SpotInspectorDialog(QDialog):
    """クリックした spot の全 status を一覧表示するデバッグ用ダイアログ。

    同一 (RT, m/z) 座標にある全 matched 行(複数候補クラス含む、kept も
    rejected も)を表示し、それぞれの:
      - lipid_class, compound, adduct, Δppm
      - is_filter_status / adduct_filter_status / rt_outlier_status /
        coherence_status / manual_status / final_status
    を確認できる。

    用途: 「なぜこの spot が rejected なのか」「kept なのにフィルタを通った
    のはどの行か」を診断する。
    """

    def __init__(self, mode: str, coord: tuple, rows_df, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Spot Status Inspector")
        self.resize(880, 380)

        lay = QVBoxLayout(self)

        # ヘッダ
        if rows_df is None or len(rows_df) == 0:
            n_rows = 0
        else:
            n_rows = len(rows_df)
        lbl = QLabel(
            f"Spot: [{mode}] RT = {coord[0]:.3f} min, m/z = {coord[1]:.4f}\n"
            f"{n_rows} matched library row(s) at this coordinate.")
        lbl.setStyleSheet(
            "font-weight: bold; padding: 6px; "
            "background: #F0F4FA; border: 1px solid #5B7AB0;")
        lbl.setWordWrap(True)
        lay.addWidget(lbl)

        # テーブル
        headers = [
            "Class", "Compound", "Adduct", "Δppm",
            "IS", "Adduct flt.", "RT outlier",
            "Coherence", "Manual", "FINAL",
        ]
        tbl = QTableWidget()
        tbl.setColumnCount(len(headers))
        tbl.setHorizontalHeaderLabels(headers)
        tbl.setEditTriggers(QTableWidget.NoEditTriggers)
        tbl.setSelectionBehavior(QTableWidget.SelectRows)
        tbl.setAlternatingRowColors(True)

        if rows_df is not None and len(rows_df) > 0:
            tbl.setRowCount(n_rows)
            for r in range(n_rows):
                row = rows_df.iloc[r]
                final = str(row.get('final_status', '-'))
                cells = [
                    str(row.get('lipid_class', '')),
                    str(row.get('compound', '')),
                    str(row.get('adduct', '')),
                    f"{float(row.get('delta_ppm', 0.0) or 0.0):+.2f}",
                    str(row.get('is_filter_status', '-')),
                    str(row.get('adduct_filter_status', '-')),
                    str(row.get('rt_outlier_status', '-')),
                    str(row.get('coherence_status', '-')),
                    str(row.get('manual_status', '-')),
                    final,
                ]
                for c, txt in enumerate(cells):
                    item = QTableWidgetItem(txt)
                    # FINAL 列に色付け
                    if c == len(headers) - 1:
                        if txt == 'kept':
                            item.setBackground(QColor('#DCF3D6'))
                            item.setForeground(QBrush(QColor('#1F5C12')))
                        elif txt == 'rejected':
                            item.setBackground(QColor('#F8DCDC'))
                            item.setForeground(QBrush(QColor('#7A2222')))
                    tbl.setItem(r, c, item)
        tbl.resizeColumnsToContents()
        lay.addWidget(tbl, stretch=1)

        # 凡例
        note = QLabel(
            "Status legend:\n"
            "  IS:  'kept' / 'rejected_window' (outside IS RT window)\n"
            "  Adduct flt.:  'rejected_adduct' (loser at conflict spot)\n"
            "  RT outlier:  'rejected_outlier' (deviates from class median)\n"
            "  Coherence:  'rejected_residual' / 'rejected_loser'\n"
            "  Manual:  'kept' / 'rejected_manual' (user override)\n"
            "  FINAL:  consolidated decision used for export "
            "(kept = visible by default).")
        note.setStyleSheet(
            "font-size: 10px; color: #555; padding: 4px;")
        note.setWordWrap(True)
        lay.addWidget(note)

        # ボタン
        btns = QHBoxLayout()
        btns.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btns.addWidget(btn_close)
        lay.addLayout(btns)


class QuantIonSelectorDialog(QDialog):
    """各脂質クラスの定量イオンモード(pos / neg / skip)選択ダイアログ。

    フィルタ + 手動帰属が完了した後、最終的な定量に使うイオンモードを
    クラス単位で決定する。両モードに matched & kept スポットがある場合の
    曖昧さを解決する。

    機能:
      - クラスごとの Pos / Neg matched 件数を表示
      - Pos / Neg / Skip のラジオで選択(初期値は件数の多い方)
      - Save JSON / Load JSON で決定を外部ファイルに保存・読込
      - Reset で全てを「件数の多い方」に戻す
      - クラス名クリックで散布図をそのクラスにフォーカス
    """

    # クラス名クリック時に発火(空文字 → 全クラス表示)
    classFocusRequested = Signal(str)

    def __init__(self, per_class: dict, initial: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Quant Ion")
        self.resize(620, 540)
        # per_class: {cls: {'pos': int, 'neg': int}}
        self._per_class = dict(per_class)
        self._initial = dict(initial or {})

        lay = QVBoxLayout(self)

        intro = QLabel(
            "For each lipid class, choose the ion mode used for "
            "quantification.\n"
            "Pos / Neg = use that mode's data,  Skip = exclude the class "
            "from export.\n"
            "Defaults: mode with more matched & kept spots (Skip if both "
            "are zero).\n"
            "Tip: click a class name to focus the scatter plots on that "
            "class.")
        intro.setStyleSheet("color:#444; padding:4px;")
        intro.setWordWrap(True)
        lay.addWidget(intro)

        # クラスごとに Radio をテーブル状に配置
        self._radio_groups: dict = {}   # cls -> (rb_pos, rb_neg, rb_skip)
        # ヘッダ
        hdr = QHBoxLayout()
        for label, w in [("Class", 80), ("Pos n", 60), ("Neg n", 60),
                         ("Quant", 280)]:
            l = QLabel(label)
            l.setStyleSheet("font-weight:bold; color:#333;")
            l.setMinimumWidth(w)
            hdr.addWidget(l)
        hdr.addStretch()
        lay.addLayout(hdr)

        # 行をスクロールエリアに入れる
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        rows_lay = QVBoxLayout(content)
        rows_lay.setSpacing(2)

        for cls in sorted(self._per_class.keys()):
            counts = self._per_class[cls]
            n_pos = int(counts.get('pos', 0))
            n_neg = int(counts.get('neg', 0))
            # 初期値の決定
            if cls in self._initial:
                initial_choice = self._initial[cls]
            else:
                if n_pos == 0 and n_neg == 0:
                    initial_choice = 'skip'
                elif n_pos >= n_neg:
                    initial_choice = 'pos'
                else:
                    initial_choice = 'neg'

            row = QHBoxLayout()
            # クラス名は QToolButton 化、クリックでフォーカス
            btn_cls = QToolButton()
            btn_cls.setText(str(cls))
            btn_cls.setAutoRaise(True)   # flat 表示で label 風
            btn_cls.setCursor(Qt.PointingHandCursor)
            btn_cls.setMinimumWidth(80)
            btn_cls.setStyleSheet(
                "QToolButton {"
                " font-weight: bold; padding: 2px 6px; text-align: left;"
                "}"
                "QToolButton:hover {"
                " background-color: #E8F0FB; color: #1A4F8F;"
                "}")
            btn_cls.setToolTip(
                f"Click to focus the scatter plots on {cls} only.")
            btn_cls.clicked.connect(
                lambda _checked=False, c=str(cls):
                    self.classFocusRequested.emit(c))
            row.addWidget(btn_cls)
            lbl_pos = QLabel(str(n_pos))
            lbl_pos.setMinimumWidth(60)
            lbl_pos.setStyleSheet(
                "color:#888;" if n_pos == 0 else "color:#222;")
            row.addWidget(lbl_pos)
            lbl_neg = QLabel(str(n_neg))
            lbl_neg.setMinimumWidth(60)
            lbl_neg.setStyleSheet(
                "color:#888;" if n_neg == 0 else "color:#222;")
            row.addWidget(lbl_neg)
            # Radio ボタン
            from PySide6.QtWidgets import QButtonGroup
            grp = QButtonGroup(self)
            rb_pos = QRadioButton("Pos")
            rb_neg = QRadioButton("Neg")
            rb_skip = QRadioButton("Skip")
            # Pos / Neg が件数 0 なら disable
            if n_pos == 0:
                rb_pos.setEnabled(False)
            if n_neg == 0:
                rb_neg.setEnabled(False)
            grp.addButton(rb_pos)
            grp.addButton(rb_neg)
            grp.addButton(rb_skip)
            if initial_choice == 'pos' and rb_pos.isEnabled():
                rb_pos.setChecked(True)
            elif initial_choice == 'neg' and rb_neg.isEnabled():
                rb_neg.setChecked(True)
            else:
                rb_skip.setChecked(True)
            row.addWidget(rb_pos)
            row.addWidget(rb_neg)
            row.addWidget(rb_skip)
            row.addStretch()
            rows_lay.addLayout(row)
            self._radio_groups[str(cls)] = (rb_pos, rb_neg, rb_skip)

        rows_lay.addStretch()
        scroll.setWidget(content)
        lay.addWidget(scroll, stretch=1)

        # ボタン
        btns = QHBoxLayout()
        # 散布図のフォーカス解除(全クラス表示)
        btn_show_all = QPushButton("Show All Classes")
        btn_show_all.setToolTip(
            "Restore the scatter plots to show all classes (undo focus).")
        btn_show_all.clicked.connect(
            lambda: self.classFocusRequested.emit(""))
        btns.addWidget(btn_show_all)
        btn_reset = QPushButton("Reset (default by count)")
        btn_reset.clicked.connect(self._on_reset)
        btns.addWidget(btn_reset)
        # ボタン文言を Session File に統一
        btn_save = QPushButton("Save LipidZoner Session File…")
        btn_save.setToolTip(
            "Save the unified LipidZoner Session File (JSON).\n"
            "Current ion selection is written. Existing sections in\n"
            "the file are preserved (merged).")
        btn_save.clicked.connect(self._on_save_json)
        btns.addWidget(btn_save)
        btn_load = QPushButton("Load LipidZoner Session File…")
        btn_load.setToolTip(
            "Load the unified LipidZoner Session File (JSON).\n"
            "Ion selection (and other available sections) is restored.")
        btn_load.clicked.connect(self._on_load_json)
        btns.addWidget(btn_load)
        btns.addStretch()
        btn_cancel = QPushButton("Cancel")
        btn_cancel.clicked.connect(self.reject)
        btns.addWidget(btn_cancel)
        btn_ok = QPushButton("OK")
        btn_ok.setDefault(True)
        btn_ok.clicked.connect(self.accept)
        btns.addWidget(btn_ok)
        lay.addLayout(btns)

    def _on_reset(self):
        """全クラスを「件数が多い方」の初期値に戻す。"""
        for cls, (rb_pos, rb_neg, rb_skip) in self._radio_groups.items():
            counts = self._per_class.get(cls, {})
            n_pos = int(counts.get('pos', 0))
            n_neg = int(counts.get('neg', 0))
            if n_pos == 0 and n_neg == 0:
                rb_skip.setChecked(True)
            elif n_pos >= n_neg:
                rb_pos.setChecked(True)
            else:
                rb_neg.setChecked(True)

    def _apply_choices(self, choices: dict):
        """choices dict を radio に適用。"""
        for cls, choice in choices.items():
            grp = self._radio_groups.get(str(cls))
            if grp is None:
                continue
            rb_pos, rb_neg, rb_skip = grp
            if choice == 'pos' and rb_pos.isEnabled():
                rb_pos.setChecked(True)
            elif choice == 'neg' and rb_neg.isEnabled():
                rb_neg.setChecked(True)
            else:
                rb_skip.setChecked(True)

    def _on_save_json(self):
        # ファイル名・タイトルを統一
        path, _ = QFileDialog.getSaveFileName(
            self, "Save LipidZoner Session File",
            "lipidzoner_session.json", "JSON (*.json)")
        if not path:
            return
        try:
            # parent chain (QuantIonSelectorDialog →
            # PreviewRTDialog → MainWindow) を辿って MainWindow の
            # save_unified_state に委譲。session + params + quant_ion +
            # filter_exemptions の完全な統合 JSON を書き込む。
            mw = None
            p = self.parent()
            while p is not None:
                if hasattr(p, 'save_unified_state'):
                    mw = p
                    break
                # _host 経由(PreviewRTDialog は parent=None 生成のため)
                next_p = getattr(p, '_host', None)
                if next_p is None:
                    try:
                        next_p = p.parent() if hasattr(p, 'parent') else None
                    except Exception:
                        next_p = None
                p = next_p
            if mw is not None:
                mw.save_unified_state(
                    path, override_quant_ion=self.result_choices())
            else:
                # MainWindow に辿り着けなければ保存を諦める
                # (host が無い状態は本来発生しない)
                QMessageBox.warning(
                    self, "Save partial",
                    "Host main window not available. Session save skipped.")
                return
            QMessageBox.information(self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))

    def _on_load_json(self):
        """v2 形式の session JSON から quant_ion_choices だけを適用。"""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load LipidZoner Session File",
            "", "JSON (*.json)")
        if not path:
            return
        try:
            data = _load_lz_v2_file(path)
            qc = data.get('quant_ion_choices')
            if not isinstance(qc, dict) or not qc:
                raise ValueError(
                    "No 'quant_ion_choices' section found in the file.")
            self._apply_choices(qc)
            QMessageBox.information(self, "Loaded", f"Loaded: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Load failed", str(e))

    def result_choices(self) -> dict:
        """選択結果を {class: 'pos' | 'neg' | 'skip'} で返す。"""
        out: dict = {}
        for cls, (rb_pos, rb_neg, rb_skip) in self._radio_groups.items():
            if rb_pos.isChecked():
                out[cls] = 'pos'
            elif rb_neg.isChecked():
                out[cls] = 'neg'
            else:
                out[cls] = 'skip'
        return out








class FormatSettingsDialog(QDialog):
    """散布図のフォーマット設定(タイトル / フォント / サイズ等)
    を編集する小ダイアログ。Pos/Neg snapshot モードでは適用範囲が限定的なため、
    OK 押下後は呼び出し側が render モードに切り替えてから設定を反映する。
    """

    def __init__(self, current_settings: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Figure Format Settings")
        self.resize(420, 460)
        self._settings = dict(DEFAULT_FORMAT_SETTINGS)
        self._settings.update(current_settings or {})

        lay = QVBoxLayout(self)
        intro = QLabel(
            "Customize figure formatting. Settings apply when you click OK.\n"
            "Note: Pos / Neg will switch to render mode (Analysis snapshot off).")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #555; padding: 4px;")
        lay.addWidget(intro)

        form = QFormLayout()

        self._le_title = QLineEdit(self._settings.get('title', ''))
        self._le_title.setPlaceholderText("(empty = auto-generated)")
        form.addRow("Title text:", self._le_title)

        self._chk_show_title = QCheckBox("Show title")
        self._chk_show_title.setChecked(
            bool(self._settings.get('show_title', True)))
        form.addRow("", self._chk_show_title)

        self._spin_title_size = QSpinBox()
        self._spin_title_size.setRange(6, 32)
        self._spin_title_size.setValue(
            int(self._settings.get('title_size', 11)))
        form.addRow("Title font size:", self._spin_title_size)

        self._spin_label_size = QSpinBox()
        self._spin_label_size.setRange(6, 24)
        self._spin_label_size.setValue(
            int(self._settings.get('label_size', 10)))
        form.addRow("Axis label size:", self._spin_label_size)

        self._spin_tick_size = QSpinBox()
        self._spin_tick_size.setRange(6, 18)
        self._spin_tick_size.setValue(
            int(self._settings.get('tick_size', 9)))
        form.addRow("Tick label size:", self._spin_tick_size)

        self._spin_legend_size = QSpinBox()
        self._spin_legend_size.setRange(6, 18)
        self._spin_legend_size.setValue(
            int(self._settings.get('legend_size', 8)))
        form.addRow("Legend font size:", self._spin_legend_size)

        self._spin_marker_size = QSpinBox()
        self._spin_marker_size.setRange(2, 120)
        self._spin_marker_size.setValue(
            int(self._settings.get('marker_size', 24)))
        form.addRow("Marker size:", self._spin_marker_size)

        self._spin_alpha = QDoubleSpinBox()
        self._spin_alpha.setRange(0.05, 1.0)
        self._spin_alpha.setSingleStep(0.05)
        self._spin_alpha.setDecimals(2)
        self._spin_alpha.setValue(
            float(self._settings.get('marker_alpha', 0.7)))
        form.addRow("Marker alpha:", self._spin_alpha)

        self._spin_fig_w = QDoubleSpinBox()
        self._spin_fig_w.setRange(3.0, 20.0)
        self._spin_fig_w.setSingleStep(0.5)
        self._spin_fig_w.setDecimals(1)
        self._spin_fig_w.setValue(
            float(self._settings.get('fig_width_in', 8.0)))
        form.addRow("Figure width (in):", self._spin_fig_w)

        self._spin_fig_h = QDoubleSpinBox()
        self._spin_fig_h.setRange(2.0, 20.0)
        self._spin_fig_h.setSingleStep(0.5)
        self._spin_fig_h.setDecimals(1)
        self._spin_fig_h.setValue(
            float(self._settings.get('fig_height_in', 6.0)))
        form.addRow("Figure height (in):", self._spin_fig_h)

        self._spin_dpi = QSpinBox()
        self._spin_dpi.setRange(50, 600)
        self._spin_dpi.setSingleStep(50)
        self._spin_dpi.setValue(int(self._settings.get('save_dpi', 150)))
        form.addRow("Save DPI:", self._spin_dpi)

        lay.addLayout(form)

        # ボタン
        btn_row = QHBoxLayout()
        btn_reset = QPushButton("Reset to defaults")
        btn_reset.clicked.connect(self._on_reset)
        btn_row.addWidget(btn_reset)
        btn_row.addStretch()
        btn_cancel = QPushButton("Cancel")
        btn_cancel.clicked.connect(self.reject)
        btn_row.addWidget(btn_cancel)
        btn_ok = QPushButton("OK")
        btn_ok.setDefault(True)
        btn_ok.clicked.connect(self.accept)
        btn_row.addWidget(btn_ok)
        lay.addLayout(btn_row)

    def _on_reset(self):
        self._le_title.setText('')
        self._chk_show_title.setChecked(True)
        self._spin_title_size.setValue(11)
        self._spin_label_size.setValue(10)
        self._spin_tick_size.setValue(9)
        self._spin_legend_size.setValue(8)
        self._spin_marker_size.setValue(24)
        self._spin_alpha.setValue(0.7)
        self._spin_fig_w.setValue(8.0)
        self._spin_fig_h.setValue(6.0)
        self._spin_dpi.setValue(150)

    def result_settings(self) -> dict:
        return {
            'title': self._le_title.text(),
            'show_title': self._chk_show_title.isChecked(),
            'title_size': self._spin_title_size.value(),
            'label_size': self._spin_label_size.value(),
            'tick_size': self._spin_tick_size.value(),
            'legend_size': self._spin_legend_size.value(),
            'marker_size': self._spin_marker_size.value(),
            'marker_alpha': self._spin_alpha.value(),
            'fig_width_in': self._spin_fig_w.value(),
            'fig_height_in': self._spin_fig_h.value(),
            'save_dpi': self._spin_dpi.value(),
        }


class JointPlotDialog(QDialog):
    """Joint Plot / Scatter Preview ダイアログ。

    モード: Pos / Neg / Merged (両モードを同一軸に重ね合わせ)
    オプション:
      - include_histograms=True : 散布図 + RT/m/z 周辺ヒストグラム
      - include_histograms=False: 散布図のみ
    描画はクラス色のアノテーションを反映(host._class_colors を使用)。
    """

    def __init__(self, host, parent=None, include_histograms: bool = True,
                 initial_mode: str = "pos", curated: bool = False):
        super().__init__(parent or host)
        self._include_hist = bool(include_histograms)
        # Curated view フラグ。host._quant_ion_choices
        # に基づきクラス×モードでフィルタする。
        self._curated = bool(curated)
        base_title = ("Joint Plot — Preview & Save" if self._include_hist
                      else "Scatter Preview — Preview & Save")
        title = (f"Curated view — {base_title}"
                 if self._curated else base_title)
        self.setWindowTitle(title)
        self.resize(900, 720)
        self._host = host

        lay = QVBoxLayout(self)

        # 上部: モード選択 + Refresh + Save
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Mode:"))
        self._radio_pos = QRadioButton("Pos")
        self._radio_neg = QRadioButton("Neg")
        self._radio_merged = QRadioButton("Merged")
        if initial_mode == 'neg':
            self._radio_neg.setChecked(True)
        elif initial_mode == 'merged':
            self._radio_merged.setChecked(True)
        else:
            self._radio_pos.setChecked(True)
        self._radio_pos.toggled.connect(self._refresh)
        self._radio_neg.toggled.connect(self._refresh)
        self._radio_merged.toggled.connect(self._refresh)
        ctrl.addWidget(self._radio_pos)
        ctrl.addWidget(self._radio_neg)
        ctrl.addWidget(self._radio_merged)
        ctrl.addStretch()
        btn_refresh = QPushButton("Refresh")
        btn_refresh.setToolTip("Re-render with the latest peak data.")
        btn_refresh.clicked.connect(self._refresh)
        ctrl.addWidget(btn_refresh)
        # 図フォーマット設定
        btn_format = QPushButton("Format…")
        btn_format.setToolTip(
            "Customize figure format (title, fonts, sizes, DPI) before saving.")
        btn_format.clicked.connect(self._open_format)
        ctrl.addWidget(btn_format)
        btn_save = QPushButton("Save Figure PNG…")
        btn_save.setToolTip("Save the scatter (without legend) as PNG.")
        btn_save.clicked.connect(self._save)
        ctrl.addWidget(btn_save)
        # 凡例を別 PNG として保存
        btn_save_legend = QPushButton("Save Legend PNG…")
        btn_save_legend.setToolTip(
            "Save the legend (class color list) as a separate PNG.")
        btn_save_legend.clicked.connect(self._save_legend)
        ctrl.addWidget(btn_save_legend)
        lay.addLayout(ctrl)

        # フォーマット設定 + render モード強制フラグ
        self._format_settings: dict = dict(DEFAULT_FORMAT_SETTINGS)
        self._format_active: bool = False

        # キャンバス
        self._fig = Figure(figsize=(8, 6), constrained_layout=False)
        self._canvas = FigureCanvas(self._fig)
        self._canvas.setMinimumHeight(520)
        lay.addWidget(self._canvas, stretch=1)
        # 凡例エントリ(順序付き class → color)
        self._legend_entries: dict = {}

        # 下部: Close
        bottom = QHBoxLayout()
        bottom.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        bottom.addWidget(btn_close)
        lay.addLayout(bottom)

        self._refresh()

    def _get_mode_data(self, mode: str):
        """指定モードの raw peak data + match_df を取得。

        class_colors は _PER_MODE_STATE_VARS に含まれていないため
        snap.get('_class_colors') は通常 None。そのため host._class_colors に
        フォールバック(active mode と同じ色を使用)。

        戻り値: (rt_arr, mz_arr, match_df, class_colors) または all-None。
        """
        host = self._host
        if mode == host.fe.ion_mode:
            rt_arr = host._rt
            mz_arr = host._mz
            match_df = host._match_df
        else:
            snap = host._mode_state.get(mode, {}) or {}
            rt_arr = snap.get('_rt')
            mz_arr = snap.get('_mz')
            match_df = snap.get('_match_df')
        # クラス色は両モード共通(host._class_colors を最終フォールバック)
        snap_cc = (host._mode_state.get(mode, {}) or {}).get('_class_colors')
        class_colors = snap_cc or (host._class_colors or {})
        return rt_arr, mz_arr, match_df, class_colors

    def _snapshot_host_axis(self, mode: str) -> bytes:
        """host の Analysis figure の対象 axis を PNG として
        スナップショットする(他軸を一時的に hide → savefig → 復元)。
        bbox_inches='tight' で hidden axis 領域を除いた範囲のみ出力する。"""
        host = self._host
        ax_pos = getattr(host, '_ax_sc_pos', None)
        ax_neg = getattr(host, '_ax_sc_neg', None)
        if mode == 'pos':
            keep, hide = ax_pos, ax_neg
        elif mode == 'neg':
            keep, hide = ax_neg, ax_pos
        else:
            return b''
        if keep is None:
            return b''
        was_visible = True
        if hide is not None:
            try:
                was_visible = hide.get_visible()
                hide.set_visible(False)
            except Exception:
                hide = None
        try:
            from io import BytesIO
            buf = BytesIO()
            host._fig.savefig(buf, format='png', dpi=150,
                              bbox_inches='tight')
            return buf.getvalue()
        except Exception:
            return b''
        finally:
            if hide is not None:
                try:
                    hide.set_visible(was_visible)
                    host._canvas.draw_idle()
                except Exception:
                    pass

    def _build_legend_from_match_df(self, mode: str) -> dict:
        """snapshot モード用の凡例エントリを host._class_colors
        から構築。 mode==merged の場合は両モードの和集合。"""
        out: dict = {}
        host = self._host
        cc = host._class_colors or {}
        modes = ('pos', 'neg') if mode == 'merged' else (mode,)
        for m in modes:
            if m == host.fe.ion_mode:
                match_df = host._match_df
            else:
                snap = host._mode_state.get(m, {}) or {}
                match_df = snap.get('_match_df')
            if match_df is None or match_df.empty:
                continue
            mm = match_df[match_df['matched']]
            if 'final_status' in mm.columns:
                mm = mm[mm['final_status'] == STATUS_KEPT]
            if mm.empty:
                continue
            for cls in mm['lipid_class'].dropna().unique():
                key = str(cls)
                if key not in out and cls in cc:
                    out[key] = cc[cls]
        return out

    def _render_scatter(self, ax, mode: str, alpha_mul: float = 1.0,
                         show_background: bool = True):
        """クラス色アノテーションを反映した散布図を ax に描画する。
        凡例エントリを self._legend_entries に蓄積。
        Curated view では host._quant_ion_choices[cls] == mode
        のクラスのみを overlay。背景の raw peaks は変更しない(視覚コンテキスト)。
        render 中の (n_shown, n_total) を self._last_render_counts
        に保存。タイトル生成で使う(Curated/Current 両方で件数を明示)。
        show_background=False で not matched 灰色背景を抑制。
        Merged view で pos/neg 双方が呼ばれて二重描画される問題対策。
        戻り値: rt_arr / mz_arr のタプル (周辺ヒストグラム用) または None。"""
        rt_arr, mz_arr, match_df, class_colors = self._get_mode_data(mode)
        if rt_arr is None or mz_arr is None or len(rt_arr) == 0:
            return None

        # 背景: 全 peak を淡いグレーで(Not matched 含む)
        # Curated view では not matched 背景を非表示にする
        # show_background=False (Merged view から)でも非表示
        if not getattr(self, '_curated', False) and show_background:
            ax.scatter(rt_arr, mz_arr, s=4, alpha=0.25 * alpha_mul,
                       color='#B4B2A9', edgecolors='none', zorder=1)

        # ライブラリ overlay: matched & kept だけクラス色で描画
        # 件数を記録(n_total = matched & kept 全件、
        # n_shown = curated filter 適用後 / 非 curated は同じ値)
        n_total = 0
        n_shown = 0
        if match_df is not None and not match_df.empty:
            m_full = match_df[match_df['matched']]
            if 'final_status' in m_full.columns:
                m_full = m_full[m_full['final_status'] == STATUS_KEPT]
            n_total = len(m_full)
            m = m_full
            # Curated view ではこのモードで採用された
            # クラスのみに絞り込む
            if self._curated and not m.empty:
                qi = (getattr(self._host, '_quant_ion_choices', None) or {})
                curated_classes = {
                    str(c) for c, choice in qi.items() if choice == mode
                }
                if curated_classes:
                    m = m[m['lipid_class'].astype(str).isin(curated_classes)]
                else:
                    m = m.iloc[0:0]  # empty
            n_shown = len(m)
            if not m.empty:
                for cls, grp in m.groupby('lipid_class'):
                    color = class_colors.get(cls, '#888888')
                    # format settings から marker size / alpha を取得
                    ms = int(self._format_settings.get('marker_size', 24))
                    ma_base = float(
                        self._format_settings.get('marker_alpha', 0.7))
                    ax.scatter(
                        grp['obs_rt'].values, grp['obs_mz'].values,
                        s=ms, color=color, alpha=ma_base * alpha_mul,
                        edgecolors='none', zorder=3, label=str(cls))
                    # 凡例エントリ蓄積(順序維持 dict, 同 class は最初の色を保持)
                    if cls not in self._legend_entries:
                        self._legend_entries[str(cls)] = color
        # 件数を保存(タイトル生成で参照)
        if not hasattr(self, '_last_render_counts'):
            self._last_render_counts = {}
        self._last_render_counts[mode] = (n_shown, n_total)
        return rt_arr, mz_arr

    def _refresh(self):
        """Pos / Neg / Merged で散布図を描画。
        Pos/Neg + Histogram OFF の場合は Analysis 図の対応 axis を
        スナップショットして WYSIWYG にする。それ以外は _render_scatter で再描画。
        """
        self._fig.clear()
        # 凡例エントリを毎回リセット
        self._legend_entries = {}
        # snapshot 用 PNG bytes をリセット
        self._snapshot_bytes = b''
        self._is_snapshot_mode = False
        if self._radio_pos.isChecked():
            target_mode = 'pos'
        elif self._radio_neg.isChecked():
            target_mode = 'neg'
        else:
            target_mode = 'merged'

        # WYSIWYG snapshot path
        # Pos/Neg かつ Histogram OFF, かつ format 設定が default なら
        # Analysis 図の対応 axis をそのままスナップショット表示する。
        # format 設定変更後は render モードに切替えて設定を反映。
        # Curated view は snapshot path を無効化(quant_ion 絞り込み
        # が host axis に反映されていないため、必ず re-render する)。
        if ((not self._include_hist) and target_mode in ('pos', 'neg')
                and not self._format_active
                and not getattr(self, '_curated', False)):
            png = self._snapshot_host_axis(target_mode)
            if png:
                self._snapshot_bytes = png
                self._is_snapshot_mode = True
                self._legend_entries = self._build_legend_from_match_df(
                    target_mode)
                try:
                    import matplotlib.image as mpimg
                    from io import BytesIO
                    img = mpimg.imread(BytesIO(png), format='png')
                    ax = self._fig.add_subplot(111)
                    ax.imshow(img)
                    ax.set_axis_off()
                    self._main_ax = ax  # 凡例 overlay を有効化
                    self._fig.tight_layout(pad=0.2)
                except Exception:
                    pass
                # snapshot 上にも凡例 overlay を描画
                self._draw_inline_legend()
                self._canvas.draw_idle()
                return

        # データ可用性チェック
        if target_mode == 'merged':
            pos_data = self._get_mode_data('pos')
            neg_data = self._get_mode_data('neg')
            has_data = ((pos_data[0] is not None and len(pos_data[0]) > 0)
                        or (neg_data[0] is not None and len(neg_data[0]) > 0))
        else:
            data = self._get_mode_data(target_mode)
            has_data = data[0] is not None and len(data[0]) > 0

        if not has_data:
            ax = self._fig.add_subplot(111)
            ax.text(0.5, 0.5,
                    f"No {target_mode} data loaded.",
                    ha='center', va='center',
                    transform=ax.transAxes, fontsize=11)
            ax.set_axis_off()
            self._canvas.draw_idle()
            return

        if self._include_hist:
            # 散布図 + 周辺ヒストグラム
            gs = self._fig.add_gridspec(
                2, 2,
                width_ratios=[5, 1], height_ratios=[1, 5],
                hspace=0.04, wspace=0.04,
                left=0.10, right=0.96, bottom=0.10, top=0.96,
            )
            ax_top = self._fig.add_subplot(gs[0, 0])
            ax_main = self._fig.add_subplot(gs[1, 0])
            ax_right = self._fig.add_subplot(gs[1, 1], sharey=ax_main)

            hist_data = []
            if target_mode == 'merged':
                # Merged では灰色背景を抑制(二重描画回避)
                d_pos = self._render_scatter(
                    ax_main, 'pos', alpha_mul=0.85, show_background=False)
                d_neg = self._render_scatter(
                    ax_main, 'neg', alpha_mul=0.85, show_background=False)
                if d_pos:
                    hist_data.append(d_pos)
                if d_neg:
                    hist_data.append(d_neg)
            else:
                d = self._render_scatter(ax_main, target_mode)
                if d:
                    hist_data.append(d)

            # format settings からフォントサイズを取得
            fs = self._format_settings
            lbl_sz = int(fs.get('label_size', 10))
            tick_sz = int(fs.get('tick_size', 9))
            title_sz = int(fs.get('title_size', 11))
            ax_main.set_xlabel('RT (min)', fontsize=lbl_sz)
            ax_main.set_ylabel('m/z', fontsize=lbl_sz)
            ax_main.tick_params(labelsize=tick_sz)
            # グリッド削除

            # ヒストグラム(全モードを連結)
            if hist_data:
                all_rt = np.concatenate([d[0] for d in hist_data])
                all_mz = np.concatenate([d[1] for d in hist_data])
                ax_top.hist(all_rt, bins=60, color='#3a6db0',
                            alpha=0.85, edgecolor='black', linewidth=0.4)
                ax_right.hist(all_mz, bins=60, color='#3a6db0',
                              alpha=0.85, edgecolor='black', linewidth=0.4,
                              orientation='horizontal')
            ax_top.tick_params(axis='x', labelbottom=False, labelsize=tick_sz)
            ax_top.set_ylabel('Count', fontsize=lbl_sz)
            ax_right.tick_params(
                axis='y', labelleft=False, labelsize=tick_sz)
            ax_right.set_xlabel('Count', fontsize=lbl_sz)

            # Joint plot のタイトルにも matched 数を含める
            if fs.get('show_title', True):
                title_txt = fs.get('title')
                if not title_txt:
                    counts = getattr(self, '_last_render_counts', None) or {}
                    if target_mode == 'merged':
                        p_shown, p_total = counts.get('pos', (0, 0))
                        n_shown, n_total = counts.get('neg', (0, 0))
                        sub = (f"{p_shown + n_shown} / {p_total + n_total} matched"
                               if (p_total + n_total) > 0 else "no data")
                        label = "Curated " if self._curated else ""
                        title_txt = (f"{label}Merged Joint Plot — {sub} "
                                     "+ marginals")
                    else:
                        shown, total = counts.get(target_mode, (0, 0))
                        if total > 0:
                            sub = (f"{shown} / {total} matched"
                                   if shown != total
                                   else f"{total} matched")
                        else:
                            sub = "no data"
                        label = "Curated " if self._curated else ""
                        title_txt = (f"{label}{target_mode.upper()} "
                                     f"Joint Plot — {sub} + marginals")
                self._fig.suptitle(title_txt, fontsize=title_sz)
            self._main_ax = ax_main
        else:
            # 散布図のみ
            ax = self._fig.add_subplot(111)
            if target_mode == 'merged':
                # Merged では灰色背景を抑制(二重描画回避)
                self._render_scatter(
                    ax, 'pos', alpha_mul=0.85, show_background=False)
                self._render_scatter(
                    ax, 'neg', alpha_mul=0.85, show_background=False)
            else:
                self._render_scatter(ax, target_mode)
            # format settings 反映
            fs = self._format_settings
            lbl_sz = int(fs.get('label_size', 10))
            tick_sz = int(fs.get('tick_size', 9))
            title_sz = int(fs.get('title_size', 11))
            ax.set_xlabel('RT (min)', fontsize=lbl_sz)
            ax.set_ylabel('m/z', fontsize=lbl_sz)
            ax.tick_params(labelsize=tick_sz)
            # グリッド削除
            # タイトルに matched 数を含める
            # Curated view は ⑦ で絞り込まれた件数 / 全 kept 件数を表示
            if fs.get('show_title', True):
                title_txt = fs.get('title')
                if not title_txt:
                    counts = getattr(self, '_last_render_counts', None) or {}
                    if target_mode == 'merged':
                        # pos + neg の合算
                        p_shown, p_total = counts.get('pos', (0, 0))
                        n_shown, n_total = counts.get('neg', (0, 0))
                        sub = (f"{p_shown + n_shown} / {p_total + n_total} matched"
                               if (p_total + n_total) > 0 else "no data")
                        label = "Curated " if self._curated else ""
                        title_txt = f"{label}Merged Scatter — {sub}"
                    else:
                        shown, total = counts.get(target_mode, (0, 0))
                        if total > 0:
                            sub = (f"{shown} / {total} matched"
                                   if shown != total
                                   else f"{total} matched")
                        else:
                            sub = "no data"
                        label = "Curated " if self._curated else ""
                        title_txt = (f"{label}{target_mode.upper()} "
                                     f"Scatter — {sub}")
                self._fig.suptitle(title_txt, fontsize=title_sz)
            self._fig.tight_layout()
            self._main_ax = ax

        # インライン凡例を main_ax 上に表示
        self._draw_inline_legend()
        self._canvas.draw_idle()

    def _draw_inline_legend(self):
        """蓄積した凡例エントリを main_ax の右上に描画。
        保存時には _save が一時的に非表示にできるよう legend artist を保持。"""
        from matplotlib.lines import Line2D
        ax = getattr(self, '_main_ax', None)
        if ax is None or not self._legend_entries:
            self._legend = None
            return
        handles = [
            Line2D([0], [0], marker='o', linestyle='',
                   markerfacecolor=col, markeredgecolor='none',
                   markersize=7, label=lbl)
            for lbl, col in self._legend_entries.items()
        ]
        self._legend = ax.legend(
            handles=handles, loc='upper right',
            frameon=True,
            fontsize=int(self._format_settings.get('legend_size', 8)),
            ncol=1,
            framealpha=0.85, borderpad=0.6,
            handletextpad=0.5, labelspacing=0.4)

    def _open_format(self):
        """Figure Format Settings ダイアログを開く。
        OK 時は設定を適用 + render モード強制 (Pos/Neg snapshot 無効)。
        Reset 時は default 設定に戻し、snapshot モードを復活させる。"""
        dlg = FormatSettingsDialog(self._format_settings, parent=self)
        if dlg.exec():
            new_settings = dlg.result_settings()
            # default と完全一致なら format_active を解除(snapshot に戻す)
            if all(new_settings.get(k) == DEFAULT_FORMAT_SETTINGS.get(k)
                   for k in DEFAULT_FORMAT_SETTINGS):
                self._format_settings = dict(DEFAULT_FORMAT_SETTINGS)
                self._format_active = False
            else:
                self._format_settings = new_settings
                self._format_active = True
            # figure サイズも反映
            try:
                w = float(self._format_settings.get('fig_width_in', 8.0))
                h = float(self._format_settings.get('fig_height_in', 6.0))
                self._fig.set_size_inches(w, h, forward=True)
            except Exception:
                pass
            self._refresh()

    def _save_legend(self):
        """凡例だけを別 PNG ファイルとして保存。"""
        from matplotlib.lines import Line2D
        if not self._legend_entries:
            QMessageBox.information(
                self, "No legend",
                "No class entries to put in the legend.\n"
                "Load library / Match first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Legend PNG",
            timestamped_filename("legend.png"), "PNG (*.png)")
        if not path:
            return
        try:
            n = len(self._legend_entries)
            # 凡例のみの図(高さは項目数に比例)
            legend_fig = Figure(figsize=(2.6, max(1.0, 0.3 * n + 0.4)))
            handles = [
                Line2D([0], [0], marker='o', linestyle='',
                       markerfacecolor=col, markeredgecolor='none',
                       markersize=8, label=lbl)
                for lbl, col in self._legend_entries.items()
            ]
            # format settings から legend font size を取得
            leg_sz = int(self._format_settings.get('legend_size', 10))
            dpi = int(self._format_settings.get('save_dpi', 200))
            legend_fig.legend(
                handles=handles, loc='center', frameon=False,
                ncol=1, fontsize=leg_sz,
                handletextpad=0.6, labelspacing=0.5)
            legend_fig.savefig(path, dpi=dpi, bbox_inches='tight')
            QMessageBox.information(self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))

    def _save(self):
        """図 PNG を保存。
        - Snapshot モード時は元 PNG bytes を直接書き出し(劣化なし)
        - それ以外は self._fig を savefig(凡例は一時非表示でクリーンな図に)
        - DPI は format_settings から取得"""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure PNG",
            timestamped_filename("scatter.png"), "PNG (*.png)")
        if not path:
            return
        try:
            if getattr(self, '_is_snapshot_mode', False) and self._snapshot_bytes:
                with open(path, 'wb') as f:
                    f.write(self._snapshot_bytes)
                QMessageBox.information(self, "Saved", f"Saved: {path}")
                return
            legend = getattr(self, '_legend', None)
            dpi = int(self._format_settings.get('save_dpi', 150))
            try:
                if legend is not None:
                    legend.set_visible(False)
                self._fig.savefig(path, dpi=dpi, bbox_inches='tight')
                QMessageBox.information(self, "Saved", f"Saved: {path}")
            finally:
                if legend is not None:
                    try:
                        legend.set_visible(True)
                        self._canvas.draw_idle()
                    except Exception:
                        pass
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))


class HeatmapPreviewDialog(QDialog):
    """Adduct heatmap (IS / Observed) のプレビュー + 保存。

    host._render_adduct_heatmap_png で PNG bytes にレンダリングし、
    QPixmap として表示。Save PNG ボタンで同じ bytes を書き出す。
    保存内容 = プレビュー内容 が保証される。
    """

    def __init__(self, host, view: dict, title: str,
                 default_filename: str, parent=None):
        super().__init__(parent or host)
        self.setWindowTitle(f"{title} — Preview & Save")
        self.resize(900, 720)
        self._host = host
        self._title = title
        self._default_filename = default_filename

        # PNG bytes にレンダリング
        from io import BytesIO
        buf = BytesIO()
        self._png_bytes: bytes = b''
        self._error: str = ''
        try:
            host._render_adduct_heatmap_png(view, buf, title=title)
            self._png_bytes = buf.getvalue()
        except Exception as exc:
            self._error = str(exc)

        lay = QVBoxLayout(self)

        # 上部: Save ボタン
        ctrl = QHBoxLayout()
        ctrl.addStretch()
        btn_save = QPushButton("Save PNG…")
        btn_save.setEnabled(bool(self._png_bytes))
        btn_save.clicked.connect(self._save)
        ctrl.addWidget(btn_save)
        lay.addLayout(ctrl)

        # プレビュー(QPixmap)
        if self._png_bytes:
            pix = QPixmap()
            pix.loadFromData(self._png_bytes, 'PNG')
            if pix.width() > 850:
                pix = pix.scaledToWidth(850, Qt.SmoothTransformation)
            lbl = QLabel()
            lbl.setPixmap(pix)
            lbl.setAlignment(Qt.AlignCenter)
            scroll = QScrollArea()
            scroll.setWidget(lbl)
            scroll.setWidgetResizable(False)
            scroll.setAlignment(Qt.AlignCenter)
            lay.addWidget(scroll, stretch=1)
        else:
            err_lbl = QLabel(
                f"Failed to render heatmap:\n{self._error}")
            err_lbl.setAlignment(Qt.AlignCenter)
            lay.addWidget(err_lbl, stretch=1)

        # 下部: Close
        bottom = QHBoxLayout()
        bottom.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        bottom.addWidget(btn_close)
        lay.addLayout(bottom)

    def _save(self):
        if not self._png_bytes:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, f"Save {self._title} PNG",
            timestamped_filename(self._default_filename), "PNG (*.png)")
        if not path:
            return
        try:
            with open(path, 'wb') as f:
                f.write(self._png_bytes)
            QMessageBox.information(self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))


class CombinedPreviewDialog(QDialog):
    """Merged (combined pos+neg) scatter のプレビュー + 保存。

    host._fig(Analysis タブの両モード散布図)を PNG にレンダリングし、
    QPixmap として表示。Save PNG ボタンで同じ PNG を保存する。
    保存内容 = プレビュー内容 が保証される(同じ bytes)。
    """

    def __init__(self, host, parent=None):
        super().__init__(parent or host)
        self.setWindowTitle("Combined (Merged) Scatter — Preview & Save")
        self.resize(900, 720)
        self._host = host

        # host._fig を PNG bytes にレンダリング
        from io import BytesIO
        buf = BytesIO()
        self._png_bytes: bytes = b''
        self._error: str = ''
        try:
            host._fig.savefig(buf, format='png', dpi=150,
                              bbox_inches='tight')
            self._png_bytes = buf.getvalue()
        except Exception as exc:
            self._error = str(exc)

        lay = QVBoxLayout(self)

        # 上部: Save / Close ボタン
        ctrl = QHBoxLayout()
        ctrl.addStretch()
        btn_save = QPushButton("Save PNG…")
        btn_save.setEnabled(bool(self._png_bytes))
        btn_save.clicked.connect(self._save)
        ctrl.addWidget(btn_save)
        lay.addLayout(ctrl)

        # プレビュー(QPixmap を QLabel に表示)
        if self._png_bytes:
            pix = QPixmap()
            pix.loadFromData(self._png_bytes, 'PNG')
            if pix.width() > 850:
                pix = pix.scaledToWidth(850, Qt.SmoothTransformation)
            lbl = QLabel()
            lbl.setPixmap(pix)
            lbl.setAlignment(Qt.AlignCenter)
            scroll = QScrollArea()
            scroll.setWidget(lbl)
            scroll.setWidgetResizable(False)
            scroll.setAlignment(Qt.AlignCenter)
            lay.addWidget(scroll, stretch=1)
        else:
            err_lbl = QLabel(
                f"Failed to render combined scatter:\n{self._error}")
            err_lbl.setAlignment(Qt.AlignCenter)
            lay.addWidget(err_lbl, stretch=1)

        # 下部: Close
        bottom = QHBoxLayout()
        bottom.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        bottom.addWidget(btn_close)
        lay.addLayout(bottom)

    def _save(self):
        if not self._png_bytes:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Combined scatter PNG",
            timestamped_filename("scatter_combined.png"), "PNG (*.png)")
        if not path:
            return
        try:
            with open(path, 'wb') as f:
                f.write(self._png_bytes)
            QMessageBox.information(self, "Saved", f"Saved: {path}")
        except Exception as e:
            QMessageBox.warning(self, "Save failed", str(e))


class ParametersDialog(QDialog):
    """全ステップのパラメータを 1 個のダイアログで編集。

    Host(PreviewRTDialog)が保持する spinbox 群を参照するラベル + 入力欄
    を表示する。値は host._spin_xxx に直接反映される。
    """

    def __init__(self, host, parent=None):
        super().__init__(parent)
        # Parameters → Advanced + タブ構造化
        self.setWindowTitle("Advanced Setting")
        self.resize(640, 540)
        self._host = host

        lay = QVBoxLayout(self)
        intro = QLabel(
            "Edit step parameters and access step-specific reports.\n"
            "Use the tabs to navigate between sections.")
        intro.setWordWrap(True)
        intro.setStyleSheet("color: #555; padding: 4px;")
        lay.addWidget(intro)

        # ── タブ構造 ──────────────────────────────────────────────
        tabs = QTabWidget()
        lay.addWidget(tabs, stretch=1)

        # === Tab 1: Parameters ===
        params_widget = QWidget()
        params_lay = QVBoxLayout(params_widget)
        form = QFormLayout()
        form.setContentsMargins(8, 8, 8, 8)
        form.setSpacing(6)

        # 各 spinbox を Host から借りてここに表示
        sections = [
            ("① Match Overlay", [
                ("ppm tol:", host._spin_ppm),
                ("Void volume cutoff:", host._spin_void_rt),
            ]),
            ("③ IS Filter", [
                ("IS tol (min):", host._spin_is_tol),
                ("Auto-pick:", host._cb_is_auto_intensity),
            ]),
            ("④ RT Outlier Filter", [
                ("MAD\u00d7:", host._spin_iqr),
            ]),
            ("⑤ Adduct Ion Filter", [
                ("ppm:", host._spin_adduct_ppm),
                ("RT (min):", host._spin_adduct_rt),
                ("int floor:", host._spin_adduct_int),
                ("vote threshold:", host._spin_adduct_vote),
            ]),
            ("⑥ Coherence Filter", [
                ("\u03c3:", host._spin_sigma),
                ("Loose attribution:", host._cb_loose_coherence),
                # 位置異性体が分離するクラスを 2 系列で扱う
                ("Bimodal classes:", host._ed_two_series),
                ("Quant export:", host._cb_sum_two_series),
            ]),
        ]

        self._reparented_widgets = []
        for sec_title, rows in sections:
            sec_lbl = QLabel(sec_title)
            sec_lbl.setStyleSheet(
                "font-weight: bold; color: #333; padding-top: 6px;")
            form.addRow(sec_lbl)
            for label_txt, widget in rows:
                # 一時的に widget をこのダイアログ配下に reparent して表示
                widget.setParent(self)
                widget.setVisible(True)
                self._reparented_widgets.append(widget)
                form.addRow(label_txt, widget)

        params_lay.addLayout(form)
        params_lay.addStretch()
        tabs.addTab(params_widget, "Parameters")

        # Report 系のタブは Utility タブに移行(削除)
        # 残すのは Parameters + Patterns(⑤)+ Pairs(⑥)のみ。

        # === ⑤ Adduct Ion Filter タブは廃止 ===
        # 同等機能(期待アダクトパターン heatmap)は
        # Utility > Heatmap > "IS adduct ion patterns" で提供されている。

        # === Tab 2: ⑥ Coherence Filter (Pairs のみ) ===
        co_widget = QWidget()
        co_lay = QVBoxLayout(co_widget)
        co_lay.addWidget(QLabel(
            "Configure target pairs for ⑥ Coherence Filter.\n"
            "(Reports moved to the Utility tab.)"))
        btn_co_pairs = QPushButton("Pairs…  (Select target pairs)")
        btn_co_pairs.clicked.connect(
            lambda: self._close_and_invoke(host._open_coherence_pairs))
        co_lay.addWidget(btn_co_pairs)
        co_lay.addStretch()
        tabs.addTab(co_widget, "⑥ Coherence Filter")

        # === Tab 4: Filter Classes ===
        self._filter_class_checks: dict[str, dict[str, QCheckBox]] = {
            'is': {}, 'adduct': {}, 'rt_outlier': {}, 'coherence': {},
        }
        fc_widget = QWidget()
        fc_lay = QVBoxLayout(fc_widget)
        fc_intro = QLabel(
            "Per-class filter application. Uncheck a lipid class to skip "
            "that filter for it (the filter's rejection flag will not be "
            "applied to rows of that class).")
        fc_intro.setWordWrap(True)
        fc_intro.setStyleSheet("color: #555; padding: 2px;")
        fc_lay.addWidget(fc_intro)
        fc_sub_tabs = QTabWidget()
        all_classes = self._collect_all_classes()
        for _fkey, _flabel in (
            ('is', "③ IS Filter"),
            ('rt_outlier', "④ RT Outlier Filter"),
            ('adduct', "⑤ Adduct Ion Filter"),
            ('coherence', "⑥ Coherence Filter"),
        ):
            sub = QWidget()
            sub_lay = QVBoxLayout(sub)
            exempt_set = host._filter_class_exemptions.get(_fkey, set())
            btn_row = QHBoxLayout()
            btn_all = QPushButton("All")
            btn_none = QPushButton("None")
            btn_all.clicked.connect(
                lambda _checked=False, k=_fkey: self._fc_set_all(k, True))
            btn_none.clicked.connect(
                lambda _checked=False, k=_fkey: self._fc_set_all(k, False))
            btn_row.addWidget(btn_all)
            btn_row.addWidget(btn_none)
            btn_row.addStretch()
            sub_lay.addLayout(btn_row)
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            inner = QWidget()
            grid = QGridLayout(inner)
            grid.setContentsMargins(4, 4, 4, 4)
            grid.setSpacing(4)
            n_cols = 3
            for i, cls in enumerate(all_classes):
                cb = QCheckBox(cls)
                cb.setChecked(cls not in exempt_set)
                self._filter_class_checks[_fkey][cls] = cb
                grid.addWidget(cb, i // n_cols, i % n_cols)
            scroll.setWidget(inner)
            sub_lay.addWidget(scroll, stretch=1)
            fc_sub_tabs.addTab(sub, _flabel)
        fc_lay.addWidget(fc_sub_tabs, stretch=1)
        tabs.addTab(fc_widget, "Filter Classes")

        # 閉じる + Reset ボタン
        btns = QHBoxLayout()
        btn_reset = QPushButton("Reset to defaults")
        btn_reset.clicked.connect(self._on_reset)
        btns.addWidget(btn_reset)
        btns.addStretch()
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.accept)
        btns.addWidget(btn_close)
        lay.addLayout(btns)

    def _close_and_invoke(self, callable_):
        """sub-dialog を開く前に Advanced を閉じる。

        二重ダイアログ(Advanced + Report 等)では行クリック時の
        散布図ハイライトが背面の host ウィンドウで隠れてしまうため、
        sub-dialog 起動時に Advanced 自体を閉じて単一ダイアログ表示にする。
        QTimer.singleShot で Advanced の close 後に host メソッドを
        呼ぶことで、widget reparent / event loop の整合を取る。"""
        self.accept()
        QTimer.singleShot(0, callable_)

    def _on_reset(self):
        self._host._utility_reset_params()

    # Filter Classes helpers
    def _collect_all_classes(self) -> list[str]:
        host = self._host
        classes: set = set()
        if host._match_df is not None and not host._match_df.empty:
            try:
                classes.update(
                    host._match_df['lipid_class'].dropna().unique().tolist())
            except Exception:
                pass
        if not classes and host._lib_path:
            try:
                lib_df = load_lipid_library(host._lib_path, host.fe.ion_mode)
                classes.update(
                    lib_df['lipid_class'].dropna().unique().tolist())
            except Exception:
                pass
        return sorted(classes)

    def _fc_set_all(self, key: str, state: bool) -> None:
        d = self._filter_class_checks.get(key) or {}
        for cb in d.values():
            cb.setChecked(state)

    def _save_filter_class_exemptions(self) -> None:
        try:
            host = self._host
            for key, d in self._filter_class_checks.items():
                exempt = {cls for cls, cb in d.items() if not cb.isChecked()}
                host._filter_class_exemptions[key] = exempt
        except Exception as e:
            log.warning(f"[Filter Classes] save failed: {e}")

    def _return_widgets_to_host(self):
        """借用した spinbox を Host の _params_holder に返却。
        accept/reject/X いずれの経路でも呼ばれるように done() に集約。
        冪等(複数回呼ばれても害無し)。"""
        try:
            holder = getattr(self._host, '_params_holder', None)
            if holder is None:
                return
            for w in list(self._reparented_widgets):
                try:
                    # 既に削除済みの場合は AttributeError や RuntimeError
                    w.setParent(holder)
                    w.setVisible(False)
                except Exception:
                    pass
            self._reparented_widgets = []
        except Exception:
            pass

    def done(self, result):
        """accept()/reject() のいずれでも widget 返却を実行。
        closeEvent は X ボタン経由でしか発火しないため、ここに集約する。"""
        # Filter Classes exemption set を host に保存
        try:
            self._save_filter_class_exemptions()
        except Exception:
            pass
        self._return_widgets_to_host()
        super().done(result)

    def closeEvent(self, ev):
        """保険として残す(done() で先に返却済みなら no-op)。"""
        try:
            self._save_filter_class_exemptions()
        except Exception:
            pass
        self._return_widgets_to_host()
        super().closeEvent(ev)


class _EICPanelMixin:
    """EIC(抽出イオンクロマトグラム)パネルの共通部品。

    fix26 で CandidatePickerDialog に入れたものを、散布図の右クリックから
    開く SpotEICDialog と共用するために括り出した。丸ごとコピーすると
    FileEntry / StubFileEntry の二重定義と同じ轍を踏むため。

    使う側に用意しておくもの:
      self._raw_paths        : サンプル順の mzML 実パス(空なら EIC 不可)
      self._sample_col_names : サンプル名(表示用)
      self._theor_mz         : EIC を引く m/z
      self._eic_df / _eic_proc / _eic_timer / _eic_out_csv /
      self._eic_t0 / _eic_sample_name : 状態(いずれも None / 0 / "" 初期化)

    差し替えるフック:
      _eic_markers()            -> [(rt, ラベル), ...] 縦線を引く位置
      _eic_selected_marker()    -> 強調する marker の index(無ければ -1)
      _eic_on_marker_clicked(i) -> 図をクリックしたときの処理

    非ブロッキングにしてあるのは、mzML 1 本の EIC 抽出に 15〜30 秒
    かかるため。subprocess.run で待つとその間 UI が固まる。
    """

    # ── フック(既定は「マーカー無し・選択なし」)──────────────────
    def _eic_markers(self) -> list:
        return []

    def _eic_selected_marker(self) -> int:
        return -1

    def _eic_on_marker_clicked(self, index: int):
        self._draw_eic()

    # ── EIC パネル ────────────────────────────────────────
    #  候補どうしの m/z 差は非常に小さい(検証データの TG 48:1 D7 で
    #  0.0004 Da)ので、理論 m/z ± EIC_MZ_TOL の EIC を 1 本引けば
    #  全候補が同じトレースに写る。候補ごとに引き直す必要はない。
    EIC_MZ_TOL_DA = 0.005          # ① タブ Step 2 の EIC 既定と同じ
    # status 行の末尾に出す案内。選ぶものが無い画面では空にする。
    EIC_CLICK_HINT = "click a peak to select it"

    EIC_TIMEOUT_S = 180.0

    EIC_PAD_MIN = 0.30             # 候補 RT の外側にとる余白

    def _build_eic_panel(self, lay):
        box = QGroupBox("Chromatogram (EIC)")
        v = QVBoxLayout(box)
        v.setContentsMargins(6, 6, 6, 4)

        self._eic_status = QLabel("")
        self._eic_status.setStyleSheet("color:#555; font-size:11px;")
        self._eic_status.setWordWrap(True)
        v.addWidget(self._eic_status)

        self._eic_fig = Figure(figsize=(8.0, 2.8))
        self._eic_canvas = FigureCanvas(self._eic_fig)
        self._eic_ax = self._eic_fig.add_subplot(111)
        self._eic_ax.set_axis_off()
        self._eic_canvas.setMinimumHeight(230)
        v.addWidget(self._eic_canvas, stretch=1)

        row = QHBoxLayout()
        self._eic_cb_full = QCheckBox("Full RT range")
        self._eic_cb_full.setToolTip(
            "By default only the region around the candidate RT is shown.\n"
            "Check this to show the full RT range.")
        self._eic_cb_full.toggled.connect(lambda *_: self._draw_eic())
        row.addWidget(self._eic_cb_full)
        self._eic_btn = QPushButton("Extract EIC")
        self._eic_btn.clicked.connect(self._on_eic_button)
        row.addWidget(self._eic_btn)
        row.addStretch()
        v.addLayout(row)

        lay.addWidget(box, stretch=1)
        self._eic_canvas.mpl_connect('button_press_event', self._on_eic_click)

        if not self._raw_paths:
            self._eic_status.setText(
                "No raw data (mzML) was found, so no EIC can be shown. "
                "It happens with an imported external table, or when the mzML "
                "files were moved. Judge from RT, intensity and the sample profile.")
            self._eic_btn.setEnabled(False)

    def _pick_eic_sample_index(self) -> int:
        """EIC を出すサンプルを選ぶ。

        「どれかの候補が最も強く出ているサンプル」を採る。弱いサンプルを
        引くとピーク形状が読めず、選択の判断材料にならないため。
        """
        best_j, best_v = 0, -1.0
        for c in self._candidates:
            for j, v in enumerate(c.get('sample_intensities') or []):
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    continue
                if fv > best_v:
                    best_v, best_j = fv, j
        return best_j

    def _on_eic_button(self):
        if self._eic_proc is not None:
            self._stop_eic("Cancelled.")
        else:
            self._start_eic()

    def _start_eic(self):
        if self._eic_proc is not None or not self._raw_paths:
            return
        worker = resolve_align_worker()
        if worker is None:
            self._eic_status.setText(
                f"{ALIGN_WORKER_NAME} not found, so no EIC can be shown.")
            self._eic_btn.setEnabled(False)
            return
        j = self._pick_eic_sample_index()
        if j >= len(self._raw_paths):
            j = 0
        path = self._raw_paths[j]
        self._eic_sample_name = (self._sample_col_names[j]
                                 if j < len(self._sample_col_names)
                                 else Path(path).stem)
        try:
            with _tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as tmp:
                self._eic_out_csv = tmp.name
            cmd = [sys.executable, str(worker), "--eic",
                   "--mzml", str(path), "--out", self._eic_out_csv,
                   "--eic-mz", f"{self._theor_mz:.6f}",
                   "--eic-mz-tol", f"{self.EIC_MZ_TOL_DA:.6f}"]
            self._eic_proc = _subprocess.Popen(
                cmd, stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
                text=True)
        except Exception as e:
            self._eic_status.setText(f"Could not start EIC extraction: {e}")
            self._eic_proc = None
            return
        self._eic_t0 = _time.time()
        self._eic_btn.setText("Cancel")
        self._eic_status.setText(
            f"Extracting EIC from {self._eic_sample_name}… (15-30 sec per mzML)")
        self._eic_timer = QTimer(self)
        self._eic_timer.timeout.connect(self._poll_eic)
        self._eic_timer.start(300)

    def _poll_eic(self):
        p = self._eic_proc
        if p is None:
            return
        el = _time.time() - self._eic_t0
        if p.poll() is None:
            if el > self.EIC_TIMEOUT_S:
                self._stop_eic(f"EIC extraction exceeded {self.EIC_TIMEOUT_S:.0f} sec "
                               "and was aborted.")
                return
            self._eic_status.setText(
                f"Extracting EIC from {self._eic_sample_name}… {el:.0f} sec")
            return
        rc = p.returncode
        err = ""
        try:
            _o, err = p.communicate(timeout=5)
        except Exception:
            pass
        self._eic_timer.stop()
        self._eic_timer = None
        self._eic_proc = None
        self._eic_btn.setText("Re-extract EIC")
        if rc != 0:
            self._eic_status.setText(
                f"EIC extraction failed (exit={rc}): {(err or '')[-200:]}")
            self._cleanup_eic_tmp()
            return
        try:
            self._eic_df = pd.read_csv(self._eic_out_csv)
        except Exception as e:
            self._eic_status.setText(f"Could not read the EIC: {e}")
            self._cleanup_eic_tmp()
            return
        self._cleanup_eic_tmp()
        self._eic_status.setText(
            f"{self._eic_sample_name}   |   m/z {self._theor_mz:.4f} "
            f"± {self.EIC_MZ_TOL_DA:.4f} Da"
            + (f"   |   {self.EIC_CLICK_HINT}" if self.EIC_CLICK_HINT else ""))
        self._draw_eic()

    def _stop_eic(self, msg: str):
        if self._eic_timer is not None:
            self._eic_timer.stop()
            self._eic_timer = None
        p, self._eic_proc = self._eic_proc, None
        if p is not None:
            try:
                p.kill()
            except Exception:
                pass
        self._cleanup_eic_tmp()
        self._eic_btn.setText("Extract EIC")
        if msg:
            self._eic_status.setText(msg)

    def _cleanup_eic_tmp(self):
        if self._eic_out_csv:
            try:
                os.unlink(self._eic_out_csv)
            except Exception:
                pass
            self._eic_out_csv = None

    def _draw_eic(self):
        ax = self._eic_ax
        ax.clear()
        if self._eic_df is None or self._eic_df.empty:
            ax.set_axis_off()
            self._eic_canvas.draw_idle()
            return
        ax.set_axis_on()
        x = self._eic_df['rt_min'].values
        y = self._eic_df['intensity'].values
        ax.plot(x, y, lw=1.0, color='#3a6fa0')
        ax.fill_between(x, 0, y, alpha=0.18, color='#3a6fa0')

        markers = self._eic_markers()
        sel = self._eic_selected_marker()
        for i, (rt, lbl) in enumerate(markers):
            rt = float(rt)
            on = (i == sel)
            ax.axvline(rt, color=('#c0392b' if on else '#7f8c8d'),
                       lw=(1.8 if on else 1.0),
                       ls=('-' if on else '--'), zorder=3)
            ax.annotate(
                str(lbl), xy=(rt, 1.0), xycoords=('data', 'axes fraction'),
                xytext=(0, -13), textcoords='offset points',
                ha='center', va='center', fontsize=8, zorder=4,
                color=('white' if on else '#333'),
                bbox=dict(boxstyle='round,pad=0.25',
                          fc=('#c0392b' if on else '#ecf0f1'), ec='none'))

        if not self._eic_cb_full.isChecked() and markers:
            rts = [float(r) for r, _l in markers]
            lo, hi = min(rts), max(rts)
            pad = max(self.EIC_PAD_MIN, (hi - lo) * 0.8)
            ax.set_xlim(lo - pad, hi + pad)
            m = ((x >= lo - pad) & (x <= hi + pad))
            if m.any() and float(y[m].max()) > 0:
                ax.set_ylim(0, float(y[m].max()) * 1.12)
        ax.set_xlabel('RT (min)', fontsize=8)
        ax.set_ylabel('Intensity', fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(alpha=0.3)
        try:
            self._eic_fig.tight_layout()
        except Exception:
            pass
        self._eic_canvas.draw_idle()

    def _on_eic_click(self, event):
        """図をクリックしたら最も近い候補を選ぶ。"""
        if event.inaxes is not self._eic_ax or event.xdata is None:
            return
        markers = self._eic_markers()
        if not markers:
            return
        rts = [float(r) for r, _l in markers]
        i = min(range(len(rts)), key=lambda k: abs(rts[k] - float(event.xdata)))
        self._eic_on_marker_clicked(i)


class CandidatePickerDialog(_EICPanelMixin, QDialog):
    """IS の ppm tol 内の全候補ピークを表示し、ユーザーに正しい候補を
    選ばせるダイアログ。

    背景:
      混雑する m/z 領域では同じ理論 m/z に複数のピークが該当しうる。
      MS-DIAL 出力の精度限界(2〜4 桁)も相まって、m/z minimum-error
      自動マッチングは誤った候補を選んでしまうリスクがある。
      IS は LipidQuant パイプライン全体の根幹となる指標のため、
      ユーザーに明示的に選択させる仕組みが必要。

    入力:
      title: ダイアログタイトル(例: "DG 33:1 D7 candidates (5)")
      theoretical_mz: IS の理論 m/z(全候補で同じ値)
      candidates: list of dict(_collect_match_candidates の戻り値)
                  各 dict に peak_idx, obs_rt, obs_mz, delta_ppm,
                  mean_intensity, sample_intensities が含まれる
      sample_col_names: サンプル列名(プロファイル列ヘッダ用)
      initial_peak_idx: 初期選択する候補の peak_idx
                        None なら先頭(m/z 誤差最小)が選ばれる
    """

    # スパークライン用の Unicode block 文字(8 段階)
    _SPARK_CHARS = "▁▂▃▄▅▆▇█"

    def __init__(
        self,
        title: str,
        theoretical_mz: float,
        candidates: list[dict],
        sample_col_names: list[str] | None = None,
        initial_peak_idx: int | None = None,
        raw_paths: "list[str] | None" = None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(title)
        # 下に EIC パネルが付いた分だけ高くする
        self.resize(880, 760)
        # 表示用に RT 昇順でソート(同 RT は m/z 誤差絶対値)
        # 元の m/z 誤差順は match_df の isomer_rank に使われるため、ここで
        # 表示用にだけソートする。peak_idx を介して selection は維持される。
        self._candidates = sorted(
            list(candidates),
            key=lambda c: (
                float(c.get('obs_rt', float('inf'))),
                abs(float(c.get('delta_ppm', 0.0))),
            ),
        )
        self._theor_mz = float(theoretical_mz)
        self._sample_col_names = sample_col_names or []
        self._selected_peak_idx: int | None = None
        self._show_sample_profile = False
        # EIC 用
        self._raw_paths = [str(p) for p in (raw_paths or [])]
        self._eic_df = None
        self._eic_proc = None
        self._eic_timer = None
        self._eic_out_csv = None
        self._eic_t0 = 0.0
        self._eic_sample_name = ""

        lay = QVBoxLayout(self)

        # ── ヘッダ説明 ─────────────────────────────────────────────
        info = QLabel(
            f"Theoretical m/z: <b>{theoretical_mz:.4f}</b>"
            f"   /   {len(candidates)} candidate(s) within ppm tolerance.<br>"
            "Pick the correct IS peak below. RT and intensity are the "
            "primary differentiators when m/z values overlap."
        )
        info.setStyleSheet("color:#444; font-size:11px; padding:4px;")
        info.setTextFormat(Qt.RichText)
        lay.addWidget(info)

        # ── オプション: Sample profile 表示トグル ───────────────────
        opts_row = QHBoxLayout()
        self._cb_profile = QCheckBox("Show sample profile (sparkline)")
        self._cb_profile.setChecked(False)
        self._cb_profile.setToolTip(
            "When checked, an extra column shows a Unicode sparkline of\n"
            "each candidate's per-sample intensity profile. Useful when\n"
            "RT and m/z alone do not differentiate candidates.")
        self._cb_profile.toggled.connect(self._on_toggle_profile)
        opts_row.addWidget(self._cb_profile)
        opts_row.addStretch()
        lay.addLayout(opts_row)

        # ── テーブル ───────────────────────────────────────────────
        # 列: RT / obs m/z / theor m/z / Δm/z / Δppm / mean intensity
        #     / [Sample profile (optional)] / Selection
        self._table = QTableWidget()
        self._build_table_columns()
        self._table.setEditTriggers(QTableWidget.NoEditTriggers)
        self._table.setSelectionBehavior(QTableWidget.SelectRows)
        self._table.setAlternatingRowColors(True)
        lay.addWidget(self._table, stretch=0)

        # ── EIC パネル─────────────────────────────────────
        self._build_eic_panel(lay)

        # ── ボタンボックス ─────────────────────────────────────────
        btns = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._on_ok)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

        # ── 行をビルド ─────────────────────────────────────────────
        self._radio_buttons: list[QRadioButton] = []
        self._populate_rows(initial_peak_idx)

        # 開いたら自動で EIC の抽出を始める(30 秒前後かかるので、
        # 表を読んでいる間に裏で走らせる)
        QTimer.singleShot(0, self._start_eic)









    # ── _EICPanelMixin のフック─────────────────────────
    def _eic_markers(self) -> list:
        """候補の RT に 1, 2, … の番号付き縦線を引く。"""
        return [(float(c.get('obs_rt', 0.0)), str(i + 1))
                for i, c in enumerate(self._candidates or [])]

    def _eic_selected_marker(self) -> int:
        return self._selected_row_index()

    def _eic_on_marker_clicked(self, index: int):
        """図のクリックで表のラジオを動かす(表と図を連動させる)。"""
        rb = (self._radio_buttons or [])[index]             if 0 <= index < len(self._radio_buttons or []) else None
        if rb is not None and not rb.isChecked():
            rb.setChecked(True)      # toggled → _draw_eic が走る
        else:
            self._draw_eic()
        self._table.selectRow(index)

    def _selected_row_index(self) -> int:
        for i, rb in enumerate(getattr(self, '_radio_buttons', []) or []):
            if rb.isChecked():
                return i
        return -1



    def closeEvent(self, ev):
        self._stop_eic("")
        super().closeEvent(ev)

    def reject(self):
        self._stop_eic("")
        super().reject()

    def _build_table_columns(self):
        """Sample profile 表示状態に応じて列構成を再構築する。"""
        if self._show_sample_profile:
            cols = ["RT (min)", "obs m/z", "theor m/z", "Δm/z", "Δppm",
                    "mean int", "Sample profile", "Select"]
        else:
            cols = ["RT (min)", "obs m/z", "theor m/z", "Δm/z", "Δppm",
                    "mean int", "Select"]
        self._table.setColumnCount(len(cols))
        self._table.setHorizontalHeaderLabels(cols)
        # m/z 関連列を等幅、Sample profile を広めに
        for i in range(len(cols)):
            self._table.setColumnWidth(i, 95)
        if self._show_sample_profile:
            # Sample profile 列だけ広く
            sp_idx = cols.index("Sample profile")
            self._table.setColumnWidth(sp_idx, 160)

    def _on_toggle_profile(self, checked: bool):
        self._show_sample_profile = bool(checked)
        # 現在の選択を保存
        prev_pidx = self._selected_peak_idx
        for i, rb in enumerate(self._radio_buttons):
            if rb.isChecked():
                prev_pidx = self._candidates[i]['peak_idx']
                break
        self._build_table_columns()
        self._populate_rows(prev_pidx)

    def _sparkline(self, values: list[float]) -> str:
        """サンプル強度ベクトルから Unicode block スパークラインを生成。
        各行内で min/max にスケール(候補ごとの相対形状を可視化)。"""
        if not values:
            return ""
        vmin, vmax = min(values), max(values)
        n = len(self._SPARK_CHARS)
        if vmax - vmin < 1e-12:
            return self._SPARK_CHARS[n // 2] * len(values)
        return ''.join(
            self._SPARK_CHARS[
                min(n - 1, int((v - vmin) / (vmax - vmin) * (n - 1)))]
            for v in values)

    def _populate_rows(self, initial_peak_idx: int | None):
        """候補リストでテーブル行を構築。initial_peak_idx で初期選択を指定。

        ラジオを QButtonGroup にまとめる。各ラジオはセルごとの
        wrap QWidget に入っているため親が別々で、Qt の autoExclusive が
        効かず**排他になっていなかった**。行 2 を選んでも行 1 のチェックが
        外れないので、_on_ok() が拾う「最初に checked な行」は常に先頭の
        まま = ユーザーの選択が無視されていた(fix16 から存在)。
        QButtonGroup は親に関係なく排他を保証する。
        """
        self._table.setRowCount(0)
        self._radio_buttons = []
        # 行を作り直すたびに group ごと作り直す(古い group に残った
        # ボタンが排他判定に混ざらないようにする)
        old_group = getattr(self, '_radio_group', None)
        if old_group is not None:
            try:
                old_group.setParent(None)
            except Exception:
                pass
        self._radio_group = QButtonGroup(self)
        self._radio_group.setExclusive(True)
        if initial_peak_idx is None and self._candidates:
            initial_peak_idx = self._candidates[0]['peak_idx']

        for r, c in enumerate(self._candidates):
            self._table.insertRow(r)
            cells = [
                f"{c['obs_rt']:.3f}",
                f"{c['obs_mz']:.4f}",
                f"{self._theor_mz:.4f}",
                f"{c['obs_mz'] - self._theor_mz:+.4f}",
                f"{c['delta_ppm']:+.2f}",
                f"{c.get('mean_intensity', 0):.0f}",
            ]
            if self._show_sample_profile:
                spark = self._sparkline(c.get('sample_intensities', []))
                cells.append(spark)
            for col_idx, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                item.setTextAlignment(Qt.AlignCenter)
                self._table.setItem(r, col_idx, item)
            # ラジオボタン列(最終列)
            rb = QRadioButton()
            self._radio_group.addButton(rb, r)      # 排他にする
            rb.setChecked(c['peak_idx'] == initial_peak_idx)
            # 表で選び直したら EIC 上の強調も動かす
            rb.toggled.connect(self._on_radio_toggled)
            self._radio_buttons.append(rb)
            wrap = QWidget()
            wlay = QHBoxLayout(wrap)
            wlay.setContentsMargins(0, 0, 0, 0)
            wlay.addWidget(rb)
            wlay.setAlignment(Qt.AlignCenter)
            self._table.setCellWidget(r, len(cells), wrap)
        self._table.resizeColumnsToContents()
        # 行数に合わせて高さを詰める(既定だと候補 2 件でも
        # 表が縦に伸びて、下の EIC が潰れる)
        try:
            h = self._table.horizontalHeader().height() + 6
            for r in range(self._table.rowCount()):
                h += self._table.rowHeight(r)
            self._table.setMaximumHeight(max(90, min(h, 320)))
        except Exception:
            pass

    def _on_radio_toggled(self, checked: bool):
        """ラジオが変わったら EIC の強調を描き直す。

        toggled は off 側と on 側で 2 回飛ぶので、on のときだけ描く。
        """
        if not checked:
            return
        try:
            self._draw_eic()
        except Exception as e:
            log.warning(f"EIC redraw failed: {e}")

    def _on_ok(self):
        self._stop_eic("")
        # 選択された候補の peak_idx を記録して accept
        # 排他になったので checked は高々 1 つ。念のため
        # 「複数 checked なら確定しない」ようにして、静かに誤った候補を
        # 返す事故が二度と起きないようにする。
        checked = [i for i, rb in enumerate(self._radio_buttons)
                   if rb.isChecked()]
        if len(checked) == 1:
            i = checked[0]
            self._selected_peak_idx = self._candidates[i]['peak_idx']
            self._selected_candidate = self._candidates[i]
        elif len(checked) > 1:
            QMessageBox.warning(
                self, "Selection ambiguous",
                f"{len(checked)} rows are selected at once.\n"
                "Select exactly one row.")
            return
        self.accept()

    def selected_candidate(self) -> dict | None:
        """OK 押下後、選択された候補 dict を返す(Cancel なら None)。"""
        if self._selected_peak_idx is None:
            return None
        for c in self._candidates:
            if c['peak_idx'] == self._selected_peak_idx:
                return c
        return None


class SpotEICDialog(_EICPanelMixin, QDialog):
    """散布図のスポット 1 点について EIC を出すダイアログ。

    要望:
      「IS 以外でもいくつかクロマトグラムを見たい場合がある。
        プロットからポップアップしたダイアログをクリックすると
        クロマトグラムを描画(最も強度の高いサンプル 1 つだけ)」

    Candidate Picker と同じ _EICPanelMixin を使うので、抽出は
    worker(Alignment_raw_worker.py --eic)の別プロセス、非ブロッキング、
    サンプルは「そのスポットが最も強く出ている 1 本」を自動選択。

    同位体は出さない。M+1 は 1.003 Da 離れるので
    ±0.005 Da の窓には元々入らない。
    """

    def __init__(self, mz: float, rt: float, label: str,
                 raw_paths, sample_col_names=None,
                 sample_intensities=None, ion_mode: str = '',
                 parent=None):
        super().__init__(parent)
        self.setWindowTitle(
            f"EIC  {label}" if label else
            f"EIC  RT={rt:.3f}  m/z={mz:.4f}")
        self.resize(880, 460)

        self._theor_mz = float(mz)
        self._spot_rt = float(rt)
        self._spot_label = str(label or '')
        self._raw_paths = [str(p) for p in (raw_paths or [])]
        self._sample_col_names = list(sample_col_names or [])
        # _pick_eic_sample_index が読む形(候補 1 件ぶんの dict)に合わせる。
        # これで「最も強度の高いサンプル」の選び方を共通化できる。
        self._candidates = [{
            'obs_rt': float(rt), 'obs_mz': float(mz),
            'sample_intensities': list(sample_intensities or []),
        }]
        self._eic_df = None
        self._eic_proc = None
        self._eic_timer = None
        self._eic_out_csv = None
        self._eic_t0 = 0.0
        self._eic_sample_name = ""

        lay = QVBoxLayout(self)
        head = QLabel(
            f"<b>{self._spot_label or '(no label)'}</b>"
            f"   /   RT <b>{rt:.3f}</b> min"
            f"   /   m/z <b>{mz:.4f}</b>"
            + (f"   /   {ion_mode}" if ion_mode else ""))
        head.setTextFormat(Qt.RichText)
        head.setStyleSheet("color:#444; font-size:11px; padding:4px;")
        head.setWordWrap(True)
        lay.addWidget(head)

        self._build_eic_panel(lay)

        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        btns.accepted.connect(self.accept)
        lay.addWidget(btns)

        # 開いたら自動で抽出を始める(15〜30 秒かかるので裏で走らせる)
        QTimer.singleShot(0, self._start_eic)

    # ── _EICPanelMixin のフック ──────────────────────────────────
    def _eic_markers(self) -> list:
        """このスポットの RT に 1 本だけ縦線を引く。"""
        return [(self._spot_rt, '')]

    def _eic_selected_marker(self) -> int:
        return 0        # 常に強調(1 本しかない)

    # 選ぶものが無いので、クリックの案内は出さない
    EIC_CLICK_HINT = ""

    def _eic_on_marker_clicked(self, index: int):
        return          # 選ぶものが無いので何もしない

    def closeEvent(self, ev):
        self._stop_eic("")
        super().closeEvent(ev)

    def reject(self):
        self._stop_eic("")
        super().reject()


# ════════════════════════════════════════════════════════════════════
#  MainWindow
# ════════════════════════════════════════════════════════════════════

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"LipidZoner  v{__version__}")
        self.resize(1000, 700)
        self._file_entries: list[FileEntry] = []
        self._file_counter  = 0
        self._output_dir: Path | None = None
        # PreviewRTDialog の設定を Ion mode 別に保持
        # {"pos": {...}, "neg": {...}}
        self._preview_settings_by_mode: dict[str, dict] = {}
        # file_tag → reserved_coords (IS優先マッチで帰属されたobs座標の集合)
        # エクスポート時にRT範囲切り出し後、該当行を除外するために使う
        self._reserved_coords_by_tag: dict[str, set] = {}
        # file_tag → {loser_class: set of (rt, mz)}
        # Coherence 帰属で敗者となったクラス・座標。エクスポート時に
        # 対象クラス = 敗者クラスとなっている座標を除外する。
        self._coherence_loser_by_tag: dict[str, dict[str, set]] = {}
        # file_tag → match_df。PreviewRTDialog から伝播され、
        # _export_one が final_status='kept' フィルタに使う。
        self._preview_match_df_by_tag: dict[str, pd.DataFrame] = {}
        # file_tag → PreviewRTDialog インスタンス
        # ダイアログを × で閉じても同じインスタンスを保持しておき、
        # 再度開く際は show() し直すことで pipeline 結果(match_df,
        # _adduct_attribution, _coherence_assignments 等)を保全する。
        self._preview_dialogs: dict[str, QDialog] = {}

        cw = QWidget()
        self.setCentralWidget(cw)
        root = QVBoxLayout(cw)

        # ── 入力ファイル ──────────────────────────────────────────────
        # Preview RT を補助ボタンと階層分離して強調
        file_box = QGroupBox("Input Files")
        fb = QHBoxLayout(file_box)
        self._cmb_files = QComboBox()
        self._cmb_files.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        fb.addWidget(self._cmb_files)
        # 補助操作グループ(Add / Remove / Inspect)
        for txt, slot in [
            ("Add …",      self._add_file),
            ("Remove",     self._remove_file),
            ("Inspect",    self._inspect),
        ]:
            b = QPushButton(txt)
            b.clicked.connect(slot)
            fb.addWidget(b)
        # 縦区切りで階層を視覚化
        sep = QFrame()
        sep.setFrameShape(QFrame.VLine)
        sep.setFrameShadow(QFrame.Sunken)
        fb.addWidget(sep)
        # メイン操作: Annotation Pipeline(青強調)
        # "Preview RT" → "Annotation Pipeline" にリネーム
        btn_preview = QPushButton("Annotation Pipeline…")
        btn_preview.setToolTip(
            "Open the Annotation Pipeline workspace\n"
            "(Match → Conflict → IS → Adduct → RT Outlier → Coherence)\n"
            "and auto-generate tasks for the current file.")
        btn_preview.setMinimumHeight(34)
        btn_preview.setStyleSheet(
            "QPushButton {"
            " font-size: 13px; font-weight: bold;"
            " padding: 4px 14px;"
            " background-color: #2A6FC9; color: white;"
            " border: 1px solid #1E5BAA; border-radius: 4px;"
            "}"
            "QPushButton:hover { background-color: #3782D9; }"
            "QPushButton:pressed { background-color: #1E5BAA; }")
        btn_preview.clicked.connect(self._preview_rt)
        fb.addWidget(btn_preview)
        root.addWidget(file_box)

        # ── 出力フォルダ ──────────────────────────────────────────────
        out_box = QGroupBox("Output Directory")
        ob = QHBoxLayout(out_box)
        self._lbl_outdir = QLabel("(default: same as input file)")
        ob.addWidget(self._lbl_outdir, 1)
        b_out = QPushButton("Choose …")
        b_out.clicked.connect(self._choose_outdir)
        ob.addWidget(b_out)
        root.addWidget(out_box)

        # ── タスクテーブル ────────────────────────────────────────────
        task_box = QGroupBox("Tasks")
        tb = QVBoxLayout(task_box)
        self._task_table = QTableWidget(0, 5)
        self._task_table.setHorizontalHeaderLabels(
            ["Class", "RT start", "RT end", "Ion Mode", "File"])
        self._task_table.horizontalHeader().setStretchLastSection(True)
        self._task_table.setContextMenuPolicy(Qt.CustomContextMenu)
        self._task_table.customContextMenuRequested.connect(
            self._task_ctx_menu)
        tb.addWidget(self._task_table)

        th = QHBoxLayout()
        for txt, slot in [
            ("Add Row",              self._add_task_row),
            ("Delete Row",           self._delete_task_row),
            ("Paste from Clipboard", self._paste_tasks),
            ("Copy All Tasks",       self._copy_all_tasks),
        ]:
            b = QPushButton(txt)
            b.clicked.connect(slot)
            th.addWidget(b)
        tb.addLayout(th)
        root.addWidget(task_box)

        # ── 実行ボタン ────────────────────────────────────────────────
        # self._btn_run に保存して実行中は無効化する
        self._btn_run = QPushButton("▶  Run Export")
        self._btn_run.setStyleSheet("font-size:14px;padding:8px;")
        self._btn_run.clicked.connect(self._run)
        root.addWidget(self._btn_run)

        self._lbl_log = QLabel("")
        root.addWidget(self._lbl_log)

        # ── Tasks UI を hide、AP Tasks タブに集約 ────
        # Output Directory / Tasks table / Run Export / Log は AP の
        # Tasks タブに移行済み(Stage 2)。MainWindow からは視覚的に
        # 隠すが、Python オブジェクトは残しておく(AP の delegate が
        # 参照するデータストアとして必要)。
        out_box.setVisible(False)
        task_box.setVisible(False)
        self._btn_run.setVisible(False)
        self._lbl_log.setVisible(False)
        _main_info_lbl = QLabel(
            "Tasks have moved to the Annotation Pipeline window.\n"
            "Click 'Annotation Pipeline' next to a file row to access them.")
        _main_info_lbl.setStyleSheet(
            "color:#666; padding:24px; font-style:italic;"
            "background-color:#f7f7f7; border:1px solid #ddd;"
            "border-radius:4px;")
        _main_info_lbl.setAlignment(Qt.AlignCenter)
        _main_info_lbl.setWordWrap(True)
        root.addWidget(_main_info_lbl)
        root.addStretch(1)

        # ── 修正: ウィンドウサイズ + 中央配置 ──────────────────
        try:
            screen = QApplication.primaryScreen()
            if screen is not None:
                geom = screen.availableGeometry()
                target_w = min(1280, max(800, int(geom.width() * 0.7)))
                # Tasks UI hide により高さ大幅縮小
                target_h = min(520, max(360, int(geom.height() * 0.45)))
                self.resize(target_w, target_h)
                self.move(
                    geom.x() + (geom.width() - target_w) // 2,
                    geom.y() + (geom.height() - target_h) // 2,
                )
        except Exception as _e:
            log.warning(f"[MainWindow position] failed: {_e}")

        # ── メニューバー ──────────────────────────────────────────────
        # メニュー文言を Session File に統一
        session_menu = self.menuBar().addMenu("Session")
        for label, shortcut, slot in [
            ("Save LipidZoner Session File…", "Ctrl+S", self._save_session),
            ("Load LipidZoner Session File…", "Ctrl+O", self._load_session),
        ]:
            act = QAction(label, self)
            act.setShortcut(QKeySequence(shortcut))
            act.triggered.connect(slot)
            session_menu.addAction(act)

        # View メニュー撤廃
        # Windows のウィンドウ操作(最大化/復元)で代替可能なため不要。

    # ── ファイル操作 ──────────────────────────────────────────────────
    def _add_file(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Open MS-DIAL export(s)",
            "", "Text files (*.txt *.tsv);;All (*)")
        for p in paths:
            # Ion mode をダイアログで尋ねる（pos/neg 2択）
            ion_mode = self._ask_ion_mode(Path(p).name)
            if ion_mode is None:
                continue  # キャンセル時はスキップ

            self._file_counter += 1
            tag = f"#{self._file_counter}"
            try:
                fe = FileEntry(p, tag, ion_mode=ion_mode)
            except Exception as e:
                QMessageBox.critical(self, "Load Error", f"{p}:\n{e}")
                continue
            self._file_entries.append(fe)
            self._cmb_files.addItem(fe.label, tag)
        if self._file_entries:
            self._cmb_files.setCurrentIndex(self._cmb_files.count() - 1)

    @property
    def _align_params(self) -> dict:
        """per-sample alignment / 同位体グルーピングのパラメータ。

        fix19 で「Peak Detection & Alignment」タブから編集できるようにする。
        それまではここが唯一の出所で、読込経路(手動 / session)の両方が
        同じ値を使う。"""
        if getattr(self, '_align_params_store', None) is None:
            self._align_params_store = dict(DEFAULT_ALIGN_PARAMS)
        return self._align_params_store

    @_align_params.setter
    def _align_params(self, value: dict):
        d = dict(DEFAULT_ALIGN_PARAMS)
        if value:
            d.update({k: v for k, v in value.items() if k in d})
        self._align_params_store = d

    def _add_per_sample_folder(self):
        """Per-sample peak list フォルダを選択し PerSampleFileEntry を作成。

        v1.1.0-alpha.6 で追加。MS-DIAL の Peak list result(per-sample TXT)が
        並んだフォルダを指定すると、内部で simple alignment + isotope grouping
        を実行して PerSampleFileEntry を生成し、_file_entries に追加する。

        既存の FileEntry と同じリストに混在可能(FileEntry-compatible API)。
        """
        folder = QFileDialog.getExistingDirectory(
            self, "Select Per-sample Peak List Folder", "")
        if not folder:
            return
        ion_mode = self._ask_ion_mode(Path(folder).name)
        if ion_mode is None:
            return  # キャンセル時はスキップ

        # 以前はここで既定値のまま無言でアライメントが走っていた。
        # 実行前に必ず使用パラメータを提示し、詳細な調整は
        # 「Peak Detection & Alignment」タブへ誘導する。
        _p = self._align_params
        _msg = (
            f"Alignment and isotope grouping will run with these parameters:\n\n"
            f"  alignment  m/z tol : {_p['align_mz_tol']} Da"
            + (f" + {_p['align_mz_ppm']} ppm" if _p['align_mz_ppm'] else "")
            + f"\n"
            f"  alignment  RT tol  : {_p['align_rt_tol']} min\n"
            f"  refine iterations  : {_p['align_n_refine']}\n"
            f"  intensity column   : {_p['height_col']}\n"
            f"  RT drift correction: {_p['rt_drift_correct']}\n"
            f"  isotope grouping   : {_p['iso_enabled']}"
            + (f"  (mz {_p['iso_mz_tol']} / rt {_p['iso_rt_tol']} / "
               f"max M+{_p['iso_max_iso']} / ratio tol {_p['iso_ratio_tol']})"
               if _p['iso_enabled'] else "")
            + "\n\n"
              "To change them, cancel and run it from the\n"
              f"'{TAB_LABELS['align']}' tab.")
        if QMessageBox.question(
                self, "Confirm alignment parameters", _msg,
                QMessageBox.Ok | QMessageBox.Cancel,
                QMessageBox.Ok) != QMessageBox.Ok:
            return

        self._file_counter += 1
        tag = f"#{self._file_counter}"
        try:
            fe = create_per_sample_file_entry(
                folder=folder, ion_mode=ion_mode, tag=tag,
                params=self._align_params,
            )
        except Exception as e:
            QMessageBox.critical(
                self, "Per-sample Load Error",
                f"Failed to load per-sample folder:\n{folder}\n\n{e}")
            self._file_counter -= 1
            return
        self._file_entries.append(fe)
        self._cmb_files.addItem(fe.label, tag)
        if self._file_entries:
            self._cmb_files.setCurrentIndex(self._cmb_files.count() - 1)

    def _ask_ion_mode(self, filename: str) -> str | None:
        """Ion mode 選択ダイアログ。pos/neg を返す、キャンセル時は None。"""
        from PySide6.QtWidgets import QInputDialog
        mode, ok = QInputDialog.getItem(
            self, "Select Ion Mode",
            f"Select ion mode for:\n{filename}",
            ["pos", "neg"], 0, False)
        return mode if ok else None

    def _current_fe(self) -> FileEntry | None:
        tag = self._cmb_files.currentData()
        return next((fe for fe in self._file_entries if fe.tag == tag), None)

    def _remove_file(self):
        fe = self._current_fe()
        if fe:
            # 該当ファイルの Preview ダイアログがあれば破棄
            dlg = self._preview_dialogs.pop(fe.tag, None)
            if dlg is not None:
                try:
                    dlg.close()
                    dlg.deleteLater()
                except Exception:
                    pass
            self._file_entries.remove(fe)
            self._cmb_files.removeItem(self._cmb_files.currentIndex())

    def _inspect(self):
        fe = self._current_fe()
        if fe:
            InspectDialog(fe, self).exec()

    def _preview_rt(self):
        fe = self._current_fe()
        if fe is None:
            return
        # 既存ダイアログがあれば再利用(state を保全)
        existing = self._preview_dialogs.get(fe.tag)
        if existing is not None:
            try:
                # まだ生きていれば show() して前面に
                existing.show()
                existing.raise_()
                existing.activateWindow()
                # AP 表示時は MainWindow を hide
                self.hide()
                return
            except RuntimeError:
                # Qt 側で破棄済みなら辞書からも消す
                self._preview_dialogs.pop(fe.tag, None)
        # 新規作成
        mode_settings = self._preview_settings_by_mode.get(fe.ion_mode)
        dlg = PreviewRTDialog(
            fe, self._file_entries,
            initial_settings=mode_settings,
            all_initial_settings=self._preview_settings_by_mode,
            parent=self)
        dlg.task_ready.connect(self._insert_task)
        dlg.settings_saved.connect(
            lambda settings, m=fe.ion_mode:
                self._on_preview_settings_saved(m, settings))
        dlg.reserved_updated.connect(self._on_reserved_updated)
        dlg.coherence_updated.connect(self._on_coherence_updated)
        dlg.match_df_updated.connect(self._on_match_df_updated)
        # 破棄せず保持。close 時も destroy しないよう
        # WA_DeleteOnClose は明示的に false にしておく(デフォルトでも false)。
        dlg.setAttribute(Qt.WA_DeleteOnClose, False)
        self._preview_dialogs[fe.tag] = dlg
        dlg.show()
        dlg.raise_()
        dlg.activateWindow()
        # AP 表示時は MainWindow を hide
        self.hide()
        # Load Session 直後 + 保存時に match_run=True
        # だった場合は、自動的に Run All を起動して状態を再現
        # 同時に Quant Ion 選択も新規ダイアログへ伝播
        try:
            fe_exempt = getattr(self, '_loaded_session_filter_exempt', None)
            qi = getattr(self, '_loaded_session_quant_ion', None)
            if hasattr(dlg, 'auto_replay_after_load'):
                dlg.auto_replay_after_load(
                    filter_exemptions=fe_exempt,
                    quant_ion_choices=qi)
        except Exception as e:
            log.warning(f"[_open_preview] auto-replay trigger failed: {e}")

    def _on_preview_settings_saved(self, ion_mode: str, settings: dict):
        """PreviewRTDialog が閉じられたときに呼ばれ、設定を Ion mode 別に保持する"""
        self._preview_settings_by_mode[ion_mode] = settings

    def _on_reserved_updated(self, file_tag: str, reserved_coords: set):
        """PreviewRTDialog から reserved_coords の更新通知を受け取る"""
        if reserved_coords:
            self._reserved_coords_by_tag[file_tag] = set(reserved_coords)
        elif file_tag in self._reserved_coords_by_tag:
            # 空集合に更新された場合は削除（メモリ節約）
            del self._reserved_coords_by_tag[file_tag]

    def _on_coherence_updated(
        self, file_tag: str, loser_coords_by_class: dict,
    ):
        """PreviewRTDialog から Coherence 帰属(敗者クラス座標)を受け取る"""
        if loser_coords_by_class:
            # deepcopy 不要(受信 dict を保持するだけ)
            self._coherence_loser_by_tag[file_tag] = {
                k: set(v) for k, v in loser_coords_by_class.items()}
        elif file_tag in self._coherence_loser_by_tag:
            del self._coherence_loser_by_tag[file_tag]

    def _on_match_df_updated(self, file_tag: str, match_df):
        """PreviewRTDialog から match_df を受け取り、_export_one で
        final_status='kept' フィルタとして使う。"""
        if match_df is not None and isinstance(match_df, pd.DataFrame):
            self._preview_match_df_by_tag[file_tag] = match_df.copy()
        elif file_tag in self._preview_match_df_by_tag:
            del self._preview_match_df_by_tag[file_tag]

    def _choose_outdir(self):
        d = QFileDialog.getExistingDirectory(
            self, "Choose Output Directory")
        if d:
            self._output_dir = Path(d)
            self._lbl_outdir.setText(str(self._output_dir))

    # ── セッション保存 / 読み込み ─────────────────────────────────────

    def _collect_tasks(self) -> list[dict]:
        tasks = []
        for row in range(self._task_table.rowCount()):
            t = self._read_task(row)
            if t:
                tasks.append(t)
        return tasks

    # どの Save 入り口からも統合 JSON を書き込めるよう
    # MainWindow に共通の save_unified_state を持たせる。
    def _gather_state_from_dialogs(self):
        """開いている preview dialog から params / quant_ion /
        filter_exemptions を集める(初出のものを使う)。"""
        params = None
        quant_ion = None
        filter_exemptions = None
        for dlg in (self._preview_dialogs or {}).values():
            if dlg is None:
                continue
            try:
                if params is None and hasattr(dlg, '_gather_current_params'):
                    params = dlg._gather_current_params()
            except Exception:
                pass
            try:
                if quant_ion is None:
                    qc = getattr(dlg, '_quant_ion_choices', None)
                    if qc:
                        quant_ion = dict(qc)
            except Exception:
                pass
            try:
                if filter_exemptions is None:
                    fe_dict = getattr(dlg, '_filter_class_exemptions', None)
                    if fe_dict:
                        fe_serial = {
                            k: sorted(v) for k, v in fe_dict.items() if v}
                        if fe_serial:
                            filter_exemptions = fe_serial
            except Exception:
                pass
        return params, quant_ion, filter_exemptions

    def save_unified_state(self, path: str,
                           override_params=None,
                           override_quant_ion=None,
                           override_filter_exemptions=None) -> None:
        """v2 形式で session を保存する共通メソッド。
        過去の merge 挙動は廃止(format_version 2 で完全上書き)。"""
        data = _build_lz_v2_dict(
            self,
            override_params=override_params,
            override_quant_ion=override_quant_ion,
            override_filter_exemptions=override_filter_exemptions)
        # アライメント結果をサイドカーに書き出す。
        # これが無いと raw mzML から作ったエントリを復元できない
        # (フォルダに TXT が無いので再アライメントできず、黙って消える)。
        try:
            n_side = write_session_sidecar(self, path, data)
            if n_side:
                log.info(f"[sidecar] written: {n_side} entries -> "
                      f"{session_sidecar_dir(path)}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            log.warning(f"sidecar write failed: {e}")
        Path(path).write_text(
            json.dumps(data, ensure_ascii=False, indent=2),
            encoding="utf-8")

    def _save_session(self):
        # ファイル名・タイトルを統一して merge 挙動を明示
        path, _ = QFileDialog.getSaveFileName(
            self, "Save LipidZoner Session File",
            "lipidzoner_session.json",
            "JSON files (*.json);;All (*)")
        if not path:
            return
        if not path.endswith(".json"):
            path += ".json"

        # 共通の save_unified_state を呼ぶだけ。
        try:
            self.save_unified_state(path)
            self._lbl_log.setText(f"Session saved: {path}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", str(e))

    def _load_session(self):
        """v2 形式の Session ファイルを読込。
        旧 v1 形式は受け付けない(format_version チェック)。"""
        path, _ = QFileDialog.getOpenFileName(
            self, "Load LipidZoner Session File",
            "", "JSON files (*.json);;All (*)")
        if not path:
            return
        try:
            data = _load_lz_v2_file(path)
        except Exception as e:
            QMessageBox.critical(self, "Load Error", str(e))
            return

        # smart fallback: library_path が見つからなければ JSON と同じ
        # ディレクトリの basename を試す(配布パッケージ用)
        json_dir = Path(path).parent
        lib_path_str = data.get('library_path')
        if lib_path_str:
            lp = Path(lib_path_str)
            if not lp.exists():
                alt = json_dir / lp.name
                if alt.exists():
                    data['library_path'] = str(alt)

        # ── 直接トップレベルから読み出し ──────────────────────────
        params = data.get('parameters') or {}
        quant_ion = data.get('quant_ion_choices') or {}
        filter_exempt = data.get('filter_exemptions') or {}

        # 次にダイアログが開かれた時の auto-replay でも使う
        self._loaded_session_params = params or None
        self._loaded_session_quant_ion = quant_ion or None
        self._loaded_session_filter_exempt = filter_exempt or None

        # ── alignment セクションを復元 ─────────────────────
        # ここで host._align_params を先に更新しておく。直後のファイル
        # 再ロードで、per-file の条件が無いエントリはこの値を使う。
        align_sec = data.get('alignment') or {}
        if align_sec.get('params'):
            try:
                self._align_params = align_sec['params']
            except Exception as e:
                log.warning(f"[_load_session] align params restore failed: {e}")
        self._loaded_session_alignment = align_sec or None

        # ── 開いている AP に params / quant_ion / exemptions を適用 ──
        per_mode_settings = _v2_to_per_mode_settings(data)
        try:
            for dlg in (self._preview_dialogs or {}).values():
                if dlg is None:
                    continue
                if params and hasattr(dlg, '_apply_params'):
                    try:
                        dlg._apply_params(params)
                    except Exception:
                        pass
                # ① タブのパラメータも復元する
                if align_sec and hasattr(dlg, '_al_apply_session_section'):
                    try:
                        dlg._al_apply_session_section(align_sec)
                    except Exception as e:
                        log.warning(f"[_load_session] ① tab restore failed: {e}")
                if quant_ion:
                    try:
                        dlg._quant_ion_choices = dict(quant_ion)
                    except Exception:
                        pass
                if filter_exempt:
                    try:
                        dlg._filter_class_exemptions = {
                            k: set(v) for k, v in filter_exempt.items()}
                    except Exception:
                        pass
        except Exception as e:
            log.warning(f"[_load_session] section apply failed: {e}")

        # ── ファイル再ロード ──────────────────────────────────────
        self._file_entries.clear()
        self._cmb_files.clear()
        self._file_counter = 0
        skipped = []
        json_dir = Path(path).parent  # smart fallback 用: JSON と同じディレクトリ
        # 極性ごとに 1 エントリだけにする。旧 session に重複が
        # あっても、ここで最新の 1 件へ寄せる。
        _files, _dup_dropped = dedupe_entries_by_polarity(
            data.get("files", []))
        for entry in _files:
            tag = entry.get("tag", f"#{self._file_counter + 1}")
            try:
                num = int(tag.lstrip("#"))
                self._file_counter = max(self._file_counter, num)
            except ValueError:
                pass

            # ── まずサイドカーを試す ───────────────────────
            # 保存時のアライメント結果そのものを持っているので、
            # raw mzML 由来でも復元でき、再アライメントも要らない。
            # 元データフォルダが移動・削除されていても復元できる。
            fe = None
            if entry.get("per_sample") and entry.get("aligned_data"):
                try:
                    fe = read_session_sidecar(path, data, entry, tag)
                except Exception as e:
                    # 失敗しても黙って従来経路へ落ちるのは危険
                    # (raw mzML 由来だとフォルダに TXT が無く消える)ので、
                    # 理由を skipped に残してユーザーに見せる。
                    log.warning(f"[sidecar] load failed for {tag}: {e}")
                    skipped.append(
                        f"{tag}: could not read the sidecar ({e})")
                    fe = None
            if fe is not None:
                self._file_entries.append(fe)
                self._cmb_files.addItem(fe.label, tag)
                continue

            # ── 従来経路: フォルダ / ファイルから読み直す ─────────
            p = Path(entry["path"])
            if not p.exists():
                # smart fallback: JSON と同じディレクトリで basename を試す
                # (配布パッケージ用。JSON + .txt を同じフォルダに置けば動く)
                alt = json_dir / p.name
                if alt.exists():
                    p = alt
                else:
                    skipped.append(str(entry["path"]))
                    continue
            try:
                if entry.get("per_sample"):
                    # per-sample フォルダから PerSampleFileEntry を復元
                    # そのエントリを保存したときの条件を優先して使う。
                    #   無い場合(fix20 以前に保存した session)は現在の設定に
                    #   フォールバックする。
                    _fp = entry.get("align_params") or self._align_params
                    fe = create_per_sample_file_entry(
                        folder=str(p), ion_mode=entry.get("ion_mode", "pos"),
                        tag=tag, params=_fp)
                else:
                    fe = FileEntry(p, tag,
                                   ion_mode=entry.get("ion_mode", "pos"))
            except Exception as e:
                skipped.append(f"{p} ({e})")
                continue
            # 生データのパスも戻す(Candidate Picker の EIC 用)
            if entry.get("raw_paths"):
                try:
                    fe.raw_paths = list(entry["raw_paths"])
                except Exception:
                    pass
            self._file_entries.append(fe)
            self._cmb_files.addItem(fe.label, tag)
        if self._file_entries:
            self._cmb_files.setCurrentIndex(self._cmb_files.count() - 1)
        # 重複を落としたことは必ず知らせる。黙って捨てると
        # 「どちらのデータで解析しているのか」が分からなくなる。
        if _dup_dropped:
            log.info("[session] several entries shared the same polarity:\n  "
                  + "\n  ".join(_dup_dropped))
            QMessageBox.information(
                self, "Duplicate data reduced to one entry",
                "Annotation can hold only one per-sample dataset per\n"
                "polarity, so entries sharing a polarity were reduced to\n"
                "the most recent one:\n\n  " + "\n  ".join(_dup_dropped)
                + "\n\nTo use the one that was dropped, re-enter its conditions\n"
                  "in the ① tab and Send it again.")

        # ── 出力フォルダ ──────────────────────────────────────────
        od = data.get("output_dir")
        if od and Path(od).exists():
            self._output_dir = Path(od)
            self._lbl_outdir.setText(str(self._output_dir))
        else:
            self._output_dir = None
            self._lbl_outdir.setText("(default: same as input file)")

        # ── Preview 設定(per-mode)──────────────────────────────
        self._preview_settings_by_mode = per_mode_settings

        # ── reserved_coords 復元 ──────────────────────────────────
        self._reserved_coords_by_tag = {}
        for tag, coord_list in (data.get("reserved_coords_by_tag") or {}).items():
            self._reserved_coords_by_tag[tag] = {
                (float(c[0]), float(c[1])) for c in coord_list
            }

        # ── タスク復元 ────────────────────────────────────────────
        self._task_table.setRowCount(0)
        for t in data.get("tasks", []) or []:
            self._insert_task(t)

        # ── 結果通知 ──────────────────────────────────────────────
        msg = f"Session loaded: {Path(path).name}"
        if skipped:
            detail = "\n".join(skipped)
            QMessageBox.warning(
                self, "Files Not Found",
                f"The following files could not be loaded:\n\n{detail}")
        self._lbl_log.setText(
            msg + (f"  ({len(skipped)} file(s) skipped)" if skipped else ""))

    def _load_session_settings_only(self) -> bool:
        """Settings-only ロード: parameters + filter_exemptions のみ別セッション
        ファイルから適用する。データファイル / output_dir / tasks /
        reserved_coords / preview_settings / quant_ion_choices には触れない。

        Returns
        -------
        bool
            True: ファイル読込・適用成功 / False: キャンセル or エラー
        """
        path, _ = QFileDialog.getOpenFileName(
            self, "Load LipidZoner Session File (Settings only)",
            "", "JSON files (*.json);;All (*)")
        if not path:
            return False
        try:
            data = _load_lz_v2_file(path)
        except Exception as e:
            QMessageBox.critical(self, "Load Error", str(e))
            return False

        params = data.get('parameters') or {}
        filter_exempt = data.get('filter_exemptions') or {}
        # alignment / 検出パラメータも「設定」として扱う
        align_sec = data.get('alignment') or {}

        if not params and not filter_exempt and not align_sec:
            QMessageBox.warning(
                self, "Nothing to apply",
                "The session file does not contain parameters, filter\n"
                "exemptions, or alignment settings to apply.")
            return False

        if align_sec.get('params'):
            try:
                self._align_params = align_sec['params']
            except Exception as e:
                log.warning(f"[settings_only] align params failed: {e}")
        for _dlg in (self._preview_dialogs or {}).values():
            if _dlg is None or not hasattr(_dlg, '_al_apply_session_section'):
                continue
            try:
                _dlg._al_apply_session_section(align_sec)
            except Exception as e:
                log.warning(f"[settings_only] ① tab restore failed: {e}")

        # 一時保存(AP の auto_replay や次回 dialog 起動時の参照用)
        self._loaded_session_params = params or None
        self._loaded_session_filter_exempt = filter_exempt or None
        # quant_ion は新データセットに対しては引き継がない(feature が異なる)
        self._loaded_session_quant_ion = None

        # 開いている AP 全てに適用
        try:
            for dlg in (self._preview_dialogs or {}).values():
                if dlg is None:
                    continue
                if params and hasattr(dlg, '_apply_params'):
                    try:
                        dlg._apply_params(params)
                    except Exception as _e:
                        log.warning(f"[_load_session_settings_only apply_params] failed: {_e}")
                if filter_exempt:
                    try:
                        dlg._filter_class_exemptions = {
                            k: set(v) for k, v in filter_exempt.items()}
                    except Exception as _e:
                        log.warning(f"[_load_session_settings_only fe_exempt] failed: {_e}")
        except Exception as e:
            log.warning(f"[_load_session_settings_only] section apply failed: {e}")

        self._lbl_log.setText(
            f"Session settings applied (settings only): {Path(path).name}")
        return True

    # ── タスクテーブル操作 ────────────────────────────────────────────
    def _add_task_row(self):
        r = self._task_table.rowCount()
        self._task_table.insertRow(r)
        self._task_table.setCellWidget(
            r, 3, self._ion_combo())
        self._task_table.setCellWidget(
            r, 4, self._file_combo())

    def _ion_combo(self, current="pos") -> QComboBox:
        cmb = QComboBox()
        cmb.addItems(["pos", "neg"])
        cmb.setCurrentText(current)
        return cmb

    def _file_combo(self, current_tag: str = "") -> QComboBox:
        cmb = QComboBox()
        for fe in self._file_entries:
            cmb.addItem(fe.label, fe.tag)
        idx = cmb.findData(current_tag)
        if idx >= 0:
            cmb.setCurrentIndex(idx)
        return cmb

    def _delete_task_row(self):
        rows = sorted(
            set(idx.row() for idx in
                self._task_table.selectedIndexes()),
            reverse=True)
        for r in rows:
            self._task_table.removeRow(r)

    def _insert_task(self, d: dict):
        r = self._task_table.rowCount()
        self._task_table.insertRow(r)
        self._task_table.setItem(r, 0, QTableWidgetItem(d["cls"]))
        self._task_table.setItem(r, 1, QTableWidgetItem(str(d["rt_start"])))
        self._task_table.setItem(r, 2, QTableWidgetItem(str(d["rt_end"])))
        self._task_table.setCellWidget(r, 3, self._ion_combo(d["ion"]))
        self._task_table.setCellWidget(r, 4, self._file_combo(d["file_tag"]))

    def _paste_tasks(self):
        txt = QApplication.clipboard().text()
        for line in txt.strip().splitlines():
            parts = re.split(r"[\t,]", line)
            if len(parts) < 4:
                continue
            self._insert_task(dict(
                cls=parts[0].strip(),
                rt_start=float(parts[1].strip()),
                rt_end=float(parts[2].strip()),
                ion=parts[3].strip(),
                file_tag=(parts[4].strip() if len(parts) > 4
                          else (self._file_entries[0].tag
                                if self._file_entries else "")),
            ))

    def _task_ctx_menu(self, pos):
        menu = QMenu(self)
        menu.addAction("Add Row",              self._add_task_row)
        menu.addAction("Delete Selected",      self._delete_task_row)
        menu.addAction("Paste from Clipboard", self._paste_tasks)
        menu.addSeparator()
        menu.addAction("Copy All Tasks",       self._copy_all_tasks)
        menu.exec(self._task_table.viewport().mapToGlobal(pos))

    def _copy_all_tasks(self):
        if self._task_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "No tasks to copy.")
            return
        lines = ["\t".join(
            self._task_table.horizontalHeaderItem(c).text()
            for c in range(self._task_table.columnCount()))]
        for row in range(self._task_table.rowCount()):
            cells = []
            for col in range(self._task_table.columnCount()):
                w = self._task_table.cellWidget(row, col)
                if isinstance(w, QComboBox):
                    cells.append(w.currentText())
                else:
                    item = self._task_table.item(row, col)
                    cells.append(item.text() if item else "")
            lines.append("\t".join(cells))
        QApplication.clipboard().setText("\n".join(lines))
        QMessageBox.information(
            self, "Success",
            f"Copied {self._task_table.rowCount()} task(s) to clipboard.")


    # ── エクスポート ──────────────────────────────────────────────────
    def _read_task(self, row: int) -> dict | None:
        items = [self._task_table.item(row, c) for c in range(3)]
        if not all(items):
            return None
        ion_w  = self._task_table.cellWidget(row, 3)
        file_w = self._task_table.cellWidget(row, 4)
        return dict(
            cls=items[0].text().strip(),
            rt_start=float(items[1].text()),
            rt_end=float(items[2].text()),

            ion=ion_w.currentText()  if ion_w  else "pos",
            file_tag=file_w.currentData() if file_w else None,
        )

    def _fe_by_tag(self, tag: str) -> FileEntry | None:
        return next((fe for fe in self._file_entries
                     if fe.tag == tag), None)

    def _run(self):
        n = self._task_table.rowCount()
        if n == 0:
            QMessageBox.warning(self, "No Tasks",
                                "Add at least one task row.")
            return
        # 進捗表示 + 連打防止
        original_label = "▶  Run Export"
        try:
            self._btn_run.setEnabled(False)
            self._btn_run.setText("Running…")
            QApplication.processEvents()
        except Exception:
            pass
        exported = 0
        errors: list = []
        for row in range(n):
            task = self._read_task(row)
            if not task:
                continue
            fe = self._fe_by_tag(task["file_tag"])
            if not fe:
                self._lbl_log.setText(f"Row {row+1}: file not found")
                errors.append(f"Row {row+1}: file not found")
                continue
            # 進捗ラベルを更新して processEvents で UI に反映
            try:
                self._btn_run.setText(f"Running… ({row+1}/{n})")
                self._lbl_log.setText(
                    f"Exporting {row+1}/{n}: "
                    f"{task.get('cls', '?')} ({task.get('ion', '?')})…")
                QApplication.processEvents()
            except Exception:
                pass
            try:
                self._export_one(fe, task)
                exported += 1
            except Exception as e:
                errors.append(f"Row {row+1}: {e}")
        # 後処理: ボタン復帰 + 完了ログ + ポップアップ通知
        try:
            self._btn_run.setEnabled(True)
            self._btn_run.setText(original_label)
        except Exception:
            pass
        summary = f"Done – {exported}/{n} tasks exported."
        if errors:
            summary += f"  ({len(errors)} error(s))"
        try:
            self._lbl_log.setText(summary)
        except Exception:
            pass
        # 完了ポップアップで連打防止 + 結果明示
        if errors:
            err_detail = "\n".join(errors[:10])
            if len(errors) > 10:
                err_detail += f"\n  … and {len(errors) - 10} more"
            QMessageBox.warning(
                self, "Export finished with errors",
                f"{summary}\n\n{err_detail}")
        else:
            QMessageBox.information(
                self, "Export complete",
                f"{summary}\nOutput directory:\n"
                f"{self._output_dir or '(same as input file)'}")

    def _export_one(self, fe: FileEntry, task: dict):
        """単一タスクのエクスポート(LipidQuant 入力形式の TXT を出力)。

          - PreviewRTDialog で生成された match_df があり status 列を持つ場合、
            (lipid_class == task["cls"]) かつ (final_status == 'kept') かつ
            RT が [rt_start, rt_end] に入る match のピークインデックスのみ出力。
          - match_df が無い・古い場合はレガシー挙動にフォールバック
            (RT 範囲 + reserved + Coherence loser の除外ロジック)。

        除外ロジック:
          1. RT 範囲外 → 除外
          2. _reserved_coords_by_tag(IS 帰属座標)→ 除外
          3. _coherence_loser_by_tag[fe.tag][task["cls"]]
             (Coherence 帰属で敗者になった座標)→ 除外
        """
        rt     = fe.rt_array()
        mz_raw = fe.mz_array()
        if rt is None or mz_raw is None:
            raise ValueError("RT or m/z column not found.")

        # ── match_df ベースのフィルタ ────────────────────
        used_v53_path = False
        excluded_v53 = {'rejected_window':0, 'rejected_outlier':0,
                        'rejected_residual':0, 'rejected_loser':0,
                        'rejected_manual':0}
        match_df = getattr(self, '_preview_match_df_by_tag', {}).get(fe.tag)
        if (match_df is not None
                and isinstance(match_df, pd.DataFrame)
                and not match_df.empty
                and 'final_status' in match_df.columns):
            used_v53_path = True
            cls = task['cls']
            sub_mdf = match_df[
                (match_df['lipid_class'] == cls) &
                (match_df.get('matched', False) == True) &
                (match_df['obs_rt'] >= task['rt_start']) &
                (match_df['obs_rt'] <= task['rt_end'])
            ]
            for col in ('is_filter_status', 'rt_outlier_status',
                        'coherence_status', 'manual_status'):
                if col in sub_mdf.columns:
                    for st in excluded_v53.keys():
                        excluded_v53[st] += int((sub_mdf[col] == st).sum())
            kept = sub_mdf[sub_mdf['final_status'] == STATUS_KEPT]
            if kept.empty:
                raise ValueError(
                    f"No kept matches for {task['cls']} in RT "
                    f"[{task['rt_start']}, {task['rt_end']}].")
            idx = np.array(sorted(set(kept['matched_idx'].astype(int).tolist())))
            excluded_count = 0  # path doesn't double-count
            coh_excluded = 0
        else:
            # フォールバック
            mask = (rt >= task["rt_start"]) & (rt <= task["rt_end"])
            idx  = np.where(mask)[0]
            if len(idx) == 0:
                raise ValueError(
                    f"No rows in RT [{task['rt_start']}, {task['rt_end']}].")

            reserved = self._reserved_coords_by_tag.get(fe.tag, set())
            excluded_count = 0
            if reserved:
                keep_mask = np.array([
                    (round(float(rt[i]), 6), round(float(mz_raw[i]), 6))
                    not in reserved
                    for i in idx
                ])
                excluded_count = int((~keep_mask).sum())
                idx = idx[keep_mask]
                if len(idx) == 0:
                    raise ValueError(
                        f"All rows in RT [{task['rt_start']}, {task['rt_end']}] "
                        f"were excluded as reserved for IS-based classes.")

            coh_excluded = 0
            coh_loser_set = (self._coherence_loser_by_tag
                             .get(fe.tag, {})
                             .get(task['cls'], set()))
            if coh_loser_set:
                keep_mask = np.array([
                    (round(float(rt[i]), 6), round(float(mz_raw[i]), 6))
                    not in coh_loser_set
                    for i in idx
                ])
                coh_excluded = int((~keep_mask).sum())
                idx = idx[keep_mask]
                if len(idx) == 0:
                    raise ValueError(
                        f"All rows in RT [{task['rt_start']}, {task['rt_end']}] "
                        f"were excluded by Coherence attribution as losers.")

        # ── 出力ファイル生成 ─────────────────────────────────────
        sub     = fe.df.iloc[idx].reset_index(drop=True)
        mz_sub  = make_unique_mz(mz_raw[idx])
        s_cols  = fe.sample_columns()
        now     = datetime.datetime.now()

        header = "\t".join(
            ["ID", "Ret. Time", "m/z",
             "Biotransformations/Adducts", "Included", "Saturated"]
            + s_cols)
        info   = f"LipidZoner v{__version__}\t{now:%Y-%m-%d %H:%M:%S}"

        # Ret. Time に実測 RT を書く。LipidQuant LoadRawFiles.m は
        #   この列を読むが保存も使用もしない(ローカル変数 rt のまま捨てる)
        #   ので、LipidQuant 側の変更は要らない。実値にしておくと TXT 単体で
        #   「どの feature だったか」「同じ m/z で RT 違いの行があるか」が
        #   追える。空欄にすると str2num("") が [] になって落ちるので不可。
        rows = []
        for i in range(len(sub)):
            try:
                _rt_v = float(rt[idx[i]])
            except (IndexError, TypeError, ValueError):
                _rt_v = 0.0
            if not np.isfinite(_rt_v):
                _rt_v = 0.0
            vals = ["", f"{_rt_v:.4f}", f"{mz_sub[i]:.6f}", "", "Yes", "No"]
            vals += [str(sub[sc].iloc[i]) for sc in s_cols]
            rows.append("\t".join(vals))

        content = "\n".join(
            ["Feature table", "", info, "", header] + rows)

        fname = (f"{task['cls']}_"
                 f"{task['rt_start']:.2f}-{task['rt_end']:.2f}_"
                 f"{task['ion']}_"
                 f"{now:%Y%m%d}_{now:%H%M}_{now:%S}.txt")
        outdir  = self._output_dir or fe.path.parent
        outpath = outdir / fname
        outpath.write_text(content, encoding="utf-8")

        # ログ
        log_msg = f"Exported: {outpath}  ({len(idx)} kept)"
        if used_v53_path:
            details = []
            for k, v in excluded_v53.items():
                if v:
                    details.append(f"{v} {k}")
            if details:
                log_msg += "  [excluded: " + ", ".join(details) + "]"
        else:
            if excluded_count:
                log_msg += f"  ({excluded_count} reserved excluded)"
            if coh_excluded:
                log_msg += f"  ({coh_excluded} Coherence-attributed loser excluded)"
        self._lbl_log.setText(log_msg)


# ════════════════════════════════════════════════════════════════════
#  エントリーポイント
# ════════════════════════════════════════════════════════════════════

# ════════════════════════════════════════════════════════════════════
#  JSON Session v2 helpers
# ════════════════════════════════════════════════════════════════════

LZ_SESSION_FORMAT_VERSION = 2


# ════════════════════════════════════════════════════════════════════
#  session サイドカー
#  ----
#  PerSampleFileEntry のアライメント結果を session JSON の隣に置く。
#  JSON はテキストのまま(差分が読める)にしておきたいので、DataFrame は
#  別ファイルにする。
#
#  形式は pandas の pickle。CSV だと dtype が往復で変わり
#  (iso_group_id が float 化する等)、「保存した時と同じ状態に戻す」
#  という目的を満たせない。parquet は pyarrow が要るので、依存を
#  増やさない pickle を選んだ。読めなかった場合は従来の
#  「フォルダから再アライメント」に落ちるので、壊れても致命傷には
#  ならない。
# ════════════════════════════════════════════════════════════════════

SESSION_DATA_SUFFIX = "_data"


def session_sidecar_dir(session_path) -> Path:
    """session JSON に対応するサイドカーディレクトリのパス。"""
    p = Path(session_path)
    return p.parent / (p.stem + SESSION_DATA_SUFFIX)


def _sidecar_safe_name(tag: str) -> str:
    """'#1' のようなタグをファイル名に使える形にする。"""
    out = ''.join(ch if (ch.isalnum() or ch in '-_') else '_'
                  for ch in str(tag))
    return out or 'entry'


def _sidecar_plain_dtypes(df: "pd.DataFrame") -> "pd.DataFrame":
    """pickle に入れる前に、素の numpy dtype へ落とす。

    pandas 3.x は文字列列の既定 dtype を pyarrow バックエンドの `str` に
    する。pickle は dtype ごと保存するので、そのまま書くと**読む側に
    pyarrow が必須**になる。pyarrow の無い環境では unpickle が例外に
    なり、従来経路(フォルダから TXT を読み直す)へ落ちてしまう。
    raw mzML 由来のエントリはそこに TXT が無いので、fix24 で直した
    「エントリが黙って消える」に逆戻りする。

    実際、検証用に保存した session を pyarrow 無しの環境で
    開いたら、この経路で 2 エントリとも復元されなかった。

    float64 / int64 / bool は pyarrow を要求しないのでそのまま残す。
    """
    out = df.copy()
    # 列名・行ラベルの Index も pandas 3 では dtype='str'(pyarrow 由来)に
    # なる。中身の dtype だけ直しても、ここが残ると pickle は pyarrow を
    # 参照し続ける(実際そうなっていた)。
    for _attr in ('columns', 'index'):
        _idx = getattr(out, _attr)
        if str(_idx.dtype) in ('str', 'string') or 'arrow' in str(_idx.dtype).lower():
            setattr(out, _attr, pd.Index(list(_idx), dtype=object))
    for col in out.columns:
        dt = out[col].dtype
        name = str(dt)
        if (name in ('str', 'string')
                or 'arrow' in name.lower()
                or getattr(dt, 'storage', None) is not None
                or (hasattr(dt, 'na_value') and name not in (
                    'float64', 'int64', 'bool'))):
            try:
                out[col] = out[col].astype(object).where(
                    out[col].notna(), None)
            except Exception:
                out[col] = out[col].astype(object)
    return out


def write_session_sidecar(host, session_path, data: dict) -> int:
    """PerSampleFileEntry の中身をサイドカーに書き出し、data に印を付ける。

    Returns
    -------
    int : 書き出したエントリ数
    """
    entries = [fe for fe in getattr(host, '_file_entries', [])
               if isinstance(fe, PerSampleFileEntry)]
    if not entries:
        return 0
    d = session_sidecar_dir(session_path)
    d.mkdir(parents=True, exist_ok=True)
    written = {}
    for fe in entries:
        fn = f"{_sidecar_safe_name(fe.tag)}_aligned.pkl"
        payload = {
            'format': 1,
            'app_version': __version__,
            # pyarrow バックエンドの dtype を落としてから書く
            'aligned_df': _sidecar_plain_dtypes(fe._aligned_df),
            'sample_names': list(getattr(fe, '_sample_names', []) or []),
            'ion_mode': fe.ion_mode,
            'tag': fe.tag,
            'source_folder': getattr(fe, 'source_folder', None),
            'raw_paths': list(getattr(fe, 'raw_paths', None) or []),
            'align_params': dict(getattr(fe, 'align_params', None) or {}),
            'align_info': dict(getattr(fe, 'align_info', None) or {}),
        }
        pd.to_pickle(payload, d / fn)
        written[fe.tag] = fn
    for ent in data.get('files', []):
        fn = written.get(ent.get('tag'))
        if fn:
            ent['aligned_data'] = fn
    data['session_data_dir'] = d.name
    return len(written)


def dedupe_entries_by_polarity(files: list) -> tuple[list, list]:
    """極性ごとに 1 エントリだけ残す。

    Annotation は極性ごとに 1 つの per-sample データしか使えない
    (_load_mode_state は極性ごとに最初の一致を採る)。fix30 以前は
    Send のたびに append していたので、同じ極性のエントリが増えることが
    あった。読み込み時にここで 1 件へ寄せる。

    **最後のエントリ**(= 最後に Send した最新の作業)を採用する。
    Send 側の「後から来たもので置き換える」規則と揃えてある。

    戻り値: (残すエントリのリスト, 落としたものの説明のリスト)
    """
    if not files:
        return [], []
    last_of: dict = {}
    for ent in files:
        if not isinstance(ent, dict):
            continue
        pol = ent.get('ion_mode') or 'pos'
        last_of[pol] = ent

    def _desc(ent) -> str:
        ai = ent.get('align_info') or {}
        dp = ai.get('detect_params') or {}
        bits = [str(ent.get('tag', '?'))]
        if ai.get('n_features') is not None:
            bits.append(f"{ai['n_features']:,} features")
        if dp.get('noise_threshold_int') is not None:
            bits.append(f"noise {dp['noise_threshold_int']:g}")
        if dp.get('mass_error_ppm') is not None:
            bits.append(f"{dp['mass_error_ppm']:g} ppm")
        return "  ".join(bits)

    kept, dropped = [], []
    for ent in files:
        if not isinstance(ent, dict):
            continue
        pol = ent.get('ion_mode') or 'pos'
        if ent is last_of.get(pol):
            kept.append(ent)
        else:
            dropped.append(
                f"[{pol}] {_desc(ent)}  →  kept: {_desc(last_of[pol])}")
    return kept, dropped


def read_session_sidecar(session_path, data: dict, entry: dict, tag: str):
    """サイドカーから PerSampleFileEntry を組み直す。

    見つからない / 読めない場合は None を返す(呼び出し側で従来経路へ)。
    """
    fn = entry.get('aligned_data')
    if not fn:
        return None
    base = Path(session_path).parent
    dirname = data.get('session_data_dir') or session_sidecar_dir(session_path).name
    p = base / dirname / fn
    if not p.exists():
        # session フォルダごと移動された場合に備えて JSON と同階層も見る
        alt = base / fn
        if not alt.exists():
            return None
        p = alt
    payload = pd.read_pickle(p)
    if not isinstance(payload, dict) or 'aligned_df' not in payload:
        raise ValueError(f"unexpected sidecar payload: {p}")
    fe = PerSampleFileEntry(
        aligned_df=payload['aligned_df'],
        sample_names=payload.get('sample_names') or [],
        ion_mode=entry.get('ion_mode') or payload.get('ion_mode') or 'pos',
        tag=tag,
        source_folder=(entry.get('path')
                       or payload.get('source_folder') or None),
    )
    fe.align_params = dict(entry.get('align_params')
                           or payload.get('align_params') or {})
    fe.align_info = dict(entry.get('align_info')
                         or payload.get('align_info') or {})
    fe.raw_paths = list(entry.get('raw_paths')
                        or payload.get('raw_paths') or []) or None
    return fe


def _build_lz_v2_dict(host, override_params=None, override_quant_ion=None,
                      override_filter_exemptions=None) -> dict:
    """v2 形式の JSON dict を MainWindow + 開いている AP ダイアログから構築。

    Args:
        host: MainWindow instance
        override_params: dict | None — global parameters (ppm/mad/is/adduct/sigma)
        override_quant_ion: dict | None — class → 'pos'/'neg'/'skip'
        override_filter_exemptions: dict | None — filter → [class names]

    Returns:
        v2 形式の dict(json.dumps 可能)
    """
    import datetime as _dt

    # 開いている AP ダイアログから per-mode state を refresh
    try:
        for dlg in (host._preview_dialogs or {}).values():
            if dlg is None or not hasattr(dlg, '_current_preview_settings_for_mode'):
                continue
            active_mode = (getattr(dlg, '_active_mode', None)
                           or getattr(dlg.fe, 'ion_mode', None))
            modes_to_refresh = set()
            if active_mode:
                modes_to_refresh.add(active_mode)
            modes_to_refresh.update(
                (getattr(dlg, '_mode_state', {}) or {}).keys())
            for m in modes_to_refresh:
                try:
                    s = dlg._current_preview_settings_for_mode(m)
                    if s:
                        host._preview_settings_by_mode[m] = s
                except Exception as e:
                    log.warning(f"[_build_lz_v2_dict] mode {m} refresh failed: {e}")
    except Exception as e:
        log.warning(f"[_build_lz_v2_dict] preview refresh outer error: {e}")

    # Params / quant_ion / filter_exemptions を集める(override 優先)
    if (override_params is None or override_quant_ion is None
            or override_filter_exemptions is None):
        params_d, qi_d, fe_d = host._gather_state_from_dialogs()
        if override_params is None:
            override_params = params_d
        if override_quant_ion is None:
            override_quant_ion = qi_d
        if override_filter_exemptions is None:
            override_filter_exemptions = fe_d

    # global library_path: 開いている AP の self._lib_path から取得
    library_path = None
    for dlg in (host._preview_dialogs or {}).values():
        if dlg is None:
            continue
        lp = getattr(dlg, '_lib_path', None)
        if lp:
            library_path = str(lp)
            break
    # フォールバック: per-mode 設定の library_path
    if library_path is None:
        for m, s in (host._preview_settings_by_mode or {}).items():
            if isinstance(s, dict) and s.get('library_path'):
                library_path = s['library_path']
                break

    # per-mode: class_ref_rt / is_filter_choices / manual_winner_coords /
    # manual_overrides / pipeline_flags を抽出
    modes_dict = {}
    for m, s in (host._preview_settings_by_mode or {}).items():
        if not isinstance(s, dict):
            continue
        modes_dict[m] = {
            'class_ref_rt': s.get('class_ref_rt') or {},
            'is_filter_choices': s.get('is_filter_choices') or {},
            'manual_winner_coords': s.get('manual_winner_coords') or [],
            'manual_overrides': s.get('manual_overrides') or [],
            'pipeline_flags': s.get('pipeline_state') or {},
        }

    return {
        'lipidzoner_session': {
            'format_version': LZ_SESSION_FORMAT_VERSION,
            'app_version': __version__,
            'saved_at': _dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        },
        'files': [
            # per-sample は実フォルダパス + per_sample フラグを保存し、
            # 読込時に create_per_sample_file_entry で再ロードできるようにする。
            # そのエントリを作ったときのアライメント条件も一緒に持つ。
            #   これが無いと、session を読み直しても常に既定値で再アライメント
            #   されてしまい、同じ結果にならなかった。
            ({'tag': fe.tag, 'ion_mode': fe.ion_mode, 'per_sample': True,
              'path': str(getattr(fe, 'source_folder', None) or fe.path),
              # Candidate Picker の EIC 用。生データが無い経路では空。
              'raw_paths': list(getattr(fe, 'raw_paths', None) or []),
              'align_params': dict(getattr(fe, 'align_params', None)
                                   or DEFAULT_ALIGN_PARAMS),
              'align_info': {k: v for k, v in
                             (getattr(fe, 'align_info', None) or {}).items()
                             if k in ('n_samples', 'n_peaks_in', 'n_features',
                                      'n_monoisotopic', 'input_mode',
                                      'detect_params')}}
             if isinstance(fe, PerSampleFileEntry)
             else {'tag': fe.tag, 'path': str(fe.path),
                   'ion_mode': fe.ion_mode})
            for fe in host._file_entries
        ],
        # ① タブの現在値(次に Run するときの既定として復元する)
        'alignment': _collect_alignment_section(host),
        'output_dir': (str(host._output_dir)
                       if host._output_dir else None),
        'library_path': library_path,
        'parameters': override_params,
        'modes': modes_dict,
        'quant_ion_choices': override_quant_ion,
        'filter_exemptions': override_filter_exemptions,
        'tasks': host._collect_tasks(),
        'reserved_coords_by_tag': {
            tag: [list(coord) for coord in coords]
            for tag, coords in host._reserved_coords_by_tag.items()
        },
    }


def _collect_alignment_section(host) -> dict:
    """session に載せる 'alignment' セクションを組み立てる。

    開いている AP(① タブ)があればその UI 値を、無ければ
    host._align_params と既定の検出パラメータを使う。
    """
    for dlg in (getattr(host, '_preview_dialogs', None) or {}).values():
        if dlg is None:
            continue
        if hasattr(dlg, '_al_collect_session_section'):
            try:
                return dlg._al_collect_session_section()
            except Exception as e:
                log.warning(f"[_collect_alignment_section] {e}")
    return {
        'params': dict(getattr(host, '_align_params', None)
                       or DEFAULT_ALIGN_PARAMS),
        'detect': dict(DEFAULT_DETECT_PARAMS),
        'input_mode': 'txt',
        'ion_mode': 'pos',
    }


def _load_lz_v2_file(path: str) -> dict:
    """v2 JSON ファイルを読み込んで dict を返す。

    v2 形式でない(または lipidzoner_session キーがない)場合は ValueError を投げる。
    古いバージョン(v1: _lipidzoner_unified)は受け付けない(意図的に破壊)。
    """
    import json as _j
    with open(path, 'r', encoding='utf-8') as f:
        data = _j.load(f)
    if not isinstance(data, dict):
        raise ValueError("JSON root must be a dict")
    meta = data.get('lipidzoner_session')
    if not isinstance(meta, dict):
        raise ValueError(
            "Not a v2 LipidZoner session file (missing 'lipidzoner_session' "
            "key). Old v1 format is not supported.")
    fv = meta.get('format_version')
    if fv != LZ_SESSION_FORMAT_VERSION:
        raise ValueError(
            f"Unsupported format_version: {fv}. Expected "
            f"{LZ_SESSION_FORMAT_VERSION}.")
    return data


def _v2_to_per_mode_settings(data: dict) -> dict:
    """v2 dict から per-mode settings dict を構築(内部 API 互換性のため)。

    返値: {mode: {library_path, ion_mode, ppm_tol, mad_k, is_tol,
                  class_ref_rt, is_filter_choices, pipeline_state,
                  manual_winner_coords, manual_overrides}}
    これは PreviewRTDialog._restore_settings が期待する形式。
    """
    library_path = data.get('library_path')
    params = data.get('parameters') or {}
    modes = data.get('modes') or {}
    result = {}
    for m, md in modes.items():
        if not isinstance(md, dict):
            continue
        result[m] = {
            'library_path': library_path,
            'ion_mode': m,
            'ppm_tol': params.get('ppm_tol'),
            'mad_k': params.get('mad_k'),
            'is_tol': params.get('is_tol'),
            'class_ref_rt': md.get('class_ref_rt') or {},
            'is_filter_choices': md.get('is_filter_choices') or {},
            'pipeline_state': md.get('pipeline_flags') or {},
            'manual_winner_coords': md.get('manual_winner_coords') or [],
            'manual_overrides': md.get('manual_overrides') or [],
        }
    return result


def main():
    import os
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    # Windows のタスクバーでアイコンを正しく表示するための設定
    # (Python 実行ファイルではなく、独立したアプリとして認識させる)
    if sys.platform == "win32":
        try:
            from ctypes import windll
            windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                "LipidZoner.1.1.0"
            )
        except Exception:
            pass

    # 起動時バリデーション
    for label, fn in (
        ("CURATED_ADDUCTS", _validate_curated_adducts),
        ("CURATED_ADDUCT_FINGERPRINTS", _validate_curated_adduct_fingerprints),
    ):
        wmsgs = fn()
        if wmsgs:
            log.info(f"[{label}] warnings at startup:")
            for w in wmsgs:
                log.info(f"  {w}")

    app = QApplication(sys.argv)
    app.setFont(QFont("Arial", UI_FONT_POINT_SIZE))

    # アプリケーション・アイコン設定
    # PyInstaller でパッケージ化された場合は sys._MEIPASS の下から、
    # 通常の Python 実行時はスクリプトと同じディレクトリから読む。
    if getattr(sys, "frozen", False):
        icon_path = os.path.join(sys._MEIPASS, "LipidZoner_icon.ico")
    else:
        icon_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "LipidZoner_icon.ico",
        )
    if os.path.exists(icon_path):
        app.setWindowIcon(QIcon(icon_path))

    w = MainWindow()
    # 空の AP を直接起動。データ選択は AP 内の [Load data…] から。
    # MainWindow は hidden のままバックエンドとして動作。
    stub = StubFileEntry()
    # parent=w は PreviewRTDialog 内で self._host に保持される。
    # 実 Qt parent は None になり、AP は top-level window として
    # 独立したタスクバーボタンを持つ。
    dlg = PreviewRTDialog(stub, [], parent=w)
    dlg.task_ready.connect(w._insert_task)
    dlg.settings_saved.connect(
        lambda settings, m=stub.ion_mode:
            w._on_preview_settings_saved(m, settings))
    dlg.reserved_updated.connect(w._on_reserved_updated)
    dlg.coherence_updated.connect(w._on_coherence_updated)
    dlg.match_df_updated.connect(w._on_match_df_updated)
    dlg.setAttribute(Qt.WA_DeleteOnClose, False)
    w._preview_dialogs[stub.tag] = dlg
    dlg.show()
    dlg.raise_()
    dlg.activateWindow()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()