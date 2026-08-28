# ABOUTME: Read-only pre-flight checks on a relaxation dataset before the pipeline commits to a run.
# ABOUTME: Labels and identity come from the pipeline's own parsers; intensities are cheap proxies.

"""Pre-flight inspection of spectra against their peak list.

A wrong or mis-registered peak list does not raise — it produces zero heights and
shoulder fits that read as noisy data, after several minutes of integration. These
checks run in seconds and answer the question first.

Two rules shape the module:

* Anything about **labels or identity** — delay values, column names, peak-list
  contents — is delegated to the code the pipeline itself uses (`DelayExtractor`,
  `NMRFileManager`). A pre-flight check that parses filenames its own way will
  eventually certify something the run cannot reproduce.
* Anything about **intensity** is a deliberate approximation: the maximum in a small
  box, not a fit. Fitting properly would cost as much as the run it is meant to
  precede.
"""

import glob
import os
import re
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

import numpy as np

SEARCH_WINDOW = 0.070          # ppm, the detection default in both dimensions
_STRONG_SNR = 20.0             # peaks this far above noise carry the intensity comparisons
_DETECT_SNR = 10.0             # a local maximum this far above noise counts as detected
_DECAYED_FRACTION = 0.10       # this many peaks under the detection threshold = decay setting in
_CAPTURE_FLOOR = 0.85


class Finding(NamedTuple):
    """One problem worth a human's attention."""
    severity: str              # 'FAIL' | 'WARN'
    check: str
    message: str


# --------------------------------------------------------------------------- spectra

def _read(path):
    import nmrglue as ng
    dic, data = ng.pipe.read(str(path))
    data = np.asarray(data)
    return (dic, data,
            ng.pipe.make_uc(dic, data, dim=0),
            ng.pipe.make_uc(dic, data, dim=1))


def estimate_noise(data: np.ndarray) -> float:
    """Median-absolute-deviation noise from the quietest corner.

    A whole-spectrum MAD reports the signal on a crowded spectrum, so the corners
    are measured separately and the smallest is taken.
    """
    corners = [data[:80, :150], data[-80:, :150], data[:80, -150:], data[-80:, -150:]]
    mad = float(min(np.median(np.abs(c - np.median(c))) * 1.4826 for c in corners))
    if mad > 0:
        return mad
    # A corner with no noise at all turns every S/N into a division by zero. Floor it
    # far below any real signal so the ratios stay finite and read as "no noise".
    peak = float(np.max(np.abs(data)))
    return peak * 1e-9 if peak > 0 else 1.0


def peak_heights(path, peaks: Sequence[Tuple[str, float, float]]) -> Tuple[np.ndarray, float]:
    """Maximum in a 7x7 box at each peak, plus the noise level. Not a fit — see module docstring."""
    _, data, uy, ux = _read(path)
    heights = np.array([
        data[max(0, uy(y, 'ppm') - 3):uy(y, 'ppm') + 4,
             max(0, ux(x, 'ppm') - 3):ux(x, 'ppm') + 4].max()
        for _, x, y in peaks])
    return heights, estimate_noise(data)


def assess_registration(path, peaks: Sequence[Tuple[str, float, float]],
                        quick: bool = False) -> Tuple[float, float]:
    """Offset (d_1H, d_15N) in ppm that best aligns the peak list with the spectrum.

    Slides the whole list over a grid and keeps the shift maximising summed intensity
    at the peak positions. A list belonging to a different experiment shows up as a
    consistent offset far larger than the search window.
    """
    _, data, uy, ux = _read(path)
    ny, nx = data.shape
    span, step = (0.3, 0.03) if quick else (0.9, 0.015)
    best = (-np.inf, 0.0, 0.0)
    for dy in np.arange(-span, span + 1e-9, step):
        rows = [uy(y + dy, 'ppm') for _, _, y in peaks]
        for dx in np.arange(-0.09, 0.0901, 0.0025):
            total = 0.0
            for (_, x, _), iy in zip(peaks, rows):
                ix = ux(x + dx, 'ppm')
                if 0 <= iy < ny and 0 <= ix < nx:
                    total += data[iy, ix]
            if total > best[0]:
                best = (total, dx, dy)
    return float(best[1]), float(best[2])


def _capture(path, peaks) -> Dict:
    """How many listed peaks have a real local maximum within the search window."""
    from scipy.ndimage import maximum_filter
    _, data, uy, ux = _read(path)
    noise = estimate_noise(data)
    local = (data == maximum_filter(data, size=5)) & (data > _DETECT_SNR * noise)
    iy, ix = np.where(local)
    found_y = np.array([uy.ppm(i) for i in iy])
    found_x = np.array([ux.ppm(i) for i in ix])
    captured = sum(1 for _, x, y in peaks
                   if ((np.abs(found_x - x) <= SEARCH_WINDOW)
                       & (np.abs(found_y - y) <= SEARCH_WINDOW)).any())
    snr = np.array([
        data[max(0, uy(y, 'ppm') - 3):uy(y, 'ppm') + 4,
             max(0, ux(x, 'ppm') - 3):ux(x, 'ppm') + 4].max() / noise
        for _, x, y in peaks])
    return dict(captured=captured, n_peaks=len(peaks), noise=noise,
                median_snr=float(np.median(snr)), weak=int((snr < _DETECT_SNR).sum()))


# ------------------------------------------------------------------- peaks and delays

def load_peak_list(path) -> Dict:
    """Load with the pipeline's own reader, then judge hygiene on what it returns."""
    from lunaNMR.utils.file_manager import NMRFileManager
    df = NMRFileManager().load_peak_list(str(path))
    if df is None or df.empty:
        return dict(path=str(path), assignments=[], peaks=[], findings=[
            Finding('FAIL', 'peak_list', f'{os.path.basename(path)}: no peaks could be read')])

    names = [str(a) for a in df['Assignment']]
    peaks = list(zip(names, df['Position_X'].astype(float), df['Position_Y'].astype(float)))

    findings = []
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        findings.append(Finding('FAIL', 'peak_list',
                                f'{os.path.basename(path)}: duplicate assignments {duplicates}'))
    dummies = sum(1 for n in names if n.lower().startswith('dummy'))
    if dummies:
        findings.append(Finding('WARN', 'peak_list',
                                f'{os.path.basename(path)}: {dummies} dummy_* entries '
                                f'(relaxation fitters exclude these)'))
    return dict(path=str(path), assignments=names, peaks=peaks, findings=findings)


def read_delays(spectrum_names: Sequence[str], series_mode: str = 'time') -> Dict:
    """Parse delays with the extractor `series` uses, so the two cannot disagree."""
    from lunaNMR.utils.delay_extractor import DelayExtractor
    extractor = DelayExtractor(mode=series_mode)
    parsed = {n: extractor.extract_value(n) for n in spectrum_names}
    unparsed = sorted(n for n, v in parsed.items() if v is None)
    values = sorted({v for v in parsed.values() if v is not None})
    counts = [v for v in parsed.values() if v is not None]
    findings = []
    if unparsed and len(unparsed) != len(parsed):
        findings.append(Finding(
            'WARN', 'delay_parsing',
            f'{len(unparsed)} of {len(parsed)} spectra have no parseable delay '
            f'({", ".join(unparsed)}); their matrix columns will be named after the file stem'))
    return dict(parsed=parsed, values=values, unparsed=unparsed,
                n_parsed=len(parsed) - len(unparsed), n_unparsed=len(unparsed),
                repeats={v: counts.count(v) for v in values if counts.count(v) > 1},
                findings=findings)


# ------------------------------------------------------------------- derived checks

def _hetnoe_planes(paths, peaks) -> Optional[Dict]:
    """Decide which of two planes is saturated by measuring, never by filename."""
    if len(paths) != 2:
        return None
    h1, n1 = peak_heights(paths[0], peaks)
    h2, _ = peak_heights(paths[1], peaks)
    strong = (h1 > _STRONG_SNR * n1) & (h2 > _STRONG_SNR * n1)
    if strong.sum() < 5:
        return None
    raw = float(np.median(h2[strong] / h1[strong]))
    sat, unsat = (paths[1], paths[0]) if raw < 1 else (paths[0], paths[1])
    ratio = min(raw, 1 / raw)
    findings = []
    if not 0.6 <= ratio <= 0.95:
        findings.append(Finding('WARN', 'hetnoe',
                                f'sat/unsat = {ratio:.3f}, outside the usual 0.75-0.85'))
    return dict(saturated=os.path.basename(sat), unsaturated=os.path.basename(unsat),
                ratio=ratio, findings=findings)


def _ordered_by_delay(spectrum_names: Sequence[str]) -> List[str]:
    """Spectra in the order `series` assigns column labels: by (delay, basename)."""
    from lunaNMR.utils.delay_extractor import DelayExtractor
    return [name for name, _, _ in
            DelayExtractor(mode='time').sort_files_with_sequence(list(spectrum_names))]


def _subseries_scale(folder, spectrum_names, delays, peaks) -> Optional[Dict]:
    """Compare repeat acquisitions at shared delays.

    Grouped by delay collision rather than by filename, so `_b_` infixes, `b` prefixes
    and any other convention are all found. Which file of a colliding pair holds the
    bare column label and which gets `_2` is taken from the same sort `series` uses,
    so the ratio is unambiguously second-relative-to-first.
    """
    from lunaNMR.utils.delay_extractor import DelayExtractor
    shared = sorted(d for d, n in delays['repeats'].items() if n == 2)
    if not shared:
        return None

    ordered = DelayExtractor(mode='time').sort_files_with_sequence(list(spectrum_names))
    by_delay: Dict[float, List[str]] = {}
    for name, delay, sequence in ordered:
        by_delay.setdefault(delay, []).append(name)

    per_delay, ratios = {}, []
    for delay in shared:
        first, second = by_delay[delay][:2]
        ha, na = peak_heights(os.path.join(folder, first), peaks)
        hb, _ = peak_heights(os.path.join(folder, second), peaks)
        strong = (ha > _STRONG_SNR * na) & (hb > _STRONG_SNR * na)
        if strong.sum() < 5:
            continue
        r = float(np.median(hb[strong] / ha[strong]))
        per_delay[delay] = dict(first=first, second=second, ratio=r, n_strong=int(strong.sum()))
        ratios.append(r)
    if not ratios:
        return None

    ratio = float(np.median(ratios))
    findings = []
    if abs(ratio - 1) > 0.05:
        findings.append(Finding(
            'FAIL', 'subseries',
            f'repeat acquisitions differ by {abs(ratio - 1) * 100:.1f}% '
            f'(ratio {ratio:.3f}; true replicates give ~0.99) — the low sub-series needs '
            f'cross-normalising x{1 / ratio:.4f}. Re-measure the ratio on fitted intensities '
            f'before applying it; this one is a box-maximum estimate.'))
    return dict(shared=sorted(per_delay), per_delay=per_delay, ratio=ratio,
                scale=1.0 / ratio, findings=findings)


# ------------------------------------------------------------------------- entry points

def _spectrum_extensions():
    """What the loader accepts, so discovery here cannot drift from `series`."""
    from lunaNMR.utils.file_manager import NMRFileManager
    return tuple(NMRFileManager().supported_nmr_formats)


SPECTRUM_EXTENSIONS = _spectrum_extensions()


def find_spectra(folder) -> List[str]:
    """Spectrum files in a folder, ordered the way `series` discovers them."""
    files: List[str] = []
    for ext in SPECTRUM_EXTENSIONS:
        files.extend(glob.glob(os.path.join(str(folder), f'*.{ext}')))
    natural = lambda f: [int(t) if t.isdigit() else t.lower()
                         for t in re.split(r'(\d+)', os.path.basename(f))]
    return sorted(files, key=natural)


def check_experiment(spectra_dir, peaks_path, *, quick: bool = False,
                     sample: bool = False, series_mode: str = 'time') -> Dict:
    """Inspect one experiment folder against its own peak list.

    `sample` assesses only the first and last spectrum; the default assesses all of
    them, because the middle of a series is where decay-driven capture loss lives.
    """
    spectra_dir = str(spectra_dir)
    paths = find_spectra(spectra_dir)
    names = [os.path.basename(p) for p in paths]
    result: Dict = dict(experiment=os.path.basename(os.path.normpath(spectra_dir)),
                        folder=spectra_dir, n_spectra=len(paths), spectra=[],
                        hetnoe=None, subseries=None, findings=[])
    if not paths:
        result['findings'].append(Finding('FAIL', 'spectra', f'no spectra in {spectra_dir}'))
        return result

    peak_list = load_peak_list(peaks_path)
    result['peak_list'] = peak_list
    result['findings'].extend(peak_list['findings'])
    peaks = peak_list['peaks']
    if not peaks:
        return result

    delays = read_delays(names, series_mode)
    result['delays'] = delays
    result['findings'].extend(delays['findings'])

    # Capture loss is only diagnostic where nothing has decayed yet. At the shortest
    # delay every listed peak should still be present, so a shortfall there means the
    # list or the data is wrong; at long delay it just means the experiment worked.
    # Ties at the shortest delay are broken by the same sort `series` uses, so the
    # reference is the file that owns the bare column label rather than the `_2` one.
    ordered = _ordered_by_delay(names)
    reference = ordered[0] if ordered else names[0]
    result['reference_spectrum'] = reference

    assessed = ([paths[0], paths[-1]] if sample and len(paths) > 1
                else [paths[0]] if sample else paths)
    for path in assessed:
        name = os.path.basename(path)
        dx, dy = assess_registration(path, peaks, quick=quick)
        cap = _capture(path, peaks)
        rate = cap['captured'] / cap['n_peaks']
        decayed = cap['weak'] / cap['n_peaks'] > _DECAYED_FRACTION
        result['spectra'].append(dict(spectrum=name, dx=dx, dy=dy, capture=cap['captured'],
                                      n_peaks=cap['n_peaks'], capture_rate=rate,
                                      median_snr=cap['median_snr'], weak=cap['weak'],
                                      decayed=decayed, is_reference=name == reference))
        offset = max(abs(dx), abs(dy)) / SEARCH_WINDOW
        if offset > 0.5:
            result['findings'].append(Finding(
                'FAIL', 'registration',
                f'{name}: peak list off by ({dx:+.4f}, {dy:+.3f}) ppm = {offset:.1f}x the '
                f'{SEARCH_WINDOW} ppm window — wrong list for these spectra'))
        elif offset > 0.15:
            result['findings'].append(Finding(
                'WARN', 'registration',
                f'{name}: peak list off by ({dx:+.4f}, {dy:+.3f}) ppm '
                f'({offset:.0%} of the window) — tolerable but check'))
        # Only the reference point carries a warning. Elsewhere the same number is decay,
        # and warning on it buries the real cases — every relaxation series has a tail.
        if name == reference and rate < _CAPTURE_FLOOR:
            result['findings'].append(Finding(
                'WARN', 'capture',
                f'{name}: only {rate:.0%} of peaks captured at the reference point '
                f'(median S/N {cap["median_snr"]:.0f}, {cap["weak"]} weak) — peaks are missing '
                f'before any decay, so check the list against these spectra'))

    if delays['n_parsed'] == 0:
        hetnoe = _hetnoe_planes(paths, peaks)
        if hetnoe:
            result['hetnoe'] = hetnoe
            result['findings'].extend(hetnoe['findings'])

    subseries = _subseries_scale(spectra_dir, names, delays, peaks)
    if subseries:
        result['subseries'] = subseries
        result['findings'].extend(subseries['findings'])
    return result


def find_peak_list(folder) -> Optional[str]:
    """The peak list belonging to an experiment lives beside its spectra."""
    lists = sorted(glob.glob(os.path.join(str(folder), '*.txt')))
    return lists[0] if lists else None


def check_dataset(root, *, quick: bool = False, sample: bool = False,
                  series_mode: str = 'time') -> Dict:
    """Inspect every experiment under a dataset root, then compare them.

    The cross-experiment comparison is the part a single-experiment check cannot do:
    residue merging is an intersection, so the smallest assignment set is the ceiling
    on every downstream result.
    """
    root = str(root)
    skip = ('/analysis/', '/lunaNMR_analysis/', '/staging/', '/series/', '/run_')
    folders = sorted({os.path.dirname(p) for ext in SPECTRUM_EXTENSIONS
                      for p in glob.glob(f'{root}/**/*.{ext}', recursive=True)
                      if not any(k in p for k in skip)})

    result: Dict = dict(root=root, experiments=[], assignments=None, findings=[])
    if not folders:
        result['findings'].append(Finding('FAIL', 'spectra', f'no spectra found under {root}'))
        return result

    sets = {}
    for folder in folders:
        peaks_path = find_peak_list(folder)
        rel = os.path.relpath(folder, root)
        if not peaks_path:
            result['findings'].append(Finding('WARN', 'peak_list', f'{rel}: no peak list in this folder'))
            continue
        experiment = check_experiment(folder, peaks_path, quick=quick, sample=sample,
                                      series_mode=series_mode)
        experiment['experiment'] = rel
        result['experiments'].append(experiment)
        result['findings'].extend(experiment['findings'])
        if experiment.get('peak_list', {}).get('assignments'):
            sets[rel] = set(experiment['peak_list']['assignments'])

    if len(sets) > 1:
        common = set.intersection(*sets.values())
        union = set.union(*sets.values())
        result['assignments'] = dict(per_experiment={k: len(v) for k, v in sets.items()},
                                     common=len(common), union=len(union))
        if len(common) != len(union):
            result['findings'].append(Finding(
                'WARN', 'assignments',
                f'assignment sets differ across experiments; merging is an intersection, '
                f'so {len(common)} residues is the ceiling (union {len(union)})'))
    return result
