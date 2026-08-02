# ABOUTME: Loads 1D NMR spectra (NMRPipe, Bruker pdata) into a OneDSpectrum with a ppm axis.
# ABOUTME: Includes the segmented-MAD noise estimator that ignores peaks and negative artefacts.

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# The spectrum is split into this many segments for noise estimation and the quietest
# fraction is kept. Signal occupies a minority of a 1D spectrum, so the quiet segments
# carry noise alone; a MAD over the whole spectrum reports the signal instead (7.6x too
# high on the reference 800 MHz spectrum).
NOISE_SEGMENTS = 32
NOISE_QUIET_FRACTION = 0.25

MAD_TO_SIGMA = 1.4826


def estimate_noise(data):
    """Noise sigma from the quietest segments of the spectrum.

    Robust to both tall peaks and inverted water-suppression artefacts, because a
    segment containing either lands outside the quiet fraction and is discarded.
    """
    data = np.asarray(data, dtype=float)
    if data.size == 0:
        return 0.0

    n_segments = min(NOISE_SEGMENTS, data.size)
    mads = np.array([np.median(np.abs(seg - np.median(seg))) * MAD_TO_SIGMA
                     for seg in np.array_split(data, n_segments)])

    keep = max(1, int(len(mads) * NOISE_QUIET_FRACTION))

    return float(np.median(np.sort(mads)[:keep]))


@dataclass
class OneDSpectrum:
    """A processed 1D spectrum on a descending ppm axis."""

    data: np.ndarray
    ppm_axis: np.ndarray
    noise: float
    metadata: dict = field(default_factory=dict)
    intensity_scale: float = 1.0

    @property
    def ppm_step(self):
        """Absolute ppm per point."""
        return abs(float(self.ppm_axis[1] - self.ppm_axis[0]))

    def hz_per_point(self):
        obs = self.metadata.get('obs_mhz')
        return self.ppm_step * obs if obs else None


def from_arrays(data, ppm_axis, noise=None, metadata=None, intensity_scale=1.0):
    """Build a OneDSpectrum from arrays, orienting the axis descending."""
    data = np.asarray(data, dtype=float)
    ppm_axis = np.asarray(ppm_axis, dtype=float)

    if data.shape != ppm_axis.shape:
        raise ValueError(f"data and ppm_axis lengths differ: {data.shape} vs {ppm_axis.shape}")

    if ppm_axis.size > 1 and ppm_axis[0] < ppm_axis[-1]:
        data, ppm_axis = data[::-1], ppm_axis[::-1]

    if noise is None:
        noise = estimate_noise(data)

    return OneDSpectrum(data=data, ppm_axis=ppm_axis, noise=float(noise),
                        metadata=dict(metadata or {}), intensity_scale=float(intensity_scale))


def _read_nmrpipe(path):
    import nmrglue as ng

    dic, data = ng.pipe.read(str(path))
    if data.ndim != 1:
        raise ValueError(f"{path} is {data.ndim}D, not a 1D spectrum")

    ppm_axis = ng.pipe.make_uc(dic, data, dim=0).ppm_scale()
    metadata = {
        'format': 'nmrpipe',
        'obs_mhz': float(dic.get('FDF2OBS', 0.0)) or None,
        'sw_hz': float(dic.get('FDF2SW', 0.0)) or None,
        'nucleus': (dic.get('FDF2LABEL') or '').strip() or None,
        'source': str(path),
    }
    return np.asarray(data, dtype=float), np.asarray(ppm_axis, dtype=float), metadata


def _read_bruker(path):
    import nmrglue as ng

    dic, data = ng.bruker.read_pdata(str(path))
    if data.ndim != 1:
        raise ValueError(f"{path} is {data.ndim}D, not a 1D spectrum")

    udic = ng.bruker.guess_udic(dic, data)
    ppm_axis = ng.fileiobase.uc_from_udic(udic, dim=0).ppm_scale()

    acqus = dic.get('acqus', {})
    metadata = {
        'format': 'bruker',
        'obs_mhz': udic[0].get('obs') or None,
        'sw_hz': udic[0].get('sw') or None,
        'nucleus': udic[0].get('label') or None,
        'scans': acqus.get('NS'),
        'receiver_gain': acqus.get('RG'),
        'source': str(path),
    }
    return np.asarray(data, dtype=float), np.asarray(ppm_axis, dtype=float), metadata


def acquisition_scale(metadata):
    """Per-spectrum factor making intensities comparable across a series.

    Bruker intensities scale with scan count and receiver gain, so a series acquired
    with different NS or RG is not directly comparable. NMRPipe headers carry neither,
    so those spectra fall back to 1.0 and need a manual factor if they differ.
    """
    scans = metadata.get('scans')
    gain = metadata.get('receiver_gain')

    if not scans or not gain:
        return 1.0

    return 1.0 / (float(scans) * float(gain))


def load_spectrum(path, intensity_scale=None):
    """Load a processed 1D spectrum from NMRPipe or Bruker pdata."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"spectrum not found: {path}")

    reader = _read_bruker if path.is_dir() or path.name in ('1r', '1i') else _read_nmrpipe

    try:
        data, ppm_axis, metadata = reader(path)
    except ValueError:
        raise
    except Exception as exc:
        raise ValueError(f"could not read {path} as a 1D spectrum: {exc}") from exc

    if intensity_scale is None:
        intensity_scale = acquisition_scale(metadata)

    return from_arrays(data, ppm_axis, metadata=metadata, intensity_scale=intensity_scale)
