"""
PS2D Data Selection Module - Exact Clone of 
===============================================================================

It replaces lunaNMR's adaptive window multiplier approach with 
fixed elliptical window selection.

Key Implementation from spectrum.cpp lines 1010-1020:
1. FIXED integration radius (radF1, radF2) - never changes during fitting
2. Elliptical boundary test: radius = (F1-pos)²/radF1² + (F2-pos)²/radF2² <= 1.0
3. Include data points inside ellipse, exclude outside
4. NO adaptive windows, NO multipliers, NO linewidth-based scaling

Default values from peak.cpp lines 181-182:
- radF1 = 0.6 ppm (indirect dimension, 15N/13C)
- radF2 = 0.06 ppm (direct dimension, 1H)

Date: 2025-10-04
Version: 1.0 - EXACT_PS2D_CLONE
"""

import numpy as np
from typing import Dict, Tuple, Optional


class Ps2dDataSelector:
    """
    Exact clone of  data selection logic from spectrum.cpp.

    This class implements FIXED elliptical windows exactly as  does.
    NO adaptive behavior, NO multipliers, NO creativity - just pure cloning.
    """

    def __init__(self, spectrum_type: str = '15N-HSQC'):
        """
        Initialize with  default values.

        Parameters
        ----------
        spectrum_type : str
            Type of spectrum ('15N-HSQC', '13C-HSQC', etc.)
        """
        # radF1(0.6), radF2(0.06), radF3(0.6)
        #
        # NMRPipe convention (from spectrum.h):
        # F1 = indirect dimension (vertical) = 15N/13C in HSQC
        # F2 = direct dimension (horizontal) = 1H in HSQC
        # F3 = 3rd dimension (for 3D spectra)

        if spectrum_type in ['15N-HSQC', '13C-HSQC', '1H-15N', '1H-13C']:
            self.radF1 = 0.6   # Indirect dimension (15N or 13C) - ppm @0.6
            self.radF2 = 0.06  # Direct dimension (1H) - ppm #0.06
        else:
            # Default to 15N-HSQC values
            self.radF1 = 0.6
            self.radF2 = 0.06

        self.spectrum_type = spectrum_type

    def select_data_elliptical(self,
                               x_data: np.ndarray,
                               y_data: np.ndarray,
                               peak_x_pos: float,
                               peak_y_pos: float,
                               dimension: str = 'x') -> Dict:
        """
        Select data using FIXED elliptical window (EXACT clone of spectrum.cpp lines 1010-1020).

        This is the CORE  logic:
        ```cpp
        Doub radius = SQR((_F1-peak[k].appF1)/peak[k].radF1) +
                      SQR((_F2-peak[k].appF2)/peak[k].radF2);
        if (radius <= 1.0) {
            // Include this data point
        }
        ```

        Parameters
        ----------
        x_data : np.ndarray
            X-axis data (ppm scale)
        y_data : np.ndarray
            Y-axis intensity data
        peak_x_pos : float
            Peak position in x dimension (ppm)
        peak_y_pos : float or None
            Peak position in y dimension (ppm), only for 2D
        dimension : str
            'x' for 1D horizontal slice, 'y' for 1D vertical slice

        Returns
        -------
        dict : Selected data with keys:
            'x_selected': selected x data
            'y_selected': selected y data
            'mask': boolean mask of selected points
            'radius_used': the radF value used (for logging)
        """
        # For 1D slices, use appropriate radius
        if dimension == 'x':
            # X-dimension (1H in HSQC) uses radF2
            radius_ppm = self.radF2
            peak_pos = peak_x_pos
        else:
            # Y-dimension (15N/13C in HSQC) uses radF1
            radius_ppm = self.radF1
            peak_pos = peak_y_pos

        # For 1D: radius = ((x - peak_pos) / radF)²
        # Include if radius <= 1.0
        radius_squared = ((x_data - peak_pos) / radius_ppm) ** 2
        mask = radius_squared <= 1.0

        # Select data inside ellipse
        x_selected = x_data[mask]
        y_selected = y_data[mask]

        return {
            'x_selected': x_selected,
            'y_selected': y_selected,
            'mask': mask,
            'radius_used': radius_ppm,
            'n_points_selected': np.sum(mask),
            'window_ppm': radius_ppm  # FIXED value, not adaptive
        }

    def select_data_2d_elliptical(self,
                                   f1_data: np.ndarray,
                                   f2_data: np.ndarray,
                                   intensity: np.ndarray,
                                   peak_f1_pos: float,
                                   peak_f2_pos: float) -> Dict:
        """
        Select 2D data using FIXED elliptical window (EXACT clone of spectrum.cpp).

        ```cpp
        Doub radius = SQR((_F1-peak[k].appF1)/peak[k].radF1) +
                      SQR((_F2-peak[k].appF2)/peak[k].radF2);
        if (radius <= 1.0) {
            f1[sizecounter] = _F1;
            f2[sizecounter] = _F2;
            y[sizecounter] = intensity;
            ++sizecounter;
        }
        ```

        Parameters
        ----------
        f1_data : np.ndarray
            F1 dimension data (indirect, 15N/13C ppm)
        f2_data : np.ndarray
            F2 dimension data (direct, 1H ppm)
        intensity : np.ndarray
            Intensity values
        peak_f1_pos : float
            Peak position in F1 (ppm)
        peak_f2_pos : float
            Peak position in F2 (ppm)

        Returns
        -------
        dict : Selected data following  structure
        """
        radius_squared = ((f1_data - peak_f1_pos) / self.radF1) ** 2 + \
                        ((f2_data - peak_f2_pos) / self.radF2) ** 2

        # Include points where radius <= 1.0 (inside ellipse)
        mask = radius_squared <= 1.0

        return {
            'f1_selected': f1_data[mask],
            'f2_selected': f2_data[mask],
            'intensity_selected': intensity[mask],
            'mask': mask,
            'radF1_used': self.radF1,
            'radF2_used': self.radF2,
            'n_points_selected': np.sum(mask)
        }

    def get_window_boundaries(self, peak_pos: float, dimension: str = 'x') -> Tuple[float, float]:
        """
        Get window boundaries for plotting/display (from spectrum.cpp lines 1584-1587).

         calculation:
        ```cpp
        imin = (Int)(((peak.f1+peak.radF1)*obsF1 - firstF1)/deltaF1 + 0.5);
        imax = (Int)(((peak.f1-peak.radF1)*obsF1 - firstF1)/deltaF1 + 0.5);
        ```

        Simplified for ppm scale: [peak_pos - radius, peak_pos + radius]

        Parameters
        ----------
        peak_pos : float
            Peak position (ppm)
        dimension : str
            'x' or 'y'

        Returns
        -------
        tuple : (min_ppm, max_ppm)
        """
        if dimension == 'x':
            radius = self.radF2
        else:
            radius = self.radF1

        return (peak_pos - radius, peak_pos + radius)

    def get_integration_info(self) -> Dict:
        """
        Get current integration window settings (for logging/debugging).

        Returns
        -------
        dict : Window configuration matching 
        """
        return {
            'radF1_ppm': self.radF1,
            'radF2_ppm': self.radF2,
            'spectrum_type': self.spectrum_type,
            'source': ' lines 181-182',
            'method': 'FIXED elliptical window (no adaptation)',
            'f1_dimension': '15N/13C (indirect)',
            'f2_dimension': '1H (direct)'
        }


def select_data_2d_for_overlap_group(f1_data: np.ndarray,
                                      f2_data: np.ndarray,
                                      intensity: np.ndarray,
                                      peak_positions: list,
                                      radF1: float = 0.6,
                                      radF2: float = 0.06) -> Dict:
    """
    Select 2D data for overlap group using UNION of elliptical windows
    (EXACT clone of spectrum.cpp lines 1009-1020)

    This is the CORE function for 2D multi-peak fitting in 
    For an overlap group with N peaks, select the union of all elliptical
    windows. Each data point is included if it falls inside ANY peak's ellipse.

    Key C++ logic (spectrum.cpp:1009-1020):
    ```cpp
    for (Uint k=0; k<peak.size(); ++k) {
        Doub radius = SQR((_F1-peak[k].appF1)/peak[k].radF1) +
                      SQR((_F2-peak[k].appF2)/peak[k].radF2);
        if (radius <= 1.0) {
            f1[sizecounter] = _F1;
            f2[sizecounter] = _F2;
            y[sizecounter] = _yvec[j-jmin];
            ++sizecounter;
            break;  // CRITICAL: Only include each point once
        }
    }
    ```

    Parameters:
    -----------
    f1_data : np.ndarray
        F1 dimension data (indirect, 15N/13C ppm), can be 1D or flattened 2D
    f2_data : np.ndarray
        F2 dimension data (direct, 1H ppm), same shape as f1_data
    intensity : np.ndarray
        Intensity values, same shape as f1_data
    peak_positions : list of tuples
        List of (f1_pos, f2_pos) for each peak in overlap group
    radF1 : float
        Ellipse radius in F1 dimension (default 0.6 ppm from peak.cpp:181)
    radF2 : float
        Ellipse radius in F2 dimension (default 0.06 ppm from peak.cpp:182)

    Returns:
    --------
    dict : Selected data with keys:
        'f1_selected': F1 values inside union
        'f2_selected': F2 values inside union
        'intensity_selected': Intensity values inside union
        'mask': Boolean mask showing which points were selected
        'n_points_selected': Number of selected points
        'radF1_used': radF1 value used
        'radF2_used': radF2 value used
    """

    # Flatten inputs if needed
    f1_flat = f1_data.ravel()
    f2_flat = f2_data.ravel()
    intensity_flat = intensity.ravel()

    n_points = len(f1_flat)
    mask = np.zeros(n_points, dtype=bool)

    # For each point, check if it's inside ANY peak's ellipse
    # (spectrum.cpp:1009-1020)
    for i in range(n_points):
        point_f1 = f1_flat[i]
        point_f2 = f2_flat[i]

        # Check against all peaks in overlap group
        for peak_f1, peak_f2 in peak_positions:
            # Elliptical boundary test (spectrum.cpp:1011-1012)
            radius_squared = ((point_f1 - peak_f1) / radF1) ** 2 + \
                           ((point_f2 - peak_f2) / radF2) ** 2

            if radius_squared <= 1.0:
                # Inside this peak's ellipse → include point
                mask[i] = True
                break  # CRITICAL: "break" prevents duplicate counting (line 1017)

    # Extract selected data
    f1_selected = f1_flat[mask]
    f2_selected = f2_flat[mask]
    intensity_selected = intensity_flat[mask]

    return {
        'f1_selected': f1_selected,
        'f2_selected': f2_selected,
        'intensity_selected': intensity_selected,
        'mask': mask,
        'n_points_selected': np.sum(mask),
        'radF1_used': radF1,
        'radF2_used': radF2
    }
