# ABOUTME: Calculator for T2 relaxation times from T1 and T1rho measurements
# ABOUTME: Implements per-residue theta correction for off-resonance T1rho experiments

import numpy as np
import pandas as pd


# 15N/1H gyromagnetic ratio factor (15N frequency = 1H frequency / this factor)
N15_GAMMA_RATIO = 9.869


def calculate_per_residue_theta(
    residue_ppm: float,
    carrier_ppm: float,
    omega1_hz: float,
    theta_nominal_deg: float,
    spec_freq_mhz: float
) -> float:
    """
    Calculate the effective tilt angle for a residue based on its chemical shift.

    The spin-lock is applied at an offset (cnst30) calculated from the nominal
    theta and omega1. Residues at different chemical shifts experience different
    effective offsets and thus different tilt angles.

    Args:
        residue_ppm: Chemical shift of the residue (15N ppm)
        carrier_ppm: Carrier frequency position (15N ppm)
        omega1_hz: Spin-lock field strength in Hz (cnst27)
        theta_nominal_deg: Nominal tilt angle in degrees (cnst28)
        spec_freq_mhz: 1H spectrometer frequency in MHz (e.g., 600, 700, 800)

    Returns:
        Effective tilt angle in degrees for this residue
    """
    # Calculate 15N frequency in MHz
    n15_freq_mhz = spec_freq_mhz / N15_GAMMA_RATIO

    # Calculate the spin-lock offset (cnst30) from nominal theta and omega1
    theta_nominal_rad = np.radians(theta_nominal_deg)

    # Handle theta = 90 (on-resonance) case
    if np.isclose(theta_nominal_deg, 90.0, atol=0.1):
        # On-resonance: cnst30 = 0, all residues have their natural offset
        cnst30_hz = 0.0
    else:
        cnst30_hz = omega1_hz / np.tan(theta_nominal_rad)

    # Calculate residue offset from carrier in Hz
    residue_offset_hz = (residue_ppm - carrier_ppm) * n15_freq_mhz

    # Effective offset during spin-lock
    # The spin-lock carrier is moved to cnst30 from the original carrier
    delta_omega_hz = cnst30_hz - residue_offset_hz

    # Calculate effective theta
    # theta = arctan(omega1 / |delta_omega|)
    if np.isclose(delta_omega_hz, 0.0, atol=0.1):
        # Residue is at the spin-lock carrier position
        theta_effective_deg = 90.0
    else:
        theta_effective_rad = np.arctan2(omega1_hz, np.abs(delta_omega_hz))
        theta_effective_deg = np.degrees(theta_effective_rad)

    return theta_effective_deg


def calculate_R2_from_R1_R1rho(
    R1: float,
    R1rho: float,
    theta_deg: float
) -> float:
    """
    Calculate R2 from R1 and R1rho using the tilt angle relationship.

    R1rho = R1 * cos^2(theta) + R2 * sin^2(theta)

    Solving for R2:
    R2 = (R1rho - R1 * cos^2(theta)) / sin^2(theta)

    Args:
        R1: Longitudinal relaxation rate (s^-1)
        R1rho: Rotating-frame relaxation rate (s^-1)
        theta_deg: Tilt angle in degrees

    Returns:
        R2: Transverse relaxation rate (s^-1)
    """
    theta_rad = np.radians(theta_deg)
    sin2_theta = np.sin(theta_rad) ** 2
    cos2_theta = np.cos(theta_rad) ** 2

    # Avoid division by zero for theta near 0
    if sin2_theta < 1e-6:
        raise ValueError(f"Theta {theta_deg}° is too close to 0°, cannot calculate R2")

    R2 = (R1rho - R1 * cos2_theta) / sin2_theta

    return R2


def calculate_T2_from_T1rho(
    r1_data: pd.DataFrame,
    r1rho_data: pd.DataFrame,
    theta_deg: float,
    omega1_hz: float,
    carrier_ppm: float,
    spec_freq_mhz: float,
    peak_list: pd.DataFrame,
    chemical_shift_column: str = 'N_ppm'
) -> pd.DataFrame:
    """
    Calculate T2 from T1 and T1rho measurements with per-residue theta correction.

    Args:
        r1_data: DataFrame with columns ['residue', 'R1', 'R1_err']
        r1rho_data: DataFrame with columns ['residue', 'R1rho', 'R1rho_err']
        theta_deg: Nominal tilt angle in degrees (cnst28)
        omega1_hz: Spin-lock field strength in Hz (cnst27)
        carrier_ppm: 15N carrier position in ppm
        spec_freq_mhz: 1H spectrometer frequency in MHz
        peak_list: DataFrame with residue chemical shifts
        chemical_shift_column: Column name for 15N chemical shift in peak_list

    Returns:
        DataFrame with columns ['residue', 'R1', 'R1rho', 'theta', 'R2', 'T2', 'T2_err']
    """
    # Find common residues in both datasets
    common_residues = set(r1_data['residue']) & set(r1rho_data['residue'])

    if len(common_residues) == 0:
        raise ValueError("No common residues found between R1 and R1rho datasets")

    # Build result dataframe
    results = []
    unmatched = []

    for residue in sorted(common_residues):
        # Get R1 and R1rho for this residue
        r1_row = r1_data[r1_data['residue'] == residue].iloc[0]
        r1rho_row = r1rho_data[r1rho_data['residue'] == residue].iloc[0]

        R1 = r1_row['R1']
        R1_err = r1_row.get('R1_err', 0.0)
        R1rho = r1rho_row['R1rho']
        R1rho_err = r1rho_row.get('R1rho_err', 0.0)

        # Get chemical shift for this residue
        peak_row = peak_list[peak_list['residue'] == residue]
        if len(peak_row) == 0:
            # Try to find by residue number only
            residue_num = residue.split('.')[0] if '.' in residue else residue
            peak_row = peak_list[peak_list['residue'].str.startswith(residue_num + '.')]

        if len(peak_row) == 0:
            # No shift for this residue, so the offset from the carrier is unknown.
            # Falling back to the carrier makes the offset zero, which collapses theta
            # to the nominal angle and quietly turns the per-residue correction into no
            # correction at all. Still done -- dropping the residue loses more -- but
            # counted, so the caller can see how much of the table is uncorrected.
            residue_ppm = carrier_ppm
            unmatched.append(str(residue))
        else:
            residue_ppm = peak_row[chemical_shift_column].iloc[0]

        # Calculate per-residue theta
        theta_effective = calculate_per_residue_theta(
            residue_ppm=residue_ppm,
            carrier_ppm=carrier_ppm,
            omega1_hz=omega1_hz,
            theta_nominal_deg=theta_deg,
            spec_freq_mhz=spec_freq_mhz
        )

        # Calculate R2
        try:
            R2 = calculate_R2_from_R1_R1rho(R1, R1rho, theta_effective)
        except ValueError:
            # Skip this residue if calculation fails
            continue

        # Calculate T2 in ms
        if R2 > 0:
            T2 = 1000.0 / R2  # Convert to ms
        else:
            # Negative R2 is unphysical, skip
            continue

        # Error propagation (simplified - assumes independent errors)
        # dR2/dR1rho = 1/sin^2(theta)
        # dR2/dR1 = -cos^2(theta)/sin^2(theta)
        theta_rad = np.radians(theta_effective)
        sin2_theta = np.sin(theta_rad) ** 2
        cos2_theta = np.cos(theta_rad) ** 2

        dR2_dR1rho = 1.0 / sin2_theta
        dR2_dR1 = -cos2_theta / sin2_theta

        R2_err = np.sqrt((dR2_dR1rho * R1rho_err) ** 2 + (dR2_dR1 * R1_err) ** 2)

        # T2 error: dT2/dR2 = -1000/R2^2
        T2_err = 1000.0 * R2_err / (R2 ** 2)

        results.append({
            'residue': residue,
            'R1': R1,
            'R1_err': R1_err,
            'R1rho': R1rho,
            'R1rho_err': R1rho_err,
            'theta': theta_effective,
            'R2': R2,
            'R2_err': R2_err,
            'T2': T2,
            'T2_err': T2_err
        })

    if unmatched and len(unmatched) == len(results):
        raise ValueError(
            f"No residue in the R1/R1rho data was found in the peak list, so every "
            f"tilt angle fell back to the nominal {theta_deg} deg. The result would be "
            f"the uncorrected R2 wearing a per-residue label. Check that the peak "
            f"list's assignments match the relaxation matrices'."
        )
    table = pd.DataFrame(results)
    # attrs, not columns: this is about the table, not about any one residue.
    table.attrs['n_shift_unmatched'] = len(unmatched)
    table.attrs['shift_unmatched'] = unmatched
    return table


def r2_table_from_fits(t1_fits, t1rho_fits, peak_list, *, omega1_hz, carrier_ppm,
                       theta_deg, spec_freq_mhz, time_units='ms'):
    """Build a per-residue R2 table from fitted T1 and T1rho decays.

    The fitters return {residue: {'value': T, 'error': T_err}}; the tilt-angle
    correction wants rates and a chemical shift per residue. This is the join between
    them, so callers do not re-derive the algebra.

    Returns a DataFrame with residue, R1, R1_err, R1rho, R1rho_err, theta, R2, R2_err,
    T2 and T2_err — the shape the spectral-density step consumes.
    """
    scale = {'ms': 1000.0, 's': 1.0, 'us': 1e6}[time_units]

    def _rates(fits, value_key, err_key):
        rows = []
        for residue, rec in fits.items():
            t = float(rec['value'])
            if not np.isfinite(t) or t <= 0:
                continue
            err = float(rec.get('error', 0.0) or 0.0)
            rows.append({'residue': str(residue), value_key: scale / t,
                         err_key: scale * err / (t ** 2)})
        return pd.DataFrame(rows)

    r1 = _rates(t1_fits, 'R1', 'R1_err')
    r1rho = _rates(t1rho_fits, 'R1rho', 'R1rho_err')
    if r1.empty or r1rho.empty:
        raise ValueError('T1 and T1rho fits must both contain usable residues')

    peaks = peak_list.copy()
    peaks['residue'] = peaks['residue'].astype(str)
    return calculate_T2_from_T1rho(
        r1_data=r1, r1rho_data=r1rho, theta_deg=theta_deg, omega1_hz=omega1_hz,
        carrier_ppm=carrier_ppm, spec_freq_mhz=spec_freq_mhz, peak_list=peaks,
        chemical_shift_column='N_ppm')
