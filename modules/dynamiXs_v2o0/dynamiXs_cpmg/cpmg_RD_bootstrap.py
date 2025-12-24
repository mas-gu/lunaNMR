import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

# Constants
T_CPMG = 0.015  # Total CPMG time per block (d21, in seconds)
P30 = 120e-6    # 180° pulse length (in seconds)
FIELD = 700     # MHz (for 15N Larmor = field * 0.1013)
N_BOOTSTRAP = 10  # Number of bootstrap iterations for fit parameters
NOISE = 45475  # User-specified noise value (standard deviation) for all delays

def read_peak_list(filename):
    df = pd.read_csv(filename, header=None)
    vclist = df.iloc[0, 1:].astype(int).to_numpy()
    residue_names = df.iloc[1:, 0].to_numpy()
    intensities = df.iloc[1:, 1:].astype(float).to_numpy()
    return residue_names, vclist, intensities

def calculate_nu_cpmg(vclist, T_CPMG=0.015, P30=120e-6):
    nu_cpmg = []
    for l0 in vclist:
        if l0 == 0:
            nu_cpmg.append(0)
        else:
            tau_cp = (T_CPMG / l0 - P30) / 2
            nu_cpmg.append(1 / (2 * tau_cp) if tau_cp > 0 else 1000)
    return np.array(nu_cpmg)

def cpmg_model(nu_cpmg, R2_0, k_ex, p_a, delta_omega_ppm):
    delta_omega_hz = delta_omega_ppm * FIELD * 0.1013
    p_b = 1 - p_a
    phi = p_a * p_b * np.square(delta_omega_hz * 2 * np.pi)
    R_ex = phi / k_ex / (1 + np.square(2 * np.pi * nu_cpmg / k_ex))
    return R2_0 + R_ex

def calculate_r2_eff(intensities, I_ref, T_total=2 * T_CPMG):
    r2_eff = np.full_like(intensities, np.nan, dtype=float)
    n_residues, n_points = intensities.shape
    for i in range(n_residues):
        for j in range(n_points):
            if intensities[i, j] > 0 and I_ref[i] > 0:
                r2_eff[i, j] = -np.log(intensities[i, j] / I_ref[i]) / T_total
    return r2_eff

def calculate_r2_eff_error(intensities, I_ref, noise, T_total=2 * T_CPMG):
    r2_eff_error = np.full_like(intensities, np.nan, dtype=float)
    n_residues, n_points = intensities.shape
    for i in range(n_residues):
        for j in range(n_points):
            if intensities[i, j] > 0 and I_ref[i] > 0:
                # Error propagation: σ(R2,eff) = (1/T_total) * sqrt((σ_I/I)^2 + (σ_Iref/I_ref)^2)
                rel_error_I = noise / intensities[i, j]
                rel_error_Iref = noise / I_ref[i]
                r2_eff_error[i, j] = (1 / T_total) * np.sqrt(rel_error_I**2 + rel_error_Iref**2)
    return r2_eff_error

def fit_dispersion(r2_eff, nu_cpmg):
    valid = (~np.isnan(r2_eff)) & (nu_cpmg != 0)  # Exclude nu_cpmg = 0
    if np.sum(valid) < 5:
        return None, None
    initial_guess = [10, 1000, 0.9, 1.0]
    try:
        popt, pcov = curve_fit(
            cpmg_model,
            nu_cpmg[valid],
            r2_eff[valid],
            p0=initial_guess,
            bounds=([0, 100, 0.5, 0], [50, 10000, 0.999, 5])
        )
        return popt, np.sqrt(np.diag(pcov))
    except RuntimeError:
        return None, None

def bootstrap_fit(intensities, I_ref, noise, nu_cpmg, residue_names):
    n_residues, n_points = intensities.shape
    bootstrap_params = np.zeros((N_BOOTSTRAP, n_residues, 4))  # Store R2_0, k_ex, p_a, delta_omega_ppm
    for i in range(n_residues):
        print(f"Starting bootstrap resampling for residue {residue_names[i]}")
        valid_indices = np.where((intensities[i] > 0) & (~np.isnan(intensities[i])))[0]
        if len(valid_indices) == 0:
            print(f"  Residue {residue_names[i]}: No valid intensities for resampling")
            continue
        for n in range(N_BOOTSTRAP):
            if n % 10 == 0:
                print(f"  Residue {residue_names[i]}: Iteration {n}/{N_BOOTSTRAP}")
            # Resample intensities with Gaussian noise
            resampled_intensities = np.copy(intensities[i])
            resampled_I_ref = I_ref[i]
            for j in valid_indices:
                resampled_intensities[j] += np.random.normal(0, noise)
            if I_ref[i] > 0:
                resampled_I_ref += np.random.normal(0, noise)
            # Ensure intensities remain positive
            resampled_intensities = np.maximum(resampled_intensities, 1e-10)
            resampled_I_ref = max(resampled_I_ref, 1e-10)
            # Calculate R2,eff for resampled data
            r2_eff_resampled = calculate_r2_eff(resampled_intensities[np.newaxis, :], np.array([resampled_I_ref]))[0]
            # Fit the resampled R2,eff
            popt, _ = fit_dispersion(r2_eff_resampled, nu_cpmg)
            if popt is not None:
                bootstrap_params[n, i] = popt
            else:
                bootstrap_params[n, i] = [np.nan, np.nan, np.nan, np.nan]
        print(f"Completed bootstrap resampling for residue {residue_names[i]}")
    # Calculate mean and standard deviation of parameters
    param_mean = np.nanmean(bootstrap_params, axis=0)
    param_std = np.nanstd(bootstrap_params, axis=0, ddof=1)
    return param_mean, param_std

def main():
    filename = "A1_WT_CPMG_15ms.csv"
    noise = NOISE  # User-specified noise value
    residue_names, vclist, intensities = read_peak_list(filename)
    vclist = vclist[:intensities.shape[1]]
    nu_cpmg = calculate_nu_cpmg(vclist)

    # Find reference column where vclist == 0
    try:
        ref_index = np.where(vclist == 0)[0][0]
    except IndexError:
        raise ValueError("No reference point (vclist == 0) found. Ensure your vclist contains a 0 entry.")

    I_ref = intensities[:, ref_index]
    r2_eff = calculate_r2_eff(intensities, I_ref)
    r2_eff_error = calculate_r2_eff_error(intensities, I_ref, noise)

    # Bootstrap for fit parameters
    param_mean, param_std = bootstrap_fit(intensities, I_ref, noise, nu_cpmg, residue_names)

    with open("cpmg_results.txt", "w") as f:
        f.write("Residue R2_0 k_ex p_a delta_omega_ppm R2_0_err k_ex_err p_a_err delta_omega_err\n")

        n_residues = len(residue_names)
        n_plots_per_fig = 20
        n_figs = int(np.ceil(n_residues / n_plots_per_fig))

        for fig_idx in range(n_figs):
            fig, axes = plt.subplots(5, 4, figsize=(20, 25), sharex=True)
            axes = axes.flatten()
            start_idx = fig_idx * n_plots_per_fig
            end_idx = min((fig_idx + 1) * n_plots_per_fig, n_residues)

            for ax_idx, i in enumerate(range(start_idx, end_idx)):
                residue = residue_names[i]
                ax = axes[ax_idx]
                # Plot points with error bars where nu_cpmg != 0
                mask = nu_cpmg != 0
                ax.errorbar(nu_cpmg[mask], r2_eff[i][mask], yerr=r2_eff_error[i][mask],
                            fmt='o', color='blue', label='Data', capsize=3)

                if not np.any(np.isnan(param_mean[i])):
                    f.write(f"{residue} {param_mean[i][0]:.2f} {param_mean[i][1]:.2f} {param_mean[i][2]:.4f} {param_mean[i][3]:.4f} "
                            f"{param_std[i][0]:.2f} {param_std[i][1]:.2f} {param_std[i][2]:.4f} {param_std[i][3]:.4f}\n")
                    nu_cpmg_fine = np.linspace(0, max(nu_cpmg)*1.1, 100)
                    r2_fit = cpmg_model(nu_cpmg_fine, *param_mean[i])
                    ax.plot(nu_cpmg_fine, r2_fit, 'r-', label='Fit')
                    textstr = f"k_ex = {param_mean[i][1]:.1f} ({param_std[i][1]:.1f}) s⁻¹\nΔω = {param_mean[i][3]:.2f} ({param_std[i][3]:.2f}) ppm"
                    ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
                            fontsize=10, verticalalignment='top', horizontalalignment='right',
                            bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
                else:
                    f.write(f"{residue} Fit failed\n")

                ax.set_title(f"Residue {residue}")
                ax.set_xlabel("ν_CPMG (Hz)")
                ax.set_ylabel("R2,eff (s⁻¹)")
                ax.grid(True)
                ax.legend()

            # Hide unused axes
            for j in range(ax_idx + 1, len(axes)):
                fig.delaxes(axes[j])

            plt.tight_layout()
            fig.savefig(f"testCPMG_dispersion_page_{fig_idx + 1}.pdf", format='pdf')
            plt.close(fig)

if __name__ == "__main__":
    main()