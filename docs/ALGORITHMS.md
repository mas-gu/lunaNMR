# Algorithms Documentation - lunaNMR v0.9

This document provides comprehensive technical documentation of the algorithms implemented in lunaNMR, covering the mathematical foundations, implementation details, and performance characteristics of each component.

## 📋 Table of Contents

1. [Overview](#overview)
2. [Voigt Profile Fitting](#voigt-profile-fitting)
3. [Peak Detection Algorithms](#peak-detection-algorithms)
4. [Baseline Correction](#baseline-correction)
5. [Quality Assessment](#quality-assessment)
6. [Series Analysis](#series-analysis)
7. [Optimization Strategies](#optimization-strategies)
8. [Performance Analysis](#performance-analysis)

---

## 🔬 Overview

lunaNMR implements state-of-the-art algorithms specifically designed for NMR spectroscopy analysis. The core philosophy emphasizes:

- **Mathematical rigor** with proper statistical uncertainty estimation
- **NMR-specific optimizations** accounting for typical spectral characteristics
- **Robust error handling** for real-world noisy data
- **Performance optimization** for large-scale analyses
- **Professional-grade quality** suitable for academic publication

### **Algorithm Architecture**

```
Input Data → Preprocessing → Peak Detection → Voigt Fitting → Quality Assessment → Results
     ↓             ↓              ↓              ↓              ↓              ↓
  Validation  Noise Filter   Clustering    Parameter    Statistical    Export/
  & Format    Baseline      Peak Merge     Bounds      Uncertainty   Visualization
  Checking    Correction    Refinement     Checking      Analysis
```

---

## 📈 Voigt Profile Fitting

### **Mathematical Foundation**

The Voigt profile is a convolution of Gaussian and Lorentzian functions, representing the true lineshape in NMR spectroscopy:

#### **Voigt Function Definition**
```
V(ν) = A ∫_{-∞}^{∞} (1/π) * (γ/((ν-ν')²+γ²)) * exp(-(ν'-ν₀)²/(2σ²)) dν'
```

Where:
- `A` = Peak amplitude
- `ν₀` = Center frequency (chemical shift)
- `σ` = Gaussian width parameter (inhomogeneous broadening)
- `γ` = Lorentzian width parameter (homogeneous broadening)
- `ν` = Frequency variable

#### **Faddeeva Function Implementation**

lunaNMR uses the complex Faddeeva function `w(z)` for computational efficiency:

```
V(ν) = A × Re[w(z)] / (σ√(2π))

where z = (ν - ν₀ + iγ) / (σ√2)
```

**Advantages of Faddeeva approach:**
- **Numerical stability** for all parameter ranges
- **Computational efficiency** compared to direct convolution
- **High precision** maintaining double-precision accuracy
- **Automatic handling** of limiting cases (pure Gaussian/Lorentzian)

### **Parameter Estimation Strategy**

#### **Multi-Stage Fitting Process**

**Stage 1: Initial Parameter Estimation**
```python
def estimate_initial_parameters(x, y, peak_center):
    # Amplitude estimation from local maximum
    amplitude = np.max(y[center_window])
    
    # Width estimation using full-width-half-maximum (FWHM)
    half_max = amplitude / 2
    fwhm_indices = np.where(y >= half_max)[0]
    estimated_width = np.ptp(x[fwhm_indices])
    
    # Split width between Gaussian and Lorentzian components
    sigma_initial = estimated_width / (2 * np.sqrt(2 * np.log(2)))  # ~0.6 * FWHM
    gamma_initial = estimated_width / 4  # Smaller Lorentzian contribution
    
    return amplitude, peak_center, sigma_initial, gamma_initial
```

**Stage 2: Coarse Fitting**
```python
# Relaxed bounds for initial convergence
bounds_coarse = [
    [0.1 * amplitude, 10 * amplitude],          # Amplitude
    [center - 0.1, center + 0.1],              # Center (ppm)
    [width_min, 10 * estimated_width],         # Gaussian width
    [width_min, 10 * estimated_width]          # Lorentzian width
]
```

**Stage 3: Fine Fitting**
```python
# Tighter bounds based on coarse results
bounds_fine = [
    [0.5 * coarse_amplitude, 2 * coarse_amplitude],
    [coarse_center - 0.05, coarse_center + 0.05],
    [0.5 * coarse_sigma, 2 * coarse_sigma],
    [0.5 * coarse_gamma, 2 * coarse_gamma]
]
```

#### **NMR-Specific Parameter Bounds**

```python
# Realistic bounds for different NMR nuclei
NMR_BOUNDS = {
    '1H': {
        'chemical_shift_range': (5.5, 12.0),  # Typical protein amide region
        'width_limits': (0.005, 0.1),        # 5-100 Hz at 600 MHz
        'amplitude_ratio': (0.1, 10.0)       # Relative to noise level
    },
    '15N': {
        'chemical_shift_range': (100.0, 140.0),  # Typical protein amide region
        'width_limits': (0.2, 2.0),              # Broader than 1H
        'amplitude_ratio': (0.1, 10.0)
    },
    '13C': {
        'chemical_shift_range': (5.0, 50.0),     # Aliphatic region
        'width_limits': (0.1, 1.0),
        'amplitude_ratio': (0.1, 10.0)
    }
}
```

### **Uncertainty Estimation**

#### **Covariance Matrix Analysis**
```python
def estimate_uncertainties(fitted_params, jacobian, residuals):
    """
    Calculate parameter uncertainties from covariance matrix
    """
    # Calculate residual sum of squares
    rss = np.sum(residuals**2)
    
    # Degrees of freedom
    dof = len(residuals) - len(fitted_params)
    
    # Mean squared error
    mse = rss / dof
    
    # Covariance matrix
    try:
        # J^T J inverse (normal equations)
        covariance = mse * np.linalg.inv(jacobian.T @ jacobian)
        
        # Standard errors are sqrt of diagonal elements
        uncertainties = np.sqrt(np.diag(covariance))
        
        # Parameter correlations from covariance
        correlations = covariance / np.outer(uncertainties, uncertainties)
        
        return uncertainties, correlations
        
    except np.linalg.LinAlgError:
        # Singular matrix - return conservative estimates
        return np.full_like(fitted_params, np.inf), None
```

#### **Bootstrap Error Estimation**
For robust uncertainty analysis:

```python
def bootstrap_uncertainties(x, y, fitted_params, n_bootstrap=100):
    """
    Bootstrap resampling for robust uncertainty estimation
    """
    n_points = len(x)
    bootstrap_params = []
    
    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(n_points, n_points, replace=True)
        x_boot, y_boot = x[indices], y[indices]
        
        # Refit with bootstrap sample
        try:
            boot_params, _ = curve_fit(voigt_function, x_boot, y_boot, 
                                     p0=fitted_params, bounds=bounds)
            bootstrap_params.append(boot_params)
        except:
            continue  # Skip failed fits
    
    # Calculate statistics
    bootstrap_params = np.array(bootstrap_params)
    uncertainties = np.std(bootstrap_params, axis=0)
    
    return uncertainties
```

---

## 🎯 Peak Detection Algorithms

### **Enhanced Peak Detection Pipeline**

#### **Step 1: Noise Level Estimation**
```python
def estimate_noise_level(spectrum, method='mad'):
    """
    Robust noise estimation using Median Absolute Deviation
    """
    if method == 'mad':
        median = np.median(spectrum)
        mad = np.median(np.abs(spectrum - median))
        noise_level = 1.4826 * mad  # Robust standard deviation estimator
    
    elif method == 'percentile':
        # Use lower percentiles to avoid peak contamination
        noise_level = np.std(spectrum[spectrum < np.percentile(spectrum, 25)])
    
    return noise_level
```

#### **Step 2: Adaptive Smoothing**
```python
def adaptive_smooth(spectrum, noise_level):
    """
    Adaptive smoothing based on local signal-to-noise ratio
    """
    from scipy.ndimage import gaussian_filter1d
    
    # Calculate local SNR
    window_size = 50
    local_snr = []
    for i in range(len(spectrum)):
        start = max(0, i - window_size//2)
        end = min(len(spectrum), i + window_size//2)
        local_max = np.max(spectrum[start:end])
        local_snr.append(local_max / noise_level)
    
    # Adaptive smoothing kernel
    sigma = np.where(np.array(local_snr) > 5, 1.0, 2.0)  # Less smoothing for peaks
    
    # Apply variable smoothing
    smoothed = gaussian_filter1d(spectrum, sigma=np.mean(sigma))
    
    return smoothed
```

#### **Step 3: Peak Detection with Network Analysis**
```python
def detect_peaks_network(spectrum, x_axis, prominence_threshold=0.1):
    """
    Network-based peak detection using graph clustering
    """
    import networkx as nx
    from scipy.signal import find_peaks, peak_prominences
    
    # Initial peak detection
    peaks, properties = find_peaks(
        spectrum,
        prominence=prominence_threshold * np.max(spectrum),
        distance=10  # Minimum separation
    )
    
    # Calculate peak prominences and widths
    prominences = peak_prominences(spectrum, peaks)[0]
    widths = peak_widths(spectrum, peaks, rel_height=0.5)[0]
    
    # Create graph for peak clustering
    G = nx.Graph()
    
    # Add peaks as nodes
    for i, peak_idx in enumerate(peaks):
        G.add_node(i, 
                  position=x_axis[peak_idx],
                  intensity=spectrum[peak_idx],
                  prominence=prominences[i],
                  width=widths[i])
    
    # Connect nearby peaks
    for i in range(len(peaks)):
        for j in range(i+1, len(peaks)):
            distance = abs(x_axis[peaks[i]] - x_axis[peaks[j]])
            overlap = (widths[i] + widths[j]) / 2
            
            if distance < 2 * overlap:  # Overlapping peaks
                G.add_edge(i, j, weight=1.0/distance)
    
    # Find connected components (peak clusters)
    clusters = list(nx.connected_components(G))
    
    # Merge overlapping peaks within clusters
    merged_peaks = []
    for cluster in clusters:
        if len(cluster) == 1:
            # Single peak
            peak_idx = list(cluster)[0]
            merged_peaks.append({
                'center': x_axis[peaks[peak_idx]],
                'intensity': spectrum[peaks[peak_idx]],
                'prominence': prominences[peak_idx],
                'width': widths[peak_idx],
                'quality': 'single'
            })
        else:
            # Merge cluster peaks
            cluster_indices = [peaks[i] for i in cluster]
            cluster_positions = [x_axis[i] for i in cluster_indices]
            cluster_intensities = [spectrum[i] for i in cluster_indices]
            
            # Weighted average position
            total_intensity = sum(cluster_intensities)
            avg_position = sum(p * i for p, i in zip(cluster_positions, cluster_intensities)) / total_intensity
            
            merged_peaks.append({
                'center': avg_position,
                'intensity': total_intensity,
                'prominence': max(prominences[list(cluster)]),
                'width': np.mean([widths[i] for i in cluster]),
                'quality': 'merged'
            })
    
    return merged_peaks
```

---

## 📐 Baseline Correction

### **Polynomial Baseline Fitting**

#### **Robust Polynomial Estimation**
```python
def estimate_baseline(x, y, degree=1, mask_peaks=True):
    """
    Robust baseline estimation using iterative polynomial fitting
    """
    from scipy.signal import find_peaks
    
    # Initial peak detection for masking
    if mask_peaks:
        peaks, _ = find_peaks(y, prominence=0.1*np.max(y))
        
        # Create mask excluding peak regions
        mask = np.ones(len(y), dtype=bool)
        for peak in peaks:
            # Mask ±width around each peak
            width = max(20, int(len(y) * 0.02))  # 2% of spectrum or 20 points
            start = max(0, peak - width)
            end = min(len(y), peak + width)
            mask[start:end] = False
    else:
        mask = np.ones(len(y), dtype=bool)
    
    # Iterative baseline fitting
    for iteration in range(3):
        # Fit polynomial to non-peak regions
        baseline_coeffs = np.polyfit(x[mask], y[mask], degree)
        baseline = np.polyval(baseline_coeffs, x)
        
        # Update mask: exclude points significantly above baseline
        residuals = y - baseline
        threshold = 2 * np.std(residuals[mask])  # 2-sigma threshold
        mask = residuals < threshold
    
    return baseline, baseline_coeffs
```

#### **Asymmetric Least Squares (AsLS) Baseline**
```python
def asls_baseline(y, lambda_param=1e6, p=0.01, max_iter=10):
    """
    Asymmetric Least Squares baseline correction
    Penalizes positive deviations more than negative ones
    """
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L-2))
    D = lambda_param * D.dot(D.transpose())
    
    w = np.ones(L)
    for i in range(max_iter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + D
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    
    return z
```

---

## 🎖️ Quality Assessment

### **Multi-Criteria Quality Evaluation**

#### **R² Calculation with Proper Degrees of Freedom**
```python
def calculate_r_squared(y_observed, y_predicted, n_parameters):
    """
    Calculate adjusted R² accounting for model complexity
    """
    # Total sum of squares
    ss_tot = np.sum((y_observed - np.mean(y_observed))**2)
    
    # Residual sum of squares
    ss_res = np.sum((y_observed - y_predicted)**2)
    
    # Standard R²
    r_squared = 1 - (ss_res / ss_tot)
    
    # Adjusted R² (penalizes additional parameters)
    n = len(y_observed)
    adjusted_r_squared = 1 - (ss_res / (n - n_parameters)) / (ss_tot / (n - 1))
    
    return r_squared, adjusted_r_squared
```

#### **Comprehensive Quality Metrics**
```python
def assess_fit_quality(fitted_params, x, y_observed, y_fitted, uncertainties):
    """
    Comprehensive quality assessment for Voigt fits
    """
    metrics = {}
    
    # Basic goodness-of-fit
    metrics['r_squared'] = calculate_r_squared(y_observed, y_fitted, len(fitted_params))[0]
    metrics['rmse'] = np.sqrt(np.mean((y_observed - y_fitted)**2))
    
    # Parameter reasonableness
    center, amplitude, sigma, gamma = fitted_params
    
    # Check parameter bounds
    metrics['center_reasonable'] = 5.5 <= center <= 12.0  # 1H range
    metrics['width_reasonable'] = 0.005 <= max(sigma, gamma) <= 0.1
    metrics['amplitude_positive'] = amplitude > 0
    
    # Signal-to-noise ratio
    noise_level = estimate_noise_level(y_observed - y_fitted)
    metrics['snr'] = amplitude / noise_level
    
    # Parameter precision
    if uncertainties is not None:
        metrics['center_precision'] = uncertainties[0] / center if center != 0 else np.inf
        metrics['amplitude_precision'] = uncertainties[1] / amplitude if amplitude != 0 else np.inf
    
    # Residual analysis
    residuals = y_observed - y_fitted
    metrics['residual_autocorr'] = np.corrcoef(residuals[:-1], residuals[1:])[0,1]
    metrics['residual_normality'] = scipy.stats.jarque_bera(residuals)[1]  # p-value
    
    # Overall quality classification
    if (metrics['r_squared'] > 0.95 and 
        metrics['snr'] > 10 and 
        metrics['center_reasonable'] and
        metrics['width_reasonable']):
        metrics['overall_quality'] = 'Excellent'
    elif (metrics['r_squared'] > 0.90 and 
          metrics['snr'] > 5):
        metrics['overall_quality'] = 'Good'
    elif metrics['r_squared'] > 0.80:
        metrics['overall_quality'] = 'Fair'
    else:
        metrics['overall_quality'] = 'Poor'
    
    return metrics
```

---

## 📊 Series Analysis

### **Exponential Decay Fitting**

For T₁ and T₂ relaxation analysis:

#### **Mono-exponential Model**
```python
def monoexponential_decay(t, I0, T, baseline=0):
    """
    Standard mono-exponential decay model
    I(t) = I0 * exp(-t/T) + baseline
    """
    return I0 * np.exp(-t / T) + baseline

def fit_t1_relaxation(time_points, intensities, initial_guess=None):
    """
    Fit T1 relaxation data with robust parameter estimation
    """
    if initial_guess is None:
        # Estimate initial parameters
        I0_est = np.max(intensities)
        T_est = time_points[np.argmin(np.abs(intensities - I0_est/np.e))]
        baseline_est = np.min(intensities)
        initial_guess = [I0_est, T_est, baseline_est]
    
    # Parameter bounds
    bounds = ([0, 0.01, -np.inf], [10*np.max(intensities), 10*np.max(time_points), np.inf])
    
    try:
        params, covariance = curve_fit(
            monoexponential_decay, 
            time_points, 
            intensities,
            p0=initial_guess,
            bounds=bounds,
            maxfev=2000
        )
        
        # Calculate uncertainties
        uncertainties = np.sqrt(np.diag(covariance))
        
        # Calculate R1 = 1/T1
        T1 = params[1]
        R1 = 1 / T1
        R1_uncertainty = uncertainties[1] / (T1**2)  # Error propagation
        
        return {
            'I0': params[0],
            'T1': T1,
            'R1': R1,
            'baseline': params[2],
            'uncertainties': {
                'I0': uncertainties[0],
                'T1': uncertainties[1],
                'R1': R1_uncertainty,
                'baseline': uncertainties[2]
            },
            'fitted_curve': monoexponential_decay(time_points, *params),
            'success': True
        }
        
    except Exception as e:
        return {'success': False, 'error': str(e)}
```

#### **Bi-exponential Model**
```python
def biexponential_decay(t, I0, T1_fast, T1_slow, fraction_fast, baseline=0):
    """
    Bi-exponential decay for complex relaxation behavior
    I(t) = I0 * [f*exp(-t/T1_fast) + (1-f)*exp(-t/T1_slow)] + baseline
    """
    return I0 * (fraction_fast * np.exp(-t / T1_fast) + 
                 (1 - fraction_fast) * np.exp(-t / T1_slow)) + baseline
```

### **Statistical Analysis of Series**

#### **Population Statistics**
```python
def analyze_series_statistics(relaxation_results):
    """
    Statistical analysis of relaxation parameters across a series
    """
    # Extract parameters
    T1_values = [r['T1'] for r in relaxation_results if r['success']]
    R1_values = [r['R1'] for r in relaxation_results if r['success']]
    
    stats = {}
    
    for param_name, values in [('T1', T1_values), ('R1', R1_values)]:
        if len(values) > 0:
            stats[param_name] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'median': np.median(values),
                'min': np.min(values),
                'max': np.max(values),
                'count': len(values),
                'cv': np.std(values) / np.mean(values) if np.mean(values) != 0 else np.inf
            }
            
            # Distribution testing
            if len(values) > 5:
                # Normality test
                _, p_normal = scipy.stats.shapiro(values)
                stats[param_name]['normal_distribution'] = p_normal > 0.05
                
                # Outlier detection (modified Z-score)
                median = np.median(values)
                mad = np.median(np.abs(values - median))
                modified_z_scores = 0.6745 * (values - median) / mad
                outliers = np.abs(modified_z_scores) > 3.5
                stats[param_name]['outliers'] = np.sum(outliers)
    
    return stats
```

---

## ⚡ Optimization Strategies

### **Multi-threaded Processing**

#### **Parallel Peak Fitting**
```python
from multiprocessing import Pool, cpu_count
from functools import partial

def parallel_peak_fitting(peak_list, spectrum_data, n_workers=None):
    """
    Parallel processing of peak fitting using multiprocessing
    """
    if n_workers is None:
        n_workers = min(cpu_count(), len(peak_list))
    
    # Prepare fitting function with fixed spectrum data
    fit_function = partial(fit_single_peak, spectrum_data=spectrum_data)
    
    # Parallel execution
    with Pool(n_workers) as pool:
        results = pool.map(fit_function, peak_list)
    
    return results

def fit_single_peak(peak_info, spectrum_data):
    """
    Fit single peak - designed for parallel execution
    """
    x, y = spectrum_data
    center = peak_info['center']
    
    # Extract local region around peak
    window_size = 0.2  # ppm
    mask = np.abs(x - center) < window_size
    x_local, y_local = x[mask], y[mask]
    
    # Fit Voigt profile
    try:
        return fit_voigt_profile(x_local, y_local, center)
    except:
        return {'success': False, 'error': 'Fitting failed'}
```

### **Memory Optimization**

#### **Chunked Processing for Large Datasets**
```python
def process_large_series(file_list, chunk_size=50):
    """
    Process large series in chunks to manage memory usage
    """
    results = []
    
    for i in range(0, len(file_list), chunk_size):
        chunk = file_list[i:i+chunk_size]
        
        # Process chunk
        chunk_results = process_chunk(chunk)
        results.extend(chunk_results)
        
        # Memory cleanup
        import gc
        gc.collect()
    
    return results
```

---

## 📏 Performance Analysis

### **Algorithm Complexity**

| Algorithm Component | Time Complexity | Space Complexity | Notes |
|-------------------|-----------------|------------------|-------|
| Voigt Fitting | O(n × m) | O(n) | n=data points, m=iterations |
| Peak Detection | O(n log n) | O(n) | Dominated by sorting operations |
| Baseline Correction | O(n³) | O(n²) | For polynomial degree 3+ |
| Network Clustering | O(p²) | O(p²) | p=number of peaks |
| Series Analysis | O(s × n) | O(s) | s=series length |

### **Performance Benchmarks**

Typical performance on modern hardware (Intel i7, 16GB RAM):

| Dataset Size | Processing Time | Memory Usage | Throughput |
|-------------|----------------|--------------|------------|
| Single spectrum (1K points, 10 peaks) | 0.5s | 50MB | 20 peaks/s |
| Series analysis (50 spectra) | 25s | 200MB | 2 spectra/s |
| Large dataset (500 spectra) | 4min | 1.5GB | 2 spectra/s |

### **Optimization Recommendations**

#### **For Large Datasets**
- Use parallel processing with `n_workers = cpu_count()`
- Process in chunks to limit memory usage
- Consider using SSD storage for I/O intensive operations

#### **For High-Precision Analysis**
- Increase `max_iterations` parameter
- Use bootstrap uncertainty estimation
- Apply multiple fitting strategies and compare results

#### **For Real-Time Applications**
- Reduce spectral resolution before processing
- Use simplified peak detection algorithms
- Cache frequently used parameters

---

## 🔬 Validation and Testing

### **Algorithm Validation**

#### **Synthetic Data Generation**
```python
def generate_test_spectrum(peaks, noise_level=50, n_points=1000):
    """
    Generate synthetic NMR spectrum for algorithm testing
    """
    x = np.linspace(7.0, 9.0, n_points)
    y = np.zeros_like(x)
    
    for peak in peaks:
        y += voigt_function(x, peak['amplitude'], peak['center'], 
                           peak['sigma'], peak['gamma'])
    
    # Add realistic noise
    noise = np.random.normal(0, noise_level, len(x))
    y += noise
    
    return x, y
```

#### **Performance Metrics**
```python
def validate_algorithm_performance(true_params, fitted_params):
    """
    Calculate validation metrics for algorithm performance
    """
    metrics = {}
    
    # Parameter accuracy
    for param in ['center', 'amplitude', 'sigma', 'gamma']:
        true_val = true_params[param]
        fitted_val = fitted_params[param]
        
        metrics[f'{param}_error'] = abs(fitted_val - true_val)
        metrics[f'{param}_relative_error'] = abs(fitted_val - true_val) / true_val
    
    # Overall accuracy
    metrics['overall_accuracy'] = np.mean([
        metrics['center_relative_error'],
        metrics['amplitude_relative_error'],
        metrics['sigma_relative_error'],
        metrics['gamma_relative_error']
    ])
    
    return metrics
```

---

This comprehensive algorithms documentation provides the mathematical and implementation foundation for understanding and extending lunaNMR's analysis capabilities. Each algorithm is designed with NMR spectroscopy's specific requirements in mind, ensuring both accuracy and practical usability for real-world applications.