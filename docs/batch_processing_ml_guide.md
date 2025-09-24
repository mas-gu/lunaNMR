# Batch Processing and ML Data Collection Guide

## Overview

The LunaNMR Batch Processor provides automated analysis of multiple NMR spectra with integrated Machine Learning (ML) training data collection. This system is designed to process large datasets efficiently while maintaining the same high-quality Voigt fitting as the interactive GUI.

## Key Features

- **Automated Spectrum Processing**: Process multiple NMR files with minimal user intervention
- **Integrated ML Data Collection**: Automatically collect training data during processing
- **Nucleus-Specific Optimization**: Tailored parameters for different nucleus types (¹H, ¹³C-¹H, ¹⁵N-¹H)
- **Configurable Presets**: Multiple processing strategies for different use cases
- **Parallel Processing**: Multi-core support for enhanced performance
- **Quality Control**: Automated quality assessment and filtering
- **Comprehensive Logging**: Detailed progress tracking and error reporting

## Architecture

The batch processing system consists of several key components:

```
BatchProcessor
├── ConfigManager (presets and parameter validation)
├── SingleSpectrumProcessor (individual spectrum processing)
├── Enhanced ML Data Collector (training data generation)
└── Parameter Synchronization (GUI compatibility)
```

## Processing Workflow

1. **File Discovery**: Automatically detect NMR files in specified directories
2. **Nucleus Detection**: Identify spectrum type from filename or characteristics
3. **Parameter Selection**: Apply nucleus-specific or user-defined parameters
4. **Peak Detection**: Use S/N-based detection with adaptive thresholds
5. **Voigt Fitting**: Apply enhanced Voigt profile fitting (same as GUI)
6. **ML Data Collection**: Automatically extract training features
7. **Quality Assessment**: Filter results based on fitting quality
8. **Result Compilation**: Generate comprehensive processing reports

## Parameter System

### Built-in Parameters

The batch processor includes optimized default parameters located in `batch_processor.py`:

```python
config = {
    'sn_thresholds': {
        '15N1H': 2.5,        # 2D HSQC spectra
        '13C1H': 2.0,        # Carbon-proton correlations
        '1H': 2.0,           # 1D proton spectra
        'default': 2.2       # Generic fallback
    },
    'expected_peaks': {
        '15N1H': 100,        # Typical protein backbone
        '13C1H': 80,         # Carbon correlations
        '1H': 60,            # Proton peaks
        'default': 150       # Conservative estimate
    },
    'quality_threshold': 0.6,        # Minimum R² for acceptance
    'max_fitting_attempts': 3,       # Retry limit for failed fits
    'skip_on_error': True,           # Continue on individual failures
    'detailed_logging': True,        # Comprehensive log output
    'progress_interval': 5           # Progress reporting frequency
}
```

### Parameter Descriptions

#### S/N Thresholds
- **Purpose**: Minimum signal-to-noise ratio for peak detection
- **15N1H (2.5)**: Higher threshold for cleaner 2D HSQC spectra
- **13C1H (2.0)**: Moderate threshold for carbon-proton correlations
- **1H (2.0)**: Lower threshold for potentially crowded 1D spectra
- **Range**: 1.0 - 5.0 (typical values)

#### Expected Peaks
- **Purpose**: Expected number of peaks for optimization
- **15N1H (100)**: Typical protein backbone resonances
- **13C1H (80)**: Carbon-proton correlations
- **1H (60)**: 1D proton spectrum complexity
- **Usage**: Guides detection sensitivity and processing time

#### Quality Threshold
- **Purpose**: Minimum R² value for peak acceptance
- **Default (0.6)**: Balanced between quality and quantity
- **Range**: 0.3 (permissive) to 0.9 (strict)
- **Impact**: Higher values = fewer but higher quality peaks

#### Max Fitting Attempts
- **Purpose**: Number of retry attempts for difficult peaks
- **Default (3)**: Balance between thoroughness and speed
- **Range**: 1 (fast) to 5 (thorough)

## Preset Configurations

### Available Presets

#### 1. **15N1H Preset**
```json
{
  "description": "Optimized for 15N-1H HSQC spectra",
  "sn_threshold": 2.5,
  "expected_peaks": 60,
  "optimal_range": [30, 80],
  "quality_threshold": 0.8
}
```
- **Use Case**: Protein backbone assignment, structural studies
- **Characteristics**: Higher S/N threshold, moderate peak count

#### 2. **13C1H Preset**
```json
{
  "description": "Optimized for 13C-1H HMQC/HSQC spectra",
  "sn_threshold": 2.0,
  "expected_peaks": 40,
  "optimal_range": [20, 60],
  "quality_threshold": 0.8
}
```
- **Use Case**: Carbon-proton correlation analysis
- **Characteristics**: Moderate settings for diverse carbon environments

#### 3. **1H Preset**
```json
{
  "description": "Optimized for 1H spectra",
  "sn_threshold": 2.0,
  "expected_peaks": 30,
  "optimal_range": [15, 40],
  "quality_threshold": 0.8
}
```
- **Use Case**: 1D proton spectrum analysis
- **Characteristics**: Lower expected peaks, handles overlap well

#### 4. **Conservative Preset**
```json
{
  "description": "High-quality data only",
  "sn_thresholds": {
    "15N1H": 3.0,
    "13C1H": 2.5,
    "1H": 2.5
  },
  "quality_threshold": 0.9,
  "optimization_strategy": "conservative"
}
```
- **Use Case**: Publication-quality datasets, method validation
- **Characteristics**: Strict quality control, fewer but excellent peaks

#### 5. **Aggressive Preset**
```json
{
  "description": "Maximum data collection",
  "sn_thresholds": {
    "15N1H": 1.8,
    "13C1H": 1.5,
    "1H": 1.5
  },
  "expected_peaks": {
    "15N1H": 200,
    "13C1H": 150,
    "1H": 100
  },
  "quality_threshold": 0.6
}
```
- **Use Case**: ML training data collection, exploratory analysis
- **Characteristics**: Lower thresholds, maximum peak collection

#### 6. **High-Throughput Preset**
```json
{
  "description": "Speed-optimized processing",
  "max_fitting_attempts": 1,
  "skip_on_error": true,
  "enable_auto_optimization": false,
  "detailed_logging": false
}
```
- **Use Case**: Large dataset processing, screening applications
- **Characteristics**: Minimal processing time, basic quality control

## ML Data Collection

### Feature Extraction

The ML data collector automatically extracts comprehensive features during processing:

#### Spectral Features (27 features)
- Peak height, baseline level, dynamic range
- Signal-to-noise ratio, spectral range
- Peak width measurements (FWHM, moments, percentiles)
- Peak asymmetry, tailing, fronting characteristics
- Baseline properties (drift, slope, curvature)
- Noise characteristics (variance, correlation)
- Intensity statistics (mean, std, skewness, kurtosis)

#### Chemical Features (11 features)
- Nucleus type, chemical shift center
- Chemical region classification
- Shift range and multiplicity estimates
- Coupling environment analysis
- Relaxation and CSA estimates

#### Context Features (11 features)
- Detection confidence, peak complexity
- Nearby peak count, overlap severity
- Noise level, interference sources
- Peak isolation score, baseline stability

#### Optimization Features (14 features)
- Initial parameters, parameter bounds
- Convergence iterations, optimization method
- Parameter uncertainties, gradient information
- Optimization difficulty assessment

#### Target Parameters (6+ features)
- Voigt parameters: amplitude, center, sigma, gamma
- Baseline correction, parameter uncertainties
- Quality metrics: R², RMSE, goodness-of-fit tests

### Data Quality Control

The ML collector implements automatic quality filtering:

```python
quality_metrics = {
    'r_squared': 0.85,           # Fit quality
    'parameter_validity': True,   # Physical constraints
    'convergence_success': True,  # Optimization success
    'residual_analysis': 'good'   # Systematic error check
}
```

## Usage Examples

### Basic Usage

```python
from lunaNMR.batch_processing.batch_processor import BatchProcessor

# Initialize with default settings
processor = BatchProcessor()

# Process all files in a directory
results = processor.process_folder(
    folder_path="data/nmr_spectra",
    nucleus_type="15N1H",  # Optional: auto-detect if None
    auto_optimize=True     # Enable parameter optimization
)

print(f"Processed {results['processed_files']} files")
print(f"ML samples collected: {results['total_ml_samples']}")
```

### Using Presets

```python
from lunaNMR.batch_processing.config_manager import ConfigManager
from lunaNMR.batch_processing.batch_processor import BatchProcessor

# Load a specific preset
config_manager = ConfigManager()
config = config_manager.get_preset_config('aggressive')

# Initialize processor with preset
processor = BatchProcessor(config)

# Process with aggressive settings
results = processor.process_folder("data/spectra")
```

### Custom Configuration

```python
# Create custom configuration
custom_config = {
    'sn_thresholds': {'default': 2.8},
    'quality_threshold': 0.75,
    'max_fitting_attempts': 2,
    'detailed_logging': True
}

processor = BatchProcessor(custom_config)
results = processor.process_folder("data/custom_analysis")
```

### Configuration Files

#### Creating Preset Files

```python
from lunaNMR.batch_processing.config_manager import ConfigManager

cm = ConfigManager()

# Create all preset configuration files
cm.create_preset_config_files("config_presets/")

# Create example configuration with comments
cm.create_example_config("my_config.json")
```

#### Using Configuration Files

```python
# Load configuration from file
config = cm.load_config("my_custom_config.json")
processor = BatchProcessor(config)
```

### Example Configuration File

```json
{
  "_comment": "Custom configuration for protein NMR analysis",
  "processing": {
    "sn_thresholds": {
      "15N1H": 2.8,
      "13C1H": 2.2,
      "default": 2.5
    },
    "expected_peaks": {
      "15N1H": 80,
      "default": 100
    },
    "quality_threshold": 0.75,
    "max_fitting_attempts": 2
  },
  "file_handling": {
    "extensions": [".ft2", ".pipe"],
    "exclude_patterns": ["*backup*", "*temp*"]
  },
  "optimization": {
    "enable_auto_optimization": true,
    "optimization_strategy": "balanced"
  },
  "logging": {
    "level": "INFO",
    "detailed_logging": true
  }
}
```

## Performance Considerations

### Parallel Processing

The batch processor supports multi-core parallel processing:

- **Threading**: Used for I/O-intensive operations
- **Multiprocessing**: Available for CPU-intensive fitting
- **Automatic Scaling**: Worker count adapts to available cores
- **Memory Management**: Efficient handling of large datasets

### Optimization Strategies

1. **Conservative**: Thorough processing, highest quality
   - Use case: Publication datasets, method validation
   - Trade-off: Slower processing, fewer peaks

2. **Balanced**: Good compromise between speed and quality
   - Use case: General research applications
   - Trade-off: Moderate speed and quality

3. **Aggressive**: Maximum data collection
   - Use case: ML training, exploratory analysis
   - Trade-off: Faster processing, more noise

## Output and Results

### Processing Summary

```
============================================================
BATCH PROCESSING SUMMARY
============================================================
Total files found: 25
Successfully processed: 24
Failed processing: 1
Total peaks detected: 3750
Total peaks fitted: 3245
ML training samples collected: 3245
Processing time: 1847.3 seconds

Success rate: 96.0%
ML training data generation: SUCCESSFUL
Ready for Phase 2 ML model development
============================================================
```

### ML Training Data Files

The system generates several output files:

1. **enhanced_training_batch_YYYYMMDD_HHMMSS.pkl**
   - Binary training data with all features
   - Used for ML model development

2. **enhanced_training_batch_YYYYMMDD_HHMMSS.json**
   - Metadata and collection statistics
   - Human-readable processing summary

3. **batch_processing_YYYYMMDD_HHMMSS.log**
   - Detailed processing log
   - Error reports and diagnostics

### Data Analysis Tools

Use the provided analyzer to inspect collected data:

```bash
cd ml_training_data/
python analyze_ml_data.py  # Analyze latest file
python analyze_ml_data.py my_data.pkl --detailed
```

## Troubleshooting

### Common Issues

1. **No ML Data Collected**
   - Check quality thresholds (try lower values)
   - Verify peak detection sensitivity
   - Review processing logs for errors

2. **Low Processing Success Rate**
   - Reduce S/N thresholds
   - Increase max_fitting_attempts
   - Check file format compatibility

3. **Performance Issues**
   - Disable detailed logging for speed
   - Reduce expected_peaks for large datasets
   - Use high_throughput preset

### Debug Mode

Enable comprehensive debugging:

```python
import logging
logging.basicConfig(level=logging.DEBUG)

processor = BatchProcessor({
    'detailed_logging': True,
    'progress_interval': 1
})
```

## Integration with GUI

The batch processor maintains full compatibility with the GUI system:

- **Parameter Synchronization**: Uses same NMRParameterManager
- **Processing Methods**: Identical enhanced_peak_fitting approach
- **Quality Standards**: Consistent R² and validation criteria
- **ML Collection**: Same Enhanced ML Training Data Collector v2.0

This ensures seamless transition between interactive analysis and batch processing workflows.

## Future Enhancements

### Phase 2 Developments

1. **Intelligent Parameter Optimization**
   - Automatic S/N threshold tuning
   - Adaptive quality thresholds
   - Spectrum-specific parameter learning

2. **Advanced ML Features**
   - Real-time model training
   - Active learning for difficult spectra
   - Uncertainty quantification

3. **Enhanced Parallel Processing**
   - GPU acceleration for fitting
   - Distributed processing capabilities
   - Cloud computing integration

## References

- Enhanced Voigt Fitting Documentation
- NMR Parameter Manager Guide
- ML Training Data Collection Specifications
- GUI Processing Workflow Documentation

---

*This guide covers the complete batch processing and ML data collection system in LunaNMR. For additional support or feature requests, please refer to the project documentation or contact the development team.*