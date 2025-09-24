#!/usr/bin/env python3
"""
Memory Stability Test for Parallel Voigt Processing

This script tests the memory management fixes in ParallelVoigtProcessor
to ensure no memory leaks occur during parallel processing operations.

Author: Claude Code
Date: 2025-09-19
"""

import sys
import os
import time
import numpy as np
import pandas as pd
from pathlib import Path

# Add lunaNMR to path
lunaNMR_path = Path(__file__).parent.parent
sys.path.insert(0, str(lunaNMR_path))

try:
    from validation.memory_monitoring import MemoryMonitor
    from core.parallel_voigt_processor import ParallelVoigtProcessor
    from core.core_integrator import EnhancedVoigtIntegrator
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("This test requires the lunaNMR package structure")
    sys.exit(1)


class MockEnhancedVoigtFitter:
    """
    Mock fitter for testing memory management without full NMR data
    """

    def __init__(self, data_size_mb=50):
        """Create mock fitter with specified data size"""
        # Create mock NMR data
        data_points = int((data_size_mb * 1024 * 1024) / 8)  # 8 bytes per float64
        self.nmr_data = np.random.random((data_points,)).astype(np.float64)

        # Create mock PPM axes
        self.ppm_x_axis = np.linspace(0, 10, data_points)
        self.ppm_y_axis = np.linspace(0, 10, 1000)

        # Mock parameters
        self.baseline_lambda = 1e6
        self.baseline_p = 0.001
        self.max_iterations = 1000
        self.min_r_squared = 0.8

        print(f"🔬 Mock fitter created with {data_size_mb}MB data ({data_points} points)")

    def enhanced_peak_fitting(self, peak_x, peak_y, assignment):
        """Mock fitting that simulates processing time"""
        time.sleep(0.1)  # Simulate processing time
        return {
            'assignment': assignment,
            'position_x': peak_x,
            'position_y': peak_y,
            'avg_r_squared': 0.95,
            'fit_quality': 'good',
            'mock_result': True
        }


def create_test_peak_list(num_peaks=20):
    """Create a test peak list for memory testing"""
    peak_data = []
    for i in range(num_peaks):
        peak_data.append({
            'Assignment': f'Peak_{i+1}',
            'Position_X': np.random.uniform(1, 9),
            'Position_Y': np.random.uniform(1, 9),
            'Height': np.random.uniform(0.1, 1.0),
            'Linewidth': np.random.uniform(0.01, 0.1)
        })

    return pd.DataFrame(peak_data)


def test_single_execution_memory():
    """Test memory usage for a single parallel execution"""
    print("\n🧪 Test 1: Single Execution Memory Management")
    print("=" * 60)

    monitor = MemoryMonitor(sample_interval=0.5)
    monitor.start_monitoring()

    try:
        # Get baseline memory
        time.sleep(1)
        baseline = monitor.get_current_memory_usage()
        baseline_rss = baseline['rss_mb'] if baseline else 0

        print(f"📊 Baseline memory: {baseline_rss:.1f}MB")

        # Create mock fitter and processor
        mock_fitter = MockEnhancedVoigtFitter(data_size_mb=100)
        processor = ParallelVoigtProcessor(mock_fitter, max_workers=4)

        # Create test peak list
        peak_list = create_test_peak_list(num_peaks=10)
        print(f"📋 Created peak list with {len(peak_list)} peaks")

        # Run parallel processing
        print("⚡ Starting parallel processing...")
        results = processor.fit_all_peaks_parallel(peak_list)

        print(f"✅ Processing completed, {len(results)} results")

        # Allow time for cleanup
        time.sleep(2)

        # Force cleanup test
        freed_memory = monitor.force_cleanup_test()

        # Final memory check
        final_usage = monitor.get_current_memory_usage()
        final_rss = final_usage['rss_mb'] if final_usage else 0

        print(f"📊 Final memory: {final_rss:.1f}MB")
        print(f"📈 Memory growth: {final_rss - baseline_rss:.1f}MB")

        # Check for leaks
        leak_check = monitor.check_memory_leaks(baseline_rss)
        if leak_check.get('has_potential_leaks', False):
            print("❌ Potential memory leaks detected:")
            for indicator in leak_check['leak_indicators']:
                print(f"   - {indicator}")
            return False
        else:
            print("✅ No significant memory leaks detected")
            return True

    finally:
        monitor.stop_monitoring()


def test_repeated_execution_memory():
    """Test memory stability over repeated executions"""
    print("\n🧪 Test 2: Repeated Execution Memory Stability")
    print("=" * 60)

    monitor = MemoryMonitor(sample_interval=0.5)
    monitor.start_monitoring()

    try:
        # Get baseline
        time.sleep(1)
        baseline = monitor.get_current_memory_usage()
        baseline_rss = baseline['rss_mb'] if baseline else 0

        print(f"📊 Baseline memory: {baseline_rss:.1f}MB")

        memory_per_run = []

        for run_num in range(5):
            print(f"\n🔄 Run {run_num + 1}/5")

            # Create fresh instances each time
            mock_fitter = MockEnhancedVoigtFitter(data_size_mb=75)
            processor = ParallelVoigtProcessor(mock_fitter, max_workers=3)

            # Create test peak list
            peak_list = create_test_peak_list(num_peaks=8)

            # Run processing
            results = processor.fit_all_peaks_parallel(peak_list)
            print(f"   Results: {len(results)}")

            # Check memory after each run
            current_usage = monitor.get_current_memory_usage()
            current_rss = current_usage['rss_mb'] if current_usage else 0
            memory_per_run.append(current_rss)

            print(f"   Memory: {current_rss:.1f}MB (+{current_rss - baseline_rss:.1f}MB)")

            # Small delay between runs
            time.sleep(1)

        # Analyze memory growth pattern
        print(f"\n📈 Memory Growth Analysis:")
        for i, memory in enumerate(memory_per_run):
            growth = memory - baseline_rss
            print(f"   Run {i+1}: {memory:.1f}MB (+{growth:.1f}MB)")

        # Check if memory is stabilizing (last run shouldn't be much higher than first)
        first_run_growth = memory_per_run[0] - baseline_rss
        last_run_growth = memory_per_run[-1] - baseline_rss
        growth_increase = last_run_growth - first_run_growth

        print(f"\n📊 Growth Analysis:")
        print(f"   First run growth: {first_run_growth:.1f}MB")
        print(f"   Last run growth:  {last_run_growth:.1f}MB")
        print(f"   Additional growth: {growth_increase:.1f}MB")

        if growth_increase > 50:  # More than 50MB additional growth
            print("❌ Significant memory accumulation detected")
            return False
        else:
            print("✅ Memory usage appears stable")
            return True

    finally:
        monitor.stop_monitoring()


def test_exception_handling_memory():
    """Test memory cleanup when exceptions occur"""
    print("\n🧪 Test 3: Exception Handling Memory Cleanup")
    print("=" * 60)

    monitor = MemoryMonitor(sample_interval=0.5)
    monitor.start_monitoring()

    try:
        baseline = monitor.get_current_memory_usage()
        baseline_rss = baseline['rss_mb'] if baseline else 0

        print(f"📊 Baseline memory: {baseline_rss:.1f}MB")

        # Test with intentional failure
        class FailingMockFitter(MockEnhancedVoigtFitter):
            def enhanced_peak_fitting(self, peak_x, peak_y, assignment):
                if assignment == 'Peak_3':  # Fail on third peak
                    raise ValueError("Intentional test failure")
                return super().enhanced_peak_fitting(peak_x, peak_y, assignment)

        mock_fitter = FailingMockFitter(data_size_mb=80)
        processor = ParallelVoigtProcessor(mock_fitter, max_workers=4)

        peak_list = create_test_peak_list(num_peaks=6)

        print("⚡ Starting processing with intentional failures...")
        try:
            results = processor.fit_all_peaks_parallel(peak_list)
            print(f"✅ Processing completed despite failures, {len(results)} results")
        except Exception as e:
            print(f"❌ Processing failed as expected: {e}")

        # Check memory after failure
        time.sleep(2)
        monitor.force_cleanup_test()

        final_usage = monitor.get_current_memory_usage()
        final_rss = final_usage['rss_mb'] if final_usage else 0
        growth = final_rss - baseline_rss

        print(f"📊 Final memory: {final_rss:.1f}MB (+{growth:.1f}MB)")

        if growth > 100:  # More than 100MB growth after failure
            print("❌ Excessive memory retention after failure")
            return False
        else:
            print("✅ Memory properly cleaned up after failure")
            return True

    finally:
        monitor.stop_monitoring()


def main():
    """Run all memory stability tests"""
    print("🚀 LunaNMR Parallel Processing Memory Stability Tests")
    print("=" * 70)

    test_results = []

    # Run all tests
    test_results.append(("Single Execution", test_single_execution_memory()))
    test_results.append(("Repeated Execution", test_repeated_execution_memory()))
    test_results.append(("Exception Handling", test_exception_handling_memory()))

    # Summary
    print("\n📋 Test Summary")
    print("=" * 50)

    all_passed = True
    for test_name, passed in test_results:
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"{test_name:20s}: {status}")
        if not passed:
            all_passed = False

    print("\n" + "=" * 50)
    if all_passed:
        print("🎉 All memory stability tests PASSED!")
        print("✅ ParallelVoigtProcessor memory management is working correctly")
    else:
        print("❌ Some memory stability tests FAILED!")
        print("⚠️ Memory leaks or management issues detected")

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)