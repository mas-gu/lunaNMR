#!/usr/bin/env python3
"""
Memory Monitoring Utilities for Parallel Processing

This module provides tools to monitor memory usage and verify
proper cleanup of shared memory blocks during parallel processing.

Author: Claude Code
Date: 2025-09-19
"""

import psutil
import os
import time
import threading
import gc


class MemoryMonitor:
    """
    Monitor memory usage and shared memory blocks during parallel processing.
    """

    def __init__(self, sample_interval=1.0):
        """
        Initialize memory monitor

        Args:
            sample_interval: Time between memory samples in seconds
        """
        self.sample_interval = sample_interval
        self.monitoring = False
        self.memory_samples = []
        self.shared_memory_blocks = []
        self.monitor_thread = None

    def start_monitoring(self):
        """Start memory monitoring in background thread"""
        if self.monitoring:
            return

        self.monitoring = True
        self.memory_samples = []
        self.shared_memory_blocks = []

        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()

        print("🔍 Memory monitoring started")

    def stop_monitoring(self):
        """Stop memory monitoring"""
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join(timeout=5)

        print("⏹️ Memory monitoring stopped")

    def _monitor_loop(self):
        """Main monitoring loop"""
        process = psutil.Process()

        while self.monitoring:
            try:
                # Get memory info
                memory_info = process.memory_info()
                memory_percent = process.memory_percent()

                # Get system memory
                system_memory = psutil.virtual_memory()

                # Count shared memory blocks
                shared_memory_count = self._count_shared_memory_blocks()

                sample = {
                    'timestamp': time.time(),
                    'rss_mb': memory_info.rss / 1024 / 1024,
                    'vms_mb': memory_info.vms / 1024 / 1024,
                    'memory_percent': memory_percent,
                    'system_available_mb': system_memory.available / 1024 / 1024,
                    'system_used_percent': system_memory.percent,
                    'shared_memory_blocks': shared_memory_count
                }

                self.memory_samples.append(sample)

                # Keep only recent samples (last 5 minutes)
                cutoff_time = time.time() - 300
                self.memory_samples = [s for s in self.memory_samples
                                     if s['timestamp'] > cutoff_time]

                time.sleep(self.sample_interval)

            except Exception as e:
                print(f"⚠️ Memory monitoring error: {e}")
                time.sleep(self.sample_interval)

    def _count_shared_memory_blocks(self):
        """Count system shared memory blocks"""
        try:
            # On Linux/macOS, shared memory is in /dev/shm
            if os.path.exists('/dev/shm'):
                return len(os.listdir('/dev/shm'))

            # On macOS, also check /tmp for shared memory
            shm_count = 0
            if os.path.exists('/tmp'):
                for item in os.listdir('/tmp'):
                    if item.startswith('psm_') or item.startswith('shm_'):
                        shm_count += 1
            return shm_count

        except Exception:
            return 0

    def get_current_memory_usage(self):
        """Get current memory usage snapshot"""
        if not self.memory_samples:
            return None

        return self.memory_samples[-1]

    def get_memory_statistics(self):
        """Get memory usage statistics"""
        if not self.memory_samples:
            return None

        rss_values = [s['rss_mb'] for s in self.memory_samples]
        shm_values = [s['shared_memory_blocks'] for s in self.memory_samples]

        return {
            'rss_min_mb': min(rss_values),
            'rss_max_mb': max(rss_values),
            'rss_avg_mb': sum(rss_values) / len(rss_values),
            'rss_current_mb': rss_values[-1],
            'shared_memory_min': min(shm_values),
            'shared_memory_max': max(shm_values),
            'shared_memory_current': shm_values[-1],
            'sample_count': len(self.memory_samples),
            'monitoring_duration': self.memory_samples[-1]['timestamp'] - self.memory_samples[0]['timestamp']
        }

    def check_memory_leaks(self, baseline_rss_mb=None):
        """
        Check for potential memory leaks

        Args:
            baseline_rss_mb: Expected baseline RSS memory in MB

        Returns:
            Dict with leak detection results
        """
        if len(self.memory_samples) < 10:
            return {"error": "Insufficient samples for leak detection"}

        stats = self.get_memory_statistics()
        current_usage = self.get_current_memory_usage()

        # Calculate memory growth trend
        first_half = self.memory_samples[:len(self.memory_samples)//2]
        second_half = self.memory_samples[len(self.memory_samples)//2:]

        first_half_avg = sum(s['rss_mb'] for s in first_half) / len(first_half)
        second_half_avg = sum(s['rss_mb'] for s in second_half) / len(second_half)

        memory_growth = second_half_avg - first_half_avg
        growth_percent = (memory_growth / first_half_avg) * 100 if first_half_avg > 0 else 0

        # Check for shared memory block accumulation
        first_shm_avg = sum(s['shared_memory_blocks'] for s in first_half) / len(first_half)
        second_shm_avg = sum(s['shared_memory_blocks'] for s in second_half) / len(second_half)
        shm_growth = second_shm_avg - first_shm_avg

        leak_indicators = []

        # Memory growth indicators
        if growth_percent > 10:
            leak_indicators.append(f"High memory growth: {growth_percent:.1f}%")

        if memory_growth > 100:  # More than 100MB growth
            leak_indicators.append(f"Large memory increase: {memory_growth:.1f}MB")

        # Shared memory indicators
        if shm_growth > 5:
            leak_indicators.append(f"Shared memory block accumulation: +{shm_growth:.1f} blocks")

        if current_usage['shared_memory_blocks'] > 50:
            leak_indicators.append(f"High shared memory block count: {current_usage['shared_memory_blocks']}")

        # Baseline comparison
        if baseline_rss_mb and current_usage['rss_mb'] > baseline_rss_mb * 2:
            leak_indicators.append(f"Memory usage doubled from baseline: {current_usage['rss_mb']:.1f}MB vs {baseline_rss_mb:.1f}MB")

        return {
            'memory_growth_mb': memory_growth,
            'memory_growth_percent': growth_percent,
            'shared_memory_growth': shm_growth,
            'leak_indicators': leak_indicators,
            'current_stats': stats,
            'has_potential_leaks': len(leak_indicators) > 0
        }

    def force_cleanup_test(self):
        """
        Force garbage collection and test cleanup effectiveness
        """
        print("🧹 Forcing garbage collection...")

        before_stats = self.get_current_memory_usage()

        # Force garbage collection
        gc.collect()
        time.sleep(1)  # Allow time for cleanup

        after_stats = self.get_current_memory_usage()

        if before_stats and after_stats:
            memory_freed = before_stats['rss_mb'] - after_stats['rss_mb']
            print(f"   Memory freed: {memory_freed:.1f}MB")
            return memory_freed

        return 0

    def print_memory_report(self):
        """Print comprehensive memory report"""
        stats = self.get_memory_statistics()
        if not stats:
            print("❌ No memory statistics available")
            return

        print("\n📊 Memory Usage Report")
        print("=" * 50)
        print("RSS Memory:")
        print(f"  Current: {stats['rss_current_mb']:.1f}MB")
        print(f"  Average: {stats['rss_avg_mb']:.1f}MB")
        print(f"  Peak:    {stats['rss_max_mb']:.1f}MB")
        print(f"  Minimum: {stats['rss_min_mb']:.1f}MB")

        print("\nShared Memory Blocks:")
        print(f"  Current: {stats['shared_memory_current']}")
        print(f"  Peak:    {stats['shared_memory_max']}")
        print(f"  Minimum: {stats['shared_memory_min']}")

        print("\nMonitoring:")
        print(f"  Duration: {stats['monitoring_duration']:.1f}s")
        print(f"  Samples:  {stats['sample_count']}")

        # Check for leaks
        leak_check = self.check_memory_leaks()
        if leak_check.get('has_potential_leaks', False):
            print("\n⚠️ Potential Memory Issues:")
            for indicator in leak_check['leak_indicators']:
                print(f"  - {indicator}")
        else:
            print("\n✅ No significant memory leaks detected")


def test_parallel_processor_memory():
    """
    Test function to verify ParallelVoigtProcessor memory management
    """
    print("🧪 Testing ParallelVoigtProcessor memory management")

    # This would be called from a test script that imports the actual classes
    monitor = MemoryMonitor(sample_interval=0.5)
    monitor.start_monitoring()

    try:
        # Simulate some processing time
        print("⏳ Simulating parallel processing...")
        time.sleep(10)

        # Test cleanup
        monitor.force_cleanup_test()
        time.sleep(2)

        # Generate report
        monitor.print_memory_report()

    finally:
        monitor.stop_monitoring()


if __name__ == "__main__":
    test_parallel_processor_memory()