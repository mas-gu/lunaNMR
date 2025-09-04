#!/usr/bin/env python3
"""
Simple Pattern Matching for NMR Peak Assignment
Logic: Detect → Filter → Graph → Pattern Match → Assign (coordinates never change)

Author: Guillaume Mas
Date: 2025
"""

import numpy as np
import networkx as nx
from scipy.spatial.distance import cdist
import itertools
import sys
import time

class SimplePatternMatcher:
    """Simple pattern matcher that preserves detected coordinates with spatial locality optimization"""

    def __init__(self, proximity_radius_x=0.2, proximity_radius_y=5.0, isolated_peak_threshold=0.1, enable_fast_isolated=True, assignment_only=False):
        #self.debug = True
        self.verbose = True

        # Spatial locality parameters (NEW)
        self.proximity_radius_x = proximity_radius_x      # 1H search window (ppm)
        self.proximity_radius_y = proximity_radius_y      # 15N/13C search window (ppm)
        self.isolated_peak_threshold = isolated_peak_threshold  # Max distance for direct assignment
        self.enable_fast_isolated = enable_fast_isolated  # Toggle for isolated peak optimization
        self.assignment_only = assignment_only            # Only transfer assignments, no new peak discovery

    def assign_peaks(self, detected_peaks, reference_peaks, max_distance=0.5):
        """
        Spatial locality pattern matching that preserves detected coordinates

        Args:
            detected_peaks: list of dicts with 'position_x', 'position_y', 'intensity'
            reference_peaks: list of dicts with 'Position_X', 'Position_Y', 'Assignment'

        Returns:
            list of assigned peaks with original detected coordinates + assignments
        """
        if self.verbose:
            print(f"🔍 Spatial Locality Pattern Matcher:")
            print(f"   Detected peaks: {len(detected_peaks)}")
            print(f"   Reference peaks: {len(reference_peaks)}")
            print(f"   Search window: {self.proximity_radius_x:.2f} × {self.proximity_radius_y:.1f} ppm")

        # CRITICAL DEBUG: Check if detected peaks have corrupted coordinates
        if len(detected_peaks) > 0:
            det_peak = detected_peaks[0]
            print(f"   First detected peak: X={det_peak.get('position_x', 'MISSING')}, Y={det_peak.get('position_y', 'MISSING')}")
        if len(reference_peaks) > 0:
            ref_peak = reference_peaks[0]
            print(f"   First reference peak: X={ref_peak.get('Position_X', 'MISSING')}, Y={ref_peak.get('Position_Y', 'MISSING')}")

        # Check if detected coordinates match reference coordinates (corruption check)
        if len(detected_peaks) > 0 and len(reference_peaks) > 0:
            det_x = detected_peaks[0].get('position_x', 0)
            ref_x = reference_peaks[0].get('Position_X', 0)
            if abs(det_x - ref_x) < 0.001:
                print(f"   COORDINATE CORRUPTION: Detected X={det_x:.4f} matches reference X={ref_x:.4f}!")
            else:
                print(f"   Coordinates different: Detected X={det_x:.4f} vs Reference X={ref_x:.4f}")

        # Step 1: Create spatial index (local candidates for each reference peak)
        if self.verbose:
            print(f"🗺️ Step 1/5: Creating spatial index...", end="", flush=True)
        local_candidates = self._create_spatial_index(detected_peaks, reference_peaks)
        if self.verbose:
            total_candidates = sum(len(candidates) for candidates in local_candidates.values())
            avg_candidates = total_candidates / len(reference_peaks) if reference_peaks else 0
            print(f" ✅ (avg {avg_candidates:.1f} candidates/ref)")

        # Step 2: Fast isolated peak assignment
        isolated_assignments = []
        if self.enable_fast_isolated:
            if self.verbose:
                print(f"⚡ Step 2/5: Fast isolated assignment...", end="", flush=True)
            isolated_assignments = self._assign_isolated_peaks(local_candidates, detected_peaks, reference_peaks)
            if self.verbose:
                print(f" ✅ {len(isolated_assignments)} isolated")

        # Step 3: Pattern matching for complex cases
        pattern_assignments = []
        if self.verbose:
            print(f"🔍 Step 3/5: Localized pattern matching...", end="", flush=True)
        pattern_assignments = self._localized_pattern_matching(local_candidates, detected_peaks, reference_peaks,
                                                             isolated_assignments, max_distance)
        if self.verbose:
            print(f" ✅ {len(pattern_assignments)} patterns")

        # Combine all assignments
        all_assignments = isolated_assignments + pattern_assignments

        # Step 4: Create final peaks (detected coordinates + reference assignments)
        if self.verbose:
            print(f"🎯 Step 4/5: Creating assignments...", end="", flush=True)
        final_peaks = []
        assigned_detected = set()

        for assignment in all_assignments:
            detected_idx = assignment['detected_idx']
            reference_idx = assignment['reference_idx']

            if detected_idx in assigned_detected:
                continue  # Already assigned

            detected_peak = detected_peaks[detected_idx]
            reference_peak = reference_peaks[reference_idx]

            # CRITICAL: Use detected coordinates, reference assignment
            final_peak = {
                'Position_X': detected_peak['position_x'],  # DETECTED coordinate
                'Position_Y': detected_peak['position_y'],  # DETECTED coordinate
                'Intensity': detected_peak['intensity'],
                'Assignment': reference_peak.get('Assignment', f'Peak_{reference_idx+1}'),
                'Volume': reference_peak.get('Volume', 0),
                'Detection_Method': assignment.get('method', 'localized_pattern_match'),
                'Match_Distance': assignment['distance'],
                'Match_Confidence': assignment['confidence']
            }

            final_peaks.append(final_peak)
            assigned_detected.add(detected_idx)

        if self.verbose:
            print(f" ✅")

        # Step 5: Add unassigned detected peaks (ONLY if not in assignment_only mode)
        unassigned_count = 0
        if not self.assignment_only:
            if self.verbose:
                print(f"➕ Step 5/5: Adding unassigned peaks...", end="", flush=True)
            for i, detected_peak in enumerate(detected_peaks):
                if i not in assigned_detected:
                    unassigned_peak = {
                        'Position_X': detected_peak['position_x'],  # DETECTED coordinate
                        'Position_Y': detected_peak['position_y'],  # DETECTED coordinate
                        'Intensity': detected_peak['intensity'],
                        'Assignment': f'Unassigned_{i+1}',
                        'Volume': 0,
                        'Detection_Method': 'detected_unassigned',
                        'Match_Distance': 0,
                        'Match_Confidence': 0
                    }
                    final_peaks.append(unassigned_peak)
                    unassigned_count += 1
        else:
            # Assignment-only mode: skip unassigned peaks
            unassigned_count = len(detected_peaks) - len(assigned_detected)
            if self.verbose:
                print(f"⚠️ Step 5/5: Assignment-only mode - skipping {unassigned_count} unassigned peaks...", end="", flush=True)

        if self.verbose:
            if not self.assignment_only:
                print(f" ✅ {unassigned_count} unassigned")
                print(f"✅ Complete: {len(final_peaks)} peaks ({len(all_assignments)} assigned: {len(isolated_assignments)} isolated + {len(pattern_assignments)} patterns, {unassigned_count} unassigned)")
            else:
                print(f" ✅ {unassigned_count} skipped")
                print(f"✅ Complete: {len(final_peaks)} peaks ({len(all_assignments)} assigned: {len(isolated_assignments)} isolated + {len(pattern_assignments)} patterns, {unassigned_count} unassigned peaks excluded)")

        return final_peaks

    def _create_spatial_index(self, detected_peaks, reference_peaks):
        """Create local candidate mapping for each reference peak (NEW - Spatial Optimization)"""
        local_candidates = {}

        for ref_idx, ref_peak in enumerate(reference_peaks):
            ref_x = ref_peak.get('Position_X', 0)
            ref_y = ref_peak.get('Position_Y', 0)
            candidates = []

            for det_idx, det_peak in enumerate(detected_peaks):
                det_x = det_peak['position_x']
                det_y = det_peak['position_y']

                # Check if detected peak is within proximity window of reference peak
                x_diff = abs(ref_x - det_x)
                y_diff = abs(ref_y - det_y)

                if x_diff <= self.proximity_radius_x and y_diff <= self.proximity_radius_y:
                    distance = np.sqrt(x_diff**2 + y_diff**2)
                    candidates.append({
                        'idx': det_idx,
                        'distance': distance
                    })

            # Sort candidates by distance (closest first)
            candidates.sort(key=lambda x: x['distance'])
            local_candidates[ref_idx] = [c['idx'] for c in candidates]

        return local_candidates

    def _assign_isolated_peaks(self, local_candidates, detected_peaks, reference_peaks):
        """Fast assignment for reference peaks with single local candidates (NEW - Performance Optimization)"""
        isolated_assignments = []

        for ref_idx, candidates in local_candidates.items():
            if len(candidates) == 1:
                det_idx = candidates[0]

                # Calculate distance
                ref_peak = reference_peaks[ref_idx]
                det_peak = detected_peaks[det_idx]

                ref_x = ref_peak.get('Position_X', 0)
                ref_y = ref_peak.get('Position_Y', 0)
                det_x = det_peak['position_x']
                det_y = det_peak['position_y']

                distance = np.sqrt((ref_x - det_x)**2 + (ref_y - det_y)**2)

                # Only assign if within threshold
                if distance <= self.isolated_peak_threshold:
                    confidence = max(0, 1.0 - distance/self.isolated_peak_threshold)
                    isolated_assignments.append({
                        'detected_idx': det_idx,
                        'reference_idx': ref_idx,
                        'distance': distance,
                        'confidence': confidence,
                        'method': 'isolated_fast'
                    })

        return isolated_assignments

    def _localized_pattern_matching(self, local_candidates, detected_peaks, reference_peaks, isolated_assignments, max_distance):
        """Pattern matching within local neighborhoods only (NEW - Locality Constraint)"""
        pattern_assignments = []

        # Track which detected peaks are already assigned by isolated matching
        assigned_detected = {assignment['detected_idx'] for assignment in isolated_assignments}
        assigned_reference = {assignment['reference_idx'] for assignment in isolated_assignments}

        # For remaining unassigned reference peaks with multiple candidates, do pattern matching
        multi_candidate_refs = []
        for ref_idx, candidates in local_candidates.items():
            if ref_idx not in assigned_reference and len(candidates) > 1:
                # Filter out already assigned detected peaks
                available_candidates = [det_idx for det_idx in candidates if det_idx not in assigned_detected]
                if available_candidates:
                    multi_candidate_refs.append({
                        'ref_idx': ref_idx,
                        'candidates': available_candidates
                    })

        # Simple greedy assignment for now (can be enhanced with graph matching later)
        for ref_info in multi_candidate_refs:
            ref_idx = ref_info['ref_idx']
            candidates = ref_info['candidates']

            ref_peak = reference_peaks[ref_idx]
            ref_x = ref_peak.get('Position_X', 0)
            ref_y = ref_peak.get('Position_Y', 0)

            # Find best candidate (closest within max_distance)
            best_candidate = None
            best_distance = float('inf')

            for det_idx in candidates:
                if det_idx in assigned_detected:
                    continue

                det_peak = detected_peaks[det_idx]
                det_x = det_peak['position_x']
                det_y = det_peak['position_y']

                distance = np.sqrt((ref_x - det_x)**2 + (ref_y - det_y)**2)

                if distance <= max_distance and distance < best_distance:
                    best_distance = distance
                    best_candidate = det_idx

            # Assign best candidate if found
            if best_candidate is not None:
                confidence = max(0, 1.0 - best_distance/max_distance)
                pattern_assignments.append({
                    'detected_idx': best_candidate,
                    'reference_idx': ref_idx,
                    'distance': best_distance,
                    'confidence': confidence,
                    'method': 'localized_pattern'
                })
                assigned_detected.add(best_candidate)

        return pattern_assignments

    def _create_graph(self, peaks, graph_type):
        """Create networkx graph from peaks"""
        G = nx.Graph()

        for i, peak in enumerate(peaks):
            if graph_type == 'detected':
                x, y = peak['position_x'], peak['position_y']
            else:  # reference
                x, y = peak.get('Position_X', 0), peak.get('Position_Y', 0)

            G.add_node(i, x=x, y=y, intensity=peak.get('intensity', peak.get('Intensity', 1)))

        # Add edges between nearby peaks
        positions = np.array([[G.nodes[i]['x'], G.nodes[i]['y']] for i in G.nodes()])
        distances = cdist(positions, positions)

        # Connect peaks within reasonable distance
        max_connect_distance = 0.25  # ppm #GM added a 2.0
        for i in range(len(peaks)):
            for j in range(i+1, len(peaks)):
                if distances[i,j] < max_connect_distance:
                    G.add_edge(i, j, distance=distances[i,j])

        return G

    def _find_pattern_matches(self, detected_graph, reference_graph, max_distance):
        """Find pattern matches between small reference networks and detected peaks"""
        assignments = []
        total_combinations = 0
        processed_combinations = 0

        # Count total combinations for progress bar
        for network_size in range(1, 5): #1, 6
            if network_size <= len(reference_graph.nodes):
                from math import factorial
                n = len(reference_graph.nodes())
                r = network_size
                if r <= n:
                    combinations = factorial(n) // (factorial(r) * factorial(n - r))
                    total_combinations += combinations

        if self.verbose and total_combinations > 50:
            print(f"\n   🔍 Trying {total_combinations} pattern combinations...")

        # Try different network sizes (1-5 peaks)
        for network_size in range(1, 5): #1, 6
            if network_size > len(reference_graph.nodes):
                continue

            if self.verbose and network_size > 1:
                print(f"   📐 Size {network_size}: ", end="", flush=True)

            # Get all reference subgraphs of this size
            size_matches = 0
            for ref_nodes in itertools.combinations(reference_graph.nodes(), network_size):
                ref_subgraph = reference_graph.subgraph(ref_nodes)

                # Try to match this pattern in detected graph
                matches = self._match_subgraph(ref_subgraph, detected_graph, max_distance)
                assignments.extend(matches)
                size_matches += len(matches)

                processed_combinations += 1
                if self.verbose and total_combinations > 50 and processed_combinations % 20 == 0:
                    progress = int(10 * processed_combinations / total_combinations)
                    bar = "█" * progress + "░" * (10 - progress)
                    print(f"\r   🔍 Progress [{bar}] {processed_combinations}/{total_combinations}", end="", flush=True)

            if self.verbose and network_size > 1:
                print(f"{size_matches} matches")

        if self.verbose and total_combinations > 50:
            print(f"\r   🔍 Progress [██████████] {processed_combinations}/{total_combinations} - Complete")

        # Remove duplicate assignments (keep best confidence)
        unique_assignments = {}
        for assignment in assignments:
            key = assignment['detected_idx']
            if key not in unique_assignments or assignment['confidence'] > unique_assignments[key]['confidence']:
                unique_assignments[key] = assignment

        return list(unique_assignments.values())

    def _match_subgraph(self, ref_subgraph, detected_graph, max_distance):
        """Match a reference subgraph pattern against detected graph"""
        matches = []
        ref_nodes = list(ref_subgraph.nodes())

        # For single peaks, just do distance matching
        if len(ref_nodes) == 1:
            ref_node = ref_nodes[0]
            ref_x = ref_subgraph.nodes[ref_node]['x']
            ref_y = ref_subgraph.nodes[ref_node]['y']

            for det_node in detected_graph.nodes():
                det_x = detected_graph.nodes[det_node]['x']
                det_y = detected_graph.nodes[det_node]['y']

                distance = np.sqrt((ref_x - det_x)**2 + (ref_y - det_y)**2)
                if distance <= max_distance:
                    confidence = max(0, 1.0 - distance/max_distance)
                    matches.append({
                        'detected_idx': det_node,
                        'reference_idx': ref_node,
                        'distance': distance,
                        'confidence': confidence
                    })

        # For multi-peak patterns, do geometric matching
        elif len(ref_nodes) > 1:
            # Try all possible mappings of detected nodes to reference nodes
            for det_nodes in itertools.combinations(detected_graph.nodes(), len(ref_nodes)):
                if self._pattern_matches(ref_subgraph, ref_nodes, detected_graph, det_nodes, max_distance):
                    # Calculate average distance and confidence
                    total_distance = 0
                    for ref_idx, det_idx in zip(ref_nodes, det_nodes):
                        ref_x = ref_subgraph.nodes[ref_idx]['x']
                        ref_y = ref_subgraph.nodes[ref_idx]['y']
                        det_x = detected_graph.nodes[det_idx]['x']
                        det_y = detected_graph.nodes[det_idx]['y']
                        total_distance += np.sqrt((ref_x - det_x)**2 + (ref_y - det_y)**2)

                    avg_distance = total_distance / len(ref_nodes)
                    confidence = max(0, 1.0 - avg_distance/max_distance)

                    # Add all pairs in this pattern match
                    for ref_idx, det_idx in zip(ref_nodes, det_nodes):
                        matches.append({
                            'detected_idx': det_idx,
                            'reference_idx': ref_idx,
                            'distance': avg_distance,
                            'confidence': confidence
                        })

        return matches

    def _pattern_matches(self, ref_subgraph, ref_nodes, detected_graph, det_nodes, max_distance):
        """Check if a detected node pattern matches reference pattern geometry"""
        # Simple geometric matching - check if distances between peaks are similar
        if len(ref_nodes) != len(det_nodes):
            return False

        # Calculate all pairwise distances in both patterns
        ref_distances = []
        det_distances = []

        for i in range(len(ref_nodes)):
            for j in range(i+1, len(ref_nodes)):
                ref_x1 = ref_subgraph.nodes[ref_nodes[i]]['x']
                ref_y1 = ref_subgraph.nodes[ref_nodes[i]]['y']
                ref_x2 = ref_subgraph.nodes[ref_nodes[j]]['x']
                ref_y2 = ref_subgraph.nodes[ref_nodes[j]]['y']
                ref_dist = np.sqrt((ref_x1-ref_x2)**2 + (ref_y1-ref_y2)**2)

                det_x1 = detected_graph.nodes[det_nodes[i]]['x']
                det_y1 = detected_graph.nodes[det_nodes[i]]['y']
                det_x2 = detected_graph.nodes[det_nodes[j]]['x']
                det_y2 = detected_graph.nodes[det_nodes[j]]['y']
                det_dist = np.sqrt((det_x1-det_x2)**2 + (det_y1-det_y2)**2)

                ref_distances.append(ref_dist)
                det_distances.append(det_dist)

        # Check if distance patterns are similar (within tolerance)
        distance_tolerance = 0.3  # ppm
        for ref_dist, det_dist in zip(ref_distances, det_distances):
            if abs(ref_dist - det_dist) > distance_tolerance:
                return False

        return True
