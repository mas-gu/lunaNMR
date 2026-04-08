# Nucleus-Adaptive Search Windows Implementation

**Date:** 2026-01-08
**Implementation:** Tight level nucleus-adaptive search windows for drift control
**Test Status:** ✅ All 9 tests passing

## Summary

Implemented nucleus-adaptive search windows that automatically scale based on nucleus type (15N vs 13C) to provide equivalent drift control relative to typical linewidths.

## Changes Made

### 1. Added Search Window Parameters to PS2DConfig

**File:** `lunaNMR/core/ps2d_config.py`

**15N-HSQC Configuration:**
```python
'search_window_f1': 0.020,  # 15N search window (ppm) - 5% of typical LW
'search_window_f2': 0.005,  # 1H search window (ppm) - 12.5% of typical LW
```

**13C-HSQC Configuration:**
```python
'search_window_f1': 0.005,  # 13C search window (ppm) - 5% of typical LW
'search_window_f2': 0.005,  # 1H search window (ppm) - 12.5% of typical LW
```

**F1 Scaling:** 4.0× (15N is 4× larger than 13C, matching typical linewidth ratio)

### 2. Updated Default Values in Main Window

**File:** `lunaNMR/gui/main_window.py:403-405`

**Before:**
```python
self.search_window_x = 0.01   # 1H dimension (ppm)
self.search_window_y = 0.04   # 15N dimension (ppm)
```

**After:**
```python
self.search_window_x = 0.005  # 1H dimension (ppm) - F2, tight level
self.search_window_y = 0.020  # 15N dimension (ppm) - F1, tight level (default 15N)
```

### 3. Enhanced Nucleus Type Change Handler

**File:** `lunaNMR/gui/main_window.py:3239-3261`

Added automatic search window updates when nucleus type changes:
```python
# Update nucleus-adaptive search windows
self.search_window_y = config.search_window_f1  # F1 (15N/13C)
self.search_window_x = config.search_window_f2  # F2 (1H)

# Update search window spinboxes
if hasattr(self, 'expert_search_x_spin'):
    self.expert_search_x_spin.setValue(self.search_window_x)
if hasattr(self, 'expert_search_y_spin'):
    self.expert_search_y_spin.setValue(self.search_window_y)

# Update integrator search windows
if hasattr(self, 'integrator') and self.integrator:
    self.integrator.set_search_window(self.search_window_x, self.search_window_y)
```

### 4. Created Comprehensive Test Suite

**File:** `tests/test_nucleus_adaptive_search_windows.py`

**Tests (9 total):**
1. ✅ 15N config has search window parameters
2. ✅ 13C config has search window parameters
3. ✅ 15N tight values (0.020/0.005 ppm)
4. ✅ 13C tight values (0.005/0.005 ppm)
5. ✅ F1 nucleus-dependent (4.0× scaling)
6. ✅ F2 nucleus-independent (same for both)
7. ✅ Search windows are ~5% F1 / ~12.5% F2 of linewidth
8. ✅ Total drift calculation (search + fit margin)
9. ✅ Global config switching updates windows

## Drift Control Comparison

### Current vs New Implementation

| Aspect | Old (Not Adaptive) | New (Tight, Adaptive) |
|--------|-------------------|---------------------|
| **15N F1 search** | 0.040 ppm | 0.020 ppm (2× tighter) |
| **13C F1 search** | 0.040 ppm ❌ | 0.005 ppm (8× tighter) ✅ |
| **Both F2 search** | 0.010 ppm | 0.005 ppm (2× tighter) |
| **15N total drift** | 0.060/0.020 ppm | 0.040/0.015 ppm |
| **13C total drift** | 0.050/0.020 ppm | 0.015/0.015 ppm |
| **Nucleus adaptive?** | ❌ NO | ✅ YES |

### Strictness Level: Tight

**Characteristics:**
- F1: 5% of typical linewidth
- F2: 12.5% of typical linewidth
- Stricter drift control than current
- Not overly restrictive
- Balanced for quality control

## Total Drift Calculation

```
Total Max Drift = search_window + fitting_margin
```

**15N-HSQC:**
- F1: 0.020 (search) + 0.020 (fit) = **0.040 ppm total**
- F2: 0.005 (search) + 0.010 (fit) = **0.015 ppm total**

**13C-HSQC:**
- F1: 0.005 (search) + 0.010 (fit) = **0.015 ppm total**
- F2: 0.005 (search) + 0.010 (fit) = **0.015 ppm total**

## User Impact

### Automatic Behavior

1. **On Application Start:** Defaults to 15N tight values (0.020/0.005 ppm)

2. **On Nucleus Switch:**
   - User selects nucleus in Expert Mode → Nucleus Type dropdown
   - Search windows automatically update to match nucleus
   - GUI spinboxes reflect new values
   - Integrator receives updated values

3. **Manual Override:**
   - User can still manually adjust in Expert Mode → Detection Parameters
   - Custom values remain until nucleus type changes again

### How to Use

**Method 1: Automatic (Recommended)**
1. Open Expert Mode
2. Select nucleus type (15N or 13C)
3. Search windows automatically adjust

**Method 2: Manual Override**
1. Open Expert Mode → Detection Parameters
2. Manually adjust "1H/15N (ppm)" spinboxes
3. Values persist until nucleus type changes

## Benefits

1. ✅ **Fixes Inconsistency:** 15N and 13C now have appropriate search windows
2. ✅ **Tighter Drift Control:** 2-8× stricter than previous implementation
3. ✅ **Nucleus-Aware:** Scales with physical linewidth differences
4. ✅ **Quality Control:** Better outlier detection in Independent mode
5. ✅ **Consistent Architecture:** Matches ps2d_config pattern for other parameters
6. ✅ **Fully Tested:** Comprehensive test suite validates behavior

## Alternative Strictness Levels

If "Tight" is too strict or too loose, alternative values are available:

### Default (Moderate)
```python
'15N': {'search_window_f1': 0.040, 'search_window_f2': 0.010}
'13C': {'search_window_f1': 0.010, 'search_window_f2': 0.010}
```

### Very Tight
```python
'15N': {'search_window_f1': 0.010, 'search_window_f2': 0.003}
'13C': {'search_window_f1': 0.0025, 'search_window_f2': 0.003}
```

### Match Cascade Limits
```python
'15N': {'search_window_f1': 0.025, 'search_window_f2': 0.005}
'13C': {'search_window_f1': 0.010, 'search_window_f2': 0.005}
```

To change strictness level, edit values in `ps2d_config.py`.

## Technical Details

### Search Window Application (Detection Stage)

**Location:** `core_integrator.py:3412-3413`

```python
# For each reference peak, find closest detected peak within search window
if x_distance <= search_window_x_ppm and y_distance <= search_window_y_ppm:
    # Match this detected peak to reference
```

**Behavior:**
- Detector finds ALL peaks in spectrum
- For each reference peak, searches within ±search_window
- If match found → uses detected position
- If NO match → creates "dummy peak" at reference position

### Relationship to Cascade Drift Limits

| Constraint | Stage | Independent Mode | Cascade Mode |
|------------|-------|------------------|--------------|
| **Search Window** | Detection | ✅ Applied | ✅ Applied |
| **Per-fit Margin** | Fitting | ✅ Applied | ✅ Applied |
| **Absolute Drift Limit** | Fitting | ❌ Not applied | ✅ Applied |

**Key Insight:** Tight search windows provide drift control in Independent mode that's conceptually similar to Cascade mode's absolute limits, but at the detection stage instead of fitting stage.

## Files Modified

1. `lunaNMR/core/ps2d_config.py` - Added search_window_f1/f2 parameters
2. `lunaNMR/gui/main_window.py` - Updated defaults and nucleus change handler
3. `tests/test_nucleus_adaptive_search_windows.py` - New comprehensive test suite
4. `NUCLEUS_ADAPTIVE_SEARCH_WINDOWS.md` - This documentation

## Verification

```bash
# Run tests
cd lunaNMR_v1o0
QT_QPA_PLATFORM=offscreen python3 -m pytest tests/test_nucleus_adaptive_search_windows.py -v

# Expected output: 9 passed
```

## Future Enhancements

Possible improvements:
1. Add strictness level selector to Expert Mode GUI
2. Add visual indicator showing current drift control level
3. Statistics on detection match rate (how many peaks found vs dummy)
4. Per-spectrum drift reports in Series Integration

---

**Implementation follows TDD:** Test written first, all tests passing ✅
