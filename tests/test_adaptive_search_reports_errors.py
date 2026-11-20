# ABOUTME: A fit that raises inside the adaptive radius search is skipped, and said so once.
# ABOUTME: Drives the real optimiser with a fit function that fails for one peak.

from lunaNMR.core.adaptive_optimizer import AdaptiveOptimizer


def _peaks(count):
    return [{'x_ppm': 8.0 + 0.01 * i, 'y_ppm': 120.0 + 0.1 * i, 'assignment': str(i)}
            for i in range(count)]


def _fit(peak, radF1, radF2):
    if peak['assignment'] == '3':
        raise ValueError('singular matrix')
    return {'r_squared': 0.95, 'pos_f2': peak['x_ppm'], 'pos_f1': peak['y_ppm']}


def test_a_fit_that_raises_is_reported_once_not_hidden(capsys):
    optimizer = AdaptiveOptimizer(noise_level=1.0)
    result = optimizer.optimize(_peaks(AdaptiveOptimizer.MIN_PEAKS_FOR_OPTIMIZATION + 4),
                                0.2, 0.02, _fit)
    assert result['success']
    printed = capsys.readouterr().out
    assert printed.count('raised and were skipped') == 1
    assert 'singular matrix' in printed


def test_a_search_with_no_error_says_nothing(capsys):
    optimizer = AdaptiveOptimizer(noise_level=1.0)
    quiet = lambda peak, radF1, radF2: {'r_squared': 0.95, 'pos_f2': peak['x_ppm'],
                                        'pos_f1': peak['y_ppm']}
    optimizer.optimize(_peaks(AdaptiveOptimizer.MIN_PEAKS_FOR_OPTIMIZATION + 4),
                       0.2, 0.02, quiet)
    assert 'raised and were skipped' not in capsys.readouterr().out
