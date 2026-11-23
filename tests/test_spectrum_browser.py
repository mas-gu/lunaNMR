# ABOUTME: The Spectrum Browser and the Spectrum Viewer it opens, on a real DNAJA1 series: one run
# ABOUTME: in the session and one opened from a project bundle, reached through both entry points.

"""
The browser lists a series' spectra with what the series measured on each: how many of
its peaks stand above the noise, the mean R-squared and S/N of its fits, and a
per-spectrum score computed from those three. The viewer shows one spectrum with its
fitted peaks. Both must read the keys a series result carries today, for a series held
in the session (arrays in every fit) and for one read back from a bundle (stored slim).
"""

import copy
import json
import math
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRA = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o5.ft", "2D_DNAJA1_HSPA8_3o0.ft")
NOT_A_SPECTRUM = "unreadable.ft"
MENU_ACTION = "Browse Individual Spectra..."

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _new_window():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    return LunaNMRMainWindow()


def _close(window):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    window.close()
    window.deleteLater()
    set_ps2d_config('15N')


@pytest.fixture(scope="module")
def session_run(app, tmp_path_factory):
    """The BatchResults of the series dialog's worker over the three spectra and a file
    that is no spectrum, as a session holds them."""
    from lunaNMR.gui.dialogs.series_integration_dialog import ProcessingWorker, SeriesIntegrationDialog
    folder = tmp_path_factory.mktemp("session_series")
    for name in SPECTRA:
        (folder / name).symlink_to(DNAJA1 / name)
    (folder / NOT_A_SPECTRUM).write_bytes(b"not a spectrum")
    window = _new_window()
    window.load_data(str(folder / SPECTRA[0]), str(DNAJA1 / "DNAJA1.txt"))
    dialog = SeriesIntegrationDialog(main_window=window)
    finished = []
    worker = ProcessingWorker([str(folder / name) for name in SPECTRA + (NOT_A_SPECTRUM,)],
                              window.integrator.peak_list.copy(), 'reference',
                              dialog._get_voigt_params(), series_name='session')
    worker.finished.connect(finished.append)
    worker.run()
    dialog.deleteLater()
    _close(window)
    assert len(finished[0].results) == 4
    return finished[0]


@pytest.fixture(scope="module")
def session_window(app, session_run):
    """A main window after the series dialog finished the run."""
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    window = _new_window()
    dialog = SeriesIntegrationDialog(main_window=window)
    dialog.current_series_name = 'session'
    dialog._on_finished(session_run)
    dialog.outcome_box.close()
    dialog.deleteLater()
    yield window
    _close(window)


@pytest.fixture(scope="module")
def bundle_window(app, dnaja1_series_bundle):
    """A main window that opened the project bundle."""
    window = _new_window()
    assert window._get_project_manager().load_project(dnaja1_series_bundle)[0]
    yield window
    _close(window)


@pytest.fixture(params=["run in the session", "opened from a bundle"])
def window(request):
    return request.getfixturevalue(
        {"run in the session": "session_window", "opened from a bundle": "bundle_window"}[request.param])


def _written(window):
    """{spectrum: its rows of the tidy CSV the window's series wrote}, for the spectra
    the series fitted."""
    from lunaNMR.utils.project_manager import series_folder
    folder = Path(series_folder(window.batch_results))
    columns = json.loads((folder / "series_metadata.json").read_text())['columns']
    tidy = pd.read_csv(folder / "series_analysis_tidy.csv",
                       dtype={'spectrum_name': str, 'assignment': str})
    rows = {c['spectrum']: tidy[tidy['spectrum_name'] == str(c['column'])] for c in columns}
    return {name: of_spectrum for name, of_spectrum in rows.items() if len(of_spectrum)}


def _trigger(window, text):
    from PySide6.QtWidgets import QMenu
    next(action for menu in window.menuBar().findChildren(QMenu)
         for action in menu.actions() if action.text() == text).trigger()


def _only_child(parent, kind):
    [child] = parent.findChildren(kind)
    return child


def _discard(dialog):
    dialog.close()
    dialog.setParent(None)
    dialog.deleteLater()


@pytest.fixture
def browser(window):
    """The Spectrum Browser, opened from the main window's menu."""
    from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
    _trigger(window, MENU_ACTION)
    dialog = _only_child(window, SpectrumBrowserDialog)
    yield dialog
    _discard(dialog)


def _listed(browser):
    """{spectrum: {column title: text shown}} of the rows the browser shows."""
    tree = browser.tree
    titles = [tree.headerItem().text(col) for col in range(tree.columnCount())]
    rows = [tree.topLevelItem(i) for i in range(tree.topLevelItemCount())]
    shown = [dict(zip(titles, (row.text(col) for col in range(tree.columnCount()))))
             for row in rows if not row.isHidden()]
    return {row['Spectrum']: row for row in shown}


def _expected_score(rows):
    """The per-spectrum score from a spectrum's tidy rows: detection 40, R-squared 40,
    S/N (full marks at 10) 20."""
    detected = (rows['detection_sigma'] >= 3).mean()
    return 40 * detected + 40 * rows['r_squared'].mean() + 20 * min(rows['snr'].mean() / 10, 1)


class TestWhatTheBrowserLists:
    def test_every_fitted_spectrum_of_the_series_is_listed(self, browser, window):
        written = _written(window)
        assert len(written) == 3
        assert set(written) <= set(_listed(browser))
        assert set(_listed(browser)) == set(window.batch_results.results)

    def test_detection_is_the_share_of_its_peaks_at_three_sigma_or_more(self, browser, window):
        for name, rows in _written(window).items():
            detected = rows['detection_sigma'] >= 3
            assert 0 < detected.sum() < len(rows)
            shown = _listed(browser)[name]
            assert shown['Detected'] == str(detected.sum())
            assert shown['Total'] == str(len(rows))
            assert shown['Detection Rate'] == f"{detected.mean():.1%}"

    def test_mean_r_squared_and_signal_to_noise_are_those_of_its_fits(self, browser, window):
        for name, rows in _written(window).items():
            shown = _listed(browser)[name]
            assert math.isclose(float(shown['Mean R²']), rows['r_squared'].mean(), abs_tol=1e-4)
            assert math.isclose(float(shown['Mean S/N']), rows['snr'].mean(), abs_tol=0.06)

    def test_the_spectrum_score_is_computed_from_those_three(self, browser, window):
        for name, rows in _written(window).items():
            score = float(_listed(browser)[name]['Spectrum score'].split()[0])
            assert math.isclose(score, _expected_score(rows), abs_tol=0.51)

    def test_no_fitted_spectrum_scores_failed(self, browser, window):
        for name in _written(window):
            assert "Failed" not in _listed(browser)[name]['Spectrum score']
            assert _listed(browser)[name]['Status'] == "✓"

    def test_the_header_names_the_series_and_counts_its_spectra(self, browser, window):
        from PySide6.QtWidgets import QLabel
        results = window.batch_results.results
        fitted = sum(1 for result in results.values() if result['status'] == 'success')
        header = (f"Series '{window.batch_results.metadata['series_name']}': "
                  f"{len(results)} spectra, {fitted} successful")
        assert any(label.text().startswith(header) for label in browser.findChildren(QLabel))


    def test_each_spectrum_name_is_shown_in_full(self, app, browser):
        browser.show()
        app.processEvents()
        assert browser.tree.sizeHintForColumn(1) <= browser.tree.columnWidth(1)

    def test_no_group_title_loses_an_ampersand_to_a_keyboard_mnemonic(self, browser):
        from PySide6.QtWidgets import QGroupBox
        for group in browser.findChildren(QGroupBox):
            assert '&' not in group.title().replace('&&', '')


class TestASpectrumThatCouldNotBeRead:
    def test_is_listed_as_failed_with_nothing_measured(self, app, session_window):
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        _trigger(session_window, MENU_ACTION)
        browser = _only_child(session_window, SpectrumBrowserDialog)
        try:
            shown = _listed(browser)[NOT_A_SPECTRUM]
            assert shown['Status'] == "✗"
            assert (shown['Detected'], shown['Total']) == ("0", "0")
            assert shown['Spectrum score'] == "0 (Failed)"
        finally:
            _discard(browser)


class TestFilteringAndSorting:
    def test_the_search_box_keeps_the_spectra_whose_name_holds_the_text(self, browser):
        browser.search_edit.setText("3o0")
        assert list(_listed(browser)) == [SPECTRA[2]]
        browser.search_edit.clear()
        assert set(SPECTRA) <= set(_listed(browser))

    def test_the_status_filter_separates_fitted_from_failed_spectra(self, app, session_window):
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        _trigger(session_window, MENU_ACTION)
        browser = _only_child(session_window, SpectrumBrowserDialog)
        try:
            browser.status_filter.setCurrentText("Failed")
            assert list(_listed(browser)) == [NOT_A_SPECTRUM]
            browser.status_filter.setCurrentText("Success")
            assert set(_listed(browser)) == set(SPECTRA)
        finally:
            _discard(browser)

    def test_the_score_filter_keeps_the_spectra_of_one_band(self, browser):
        bands = {name: row['Spectrum score'].split()[1].strip("()")
                 for name, row in _listed(browser).items()}
        assert len(set(bands.values())) > 1
        for band in set(bands.values()):
            browser.score_filter.setCurrentText(band)
            assert set(_listed(browser)) == {name for name, b in bands.items() if b == band}

    @pytest.mark.parametrize("title", ["Detected", "Detection Rate", "Mean R²", "Mean S/N",
                                       "Spectrum score"])
    def test_a_numeric_column_sorts_by_its_numbers(self, browser, title):
        from PySide6.QtCore import Qt
        tree = browser.tree
        column = [tree.headerItem().text(c) for c in range(tree.columnCount())].index(title)
        tree.sortByColumn(column, Qt.AscendingOrder)
        shown = [float(tree.topLevelItem(i).text(column).split()[0].rstrip('%'))
                 for i in range(tree.topLevelItemCount())]
        assert len(set(shown)) > 1
        assert shown == sorted(shown)


    def test_the_spectrum_column_sorts_by_name(self, browser):
        from PySide6.QtCore import Qt
        tree = browser.tree
        tree.sortByColumn(1, Qt.DescendingOrder)
        shown = [tree.topLevelItem(i).text(1) for i in range(tree.topLevelItemCount())]
        assert shown == sorted(shown, reverse=True)


    def test_refresh_lists_the_series_again_under_the_filters_in_force(self, browser):
        browser.search_edit.setText("0o5")
        _button(browser, "Refresh List").click()
        assert list(_listed(browser)) == [SPECTRA[1]]


class TestWithNoSpectrumSelected:
    @pytest.mark.parametrize("button", ["Open Spectrum Viewer", "Quick Analysis"])
    def test_a_button_asks_for_a_selection_and_opens_nothing(self, browser, monkeypatch, button):
        from PySide6.QtWidgets import QMessageBox
        from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
        said = []
        monkeypatch.setattr(QMessageBox, "warning",
                            staticmethod(lambda parent, title, text: said.append((title, text))))
        monkeypatch.setattr(QMessageBox, "information",
                            staticmethod(lambda parent, title, text: said.append((title, text))))
        _button(browser, button).click()
        assert said == [("No Selection", "Please select a spectrum from the list first.")]
        assert not browser.findChildren(SpectrumViewerDialog)


class TestQuickAnalysis:
    def test_reports_the_values_the_list_shows_for_the_selected_spectrum(self, browser, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        said = []
        monkeypatch.setattr(QMessageBox, "information",
                            staticmethod(lambda parent, title, text: said.append(text)))
        _select(browser, SPECTRA[1])
        _button(browser, "Quick Analysis").click()
        shown = _listed(browser)[SPECTRA[1]]
        [text] = said
        assert SPECTRA[1] in text
        assert f"Detected Peaks: {shown['Detected']}" in text
        assert f"Total Peaks: {shown['Total']}" in text
        assert f"Detection Rate: {shown['Detection Rate']}" in text
        assert f"Mean R²: {shown['Mean R²']}" in text
        assert f"Mean S/N: {shown['Mean S/N']}" in text
        assert f"Spectrum score: {shown['Spectrum score']}" in text


def _select(browser, spectrum):
    """Click the row of `spectrum` in the browser's list."""
    tree = browser.tree
    row = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
               if tree.topLevelItem(i).text(1) == spectrum)
    tree.setCurrentItem(row)
    tree.itemClicked.emit(row, 1)


def _button(dialog, text):
    from PySide6.QtWidgets import QPushButton
    return next(b for b in dialog.findChildren(QPushButton) if b.text() == text)


@pytest.fixture
def viewer(browser):
    """The Spectrum Viewer the browser opens for the series' second spectrum."""
    from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
    _select(browser, SPECTRA[1])
    _button(browser, "Open Spectrum Viewer").click()
    return _only_child(browser, SpectrumViewerDialog)


class TestOpeningASpectrum:
    def test_the_viewer_loads_the_spectrum_and_lists_its_fitted_peaks(self, viewer, window):
        rows = _written(window)[SPECTRA[1]]
        assert viewer.spectrum_name == SPECTRA[1]
        assert viewer.integrator is not None and viewer.integrator.nmr_data is not None
        assert viewer.peak_list.count() == len(rows)

    def test_a_double_click_on_a_row_opens_it_too(self, browser):
        from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
        tree = browser.tree
        row = next(tree.topLevelItem(i) for i in range(tree.topLevelItemCount())
                   if tree.topLevelItem(i).text(1) == SPECTRA[2])
        tree.itemDoubleClicked.emit(row, 1)
        assert _only_child(browser, SpectrumViewerDialog).spectrum_name == SPECTRA[2]


def _menu_action(window, text):
    from PySide6.QtWidgets import QMenu
    return next(action for menu in window.menuBar().findChildren(QMenu)
                for action in menu.actions() if action.text() == text)


class TestTheMainWindowsMenuEntry:
    def test_is_shown_and_enabled_beside_its_sibling(self, window):
        from PySide6.QtWidgets import QMenu
        menu = next(m for m in window.menuBar().findChildren(QMenu)
                    if MENU_ACTION in [a.text() for a in m.actions()])
        shown = [a.text() for a in menu.actions() if a.isVisible() and a.isEnabled()]
        assert MENU_ACTION in shown and "Browse Series Results..." in shown

    def test_without_a_series_it_says_so_and_opens_nothing(self, app, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        said = []
        monkeypatch.setattr(QMessageBox, "warning",
                            staticmethod(lambda parent, title, text: said.append((title, text))))
        empty = _new_window()
        try:
            assert _menu_action(empty, MENU_ACTION).isVisible()
            _menu_action(empty, MENU_ACTION).trigger()
            assert said == [("No Data", "No series results available. Run series integration first.")]
            assert not empty.findChildren(SpectrumBrowserDialog)
        finally:
            _close(empty)


# ---------- the Spectrum Viewer ----------

PEAK_ENTRY = r"^(?P<assignment>.+): \((?P<x>[-\d.]+), (?P<y>[-\d.]+)\) R²=(?P<r2>[\d.]+)$"


def _peak_entries(viewer):
    import re
    entries = [re.match(PEAK_ENTRY, viewer.peak_list.item(i).text())
               for i in range(viewer.peak_list.count())]
    assert all(entries)
    return [(m['assignment'], float(m['x']), float(m['y']), float(m['r2'])) for m in entries]


def _markers(axes):
    """[(x, y, colour)] of the peak markers drawn on the contour plot."""
    from matplotlib.collections import PathCollection
    from matplotlib.colors import to_hex
    return [(float(c.get_offsets()[0][0]), float(c.get_offsets()[0][1]), to_hex(c.get_facecolor()[0]))
            for c in axes.collections if type(c) is PathCollection]


def _index_of(window, assignment):
    rows = _written(window)[SPECTRA[1]].reset_index(drop=True)
    return int(rows.index[rows['assignment'] == assignment][0]), rows


def _clustered_fit(window):
    """The assignment of a well-fitted peak of the second spectrum fitted with a neighbour."""
    rows = window.batch_results.results[SPECTRA[1]]['integration_results']
    return next(row['assignment'] for row in rows
                if len(row.get('all_peaks') or []) > 1 and row['quality'] == 'Excellent')


def _click_peak(app, viewer, index):
    """Click entry `index` of the viewer's peak list with the mouse."""
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest
    viewer.show()
    item = viewer.peak_list.item(index)
    viewer.peak_list.scrollToItem(item)
    app.processEvents()
    QTest.mouseClick(viewer.peak_list.viewport(), Qt.LeftButton,
                     pos=viewer.peak_list.visualItemRect(item).center())


def _peak_info(viewer):
    """The selected-peak panel: its heading, and {name: value} of its lines."""
    lines = viewer.peak_info_label.text().splitlines()
    return lines[0], dict(line.split(": ", 1) for line in lines if ": " in line)


def _vertical_lines(axes):
    return [float(line.get_xdata()[0]) for line in axes.lines
            if len(line.get_xdata()) == 2 and line.get_xdata()[0] == line.get_xdata()[1]]


class TestTheViewersPeaks:
    def test_the_title_bar_carries_the_values_the_browser_lists(self, browser, viewer):
        from PySide6.QtWidgets import QLabel
        shown = _listed(browser)[SPECTRA[1]]
        title = f"Status: success | Detection: {shown['Detection Rate']} | Mean R²: {shown['Mean R²']}"
        assert title in [label.text() for label in viewer.findChildren(QLabel)]

    def test_each_entry_names_a_peak_at_its_fitted_position(self, viewer, window):
        written = _written(window)[SPECTRA[1]]
        entries = _peak_entries(viewer)
        assert len(entries) == len(written)
        for (assignment, x, y, r2), (_, row) in zip(entries, written.iterrows()):
            assert assignment == row['assignment']
            assert math.isclose(x, row['ppm_x'], abs_tol=1e-3)
            assert math.isclose(y, row['ppm_y'], abs_tol=1e-2)
            assert math.isclose(r2, row['r_squared'], abs_tol=1e-3)

    def test_the_fitted_peaks_are_marked_on_the_spectrum_where_they_were_fitted(self, viewer, window):
        written = _written(window)[SPECTRA[1]]
        markers = _markers(viewer.contour_widget.axes)
        assert len(markers) == len(written)
        for (x, y, _colour), (_, row) in zip(markers, written.iterrows()):
            assert math.isclose(x, row['ppm_x'], abs_tol=1e-4)
            assert math.isclose(y, row['ppm_y'], abs_tol=1e-4)
        ppm_x = viewer.integrator.ppm_x_axis
        low, high = sorted(viewer.contour_widget.axes.get_xlim())
        assert low >= ppm_x.min() - 0.5 and high <= ppm_x.max() + 0.5

    def test_clicking_a_peak_selects_it_and_shows_what_the_series_measured(self, app, viewer, window):
        index, rows = _index_of(window, '42GlyH')
        row = rows.loc[index]
        _click_peak(app, viewer, index)
        assert viewer.selected_peak_idx == index
        heading, info = _peak_info(viewer)
        assert heading == '42GlyH'
        assert math.isclose(float(info['Position F2 (¹H)'].split()[0]), row['ppm_x'], abs_tol=1e-4)
        assert math.isclose(float(info['Position F1'].split()[0]), row['ppm_y'], abs_tol=1e-2)
        assert math.isclose(float(info['Height']), row['height'], rel_tol=1e-2)
        assert math.isclose(float(info['Volume']), row['volume'], rel_tol=1e-2)
        assert math.isclose(float(info['R²']), row['r_squared'], abs_tol=1e-4)
        assert math.isclose(float(info['S/N']), row['snr'], abs_tol=1e-2)
        assert info['Quality'] == row['quality']

    def test_the_selected_peak_is_the_one_marked_in_red(self, app, viewer, window):
        index, rows = _index_of(window, '42GlyH')
        _click_peak(app, viewer, index)
        reds = [(x, y) for x, y, colour in _markers(viewer.contour_widget.axes) if colour == '#ff0000']
        assert len(reds) == 1
        assert math.isclose(reds[0][0], rows.loc[index, 'ppm_x'], abs_tol=1e-4)
        assert math.isclose(reds[0][1], rows.loc[index, 'ppm_y'], abs_tol=1e-4)

    def test_the_cross_sections_are_cut_through_the_selected_peak(self, app, viewer, window):
        index, rows = _index_of(window, '42GlyH')
        _click_peak(app, viewer, index)
        [along_1h] = _vertical_lines(viewer.f2_widget.axes)
        [along_indirect] = _vertical_lines(viewer.f1_widget.axes)
        assert math.isclose(along_1h, rows.loc[index, 'ppm_x'], abs_tol=1e-4)
        assert math.isclose(along_indirect, rows.loc[index, 'ppm_y'], abs_tol=1e-4)

    def test_next_and_previous_step_through_the_peaks(self, viewer, window):
        rows = _written(window)[SPECTRA[1]].reset_index(drop=True)
        _button(viewer, "Next ▶").click()
        assert (viewer.selected_peak_idx, _peak_info(viewer)[0]) == (0, rows.loc[0, 'assignment'])
        _button(viewer, "Next ▶").click()
        assert (viewer.selected_peak_idx, _peak_info(viewer)[0]) == (1, rows.loc[1, 'assignment'])
        assert viewer.peak_list.currentRow() == 1
        _button(viewer, "◀ Previous").click()
        assert (viewer.selected_peak_idx, _peak_info(viewer)[0]) == (0, rows.loc[0, 'assignment'])


class TestTheViewersFitDetail:
    def test_a_selected_peak_s_fit_is_drawn_from_its_region_and_surface(self, app, viewer, window):
        assignment = _clustered_fit(window)
        index, _rows = _index_of(window, assignment)
        _click_peak(app, viewer, index)
        drawn = viewer.voigt_3d_plotter.current_result
        assert drawn['assignment'] == assignment
        assert drawn['fitted_2d_surface'].shape == drawn['region_2d']['intensity'].shape
        titles = [ax.get_title() for ax in viewer.voigt_2d_plotter.figure.axes]
        assert f"Experimental - {assignment}" in titles

    def test_the_2d_view_states_the_r_squared_and_label_the_fit_carries(self, app, viewer, window):
        index, rows = _index_of(window, '42GlyH')
        _click_peak(app, viewer, index)
        fit = viewer.fitted_peaks[index]
        figure = viewer.voigt_2d_plotter.figure
        shown = "\n".join(text.get_text() for ax in figure.axes for text in ax.texts)
        assert f"Quality: {rows.loc[index, 'quality']}" in shown
        assert f"R² = {fit['r_squared']:.4f}" in shown

    def test_a_colour_scheme_redraws_the_3d_view_in_its_colours(self, app, viewer, window):
        from matplotlib.colors import to_rgb
        index, _rows = _index_of(window, _clustered_fit(window))
        _click_peak(app, viewer, index)
        plotter = viewer.voigt_3d_plotter

        def fit_colour():
            [fit] = [c for c in plotter.figure.axes[0].collections if c.get_label() == 'Total Fit']
            return tuple(fit.get_color()[0][:3])

        shown = viewer.colormap_3d_combo.currentText()
        assert fit_colour() == to_rgb(plotter.color_presets[shown]['total_fit']['color'])
        drawn = plotter.current_result
        viewer.colormap_3d_combo.setCurrentText('Dark')
        assert plotter.current_result is drawn
        assert fit_colour() == to_rgb(plotter.color_presets['Dark']['total_fit']['color'])

    def test_every_layer_box_changes_what_the_3d_view_draws(self, app, viewer, window):
        from PySide6.QtWidgets import QCheckBox, QGroupBox
        index, _rows = _index_of(window, _clustered_fit(window))
        _click_peak(app, viewer, index)

        def drawn():
            canvas = viewer.voigt_3d_plotter.canvas
            canvas.draw()
            return bytes(canvas.buffer_rgba())

        group = next(g for g in viewer.findChildren(QGroupBox) if g.title() == "Layer Visibility")
        boxes = {box.text(): box for box in group.findChildren(QCheckBox)}
        boxes["Individual"].setChecked(True)   # what Limit Peak limits is each peak's own surface
        without_effect = []
        for box in boxes.values():
            before = drawn()
            box.click()
            if drawn() == before:
                without_effect.append(box.text())
            box.click()
        assert not without_effect


class TestTheViewersContourSettings:
    @staticmethod
    def _contours(viewer):
        _button(viewer, "Refresh Plot").click()
        return [(tuple(c.levels), c.get_cmap().name)
                for c in viewer.contour_widget.axes.collections if hasattr(c, 'levels')]

    def test_every_setting_changes_the_contours_a_refresh_draws(self, viewer):
        from PySide6.QtWidgets import QAbstractSpinBox, QComboBox, QGroupBox, QLabel
        group = next(g for g in viewer.findChildren(QGroupBox) if g.title() == "Contour Settings")
        names = [label.text() for label in group.findChildren(QLabel)]
        controls = group.findChildren(QAbstractSpinBox) + group.findChildren(QComboBox)
        assert controls and len(names) == len(controls)
        without_effect = []
        for control in controls:
            before = self._contours(viewer)
            assert before
            if isinstance(control, QComboBox):
                control.setCurrentIndex((control.currentIndex() + 1) % control.count())
            elif control.value() == control.maximum():
                control.stepDown()
            else:
                control.stepUp()
            if self._contours(viewer) == before:
                without_effect.append(control)
        assert not without_effect


class TestExportingTheSpectrumPlot:
    def test_writes_the_contour_plot_to_the_chosen_file(self, viewer, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        target = tmp_path / "plot.png"
        monkeypatch.setattr(QFileDialog, "getSaveFileName",
                            staticmethod(lambda *a, **k: (str(target), "")))
        monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
        _button(viewer, "Export PNG").click()
        assert target.read_bytes().startswith(b"\x89PNG")


class TestASeriesWhoseSpectraHaveMoved:
    def test_the_viewer_still_lists_the_peaks_and_shows_what_was_measured(self, app, window, tmp_path):
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
        moved = copy.copy(window.batch_results)
        moved.metadata = dict(moved.metadata, data_folder=str(tmp_path))
        browser = SpectrumBrowserDialog(batch_results=moved)
        try:
            _select(browser, SPECTRA[1])
            _button(browser, "Open Spectrum Viewer").click()
            viewer = _only_child(browser, SpectrumViewerDialog)
            index, rows = _index_of(window, '42GlyH')
            assert viewer.integrator is None
            assert viewer.peak_list.count() == len(rows)
            _click_peak(app, viewer, index)
            assert _peak_info(viewer)[0] == '42GlyH'
        finally:
            _discard(browser)


# ---------- the series dialog's entry ----------

class TestTheSeriesDialogsEntry:
    @pytest.fixture
    def dialog(self, app, session_run):
        import types
        from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
        main_window = types.SimpleNamespace(current_nmr_folder=None, saved_series={})
        dialog = SeriesIntegrationDialog(main_window=main_window)
        yield dialog
        if dialog.outcome_box is not None:
            dialog.outcome_box.close()
        dialog.deleteLater()

    def test_the_button_is_in_the_dialog_and_waits_for_a_run(self, dialog):
        button = _button(dialog, "📋 Browse Spectra")
        assert button.isVisibleTo(dialog)
        assert not button.isEnabled()

    def test_after_a_run_it_opens_the_browser_on_that_run(self, dialog, session_run):
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        dialog._on_finished(session_run)
        button = _button(dialog, "📋 Browse Spectra")
        assert button.isEnabled()
        button.click()
        browser = _only_child(dialog, SpectrumBrowserDialog)
        assert set(_listed(browser)) == set(session_run.results)
        assert "Failed" not in _listed(browser)[SPECTRA[0]]['Spectrum score']
