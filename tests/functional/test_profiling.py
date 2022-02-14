import glycoproteomics
import cProfile
import os
import numpy as np

def test_peak_integration_profiling(rootdir, tmp_path):
    def main():
        top_N_peaks = 100

        rt_x_bin_size = 0.05
        mz_y_bin_size = 2.0

        x_radius_quant = rt_x_bin_size * 3.0
        y_radius_quant = mz_y_bin_size * 5.0

        x_radius_exclude = x_radius_quant * 3.0
        y_radius_exclude = y_radius_quant * 3.0

        spectrum = glycoproteomics.io.read_spectrum_file(
            os.path.join(rootdir, "tests", "data", "spectrum.txt.gz")
        )
        binned_spectrum = glycoproteomics.spectrum.bin(
            spectrum,
            rt_x_bin_size,
            mz_y_bin_size,
            np.mean
        )
        ions = glycoproteomics.spectrum.list_ions(binned_spectrum)
        ions_matrix, x_label, y_label = glycoproteomics.spectrum.to_matrix(
            binned_spectrum,
            ions
        )
        peaks = glycoproteomics.peaks.find(
            ions_matrix,
            x_label,
            y_label,
            top_N_peaks,
            x_radius_exclude,
            y_radius_exclude
        )
        peak_indicies = glycoproteomics.peaks.convert_peaks_to_indicies(
            x_label, y_label, peaks, x_radius_quant, y_radius_quant
        )
        with open(os.path.join(tmp_path, "peak_integration.tsv"), "w") as out_f:
            out_f.write("\t".join(
                [
                    "ion",
                    "peak_num",
                    "merged_rt",
                    "merged_mz",
                    "merged_height",
                    "persistence",
                    "value"
                ]) + "\n"
            )
            for ion in ions:
                ion_matrix, x_label, y_label = (
                    glycoproteomics.spectrum.to_matrix(binned_spectrum, [ion])
                )
                peak_values = glycoproteomics.peaks.integrate(
                    ion_matrix,
                    peak_indicies,
                    np.sum
                )
                for peak_idx in range(len(peak_values)):
                    out_f.write("\t".join(
                        [
                            str(ion),
                            str(peak_idx + 1),
                            str(peaks[peak_idx][0][0]),
                            str(peaks[peak_idx][0][1]),
                            str(peaks[peak_idx][1]),
                            str(peaks[peak_idx][2]),
                            str(peak_values[peak_idx])
                        ]) + "\n"
                    )

    pr = cProfile.Profile()
    pr.enable()
    main()
    pr.disable()
    pr.dump_stats(os.path.join(rootdir, "prof/peak_integration.out"))

    test_lines = open(
        os.path.join(tmp_path, "peak_integration.tsv")
    ).readlines()
    ref_lines = open(
        os.path.join(rootdir, "tests", "data", "output", "peak_integration.tsv")
    ).readlines()

    assert test_lines == ref_lines
