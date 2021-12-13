from . import persistence


def find(ion_matrix, x_label, y_label, N, x_radius, y_radius):
    peaks = persistence.imagepers.persistence(ion_matrix)
    return_peaks = []
    for peak in peaks:
        if len(return_peaks) == N:
            return return_peaks
        coord, height, persis, parent = peak
        x, y = coord
        # Check to see if the peak is too close to an already called peak
        too_close = False
        for alt_peak in return_peaks:
            alt_x, alt_y = alt_peak[0]
            if (x_label[x] - alt_x) ** 2.0 / x_radius ** 2 + (
                y_label[y] - alt_y
            ) ** 2.0 / y_radius ** 2 <= 1:
                too_close = True
        if not too_close:
            return_peaks.append(((x_label[x], y_label[y]), height, persis))
    return return_peaks


def integrate(
    ion_matrix, x_label, y_label, peaks, x_radius, y_radius, integration_function
):
    return_list = []
    for peak in peaks:
        peak_x, peak_y = peak[0]
        peak_array = []
        for i, x in enumerate(x_label):
            for j, y in enumerate(y_label):
                if (x - peak_x) ** 2 / x_radius ** 2 + (
                    y - peak_y
                ) ** 2 / y_radius ** 2 <= 1.0:
                    peak_array.append(ion_matrix[i][j])
        return_list.append(integration_function(peak_array))
    return return_list


def rt_move(peaks, rt_alignment_dict):
    return [
        ((rt_alignment_dict[peak[0][0]], peak[0][1]), peak[1], peak[2])
        for peak in peaks
    ]


def assert_no_overlap(x_label, y_label, peaks, x_radius, y_radius):
    all_bins = []
    for peak in peaks:
        peak_x, peak_y = peak[0]
        for i, x in enumerate(x_label):
            for j, y in enumerate(y_label):
                if (x - peak_x) ** 2 / x_radius ** 2 + (
                    y - peak_y
                ) ** 2 / y_radius ** 2 <= 1.0:
                    all_bins.append((i, j))
    return len(all_bins) == len(list(set(all_bins)))
