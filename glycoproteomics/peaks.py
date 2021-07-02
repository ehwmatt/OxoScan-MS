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
