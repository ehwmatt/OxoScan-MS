import numpy as np


def bin(single_spectra_dict, rt_window, mz_window, combine_function):
    return_dict = {}

    for rt, rt_dict in single_spectra_dict.items():
        for mz, mz_dict in rt_dict.items():
            for ion, ion_val in mz_dict.items():
                binned_rt = (int(rt / rt_window) + 0.5) * rt_window
                binned_mz = (int(mz / mz_window) + 0.5) * mz_window
                return_dict.setdefault(binned_rt, {}).setdefault(
                    binned_mz, {}
                ).setdefault(ion, []).append(ion_val)

    for rt, rt_dict in return_dict.items():
        for mz, mz_dict in rt_dict.items():
            for ion in mz_dict:
                mz_dict[ion] = combine_function(mz_dict[ion])

    return return_dict


def combine(multiple_spectra_dict, combine_function):
    return_dict = {}

    for spectra_name, single_spectra_dict in multiple_spectra_dict.items():
        for rt, rt_dict in single_spectra_dict.items():
            for mz, mz_dict in rt_dict.items():
                for ion, ion_val in mz_dict.items():
                    return_dict.setdefault(rt, {}).setdefault(mz, {}).setdefault(
                        ion, []
                    ).append(ion_val)

    for rt, rt_dict in return_dict.items():
        for mz, mz_dict in rt_dict.items():
            for ion in mz_dict:
                mz_dict[ion] = combine_function(mz_dict[ion])

    return return_dict


def to_matrix(dict_data, ion):
    x_data = set([])
    y_data = set([])
    for rt, rt_dict in dict_data.items():
        x_data.add(rt)
        for mz in rt_dict:
            y_data.add(mz)

    sorted_x_data = sorted(list(x_data))
    sorted_y_data = sorted(list(y_data))
    return_matrix = np.zeros((len(sorted_x_data), len(sorted_y_data)))
    for rt, rt_dict in dict_data.items():
        for mz in rt_dict:
            return_matrix[sorted_x_data.index(rt), sorted_y_data.index(mz)] = rt_dict[
                mz
            ][ion]

    return return_matrix, sorted_x_data, sorted_y_data
