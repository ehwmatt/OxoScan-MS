import numpy as np
from fastdtw import fastdtw


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


def to_matrix(single_spectra_dict, ion_list):
    x_data = set([])
    y_data = set([])
    for rt, rt_dict in single_spectra_dict.items():
        x_data.add(rt)
        for mz in rt_dict:
            y_data.add(mz)

    sorted_x_data = sorted(list(x_data))
    sorted_y_data = sorted(list(y_data))
    return_matrix = np.zeros((len(sorted_x_data), len(sorted_y_data)))
    for rt, rt_dict in single_spectra_dict.items():
        for mz in rt_dict:
            return_matrix[sorted_x_data.index(rt), sorted_y_data.index(mz)] = np.sum(
                [rt_dict[mz][ion] for ion in ion_list]
            )

    return return_matrix, sorted_x_data, sorted_y_data


def list_ions(single_spectra_dict):
    for rt, rt_dict in single_spectra_dict.items():
        for mz, mz_dict in rt_dict.items():
            return sorted(list(set(mz_dict.keys())))


def get_rt_array(matrix_spectra_dict, ions, x_data, y_data, common_y):
    return [
        np.array(
            [
                matrix_spectra_dict[ion][0][idx, [y_data.index(y) for y in common_y]]
                for ion in ions
            ]
        )
        for idx in range(len(x_data))
    ]


def align_rt(single_spectra_dict, single_spectra_ref_dict, warp_factor, normalise):
    """Puts the first spectra on the same RT axis as the reference spectra using
    Dynamic Time Warping. It is recommended to use this on binned data, as the
    MZ values have to be the same for the comparisons to be made. This can
    result in data columns (corresponding to a particular RT value) being
    duplicated (if a sample column is matched to multiple reference columns) or
    removed (if a reference column is mapped to multiple sample columns then
    only the one with the minimal distance is used).

    Outputs:
    1) Aligned spectra
    2) RT Mapping from reference RT values to sample RT values
    """
    ions = list_ions(single_spectra_dict)
    spectra_dict = {ion: to_matrix(single_spectra_dict, [ion]) for ion in ions}
    ref_dict = {ion: to_matrix(single_spectra_ref_dict, [ion]) for ion in ions}

    _discard, x_data, y_data = list(spectra_dict.values())[0]
    _discard, x_ref, y_ref = list(ref_dict.values())[0]

    common_y = list(set(y_data).intersection(set(y_ref)))
    spectra_array = get_rt_array(spectra_dict, ions, x_data, y_data, common_y)
    ref_array = get_rt_array(ref_dict, ions, x_ref, y_ref, common_y)

    dist_func = lambda a, b: np.sum(abs(a - b))

    if normalise:
        distance, path = fastdtw(
            (spectra_array - np.mean(spectra_array)) / np.std(spectra_array),
            (ref_array - np.mean(ref_array)) / np.std(ref_array),
            warp_factor,
            dist_func,
        )
    else:
        distance, path = fastdtw(spectra_array, ref_array, warp_factor, dist_func)

    path_dict = {}
    for link in path:
        spec_idx, ref_idx = link
        path_dict.setdefault(ref_idx, set([]))
        path_dict[ref_idx].add(spec_idx)

    optimal_links = {}
    for ref_idx, spec_idxs in path_dict.items():
        if len(spec_idxs) == 1:
            optimal_links[ref_idx] = list(spec_idxs)[0]
        else:
            dists = [
                dist_func(spectra_array[spec_idx], ref_array[ref_idx])
                for spec_idx in spec_idxs
            ]
            optimal_links[ref_idx] = list(spec_idxs)[dists.index(min(dists))]

    return_dict = {}
    rt_mapping = {}
    for ref_idx, spec_idx in optimal_links.items():
        return_dict[x_ref[ref_idx]] = single_spectra_dict[x_data[spec_idx]]
        rt_mapping[x_ref[ref_idx]] = x_data[spec_idx]

    return return_dict, rt_mapping
