def bin(dict_data, rt_window, mz_window, combine_function):
    return_dict = {}

    for rt, rt_dict in dict_data.items():
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
