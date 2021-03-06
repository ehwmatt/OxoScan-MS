import csv
import os
import gzip


def read_spectrum_file(file_path):
    dict_data = {}
    ion_labels = []
    if file_path.endswith(".gz"):
        infile = gzip.open(file_path, "rt")
    else:
        infile = open(file_path)
    reader = csv.reader(infile, delimiter="\t")
    header = next(reader)
    for i in range(3, len(header)):
        ion_labels.append(header[i])
    for row in reader:
        rt = float(row[0])
        mz = (float(row[1]) + float(row[2])) / 2
        ion_dict = {}
        for i in range(3, len(row)):
            ion_dict[ion_labels[i - 3]] = float(row[i])
        dict_data.setdefault(rt, {})[mz] = ion_dict
    infile.close()
    return dict_data


def read_spectra_directory(folder_path):
    return_dict = {}
    for file_path in os.listdir(folder_path):
        return_dict[file_path] = read_spectrum_file(
            os.path.join(folder_path, file_path)
        )
    return return_dict
