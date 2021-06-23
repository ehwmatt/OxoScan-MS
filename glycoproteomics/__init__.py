import csv


def read_spectra_file(filename):
    dict_data = {}
    ion_labels = []
    with open(filename, mode="r") as infile:
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
    return dict_data
