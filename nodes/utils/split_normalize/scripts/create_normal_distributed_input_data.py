import joblib as jl
import numpy as np
import re


input_data = jl.load(snakemake.input[0])

total_sequences_len = list(map(lambda tup: len(tup[0][1]), zip(*input_data)))

# compute global median
global_median = np.median(total_sequences_len)

# compute global median absolute deviation
global_mad = np.median(np.abs(np.array(total_sequences_len) - np.median(total_sequences_len)))/2

# remove outliers from original dataset
classes_median_new = list(filter(lambda x: global_median - global_mad < x < global_median + global_mad, total_sequences_len))

# compute median to seperate original data into 2 distinct datasets
global_median_new = np.median(classes_median_new)

# save updated datasets (outliers removed)
input_data_filtered = list(filter(lambda tup: global_median - global_mad < len(tup[0][1]) < global_median + global_mad, zip(*input_data)))
seq_tups_filtered = list(map(lambda tup: tup[0], input_data_filtered))
classes_filtered = list(map(lambda tup: tup[1], input_data_filtered))
input_data_global_median = (seq_tups_filtered, classes_filtered)
jl.dump(input_data_global_median, snakemake.output[0])

# based on dataset \wo outliers: get sequences < median
input_data_ds1 = list(filter(lambda tup: len(tup[0][1]) < global_median_new, zip(*input_data_global_median)))
seq_tups_ds1 = list(map(lambda tup: [re.sub("\W", "", tup[0][0]), tup[0][1]], input_data_ds1))  # clean sequence names
classes_ds1 = list(map(lambda tup: tup[1], input_data_ds1))
jl.dump((seq_tups_ds1, classes_ds1), snakemake.output[1])

# based on dataset \wo outliers: get sequences >= median
input_data_ds2 = list(filter(lambda tup: len(tup[0][1]) >= global_median_new, zip(*input_data_global_median)))
seq_tups_ds2 = list(map(lambda tup: [re.sub("\W", "", tup[0][0]), tup[0][1]], input_data_ds2))  # clean sequence names
classes_ds2 = list(map(lambda tup: tup[1], input_data_ds2))
jl.dump((seq_tups_ds2, classes_ds2), snakemake.output[2])
