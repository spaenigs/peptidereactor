import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from sklearn.externals import joblib as jl


input_data = jl.load(snakemake.input[0])
input_data_global_median = jl.load(snakemake.input[1])
input_data_ds1 = jl.load(snakemake.input[2])
input_data_ds2 = jl.load(snakemake.input[3])


def filter_by_classes(name, input_data_):
    neg_classes = list(map(
        lambda tup: len(tup[0][1]),
        filter(lambda tup: tup[1] == 0, zip(*input_data_))
    ))
    pos_classes = list(map(
        lambda tup: len(tup[0][1]),
        filter(lambda tup: tup[1] == 1, zip(*input_data_))
    ))
    print(f"For {name}, there are {str(len(neg_classes))} negative classes and {str(len(pos_classes))} positive classes. Total dataset size: {str(len(neg_classes) + len(pos_classes))}")
    return [neg_classes, pos_classes]


# create grid plot
fig = plt.figure(1)
gridspec.GridSpec(3, 4)

totalClasses = ["0 (negative)", "1 (positive)"]

raw = filter_by_classes("RAW", input_data)

# adapted from https://scientificallysound.org/2016/06/09/matplotlib-how-to-plot-subplots-of-unequal-sizes/
plt.subplot2grid((3, 4), (0, 0), colspan=2, rowspan=1)
plt.hist(raw, align="left", label=totalClasses, bins=100)
plt.xlabel("Sequence length")
plt.ylabel("Count")
plt.title("No preprocessing")

plt.subplot2grid((3, 4), (0, 2), colspan=2, rowspan=1)
plt.boxplot(raw, labels=totalClasses, flierprops={"marker": ".", "markersize": 0.5})
plt.xlabel("Class")
plt.ylabel("Sequence length")

wo_outliers = filter_by_classes("WO-OUTLIERS", input_data_global_median)

plt.subplot2grid((3, 4), (1, 0), colspan=2, rowspan=1)
plt.hist(wo_outliers, align="left", label=totalClasses, bins=100)
plt.xlabel("Sequence length")
plt.ylabel("Count")
plt.title("Median absolute deviation")

plt.subplot2grid((3, 4), (1, 2), colspan=2, rowspan=1)
plt.boxplot(wo_outliers, labels=totalClasses,
                   flierprops={"marker": ".", "markersize": 0.5})
plt.xlabel("Class")
plt.ylabel("Sequence length")
plt.title("asd", color="white")

ds1 = filter_by_classes("DS1", input_data_ds1)

plt.subplot2grid((3, 4), (2, 0), colspan=1, rowspan=1)
plt.hist(ds1, align="left", label=totalClasses, bins=100)
plt.xlabel("Sequence length")
plt.ylabel("Count")
plt.title("Dataset 1")

plt.subplot2grid((3, 4), (2, 1), colspan=1, rowspan=1)
bp = plt.boxplot(ds1, labels=totalClasses,
                   flierprops={"marker": ".", "markersize": 0.5})
plt.xlabel("Class")
plt.xticks(rotation=30)
plt.ylabel("Sequence length")

ds2 = filter_by_classes("DS2", input_data_ds2)

plt.subplot2grid((3, 4), (2, 2), colspan=1, rowspan=1)
plt.hist(ds2, align="left", label=totalClasses, bins=100)
plt.xlabel("Sequence length")
plt.ylabel("Count")
plt.title("Dataset 2")

plt.subplot2grid((3, 4), (2, 3), colspan=1, rowspan=1)
plt.boxplot(ds2, labels=totalClasses,
                   flierprops={"marker": ".", "markersize": 0.5})
plt.xlabel("Class")
plt.xticks(rotation=30)
plt.ylabel("Sequence length")

# set global plot attributes
fig.suptitle(snakemake.wildcards.dataset.upper())
fig.legend(labels=totalClasses, prop={'size': 10})
fig.tight_layout()
fig.set_size_inches(w=11,h=7)

# save figure
plt.savefig(snakemake.output[0])










# # compute global median
# global_median = np.median(neg_classes + pos_classes)
#
# # compute global median absolute deviation
# global_mad = np.median(np.abs(np.array(neg_classes + pos_classes) - np.median(neg_classes + pos_classes)))/2
#
# # remove outliers from original dataset
# neg_classes_median_new = list(filter(lambda x: global_median - global_mad < x < global_median + global_mad, neg_classes))
# pos_classes_median_new = list(filter(lambda x: global_median - global_mad < x < global_median + global_mad, pos_classes))
#
# # compute median to seperate original data into 2 distinct datasets
# global_median_new = np.median(neg_classes_median_new + pos_classes_median_new)
#
# # split negative data
# neg_classes_median_new = np.array(neg_classes_median_new)
# neg_classes_ds1 = neg_classes_median_new[neg_classes_median_new < global_median_new]
# neg_classes_ds2 = neg_classes_median_new[neg_classes_median_new >= global_median_new]
#
# # split positive data
# pos_classes_median_new = np.array(pos_classes_median_new)
# pos_classes_ds1 = pos_classes_median_new[pos_classes_median_new < global_median_new]
# pos_classes_ds2 = pos_classes_median_new[pos_classes_median_new >= global_median_new]




