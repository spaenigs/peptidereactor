import scripts.utils as utils

# localrules: compute_distance_matrix,
#             collect_distance_matrix,
#             run_clustering,
#             compute_geometric_median,
#             collect_geometric_median,
#             plot_clustering,
#             get_final_datasets

rule compute_distance_matrix:
    # """
    # Computes the distances from each dataset to another within an encoding family.
    # Every execution of this rule computes a row in the final distance matrix.
    #
    # input
    # -----
    #     A (normalized,) encoded dataset.
    #
    # output
    # ------
    #     A file with the respective distances to the other datasets. To save computation time,
    #     only the lower triangular matrix will be computed, hence other values are set to NA.
    # """
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" + \
        "{dataset}_{part}_{type}.csv",
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/normalized/" +
        "{dataset}_{part}_normalized-{normalized}.txt"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
        "{dataset}_{part}_{type}_normalized-{normalized}_vs_rest.csv"
    script:
        "scripts/compute_distance_matrix.py"


rule collect_distance_matrix:
    # """
    # Combines all rows to the final distance matrix.
    #
    # input
    # -----
    #     The files with the respective distances.
    #
    # output
    # ------
    #     The lower left triangular distance matrix.
    # """
    input:
        lambda wildcards: \
            expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/correlation/" + \
                   "{dataset}_{part}_{type}_normalized-{normalized}_vs_rest.csv",
                   dataset=wildcards.dataset,
                   part=wildcards.part,
                   encoding=wildcards.encoding,
                   normalized=wildcards.normalized,
                   type=utils.get_type(wildcards.encoding, config))
    output:
          "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
          "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    script:
        "scripts/collect_distance_matrix.py"


rule run_clustering:
    # """
    # Runs the t-SNE algorithm for a given distance matrix and reduces the dimension
    # from NxN to Nx2 components.
    #
    # input
    # -----
    #     The lower left triangular distance matrix.
    #
    # output
    # ------
    #     A file with the reduced components and meaningful subclasses for a given
    #     encoding family. The subclass is based on a regex pattern and extracted
    #     from the filenames of the encoded datasets. See scripts/utils.py for more
    #     information.
    # """
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/" + \
        "{dataset}_{part}_normalized-{normalized}_distance_matrix.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}.csv"
    run:
        from scripts.run_clustering import run_clustering
        run_clustering(utils.ENCODING_PATTERN[wildcards.encoding], # TODO replace with get_unique_types
                       str(input),
                       str(output))


rule compute_geometric_median:
    # """
    # Computes the geometric median, i.e, the median of a given set of 2D vectors
    # for each subclass.
    #
    # input
    # -----
    #     The file with the Nx2 components, i.e, a set of 2D vectors for each subclass.
    #
    # output
    # ------
    #     A file for each subclass and the results from the distances between all points
    #     within a subclass. The geometric median for a subclass is highlighted via the
    #     column 'min_gm'. Note that each of the 2D vectors can be traced back to an
    #     actual, encoded dataset.
    # """
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}.csv"
    output:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
        "{dataset}_{part}_normalized-{normalized,yes|no}_{type}_vs_rest.csv"
    script:
        "scripts/compute_geometric_median.py"


rule collect_geometric_median:
    # """
    # Collect the geometric median data for all subclasses.
    #
    # input
    # -----
    #     All files with the geometric median of a subclass.
    #
    # output:
    #     A file with all subclasses concatenated.
    # """
    input:
        lambda wildcards: \
            expand("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/geom_median/" + \
                   "{dataset}_{part}_normalized-{normalized}_{type}_vs_rest.csv",
                   dataset=wildcards.dataset,
                   part=wildcards.part,
                   normalized=wildcards.normalized,
                   encoding=wildcards.encoding,
                   type=utils.get_unique_types(wildcards.encoding))
    output:
         "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
         "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    script:
        "scripts/collect_geometric_median.py"


rule plot_clustering:
    # """
    # Plots the Nx2 components of an encoding family.
    #
    # input
    # -----
    #     The file with all subclasses concatenated.
    #
    # output
    # ------
    #     A scatter plot. All subclasses are color-coded and the respective
    #     geometric median is annotated.
    # """
    input:
        "00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/tsne/" + \
        "{dataset}_{part}_normalized-{normalized}_geometric_median.csv"
    output:
        "00_data/out/{dataset}/plots/{dataset}_{part}_{encoding}_normalized-{normalized}_tsne.svg"
    script:
        "scripts/plot_clustering.py"


rule get_final_datasets:
    # """
    # Gets the final dataset, i.e., the most representative dataset for an encoding,
    # in case the number of subclasses is 1, or for the subclass if there are more
    # than one subclasses for an encoding.
    #
    # input
    # -----
    #     The file with all subclasses concatenated.
    #
    # output
    # ------
    #     A file with a list of the final datasets. The actual datasets are silently
    #     copied to same output directory.
    # """
    input:
        lambda wildcards: utils.determine_input(wildcards, config)
    output:
        temp("00_data/out/{dataset}/{dataset}_{part}/encodings/{encoding}/csv/final/" +
             "geom_median/tsne/normalized-{normalized}/final_datasets.txt")
    script:
        "scripts/get_final_datasets.py"
