#!/usr/bin/env python

import sys
from snakemake.shell import shell
from snakemake.io import expand

DRY_RUN = sys.argv[1]

DATASETS = ["neuropeptides"]
CORES = 8

for ds in DATASETS:
    shell(f"./apps/run_pipeline -s 01_create_datasets.smk --config dataset={ds} cores={CORES} --quiet {DRY_RUN}")

NORMALIZED_DATASETS = expand("{dataset}_{id}", dataset=DATASETS, id=["ds1", "ds2"])

for ds in NORMALIZED_DATASETS:
    shell(f"./apps/run_pipeline -s 02_create_profiles.smk --config dataset={ds} cores={CORES} --quiet {DRY_RUN}")

for ds in NORMALIZED_DATASETS:
    shell(f"./apps/run_pipeline -s 03_compute_params.smk --config dataset={ds} cores={CORES} --quiet {DRY_RUN}")

for ds in NORMALIZED_DATASETS:
    shell(f"./apps/run_pipeline -s 04_create_encodings.smk --config dataset={ds} cores={CORES} --quiet {DRY_RUN}")
