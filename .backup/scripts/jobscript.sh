#!/usr/bin/env bash

# properties = {properties}

source ~/.bashrc
conda activate encoding_benchmark

. /etc/profile.d/modules.sh
module load tools/R-3.4.3

{exec_job}
