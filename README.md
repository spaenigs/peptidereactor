# PEPTIDEREACToR

## Overview

![image info](docs/images/peptidereactor.svg)

## Installation

1. Clone this repo: `git clone git@github.com:spaenigs/proteinreactor.git`
2. Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). 
3. Create conda environment: `conda env create --file apps/environment.yaml`
4. Install docker: 
    - Ubuntu: `./apps/install_docker_io`
    - Other distros: `./apps/install_docker_ce` 
5. Build image: `./apps/build_image`
6. Set required download links. See `nodes/utils/protein_structure_prediction/README.md` and
`nodes/utils/secondary_structure_profile/README.md` for further instructions.
7. Run the pipeline via `./apps/run_pipeline -s eb.smk --config cores=4 --quiet`

## Docker

### Save/load docker image

```shell script
docker save peptidereactor > peptidereactor.tar
docker load < docker/peptidereactor.tar
```

### Override entrypoint
e.g.
```shell script
docker run --entrypoint "ls" peptidereactor -l /
docker run --entrypoint "wget" peptidereactor http://raptorx.uchicago.edu/
```
or even access it interactively
```shell script
docker run -it --entrypoint "/bin/bash" peptidereactor
```

## Pipeline

### Create DAG of meta jobs

```shell script
./peptidereactor/run_pipeline -s peptidereactor.smk \ 
                    --config dataset=neuropeptides_ds3 \ 
                    --dag | dot -Tpdf > dag.pdf
```

### Split datasets, if necessary

```shell script
./peptidereactor/run_pipeline -s create_datasets.smk --config dataset=neuropeptides
```

### Encoding benchmark

Create the profile dir, if not existing yet.
```shell script
mkdir data/neuropeptides_ds1/profile 
```
Make sure, that the parameter-based encodings are set (or set them manually). See, e.g., 
`apps/iFeature/codes/ksctriad.py` for details.
```shell script
./peptidereactor/run_pipeline -s maximum_window_length.smk --config dataset=neuropeptides_ds1
```
Moreover, ngram-based encodings require a predetermined dimension. Either run the respective rule 
or set the maximum dimension manually:
```shell script
./peptidereactor/run_pipeline -s maximum_dim_size.smk --config dataset=neuropeptides_ds1
```
Finally, execute the pipeline:
```shell script
./peptidereactor/run_pipeline -s peptidereactor.smk --config dataset=neuropeptides_ds1
```

#### Run pipelines isolated

```shell script
./peptidereactor/run_pipeline -s nodes/utils/protein_structure_prediction/Snakefile \
data/neuropeptides_ds3/pdb/UniRef100_A0SIF1.pdb \  # target file
--config dataset="neuropeptides_ds3" \
         fasta_in="data/neuropeptides_ds3/seqs.fasta" \    
         classes_in="data/neuropeptides_ds3/classes.txt" \    
         download_link_in="nodes/utils/protein_structure_prediction/raptorx_download_link.txt" \
         license_key_in="nodes/utils/protein_structure_prediction/modeller_license_key.txt" \ 
         pdbs_out="data/neuropeptides_ds3/pdb/UniRef100_A0SIF1.pdb" \
         token="asd"
```

## Meta-workflow

### Add new meta-workflow

The following example demonstrates how to add a new meta-workflow, i.e., a workflow, which incorporates
sub-workflows or nodes, respectively.

1) Create the workflow file: `touch meta_workflow.smk`
2) Paste the following content into it:
    ```snakemake
    import os
    
    config["global_workdir"] = os.getcwd() + "/"  
    
    DATASET = config["dataset"]
    
    rule all:
        input:
             f"data/{DATASET}/csv/aaindex/aaindex_ANDN920101.csv",
             f"data/{DATASET}/csv/aac.csv"

    rule encoding_aaindex:
        input:
             fasta_in=f"data/{DATASET}/seqs.fasta",
             classes_in=f"data/{DATASET}/classes.txt"
        output:
             csv_out=f"data/{DATASET}/csv/aaindex/aaindex_ANDN920101.csv"
        params:
             subworkflow="aaindex",
             snakefile="nodes/encodings/aaindex/Snakefile",
             configfile="nodes/encodings/aaindex/config.yaml"
        script:
             "utils/subworkflow.py"
   
   rule encoding_aac:
        input:
             fasta_in=f"data/{DATASET}/seqs.fasta",
             classes_in=f"data/{DATASET}/classes.txt"
        output:
             csv_out=f"data/{DATASET}/csv/aac.csv"
        params:
             subworkflow="aac",
             snakefile="nodes/encodings/aac/Snakefile",
             configfile="nodes/encodings/aac/config.yaml"
        script:
             "utils/subworkflow.py"
    ```

3) In case one would like to process, e.g., `data/neuropeptides/csv/aaindex/aaindex_ANDN920101.csv` 
in a subsequent sub-workflow, __make sure to use the same file name as input__.
4) Run the meta-workflow as follows: `./apps/run_pipeline -s meta_workflow.smk --config dataset=neuropeptides`

## Nodes

### Add new node

Example: add a new node, i.e., sub-workflow, to convert 
[pdb](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format))  to 
[sdf](https://en.wikipedia.org/wiki/Chemical_table_file) files and 
find their respective, energy-minimized conformation.

1) `mkdir nodes/utils/convert_to_sdf_and_minimize`
2) `touch nodes/utils/convert_to_sdf_and_minimize/Snakefile`
3) Specify input and output via config dictionary, e.g., `config["pdbs_in"]` 
and `config["sdfs_out"]`.
4) Copy/paste into the `Snakefile` and adapt stub:

    ```snakemake
    TOKEN = config["token"]
   
    rule process_input:
       input:
            config["pdbs_in"]  # use the input of the meta rule!
       output:
            f"data/temp/{TOKEN}/intermediate.results"
       run:
            pass
   
   rule generate_output:
       input:
            f"data/temp/{TOKEN}/intermediate.results"
       output: 
            config["sdfs_out"]  # use the output of the meta rule!
       shell:
            """
            for file in {output}; do
               touch $file
            done
            """
    ```
5) Add all necessary external dependencies to `apps/environment.yaml`, in this case `openbabel`:
    ```yaml
    name: peptidereactor
    channels:
      ...
      - openbabel
    dependencies:
      ...
      - openbabel
      - pip:
        ...
    ```
    First, if necessary, remove redundant containers with 
    ```shell script
    ./peptidereactor/delete_container
    ```
    Afterwards, rebuild the `peptidereactor` docker container with 
    ```shell script
    ./peptidereactor/build_container
    ``` 
    
6) Implement the algorithm in the `Snakefile` (see #4) and call it in the meta workflow as follows:

   ```snakemake
    from modlamp.core import read_fasta
    
    rule util_convert_to_sdf_and_minimize:
   
        input:
             pdbs_in=expand("data/neuropeptides_ds2/{seq_name}.pdb",
                            seq_name=read_fasta(f"data/neuropeptides_ds2/seqs.fasta")[1])       
        output:
             sdfs_out=expand("data/neuropeptides_ds2/{seq_name}.pdb",
                            seq_name=read_fasta(f"data/neuropeptides_ds2/seqs.fasta")[1])
        params:
             subworkflow="convert_to_sdf_and_minimize",
             snakefile="nodes/utils/convert_to_sdf_and_minimize/Snakefile",
             configfile="nodes/utils/convert_to_sdf_and_minimize/config.yaml"  
        resources:
             cores=-1  # can be omitted, uses one core per default
        script:
             "utils/subworkflow.py" 
   ```
   
7) Refer to directory `nodes/*/*/` for all available nodes (and examples).
8) Make sure to use the output (`config["sdfs_out"]`) as input for another node 
(e.g., `config["sdfs_in"]`) to enable snakemake capabilities for the meta workflow!
