## Docker

### Save/load docker image

```shell script
docker save encoding_benchmark > encoding_benchmark.tar
docker load < docker/encoding_benchmark.tar
```

### Override entrypoint
e.g.
```shell script
docker run --entrypoint "ls" encoding_benchmark -l /
docker run --entrypoint "wget" encoding_benchmark http://raptorx.uchicago.edu/
```
or even access it interactively
```shell script
docker run -it --entrypoint "/bin/bash" encoding_benchmark
```

## Pipeline

### Create DAG of meta jobs

```shell script
./apps/run_pipeline -s encoding_benchmark.smk \ 
                    --config dataset=neuropeptides_ds3 \ 
                    --dag | dot -Tpdf > dag.pdf
```

### Split datasets, if necessary

```shell script
./apps/run_pipeline -s create_datasets.smk --config dataset=neuropeptides
```

### Encoding benchmark

```shell script
mkdir data/neuropeptides_ds1/profile # if not existing yet
./apps/run_pipeline -s encoding_benchmark.smk --config dataset=neuropeptides_ds1
```

#### Run pipelines isolated

```shell script
./apps/run_pipeline -s nodes/utils/protein_structure_prediction/Snakefile \
data/neuropeptides_ds3/pdb/UniRef100_A0SIF1.pdb \  # target file
--config dataset="neuropeptides_ds3" \
         fasta_in="data/neuropeptides_ds3/seqs.fasta" \    
         classes_in="data/neuropeptides_ds3/classes.txt" \    
         download_link_in="nodes/utils/protein_structure_prediction/raptorx_download_link.txt" \
         license_key_in="nodes/utils/protein_structure_prediction/modeller_license_key.txt" \ 
         pdbs_out="data/neuropeptides_ds3/pdb/UniRef100_A0SIF1.pdb" \
         token="asd"
```

## Nodes

### Add new node

Example: add a new node to convert 
[pdb]([link](https://en.wikipedia.org/wiki/Protein_Data_Bank_(file_format)) to 
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
    name: encoding_benchmark
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
    ./apps/delete_container
    ```
    Afterwards, rebuild the `encoding_benchmark` docker container with 
    ```shell script
    ./apps/build_container
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
