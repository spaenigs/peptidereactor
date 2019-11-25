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