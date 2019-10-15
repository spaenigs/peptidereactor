### Create database

```shell script
docker run -v $PWD:/snakemake \  
           -v /home/ubuntu/database/uniref90/:/snakemake/db_tmp \
           encoding_benchmark -s build_refseq.smk --config target_dir=db_tmp/  
```

### Save/load docker image

```shell script
docker save encoding_benchmark > encoding_benchmark.tar
docker load < docker/encoding_benchmark.tar
```

### Create DAG of meta jobs

```shell script
./apps/run_pipeline -s encoding_benchmark.smk --nolock \ 
                    --config dataset=neuropeptides_ds3 \ 
                    --dag | dot -Tpdf > dag.pdf
```