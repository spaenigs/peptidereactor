### Create database

```
docker run -v $PWD:/snakemake \  
           -v /home/ubuntu/database/uniref90/:/snakemake/db_tmp \
           encoding_benchmark -s build_refseq.smk --config target_dir=db_tmp/  
```

### Save/load docker image

```
docker save encoding_benchmark > encoding_benchmark.tar
docker load < docker/encoding_benchmark.tar
```