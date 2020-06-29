from snakemake import shell

db = snakemake.wildcards.database
if db in ["nr70", "nr90"]:
    target_dir = f"peptidereactor/RaptorX/databases/" + db
else:
    target_dir = f"peptidereactor/RaptorX/databases/"
shell(f"""tar -zxf {{input}} -C {target_dir};
           touch {snakemake.output[0]}""")