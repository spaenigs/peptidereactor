from snakemake import shell

db = snakemake.wildcards.database
move_files = True

if "Remain" in db:
    target = db.replace("Remain", "BC100")
    target_dir = f"peptidereactor/RaptorX/databases/{target}/"
    shell(f"mkdir -p {target_dir}")
elif "nr" in db:
    target_dir = "peptidereactor/RaptorX/databases/NR_new/"
    shell("mkdir -p peptidereactor/RaptorX/databases/NR_new")
else:
    move_files = False
    print(f"Got {db}. Nothing to do.")

if move_files:
    shell(f"""for file in `ls -1 peptidereactor/RaptorX/databases/{db}/`; do
                               mv peptidereactor/RaptorX/databases/{db}/$file {target_dir};
                         done""")

shell("touch {snakemake.output}")
