from proteinreactor.workflow_executer import WorkflowExecuter
import joblib


# bachem_profiles.smk
# bachem_encode.smk
# bachem_train.smk
# bachem_predict.smk

# 1) take input sequence and split
# 2) unpack model and get respective encodings
# 3) if necessary, compute profile
# 4) encoded splitted sequences
# 5) predict in model

rule get_encodings:
    input:
         "data/bachem/models/best_model_0.joblib"
    output:
         "data/temp/bachem/encodings.txt"
    run:
          d = joblib.load(str(input))

rule utils_sliding_windows:
    input:
         series_in=f"data/bachem/series_predict.yaml"
    output:
         fastas_out=expand(f"data/bachem_window_length_{{window_length}}/seqs.fasta", window_length=WINDOW_LENGTHS),
         classes_out=expand(f"data/bachem_window_length_{{window_length}}/classes.txt", window_length=WINDOW_LENGTHS),
    params:
         snakefile="nodes/utils/sliding_windows/sliding_windows.smk",
         configfile="nodes/utils/sliding_windows/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="bachem") as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule utils_sliding_windows_complete:
#     input:
#          series_in=f"data/bachem/series.yaml"
#     output:
#          fastas_out=expand(f"data/bachem_window_length_{{window_length}}_complete/seqs.fasta", window_length=WINDOW_LENGTHS),
#          classes_out=expand(f"data/bachem_window_length_{{window_length}}_complete/classes.yaml", window_length=WINDOW_LENGTHS),
#          classes_idx_out=expand(f"data/bachem_window_length_{{window_length}}_complete/classes.txt", window_length=WINDOW_LENGTHS),
#     params:
#          snakefile="nodes/utils/sliding_windows/sliding_windows_complete.smk",
#          configfile="nodes/utils/sliding_windows/config.yaml"
#     run:
#         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES, dataset="bachem") as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")