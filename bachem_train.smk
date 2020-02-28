from proteinreactor.workflow_executer import WorkflowExecuter

CORES = config["cores"]

# TODO split protease/bachem into separate workflows

rule all:
    input:
         # protease
         # "data/protease/machine_learning/ensemble_results_validated.csv",
         # "data/protease/machine_learning/ensemble_results_validated_tuned_hp.csv",
         # expand("data/protease/plots/best_model_cv_{tuning}.png", tuning=[0, 1]),
         # expand("data/protease/models/best_model_{tuning}.joblib", tuning=[0, 1]),
         # bachem
         # "data/bachem/machine_learning/ensemble_results_validated.csv",
         # "data/bachem/machine_learning/ensemble_results_validated_tuned_hp.csv",
         # expand("data/bachem/plots/best_model_cv_{tuning}.png", tuning=[0, 1]),
         # expand("data/bachem/models/best_model_{tuning}.joblib", tuning=[0, 1]),
         "data/bachem/machine_learning/top_encodings.csv",
         "data/bachem/machine_learning/phi_correlation/"

rule machine_learning_hold_out_datasets_bachem:
    input:
         csv_dir_in="data/bachem/csv/non_empty/filtered/total/"
         # csv_dir_in=f"data/bachem/csv/non_empty/all/"
    output:
         csv_train_dir_out=directory(f"data/bachem/csv/non_empty_train/"),
         csv_val_1_dir_out=directory(f"data/bachem/csv/non_empty_val_1/"),
         csv_val_2_dir_out=directory(f"data/bachem/csv/non_empty_val_2/"),
         csv_test_dir_out=directory(f"data/bachem/csv/non_empty_test/")
    params:
         snakefile="nodes/machine_learning/hold_out_datasets/Snakefile",
         configfile="nodes/machine_learning/hold_out_datasets/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule machine_learning_hold_out_datasets_protease:
#     input:
#          csv_dir_in=f"data/protease/csv/non_empty/all/"
#     output:
#          csv_train_dir_out=directory(f"data/protease/csv/non_empty_train/"),
#          csv_val_1_dir_out=directory(f"data/protease/csv/non_empty_val_1/"),
#          csv_val_2_dir_out=directory(f"data/protease/csv/non_empty_val_2/"),
#          csv_test_dir_out=directory(f"data/protease/csv/non_empty_test/")
#     params:
#          snakefile="nodes/machine_learning/hold_out_datasets/Snakefile",
#          configfile="nodes/machine_learning/hold_out_datasets/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

rule machine_learning_top_encodings_bachem:
    input:
         train_dir_in=f"data/bachem/csv/non_empty_train/",
         val_dir_in=f"data/bachem/csv/non_empty_val_1/"
    output:
         top_encodings_out="data/bachem/machine_learning/top_encodings.csv",
         phi_correlation_out=directory("data/bachem/machine_learning/phi_correlation/"),
    params:
         snakefile="nodes/machine_learning/top_encodings/Snakefile",
         configfile="nodes/machine_learning/top_encodings/config.yaml"
    run:
         with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
             shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule machine_learning_top_encodings_protease:
#     input:
#          train_dir_in=f"data/protease/csv/non_empty_train/",
#          val_dir_in=f"data/protease/csv/non_empty_val_1/"
#     output:
#          top_encodings_out="data/protease/machine_learning/top_encodings.csv",
#          phi_correlation_out=f"data/protease/machine_learning/phi_correlation.csv"
#     params:
#          snakefile="nodes/machine_learning/top_encodings/Snakefile",
#          configfile="nodes/machine_learning/top_encodings/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule machine_learning_best_ensemble_bachem:
#     input:
#          train_dirs_in=[f"data/bachem/csv/non_empty_train/", f"data/bachem/csv/non_empty_val_1/"],
#          val_dir_in=f"data/bachem/csv/non_empty_val_2/",
#          test_dir_in=f"data/bachem/csv/non_empty_test/",
#          phi_correlation_in=f"data/bachem/machine_learning/phi_correlation.csv"
#     output:
#          ensemble_validation_out=\
#              "data/bachem/machine_learning/ensemble_results_validated.csv",
#          ensemble_validation_tuned_hp_out=\
#              "data/bachem/machine_learning/ensemble_results_validated_tuned_hp.csv",
#          plot_cv_out=\
#               expand("data/bachem/plots/best_model_cv_{tuning}.png", tuning=[0, 1]),
#          ensemble_out=\
#               expand("data/bachem/models/best_model_{tuning}.joblib", tuning=[0, 1]),
#     params:
#          snakefile="nodes/machine_learning/best_ensemble/Snakefile",
#          configfile="nodes/machine_learning/best_ensemble/config.yaml"
#     run:
#          # TODO handle case if no well-erforming encodings have been found
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")

# rule machine_learning_best_ensemble_protease:
#     input:
#          train_dirs_in=[f"data/protease/csv/non_empty_train/", f"data/protease/csv/non_empty_val_1/"],
#          val_dir_in=f"data/protease/csv/non_empty_val_2/",
#          test_dir_in=f"data/protease/csv/non_empty_test/",
#          phi_correlation_in=f"data/protease/machine_learning/phi_correlation.csv"
#     output:
#          ensemble_validation_out=\
#              "data/protease/machine_learning/ensemble_results_validated.csv",
#          ensemble_validation_tuned_hp_out=\
#              "data/protease/machine_learning/ensemble_results_validated_tuned_hp.csv",
#          plot_cv_out=\
#               expand("data/protease/plots/best_model_cv_{tuning}.png", tuning=[0, 1]),
#          ensemble_out=\
#               expand("data/protease/models/best_model_{tuning}.joblib", tuning=[0, 1]),
#     params:
#          snakefile="nodes/machine_learning/best_ensemble/Snakefile",
#          configfile="nodes/machine_learning/best_ensemble/config.yaml"
#     run:
#          with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
#              shell(f"""{e.snakemake} -s {{params.snakefile}} --configfile {{params.configfile}}""")





