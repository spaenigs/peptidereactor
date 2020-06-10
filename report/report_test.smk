

rule all:
    input:
         expand("roc_{i}.svg", i=["hiv_protease", "ace_vaxinpad"]),
         expand("cd_{i}.svg", i=["hiv_protease", "ace_vaxinpad"]),
         expand("f_{i}.txt", i=["hiv_protease", "ace_vaxinpad"])
    output:
         "benchmark.csv"
    shell:
         "touch {output}"

rule compute_roc:
    input:
         "roc.svg"
    output:
         report("roc_{i}.svg", category="{i}", subcategory="Compute ROC", caption="report/foo.rst")
    shell:
         "cp {input[0]} {output[0]}"

rule critical_difference:
    input:
         "cd.svg"
    output:
         report("cd_{i}.svg", category="{i}", subcategory="Critical difference", caption="report/foo.rst")
    shell:
         "cp {input[0]} {output[0]}"

rule fail:
    output:
          "f_{i}.txt"
    run:
         # try:
         #     if wildcards.i == "hiv_protease":
         #        raise ValueError
         #     else:
         shell("touch {output[0]}")

         # except ValueError:
         #     print("Caught error!")
