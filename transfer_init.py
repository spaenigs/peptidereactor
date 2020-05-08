from shutil import copyfile

import sys
import re

old_path = sys.argv[1]  # "nodes/encodings/aac/__init__.py"

with open(old_path) as f:
    _main = ""
    _params = ""
    add = False
    for l in f.readlines():
        if "def rule" in l:
            hits = re.findall("def rule\((.*)\)", l)
            _params += hits[0]
        if "input:" in l:
            add = True
        if "run:" in l:
            add = False
        if add:
            _main += l.replace("            ", "")
    # print(_main)

rule_dir = old_path.split("/")[1]
rule_name = old_path.split("/")[2]

if rule_dir == "encodings":
    rule_dir = rule_dir[:-1]

scaffold = f"""import secrets


def _get_header(token):
    return f'''
rule {rule_dir}_{rule_name}_{{token}}:'''


def _get_benchmark(benchmark_out):
    return f'''
    benchmark:
        "{{benchmark_out}}"'''


def _get_main({_params}):
    return f'''
{_main.rstrip()}
    run:
        with WorkflowExecuter(dict(input), dict(output), params.configfile, cores=CORES) as e:
            shell(f\"\"\"{{{{e.snakemake}}}} -s {{{{params.snakefile}}}} --configfile {{{{params.configfile}}}}\"\"\")
'''


def rule({_params}, benchmark_dir=None):
    token = secrets.token_hex(4)
    rule = _get_header(token)
    if benchmark_dir is not None:
        benchmark_out = f"{{benchmark_dir}}{rule_dir}_{rule_name}_{{token}}.txt"
        rule += _get_benchmark(benchmark_out)
    rule += _get_main({_params})
    return rule
"""

copyfile(old_path, old_path.replace(".py", ".bckp.py"))

with open(old_path, "w") as f:
    f.write(scaffold)
    f.flush()