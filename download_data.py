from glob import glob

import subprocess
import os


def exec_cmd(cmd):
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, cwd="/home/spaenigs/")
    return process.communicate()

# dir = "/home/spaenigs/PycharmProjects/"
dir = "/media/spaenigs/4B1DB7375F3291A1/"


i = 179

for p in glob("data/*"):
    if "temp" in p:
        continue
    else:
        bash_command = f"ssh ubuntu@172.16.103.{i} " \
                       f"test -e /home/ubuntu/peptidereactor/{p}/benchmark/dataset_correlation.csv && " \
                       f"echo $?"
        output, _ = exec_cmd(bash_command)
        # only download files, if present remote and not already downloaded
        if output.decode('ascii') and not os.path.exists(f"{p}/benchmark/dataset_correlation.csv"):
            print(f"Downloading {p}...")
            # benchmark data
            exec_cmd(f"mkdir -p {dir}/peptidereactor/{p}/misc/benchmark/")
            bash_command = f"scp ubuntu@172.16.103.{i}:/home/ubuntu/peptidereactor/{p}/misc/benchmark/*  " \
                           f"{dir}peptidereactor/{p}/misc/benchmark/"
            exec_cmd(bash_command)
            bash_command = f"scp -r ubuntu@172.16.103.{i}:/home/ubuntu/peptidereactor/{p}/benchmark/  " \
                           f"{dir}/peptidereactor/{p}"
            exec_cmd(bash_command)

