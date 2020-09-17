from glob import glob

import subprocess
import os

i = 203

for p in glob("data/*"):
    if "temp" in p:
        continue
    else:
        bash_command = f"ssh ubuntu@172.16.103.{i} " \
                       f"test -e /home/ubuntu/peptidereactor/{p}/benchmark/dataset_correlation.csv && " \
                       f"echo $?"
        process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd="/home/spaenigs/")
        output, _ = process.communicate()
        # only download files, if present remote and not already downloaded
        if output.decode('ascii') and not os.path.exists(f"{p}/benchmark/dataset_correlation.csv"):
            print(f"Downloading {p}...")
            bash_command = f"scp -r ubuntu@172.16.103.{i}:/home/ubuntu/peptidereactor/{p}/benchmark/  " \
                           f"/home/spaenigs/PycharmProjects/peptidereactor/{p}"
            process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd="/home/spaenigs/")
            output, error = process.communicate()
