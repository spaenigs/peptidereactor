from glob import glob

import subprocess
import os


def exec_cmd(cmd):
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, cwd="/home/spaenigs/")
    return process.communicate()

# dir = "/home/spaenigs/PycharmProjects/"  # home-office PC
dir = "/media/spaenigs/4B1DB7375F3291A1/"  # office PC
# dir = "/home/ubuntu/"  # de.NBI VM


all = True

for user, ip in [
    ["spaenigs", "137.248.121.201"],
    ["ubuntu", "172.16.103.179"],
    ["ubuntu", "172.16.103.199"],
    ["ubuntu", "172.16.103.203"]
]:

    for p in glob("data/*"):
        if "temp" in p:
            continue
        elif all:
            print(f"Downloading {p} (complete)...")
            bash_command = f"ssh {user}@{ip} " \
                           f"test -e /home/{user}/peptidereactor/{p}/benchmark/ && " \
                           f"echo $?"
            bm_dir_exists, _ = exec_cmd(bash_command)
            bash_command = f"ssh {user}@{ip} " \
                           f"test -e /home/{user}/peptidereactor/{p}/csv/ && " \
                           f"echo $?"
            csv_dir_exists, _ = exec_cmd(bash_command)
            if bm_dir_exists and csv_dir_exists:
                bash_command = f"rm -r {dir}peptidereactor/{p}"
                exec_cmd(bash_command)
                bash_command = f"scp -r {user}@{ip}:/home/{user}/peptidereactor/{p}/  " \
                               f"{dir}/peptidereactor/data/"
                exec_cmd(bash_command)
        else:
            bash_command = f"ssh {user}@{ip} " \
                           f"test -e /home/{user}/peptidereactor/{p}/benchmark/dataset_correlation.csv && " \
                           f"echo $?"
            output, _ = exec_cmd(bash_command)
            # only download files, if present remote and not already downloaded
            if output.decode('ascii') and not os.path.exists(f"{p}/benchmark/dataset_correlation.csv"):
                print(f"Downloading {p}...")
                # benchmark data
                exec_cmd(f"mkdir -p {dir}/peptidereactor/{p}/misc/benchmark/")
                bash_command = f"scp {user}@{ip}:/home/{user}/peptidereactor/{p}/misc/benchmark/*  " \
                               f"{dir}peptidereactor/{p}/misc/benchmark/"
                exec_cmd(bash_command)
                bash_command = f"scp -r {user}@{ip}:/home/{user}/peptidereactor/{p}/benchmark/  " \
                               f"{dir}/peptidereactor/{p}"
                exec_cmd(bash_command)
