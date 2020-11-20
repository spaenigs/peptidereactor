from glob import glob

ds = "hiv_v3"

with open("jobs.txt") as f:
    jobs = sorted(l.rstrip()[:-9] for l in f.readlines())

files = [p.split("/")[-1][:-13]
         for p in sorted(glob(f"data/{ds}/misc/benchmark/*.txt"))]

res = {}

for job in jobs:
    # if job in res and job in files:
    #     res[job] += [f for f in files if f == job]
    if job not in res and job in files:
        res[job] = (len([j for j in jobs if j == job]), [f for f in files if f == job])

for k, v in res.items():
    if len(v[1]) > 1:
        print(f"### {k}")
        print(v)
        print()

print(f"--> remaining")
print(set(files).difference(jobs))


