## iFeature

See https://github.com/Superzchen/iFeature for reference. 

In order to resolve some conflicts during runtime, this version contains some enhancements:

1) Adapted paths to the `data` folder in 
    - `codes/Geary.py`
    - `codes/Moran.py`
    - `codes/NMBroto.py` 
2) `codes/TA.py` as well as `codes/ASA.py` check for the right file type, i.e., `*.spXout`.
3) The `setup.py` script to facilitate a local installation with pip.
4) A folder called `iFeature.egg-info`, which will be added by the pip installation.
5) `data/AAindex.tsv`: a reformatted version of `data/AAindex.txt`, which allows to be read 
by `pandas.read_csv()`

All changes can be verified by issuing the following commands:

```shell script
cd ~;
git clone https://github.com/Superzchen/iFeature.git
diff -r --exclude README.md iFeature/ /path/to/eb/eb/apps/iFeature/
``` 

See `../apps/Dockerfile` and `../apps/environment.yaml` for more information.  