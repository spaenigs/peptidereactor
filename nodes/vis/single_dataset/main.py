import yaml
from overview import *
from metrics import *
from roc_pr_curve import *
from similarity import *
from scatter_div import *
from crit_diff import *
from elapsed_time import *
from par_coord import *
from network import *
from dendro import *

import pandas as pd

two_charts_template = """
<!DOCTYPE html>
<html>
<head>
  <script src="https://cdn.jsdelivr.net/npm/vega@{vega_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@{vegalite_version}"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@{vegaembed_version}"></script>
</head>
<body>

<p align = "left">
    <img src="/home/spaenigs/Downloads/peptide-reactor/peptide-reactor/peptide-reactor.005.png" alt="logo" width="7%" height="7%"/>
    <hr>
</p>

<center>

<p><b>{dataset}</b></p>

</center> 

<p><u>General overview</u></p>
<center><div id="vis1"></div></center>

<p><u>Metrics</u></p>
<center><div id="vis2"></div></center>

<p><u>ROC</u></p>
<center><div id="vis3"></div></center>

<p><u>Similarity</u></p>
<center><div id="vis4"></div></center>
<center><div id="vis5"></div></center>

<p><u>Critical difference</u></p>
<center><div id="vis6"></div></center>
<center><div id="vis9"></div></center>

<p><u>Misc</u></p>
<center><div id="vis7"></div></center>
<center><div id="vis8"></div></center>

<script type="text/javascript">
    
    vegaEmbed('#vis9', {spec9}).catch(console.error);
</script>
</body>
</html>

"""


"""
    vegaEmbed('#vis1', {spec1}).catch(console.error);    
    vegaEmbed('#vis2', {spec2}).catch(console.error);
    vegaEmbed('#vis3', {spec3}).catch(console.error);
    vegaEmbed('#vis4', {spec4}).catch(console.error);
    vegaEmbed('#vis5', {spec5}).catch(console.error);
    vegaEmbed('#vis6', {spec6}).catch(console.error);  
    vegaEmbed('#vis7', {spec7}).catch(console.error);
    vegaEmbed('#vis8', {spec8}).catch(console.error);
"""

def main():

    alt.data_transformers.disable_max_rows()

    dataset = "hiv_ddi"

    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # df_mcc = pd.read_csv(f"data/{dataset}/benchmark/metrics/mcc.csv", index_col=0)
    # oc = overview_chart(df_f1, df_mcc).to_json(indent=None)
    #
    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # df_mcc = pd.read_csv(f"data/{dataset}/benchmark/metrics/mcc.csv", index_col=0)
    # mc = metrics_chart(df_f1, df_mcc).to_json(indent=None)
    #
    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # rc = roc_chart(df_f1, dataset).to_json(indent=None)
    #
    # df_div = pd.read_csv(f"data/{dataset}/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
    # df_phi = pd.read_csv(f"data/{dataset}/benchmark/similarity/seq_vs_str/phi.csv", index_col=0)
    # sc = similarity_chart(df_div, df_phi, dataset).to_json(indent=None)
    #
    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # df_div = pd.read_csv(f"data/{dataset}/benchmark/similarity/seq_vs_str/diversity.csv", index_col=0)
    # sdc = scatter_div_chart(df_f1, df_div, dataset).to_json(indent=None)
    #
    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # df_cd = pd.read_csv("data/hiv_protease/benchmark/friedman/diff_matrix.csv", index_col=0)
    # with open("data/hiv_protease/benchmark/friedman/nemenyi.yaml") as f:
    #     nm = yaml.safe_load(f)
    # cdc = dot_chart(df_f1, df_cd, nm["cd"], dataset).to_json(indent=None)
    #
    # df_time = pd.read_csv("nodes/vis/misc/benchmark_p.csv", index_col=0)
    # tc = elapsed_time_chart(df_time).to_json(indent=None)
    #
    # df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    # df_fir = pd.read_csv(f"data/{dataset}/benchmark/feature_importance.csv", index_col=0)
    # pcc = par_coord_chart(df_f1, df_fir).to_json(indent=None)

    df_f1 = pd.read_csv(f"data/{dataset}/benchmark/metrics/f1.csv", index_col=0)
    nc = dendro_chart(df_f1)

    with open('bp2.html', 'w') as f:
        f.write(two_charts_template.format(
            vega_version=alt.VEGA_VERSION,
            vegalite_version=alt.VEGALITE_VERSION,
            vegaembed_version=alt.VEGAEMBED_VERSION,
            dataset=dataset,
            # spec1=oc,
            # spec2=mc,
            # spec3=rc,
            # spec4=sc,
            # spec5=sdc,
            # spec6=cdc,
            # spec7=tc,
            # spec8=pcc,
            spec9=nc
        ))


if __name__ == '__main__':
    main()
