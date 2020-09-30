from streamlit.components.v1 import html, iframe
from glob import glob

import altair as alt
import streamlit as st

import re
import json

with open("peptidereactor-vis/resources/style.css") as f:
    st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

st.sidebar.image("http://192.168.178.30:8501/logo.png", width=200)

analysis = st.sidebar.selectbox("Choose analysis:", ["Multiple datasets", "Single datasets"])

if analysis == "Multiple datasets":
    v = glob(f"data/multiple_datasets/vis/*/")
    options = [re.findall(f"data/multiple_datasets/vis/(.*?)/", v_)[0] for v_ in v]
    option = st.sidebar.selectbox("Choose vis:", options)

    with open(f"data/multiple_datasets/vis/{option}/{option}.json") as f:
        d = json.load(f)
        if "/vega/" in d["$schema"]:
            with open("peptidereactor-vis/resources/vega_template.html") as f:
                content = f.read().replace("spec", json.dumps(d))
                html(content, width=1200, height=1200, scrolling=True)
        else:
            c2 = alt.Chart().from_dict(d)
            st.altair_chart(c2)

else:
    dataset_paths = glob("data/*/vis/")

    res = {}
    for p in dataset_paths:
        if "multiple_datasets" in p:
            continue
        dataset = re.findall("data/(.*?)/vis/", p)[0]
        v = glob(f"data/{dataset}/vis/*/")
        res[dataset] = [re.findall(f"data/{dataset}/vis/(.*?)/", v_)[0] for v_ in v]

    ds = st.sidebar.selectbox("Choose dataset:", list(res.keys()))
    option = st.sidebar.selectbox("Choose vis:", res[ds])

    with open(f"data/{ds}/vis/{option}/{option}.json") as f:
        d = json.load(f)
        if "/vega/" in d["$schema"]:
            with open("peptidereactor-vis/resources/vega_template.html") as f:
                content = f.read().replace("spec", json.dumps(d))
                html(content, width=1200, height=1200, scrolling=True)
        else:
            c2 = alt.Chart().from_dict(d)
            st.altair_chart(c2)

