from streamlit.components.v1 import html
from glob import glob

import pandas as pd
import altair as alt
import streamlit as st

import re
import json

st.set_page_config(
    page_title="PEPTIDE REACToR",
    page_icon="http://192.168.178.30:8501/favicon.png"
)

with open("peptidereactor-vis/resources/style.css") as f:
    st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

st.sidebar.image("http://127.0.0.1:8501/logo.png", width=200)

HOME, MDS, SDS = "Home", "Multiple datasets", "Single dataset"

analysis = st.sidebar.radio("Choose analysis:", [HOME, MDS, SDS])

if analysis == HOME:
    st.header("Welcome to the PEPTIDE REACToR")

    text = "A tool for <b>in-depth comparison</b> and <b>benchmarking</b> of <b>peptide encodings</b>. " \
           "All computations are <b>highly parallelized</b> and work efficiently across <b>multiple datasets and " \
           "encodings</b>."
    st.markdown(f"<div style='text-align: justify'>{text}</div>", unsafe_allow_html=True)

    st.text("")

    t2 = "Click on a dataset (left) to show tSNE-based embedding of sequences part of the positive class (right). " \
         "as well as detailed dataset information (name, size, description, and DOI) (below). Start exploring:"
    st.markdown(f"<div style='text-align: justify'>{t2}</div>", unsafe_allow_html=True)

    st.text("")

    with open(f"data/multiple_datasets/vis/home_Home_tsne/home_tsne.json") as f:
        d = json.load(f)
        c2 = alt.Chart().from_dict(d)
        st.altair_chart(c2, use_container_width=True)

    st.text("")
    st.text("")

    st.text("Copyright (c) 2020 Heiderlab")

elif analysis == MDS:
    v = glob(f"data/multiple_datasets/vis/*/")

    options = ["Overview", "Ranks", "Clustering", "Clustering (alt)", "Embedding", "Elapsed time"]
    option = st.sidebar.radio("Choose vis:", options)

    use_alt = False

    if option == options[0]:
        option = "md_overview_hm"
    elif option == options[1]:
        option = "md_ranks_hm"
    elif option == options[2]:
        option = "md_clustering_hm"
    elif option == options[3]:
        use_alt = True
        option = "md_clustering2_hm"
    elif option == options[4]:
        option = "md_tsne"
    elif option == options[5]:
        option = "md_elapsed_time"

    tmp = option.replace("clustering2", "clustering") if use_alt else option

    # options = [re.findall(f"data/multiple_datasets/vis/(.*?)/", v_)[0] for v_ in v]
    # options += ["md_clustering2_hm"]
    # option = st.sidebar.radio("Choose vis:", options, format_func=lambda x: x.replace("md_", ""))

    # tmp = option.replace("clustering2", "clustering") if "clustering2" in option else option

    with open(f"data/multiple_datasets/vis/{tmp}/{option}.json") as f:
        d = json.load(f)
        if "/vega/" in d["$schema"]:
            with open("peptidereactor-vis/resources/vega_template.html") as f:
                content = f.read().replace("spec", json.dumps(d))
                html(content, width=1200, height=1500, scrolling=True)
        else:
            c2 = alt.Chart().from_dict(d)
            st.altair_chart(c2)

elif analysis == SDS:
    dataset_paths = glob("data/*/vis/")

    res = {}
    for p in dataset_paths:
        if "multiple_datasets" in p:
            continue
        dataset = re.findall("data/(.*?)/vis/", p)[0]
        v = glob(f"data/{dataset}/vis/*/")
        res[dataset] = sorted([re.findall(f"data/{dataset}/vis/(.*?)/", v_)[0] for v_ in v])

    ds = st.sidebar.selectbox("Choose dataset:", sorted(list(res.keys())))
    option = st.sidebar.radio("Choose vis:", res[ds])

    with open(f"data/{ds}/vis/{option}/{option}.json") as f:
        d = json.load(f)
        if "/vega/" in d["$schema"]:
            with open("peptidereactor-vis/resources/vega_template.html") as f:
                content = f.read().replace("spec", json.dumps(d))
                html(content, width=1200, height=1200, scrolling=True)
        else:
            c2 = alt.Chart().from_dict(d)
            st.altair_chart(c2)
