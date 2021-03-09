from collections import OrderedDict

from streamlit.components.v1 import html
from glob import glob

import altair as alt
import streamlit as st

import re
import json


def get_options_dict(paths):
    options_dict = OrderedDict()
    for p in paths:
        name = p.split("/")[-1]
        options_dict[name.split("_")[-1]] = name
    return options_dict


def set_chart(path):
    with open(path) as f:
        d = json.load(f)
        if "/vega/" in d["$schema"]:
            with open("peptidereactor-vis/resources/vega_template.html") as f:
                content = f.read().replace("spec", json.dumps(d))
                html(content, width=1200, height=1500, scrolling=True)
        else:
            c2 = alt.Chart().from_dict(d)
            st.altair_chart(c2, use_container_width=True)


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

    url = "https://scholar.google.de/citations?user=lEVtMBMAAAAJ&hl=en"
    link_style = "text-decoration:none;"
    text = "A tool for an <b>in-depth comparison</b> and <b>benchmarking</b> of <b>peptide encodings</b>. " \
           "All computations are <b>highly parallelized</b> and work efficiently across <b>multiple datasets and " \
           f"encodings</b>. For a thorough introduction refer to Spänig <i>et al.</i> (manuscript in preparation)." \
           # f" <a href='{url}' style='{link_style}'>Spänig <i>et al.</i> (2020)</a>."
    st.markdown(f"<div style='text-align: justify'>{text}</div>", unsafe_allow_html=True)

    st.text("")

    t2 = "Click on a dataset (left) to show tSNE-based embedding of sequences part of the positive class (right) " \
         "as well as detailed dataset information (name, size, description, and DOI) (below). To download the " \
         "original dataset, click on the symbol below. Start exploring:"
    st.markdown(f"<div style='text-align: justify'>{t2}</div>", unsafe_allow_html=True)

    st.text("")

    with open(f"data/multiple_datasets/vis/home_Home_tsne/home_Home_tsne.json") as f:
        d = json.load(f)
        c2 = alt.Chart().from_dict(d)
        st.altair_chart(c2, use_container_width=True)

    st.text("")
    st.text("")
    st.text("")
    st.text("")

    st.text("Copyright (c) 2021 Heiderlab")

elif analysis == MDS:
    paths = sorted(glob("nodes/vis/mds_*"))
    options_dict = get_options_dict(paths)

    options = list(options_dict.keys())
    option = st.sidebar.radio("Choose vis:", options)

    set_chart(f"data/multiple_datasets/vis/{options_dict[option]}/{options_dict[option]}.json")

elif analysis == SDS:
    datasets = []
    for p in glob("data/*/vis/"):
        if "multiple_datasets" in p:
            continue
        dataset = re.findall("data/(.*?)/vis/", p)[0]
        datasets += [dataset]

    ds = st.sidebar.selectbox("Choose dataset:", sorted(datasets))

    paths = sorted(glob("nodes/vis/sds_*"))
    options_dict = get_options_dict(paths)

    options = list(options_dict.keys())
    option = st.sidebar.radio("Choose vis:", options)

    set_chart(f"data/{ds}/vis/{options_dict[option]}/{options_dict[option]}.json")