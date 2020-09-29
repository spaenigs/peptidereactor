from glob import glob

import altair as alt
import streamlit as st

import re
import json

with open("peptidereactor-vis/resources/style.css") as f:
    st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

st.sidebar.image("http://192.168.178.30:8501/logo.png", width=200)

dataset_paths = glob("data/*/vis/")

res = {}
for p in dataset_paths:
    dataset = re.findall("data/(.*?)/vis/", p)[0]
    v = glob(f"data/{dataset}/vis/*/")
    res[dataset] = [re.findall(f"data/{dataset}/vis/(.*?)/", v_)[0] for v_ in v]

ds = st.sidebar.selectbox("Choose dataset:", list(res.keys()))
option = st.sidebar.selectbox("Choose vis:", res[ds])

with open(f"data/{ds}/vis/{option}/{option}.json") as f:
    c2 = alt.Chart().from_dict(json.load(f))
    st.altair_chart(c2)

