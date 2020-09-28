import pandas as pd
import altair as alt
import streamlit as st
import numpy as np

import secrets

from glob import glob

# from snakemake import shell
#
# with open("server1.py", "w") as f:
#     f.write("""
# from flask import Flask, render_template
# from flask_cors import CORS
#
# app = Flask(__name__, static_folder='')
# CORS(app)
# """)
#
#     for p in glob("data/ace_vaxinpad/vis/*/*.json"):
#         f.write(f"""
# @app.route("/{p}")
# def ace_vaxinpad_{secrets.token_hex(4)}():
#     return app.send_static_file("{p}")
# """)
#
#     f.write("""
# if __name__ == '__main__':
#   app.run(debug=True)
# """)
#
#     f.flush()
#
# shell("export FLASK_APP=server1.py; flask run")

with open("peptidereactor-vis/resources/style.css") as f:
    st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

st.sidebar.image("http://192.168.178.30:8501/logo.png", width=200)

import json

ds = st.sidebar.selectbox("Choose dataset:", ["ace_vaxinpad", "hiv_protease"])

with open(f"data/{ds}/vis/overview/overview.json") as f:
    c1 = alt.Chart().from_dict(json.load(f))

with open(f"data/{ds}/vis/metrics/metrics.json") as f:
    c2 = alt.Chart().from_dict(json.load(f))

option = st.sidebar.selectbox("Choose vis:", ["overview", "metrics"])

if option == "overview":
    c = c1
elif option == "metrics":
    c = c2
else:
    c = c1

st.altair_chart(c)
