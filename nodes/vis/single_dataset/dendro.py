import numpy as np

import json


def dendro_chart(df_f1):
    def vega(values):
        json_dict = {
            "$schema": "https://vega.github.io/schema/vega/v5.json",
            "description": "An example of a radial layout for a node-link diagram of hierarchical data.",
            "width": 1720,
            "height": 1720,
            "padding": 5,
            "autosize": "none",
            "signals": [
                {
                    "name": "labels", "value": True,
                    "bind": {"input": "checkbox"}
                },
                {
                    "name": "radius", "value": 600,
                    "bind": {"input": "range", "min": 300, "max": 1000}
                },
                {
                    "name": "extent", "value": 360,
                    "bind": {"input": "range", "min": 0, "max": 360, "step": 1}
                },
                {
                    "name": "rotate", "value": 0,
                    "bind": {"input": "range", "min": 0, "max": 360, "step": 1}
                },
                {
                    "name": "layout", "value": "cluster",
                    "bind": {"input": "radio", "options": ["tidy", "cluster"]}
                },
                {
                    "name": "links", "value": "orthogonal",
                    "bind": {
                        "input": "select",
                        "options": ["line", "curve", "diagonal", "orthogonal"]
                    }
                },
                {"name": "originX", "update": "width / 2"},
                {"name": "originY", "update": "height / 2"}
            ],

            "data": [
                {
                    "name": "tree",
                    "values": values,
                    "transform": [
                        {
                            "type": "stratify",
                            "key": "id",
                            "parentKey": "parent"
                        },
                        {
                            "type": "tree",
                            "method": {"signal": "layout"},
                            "size": [1, {"signal": "radius"}],
                            "as": ["alpha", "radius", "depth", "children"]
                        },
                        {
                            "type": "formula",
                            "expr": "(rotate + extent * datum.alpha + 270) % 360",
                            "as": "angle"
                        },
                        {
                            "type": "formula",
                            "expr": "PI * datum.angle / 180",
                            "as": "radians"
                        },
                        {
                            "type": "formula",
                            "expr": "inrange(datum.angle, [90, 270])",
                            "as": "leftside"
                        },
                        {
                            "type": "formula",
                            "expr": "originX + datum.radius * cos(datum.radians)",
                            "as": "x"
                        },
                        {
                            "type": "formula",
                            "expr": "originY + datum.radius * sin(datum.radians)",
                            "as": "y"
                        }
                    ]
                },
                {
                    "name": "links",
                    "source": "tree",
                    "transform": [
                        {"type": "treelinks"},
                        {
                            "type": "linkpath",
                            "shape": {"signal": "links"}, "orient": "radial",
                            "sourceX": "source.radians", "sourceY": "source.radius",
                            "targetX": "target.radians", "targetY": "target.radius"
                        }
                    ]
                }
            ],

            "scales": [
                {
                    "name": "color",
                    "type": "linear",
                    "range": {"scheme": "viridis"},
                    "domain": {"data": "tree", "field": "color"},
                    "zero": True
                }
            ],

            "marks": [
                {
                    "type": "path",
                    "from": {"data": "links"},
                    "encode": {
                        "update": {
                            "x": {"signal": "originX"},
                            "y": {"signal": "originY"},
                            "path": {"field": "path"},
                            "stroke": {"value": "#ccc"}
                        }
                    }
                },
                {
                    "type": "symbol",
                    "from": {"data": "tree"},
                    "encode": {
                        "enter": {
                            "size": {"value": 100},
                            "stroke": {"value": "#fff"}
                        },
                        "update": {
                            "x": {"field": "x"},
                            "y": {"field": "y"},
                            "fill": {"scale": "color", "field": "color"}
                        }
                    }
                },
                {
                    "type": "text",
                    "from": {"data": "tree"},
                    "encode": {
                        "enter": {
                            "text": {"field": "name"},
                            "fontSize": {"value": 9},
                            "baseline": {"value": "middle"}
                        },
                        "update": {
                            "x": {"field": "x"},
                            "y": {"field": "y"},
                            "dx": {"signal": "(datum.leftside ? -1 : 1) * 6"},
                            "angle": {"signal": "datum.leftside ? datum.angle - 180 : datum.angle"},
                            "align": {"signal": "datum.leftside ? 'right' : 'left'"},
                            "opacity": {"signal": "labels ? 1 : 0"}
                        }
                    }
                }
            ]
        }

        return json_dict

    data = [{"id": 1, "color": 0},
            {"id": 2, "parent": 1, "color": 0},
            {"id": 3, "name": "cluster", "parent": 1, "color": 1},
            {"id": 4, "name": "cluster1", "parent": 1, "color": 1},
            {"id": 5, "name": "AgglomerativeCluster", "parent": 2, "size": 3938, "color": 2},
            {"id": 6, "name": "Agglomerative1Cluster", "parent": 2, "size": 3938, "color": 2}]

    with open("rv.json") as f:
        data = json.load(f)

    # TODO filter encodings (top 20), create df in R, json here
    # TODO check if network works for all encodings (network.py)

    df_f1.apply(np.median).to_frame("median").groupby(
        by=lambda x: "psekraac" if "lambda-corr" in x or "g-gap" in x else x[:6]).max().sort_values(by="median",
                                                                                                    ascending=False).iloc[
    :20, ]

    v = vega(data)

    return json.dumps(v)


"""
{
    "id": 218,
    "name": "Encoder",
    "parent": 216,
    "size": 4060
  },
  {
    "id": 198,
    "name": "ShapeRenderer",
    "parent": 194,
    "size": 2247
  },
"""
