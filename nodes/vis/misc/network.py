from pathos.multiprocessing import ProcessingPool as Pool

from scipy import stats
from itertools import combinations

import numpy as np
import pandas as pd

import json


def network_chart(df_cd):

    def vega(node_values, link_values):

        json_dict = {
            "$schema": "https://vega.github.io/schema/vega/v5.json",
            "description": "A node-link diagram with force-directed layout, depicting character co-occurrence in the novel Les Misérables.",
            "width": 700, "height": 500, "padding": 0, "autosize": "none",
            "signals": [
                {"name": "cx", "update": "width / 2"},
                {"name": "cy", "update": "height / 2"},
                {"name": "nodeRadius", "value": 8, "bind": {"input": "range", "min": 1, "max": 50, "step": 1}},
                {"name": "nodeCharge", "value": -30, "bind": {"input": "range", "min": -100, "max": 10, "step": 1}},
                {"name": "linkDistance", "value": 30, "bind": {"input": "range", "min": 5, "max": 100, "step": 1}},
                {"name": "static", "value": True, "bind": {"input": "checkbox"}},
                {"description": "State variable for active node fix status.",
                 "name": "fix", "value": False,
                 "on": [
                     {"events": "symbol:mouseout[!event.buttons], window:mouseup", "update": "false"},
                     {"events": "symbol:mouseover", "update": "fix || true"},
                     {"events": "[symbol:mousedown, window:mouseup] > window:mousemove!", "update": "xy()",
                      "force": True}
                 ]},
                {"description": "Graph node most recently interacted with.", "name": "node", "value": None,
                 "on": [
                     {"events": "symbol:mouseover", "update": "fix === true ? item() : node"}
                 ]},
                {"description": "Flag to restart Force simulation upon data changes.", "name": "restart",
                 "value": False,
                 "on": [
                     {"events": {"signal": "fix"}, "update": "fix && fix.length"}
                 ]}
            ], "scales": [
                {
                    "name": "color",
                    "type": "ordinal",
                    "domain": {"data": "node-data", "field": "group"},
                    "range": {"scheme": "category20c"}
                }
            ], "marks": [
                {
                    "name": "nodes",
                    "type": "symbol",
                    "zindex": 1,
                    "from": {"data": "node-data"},
                    "on": [
                        {"trigger": "fix", "modify": "node",
                         "values": "fix === true ? {fx: node.x, fy: node.y} : {fx: fix[0], fy: fix[1]}"},
                        {"trigger": "!fix", "modify": "node", "values": "{fx: null, fy: null}"}
                    ],
                    "encode": {
                        "enter": {
                            "fill": {"scale": "color", "field": "group"},
                            "stroke": {"value": "white"}
                        },
                        "update": {
                            "size": {"signal": "2 * nodeRadius * nodeRadius"},
                            "cursor": {"value": "pointer"}
                        }
                    },

                    "transform": [
                        {
                            "type": "force",
                            "iterations": 300,
                            "restart": {"signal": "restart"},
                            "static": {"signal": "static"},
                            "signal": "force",
                            "forces": [
                                {"force": "center", "x": {"signal": "cx"}, "y": {"signal": "cy"}},
                                {"force": "collide", "radius": {"signal": "nodeRadius"}},
                                {"force": "nbody", "strength": {"signal": "nodeCharge"}},
                                {"force": "link", "links": "link-data", "distance": {"signal": "linkDistance"}}
                            ]
                        }
                    ]
                },
                {
                    "type": "path",
                    "from": {"data": "link-data"},
                    "interactive": False,
                    "encode": {
                        "update": {
                            "stroke": {"value": "#ccc"},
                            "strokeWidth": {"value": 0.5}
                        }
                    },
                    "transform": [
                        {
                            "type": "linkpath",
                            "require": {"signal": "force"},
                            "shape": "line",
                            "sourceX": "datum.source.x", "sourceY": "datum.source.y",
                            "targetX": "datum.target.x", "targetY": "datum.target.y"
                        }
                    ]
                }
            ],
            "data": [{
                "name": "node-data",
                "values": node_values
            },
                {
                    "name": "link-data",
                    "values": link_values
                }
            ]
        }

        return json_dict

    node_values = [{"name": "Myriel", "group": 1, "index": 0}, {"name": "Napoleon", "group": 1, "index": 1}]
    link_values = [{"source": 1, "target": 0, "value": 10}]

    # TODO replace '.' in 'lambda.corr' in x or 'g.gap'
    groups = ["psekraac" if "lambda.corr" in x or "g.gap" in x else x[:6] for x in df_cd.index]

    node_values1 = []
    indices = {}
    for i, (e, g) in enumerate(zip(df_cd.index, groups)):
        node_values1 += [{"name": e, "group": g, "index": i}]
        indices[e] = i

    link_values1 = []
    for x, y in np.column_stack(np.triu_indices_from(df_cd.values, k=1)):
        i, c, v = df_cd.index[x], df_cd.columns[y], df_cd.iloc[x, y]
        if np.abs(v) >= 0.7:
            link_values1 += [{"source": indices[i], "target": indices[c], "value": v}]

    v = vega(node_values1, link_values1)

    return json.dumps(v)
