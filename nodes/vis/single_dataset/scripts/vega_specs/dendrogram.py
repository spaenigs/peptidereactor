def vega_dendrogram(values):
    json_dict = {
        "$schema": "https://vega.github.io/schema/vega/v5.json",
        "description": "Based on the example of a radial layout for a node-link diagram of hierarchical data.",
        "width": 1000,
        "height": 1000,
        "padding": 5,
        "autosize": "none",
        "signals": [
            {
                "name": "labels", "value": True,
                "bind": {"input": "checkbox"}
            },
            {
                "name": "radius", "value": 300,
                "bind": {"input": "range", "min": 100, "max": 1000}
            },
            {
                "name": "rotate", "value": 0,
                "bind": {"input": "range", "min": 0, "max": 360, "step": 1}
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
                        "method": "cluster",
                        "size": [1, {"signal": "radius"}],
                        "as": ["alpha", "radius", "depth", "children"]
                    },
                    {
                        "type": "formula",
                        "expr": "(rotate + 360 * datum.alpha + 270) % 360",
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
                        "shape": "orthogonal", "orient": "radial",
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
                "range": {"scheme": "tableau20"},
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