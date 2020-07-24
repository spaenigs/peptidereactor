def get_template(links):
    return f"""
<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta content="width=device-width, initial-scale=1, shrink-to-fit=no" name="viewport">

    <title>PEPTIDE REACToR</title>

    <link href="https://getbootstrap.com/docs/4.5/examples/dashboard/" rel="canonical">

    <!-- Bootstrap core CSS -->
    <link crossorigin="anonymous" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/css/bootstrap.min.css"
          integrity="sha384-9aIt2nRpC12Uk9gS9baDl411NQApFmC26EwAOH8WgZl5MYYxFfc+NcPb1dKGj7Sk" rel="stylesheet">

    <style>
        .bd-placeholder-img {{
            font-size: 1.125rem;
            text-anchor: middle;
            -webkit-user-select: none;
            -moz-user-select: none;
            -ms-user-select: none;
            user-select: none;
        }}

        @media (min-width: 768px) {{
            .bd-placeholder-img-lg {{
                font - size: 3.5rem;
            }}
        }}
    </style>

    <style>
        body {{
            font-size: .875rem;
        }}

        .feather {{
            width: 16px;
            height: 16px;
            vertical-align: text-bottom;
        }}

        /*
         * Sidebar
         */

        .sidebar {{
            position: fixed;
            top: 0;
            bottom: 0;
            left: 0;
            z-index: 100; /* Behind the navbar */
            padding: 48px 0 0; /* Height of navbar */
            box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
        }}

        @media (max-width: 767.98px) {{
            .sidebar {{
            top: 5rem;
            }}
        }}

        .sidebar-sticky {{
            position: relative;
            top: 0;
            height: calc(100vh - 48px);
            padding-top: .5rem;
            overflow-x: hidden;
            overflow-y: auto; /* Scrollable contents if viewport is shorter than content. */
        }}

        @supports ((position: -webkit-sticky) or (position: sticky)) {{
            .sidebar-sticky {{
               position: -webkit-sticky;
                position: sticky;
            }}
        }}

        .sidebar .nav-link {{
            font-weight: 500;
            color: #333;
        }}

        .sidebar .nav-link .feather {{
            margin-right: 4px;
            color: #999;
        }}

        .sidebar .nav-link.active {{
            color: #007bff;
        }}

        .sidebar .nav-link:hover .feather,
        .sidebar .nav-link.active .feather {{
            color: inherit;
        }}

        .sidebar-heading {{
            font-size: .75rem;
            text-transform: uppercase;
        }}

        /*
         * Navbar
         */

        .navbar-brand {{
            padding-top: .75rem;
            padding-bottom: .75rem;
            font-size: 1rem;
            //background-color: rgba(0, 0, 0, .25);
            box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);    
        }}

        .navbar .navbar-toggler {{
            top: .25rem;
            right: 1rem;
        }}

        .navbar .form-control {{
            padding: .75rem 1rem;
            border-width: 0;
            border-radius: 0;
        }}

        .form-control-dark {{
            color: #fff;
            background-color: rgba(255, 255, 255, .1);
            border-color: rgba(255, 255, 255, .1);
        }}

        .form-control-dark:focus {{
            border-color: transparent;
            box-shadow: 0 0 0 3px rgba(255, 255, 255, .25);
        }}

    </style>
</head>
<body>
<nav class="navbar navbar-white sticky-top bg-white flex-md-nowrap p-0 border-bottom">
    <a class="navbar-brand col-md-3 col-lg-1 mr-0 px-3" href="#">
    <img src="/home/spaenigs/Downloads/peptide-reactor/peptide-reactor/peptide-reactor.005.png" 
             alt="logo" 
             width="100%" 
             height="100%"
        />
    </a>
    <button aria-controls="sidebarMenu" aria-expanded="false" aria-label="Toggle navigation"
            class="navbar-toggler position-absolute d-md-none collapsed" data-target="#sidebarMenu" data-toggle="collapse" type="button">
        <span class="navbar-toggler-icon"></span>
    </button>
    <!--
    <input aria-label="Search" class="form-control form-control-dark w-100" placeholder="Search" type="text"> 
    <ul class="navbar-nav px-3">
        <li class="nav-item text-nowrap">
            <a class="nav-link" href="#">Sign out</a>
        </li>
    </ul>
    -->
</nav>
<div class="container-fluid"> -->
    <div class="row">
        <nav class="col-md-3 col-lg-1 d-md-block bg-white sidebar collapse" id="sidebarMenu">
            <div class="sidebar-sticky pt-5">
                <ul class="nav flex-column">
                    {links}
                </ul>
            </div>
        </nav>
        <main class="col-md-9 ml-sm-auto col-lg-11 px-md-4" role="main">
            <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3">
                <div class="tab-content" id="myTabContent">
                    <div id='vis1'>

                    </div>
                </div>
            </div>
        </main>
    </div>
</div>
<div class="tab-content" id="myTabContent">
    <div aria-labelledby="home-tab" class="tab-pane fade show active" id="home" role="tabpanel">asd</div>
    <div aria-labelledby="profile-tab" class="tab-pane fade" id="profile" role="tabpanel">ayxcyxsd</div>
    <div aria-labelledby="contact-tab" class="tab-pane fade" id="contact" role="tabpanel">.yxcc</div>
</div>
<script crossorigin="anonymous"
        integrity="sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj"
        src="https://code.jquery.com/jquery-3.5.1.slim.min.js"></script>
<script>window.jQuery || document.write('<script src="../assets/js/vendor/jquery.slim.min.js"><\/script>')</script>
<script crossorigin="anonymous"
        integrity="sha384-1CmrxMRARb6aLqgBO7yyAxTOQE2AKb9GfXnEo760AUcUmFx3ibVJJAzGytlQcNXd"
        src="https://stackpath.bootstrapcdn.com/bootstrap/4.5.0/js/bootstrap.bundle.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/feather-icons/4.9.0/feather.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/2.7.3/Chart.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@4.8.1"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
<script>
    vegaEmbed('#vis1', {{
      "$schema": "https://vega.github.io/schema/vega/v5.json",
      "description": "A heatmap showing average daily temperatures in Seattle for each hour of the day.",
      "width": 800,
      "height": 500,
      "padding": 5,
    
      "title": {{
            "text": "Seattle Annual Temperatures",
            "anchor": "middle",
            "fontSize": 16,
            "frame": "group",
            "offset": 4
        }},
    
      "signals": [
        {{
            "name": "palette", "value": "Viridis",
            "bind": {{
                "input": "select",
                "options": [
                    "Viridis",
                    "Magma",
                    "Inferno",
                    "Plasma",
                    "DarkBlue",
                    "DarkGold",
                    "DarkGreen",
                    "DarkMulti",
                    "DarkRed",
                    "LightGreyRed",
                    "LightGreyTeal",
                    "LightMulti",
                    "LightOrange",
                    "LightTealBlue",
                    "Blues",
                    "Browns",
                    "Greens",
                    "Greys",
                    "Oranges",
                    "Purples",
                    "Reds",
                    "TealBlues",
                    "Teals",
                    "WarmGreys",
                    "BlueOrange",
                    "BrownBlueGreen",
                    "PurpleGreen",
                    "PinkYellowGreen",
                    "PurpleOrange",
                    "RedBlue",
                    "RedGrey",
                    "RedYellowBlue",
                    "RedYellowGreen",
                    "BlueGreen",
                    "BluePurple",
                    "GoldGreen",
                    "GoldOrange",
                    "GoldRed",
                    "GreenBlue",
                    "OrangeRed",
                    "PurpleBlueGreen",
                    "PurpleBlue",
                    "PurpleRed",
                    "RedPurple",
                    "YellowGreenBlue",
                    "YellowGreen",
                    "YellowOrangeBrown",
                    "YellowOrangeRed"
                ]
            }}
        }},
        {{
            "name": "reverse", "value": false, "bind": {{"input": "checkbox"}}
        }}
      ],
    
      "data": [
        {{
            "name": "temperature",
            "url": "data/seattle-temps.csv",
            "format": {{"type": "csv", "parse": {{"temp": "number", "date": "date"}}}},
            "transform": [
                {{"type": "formula", "as": "hour", "expr": "hours(datum.date)"}},
                {{"type": "formula", "as": "day",
                 "expr": "datetime(year(datum.date), month(datum.date), date(datum.date))"}}
            ]
        }}
      ],
    
      "scales": [
        {{
            "name": "x",
            "type": "time",
            "domain": {{"data": "temperature", "field": "day"}},
            "range": "width"
        }},
        {{
            "name": "y",
            "type": "band",
            "domain": [
                6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
                0, 1, 2, 3, 4, 5
            ],
            "range": "height"
        }},
        {{
            "name": "color",
            "type": "linear",
            "range": {{"scheme": {{"signal": "palette"}}}},
            "domain": {{"data": "temperature", "field": "temp"}},
            "reverse": {{"signal": "reverse"}},
            "zero": false, "nice": true
        }}
      ],
    
      "axes": [
        {{"orient": "bottom", "scale": "x", "domain": false, "title": "Month", "format": "%b"}},
        {{
          "orient": "left", "scale": "y", "domain": false, "title": "Hour",
          "encode": {{
            "labels": {{
              "update": {{
                "text": {{"signal": "datum.value === 0 ? 'Midnight' : datum.value === 12 ? 'Noon' : datum.value < 12 ? datum.value + ':00 am' : (datum.value - 12) + ':00 pm'"}}
              }}
            }}
          }}
        }}
      ],
    
      "legends": [
        {{
          "fill": "color",
          "type": "gradient",
          "title": "Avg. Temp (°F)",
          "titleFontSize": 12,
          "titlePadding": 4,
          "gradientLength": {{"signal": "height - 16"}}
        }}
      ],
    
      "marks": [
        {{
          "type": "rect",
          "from": {{"data": "temperature"}},
          "encode": {{
            "enter": {{
              "x": {{"scale": "x", "field": "day"}},
              "y": {{"scale": "y", "field": "hour"}},
              "width": {{"value": 5}},
              "height": {{"scale": "y", "band": 1}},
              "tooltip": {{"signal": "timeFormat(datum.date, '%b %d %I:00 %p') + ': ' + datum.temp + '°'"}}
            }},
            "update": {{
              "fill": {{"scale": "color", "field": "temp"}}
            }}
          }}
        }}
      ]
    }})

</script>
</body>
</html>
"""

from glob import glob

import re

with open("overview.html", "w") as f:
    links = ""
    for p in glob("data/*/vis/single_dataset.html"):
        dataset = re.findall("data/(.*?)/", p)[0]
        links += f"""
        <li class="nav-item">
            <a class="nav-link" href="{p}" target="_blank">
                <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-external-link"><path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"></path><polyline points="15 3 21 3 21 9"></polyline><line x1="10" y1="14" x2="21" y2="3"></line></svg> 
                {dataset}
            </a>
        </li>
        """
        """
        <a aria-controls="home" aria-selected="true" class="nav-link" data-toggle="tab" href="{p}" target="_blank"
               id="{dataset}" role="tab">
            <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-folder"><path d="M22 19a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5l2 3h9a2 2 0 0 1 2 2z"></path></svg>
                {dataset}
            </a>
        """

    html_string = get_template(links)
    f.write(html_string)
    f.flush()



