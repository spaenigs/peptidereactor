import altair as alt

def get_template(spec_tuples, dataset):

    divs, specs = "", ""
    for i, (header, spec) in enumerate(spec_tuples, start=1):
        divs += f"""
            <h3><u>{header}</u></h3>
            <div>
                <div style="text-align: center;"><div id='vis{i}'></div></div>
            </div>
        """
        specs += f"vegaEmbed('#vis{i}', {spec}).catch(console.error);\n"

    return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <script src="https://cdn.jsdelivr.net/npm/vega@{alt.VEGA_VERSION}"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@{alt.VEGALITE_VERSION}"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@{alt.VEGAEMBED_VERSION}"></script>
    <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
    <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
    <title>PEPTIDE REACToR ({dataset})</title>
</head>
<body>
    <p align = "left">
        <img src="/home/spaenigs/Downloads/peptide-reactor/peptide-reactor/peptide-reactor.005.png" 
             alt="logo" 
             width="7%" 
             height="7%"
        />
        <hr>
    </p>    
    <div style="text-align: center;"> <p><b>{dataset}</b></p> </div>
    <div id="accordion">
    {divs}
    </div>
    <script type="text/javascript">
        {specs}
        $( function() {{
            $( "#accordion" ).accordion({{collapsible: true}});
        }});
    </script>
</body>
</html>
"""
