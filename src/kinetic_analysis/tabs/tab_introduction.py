from dash import html, dcc


def layout():
    return (
        html.Br(),
        html.P("This tool allows you to analyse translation dynamics based "
               "on SunTag system (or similar method)."),
        html.P(""),
        html.P("In this software, you can find different tabs: "),
        dcc.Markdown(children="""
        - Generate tracks
        - Choose the equation
        - Track analysis
        - Count ribosomes
        - MSD
        """),

        html.Br(),
        html.Br(),

        # Page under construction
        # html.Img(src="/assets/images/Page_Under_Construction.png",
        #          style={'height': '50%',
        #                 'width': '50%'}),
        # html.Br(),
        # html.Br(),

        # Footer
        html.Hr(style={'borderWidth': "0.3vh", "width": "25%",
                       "color": "#10D79B"}),
        html.P(["If you encounter any problems, please ",
                html.A("open an issue",
                       href="https://github.com/sophietheis/KineticAnalysis/issues"),
                " along with a detailed description of the problem.",])
    )
