from dash import html, dcc

from kineticanalysis.utils.texts import t


def layout():
    return (
        html.Br(),
        html.P(t("tab_intros.introduction")),
        html.P(""),
        html.P(t("tab_intros.introduction_tabs")),
        html.P(""),
        dcc.Markdown(children="""
        - Acquisition parameters
        - Generate tracks
        - Choose the equation
        - Track analysis
        - Count ribosomes
        - MSD
        - Combine
        """),

        html.Br(),
        html.Br(),

        # Page under construction
        # html.Img(src="/assets/images/Page_Under_Construction.png",
        #          style={'height': '50%',
        #                 'width': '50%'}),
        # html.Br(),
        # html.Br(),

    )
