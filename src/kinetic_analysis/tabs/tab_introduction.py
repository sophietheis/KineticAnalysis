from dash import html

def layout():
    return(
        html.Br(),
        html.P("In this programme, you can find different tabs: "),
        html.Li("Generate tracks"),
        html.Li("Choose the equation"),
        html.Li("Track analysis - display & all tracks"),
        html.Li("Count ribosomes"),
        html.Br(),
        html.Br(),
        html.P ("This page is in progress..."),
        html.Img(src = "/assets/images/Page_Under_Construction.png",
                 style={'height':'75%',
                        'width':'75%'})
    )