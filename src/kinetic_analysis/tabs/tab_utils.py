from dash import html, dcc
import dash_bootstrap_components as dbc


def title_h4(text: str):
    return html.H4(children=text,
                   style={"text-align": "center",
                          "color": "#10D79B"})


def color_line():
    return html.Hr(style={'borderWidth': "0.3vh", "width": "25%",
                          "color": "#10D79B"})


def button(child: str = "button", id_: str = "", classname: str = "mr-2", width: str = "150px"):
    """
    Create a Dash button with a given text and ID.

    Args:
        text_key (str): The key for the button text in the translation dictionary.
        button_id (str): The ID for the button.
        width (str): The width of the button (default is "150px").

    Returns:
        dash.html.Button: A Dash button component.
    """
    return dbc.Button(child,
                      id=id_,
                      className=classname,
                      style={"width": width}
                      )


def input_case(child: str, input_id: str, input_type: str = 'text', input_value=None):
    """
    Create a Dash input case with a given child, ID, and properties.

    Args:
        child (str): The text to display in the input case.
        input_id (str): The ID for the input.
        input_type (str): The type of the input (default is 'text').
        input_value (str): The value of the input (default is "").
        width (str): The width of the input (default is '200px').

    Returns:
        dash.html.Div: A Dash div component containing the input case.
    """
    s = html.Div([
                html.P(children=child,
                       style={"height": "auto",
                              "margin-bottom": "auto"}),
                dcc.Input(id=input_id,
                          type=input_type,
                          value=input_value,
                          style={'width': '200px'}),
            ]),
    return s


def input_case_with_tooltip(child: str,
                            input_id: str,
                            tooltip_text: str,
                            tooltip_id: str,
                            input_type: str = 'text',
                            input_value=None):
    """
    Create a Dash input case with a given child, ID, and properties.

    Args:
        child (str): The text to display in the input case.
        input_id (str): The ID for the input.
        input_type (str): The type of the input (default is 'text').
        input_value (str): The value of the input (default is "").
        width (str): The width of the input (default is '200px').
    Returns:
        dash.html.Div: A Dash div component containing the input case with a tooltip.
    """
    return html.Div([
        html.P([
            child,
            html.Span(className="fas fa-question-circle",
                      id=tooltip_id,
                      style={"cursor": "pointer",
                             "marginLeft": "5px"})
        ], style={"height": "auto",
                  "margin-bottom": "auto"}),
        dcc.Input(id=input_id,
                  type=input_type,
                  value=input_value,
                  style={'width': '200px'}),
        dbc.Tooltip(tooltip_text,
                    target=tooltip_id),
    ])
