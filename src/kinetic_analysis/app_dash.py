from dash import Dash, html, Input, Output
import dash_bootstrap_components as dbc

from kineticanalysis.tabs.tab_introduction import layout as tab0_layout

from kineticanalysis.tabs.tab_acquisition import layout as tab1_layout
from kineticanalysis.callbacks.callback_acquisition import (
    register_callbacks as tab1_callbacks)

from kineticanalysis.tabs.tab_generate_track import layout as tab2_layout
from kineticanalysis.callbacks.callback_generate_track import (
    register_callbacks as tab2_callbacks)

from kineticanalysis.tabs.tab_equation_choice import layout as tab3_layout
from kineticanalysis.callbacks.callback_equation_choice import (
    register_callbacks as tab3_callbacks)

from kineticanalysis.tabs.tab_analyse_kinetic import layout as tab4_layout
from kineticanalysis.callbacks.callback_analyse_kinetic import (
    register_callbacks as tab4_callbacks)

from kineticanalysis.tabs.tab_count_ribosome import layout as tab5_layout
from kineticanalysis.callbacks.callback_count_ribosome import (
    register_callbacks as tab5_callbacks)

from kineticanalysis.tabs.tab_MSD import layout as tab6_layout
from kineticanalysis.callbacks.callback_MSD import (
    register_callbacks as tab6_callbacks)

from kineticanalysis.tabs.tab_combine_track import layout as tab7_layout
from kineticanalysis.callbacks.callback_combine_track import (
    register_callbacks as tab7_callbacks)

from kineticanalysis.tabs.not_found_404 import layout as not_found_layout


FONT_AWESOME = "https://use.fontawesome.com/releases/v5.10.2/css/all.css"
app = Dash(__name__,
           external_stylesheets=[dbc.themes.FLATLY, FONT_AWESOME],
           )
app.server.static_folder = "assets"
app.title = "Translation dynamics app"

# Global variables to store states
# Thread safety need to be changed
app.data = {
    'directory_generation': None,
    'directory_analysis': None,
    'directory_analysis_vivo': None,
    'csv_files': [],
    'fig': None,
    'selected_file': None,
    'solver': "Exact equation",
    'csv_to_analyse': None,
}

_TAB_LAYOUTS = {
    'tab-0': tab0_layout,
    'tab-1': tab1_layout,
    'tab-2': tab2_layout,
    'tab-3': tab3_layout,
    'tab-4': tab4_layout,
    'tab-5': tab5_layout,
    'tab-6': tab6_layout,
    'tab-7': tab7_layout,
}

app.layout = dbc.Container([
    # Header
    html.Header([
        html.H1("Translation dynamic analysis app"),
        # dcc.Markdown("*Tool to estimate the kinetic parameters"
        # " of translation.*"),
    ]),

    # Body
    dbc.Tabs(children=[
        dbc.Tab(label="Introduction",
                tab_id="tab-0",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Acquisition Parameters",
                tab_id="tab-1",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Generate tracks",
                tab_id="tab-2",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Choose the equation",
                tab_id="tab-3",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Track Analysis",
                tab_id="tab-4",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Count ribosomes",
                tab_id="tab-5",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="MSD",
                tab_id="tab-6",
                activeTabClassName="fw-bold fst-italic"),
        dbc.Tab(label="Combine",
                tab_id="tab-7",
                activeTabClassName="fw-bold fst-italic"),
        ],

        id='tabs',
        active_tab='tab-0',
    ),
    html.Div(id='tabs-content',
             style={'min_height': "80vh"}),

    html.Br(),
    html.Br(),
])


@app.callback(
    Output('tabs-content', 'children'),
    Input('tabs', 'active_tab')
)
def render_content(tab):
    return _TAB_LAYOUTS.get(tab, not_found_layout)()


# Register callbacks
tab1_callbacks(app)
tab2_callbacks(app)
tab3_callbacks(app)
tab4_callbacks(app)
tab5_callbacks(app)
tab6_callbacks(app)
tab7_callbacks(app)

if __name__ == '__main__':
    # run_app()
    app.run_server(debug=False, host="0.0.0.0", port=5001)