from dash import Dash
import dash_html_components as html

from . import scatter_plot, sequence_dropdown


def create_layout(app: Dash) -> html.Div:
    return html.Div(
        className="app-div",
        children=[
            html.H1(app.title),
            html.Hr(),
            html.Div(
                className="dropdown-container",
                children=[
                    sequence_dropdown.render(app)
                ],
            ),
            scatter_plot.render(app),
            html.H5("The div after the chart.")
        ]
    )