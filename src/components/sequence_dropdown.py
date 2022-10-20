from dash import Dash
import dash_core_components as dcc
import dash_html_components as html

from src.sample_data import SAMPLES

from . import ids


def render(app: Dash) -> html.Div:
    return html.Div(
        children=[
            dcc.Dropdown(
                id=ids.SEQUENCE_DROPDOWN,
                options=[{"label": k, "value": v} for k, v in SAMPLES.items()],
                placeholder="Please select sample protein sequence from the Alberti et al. dataset below.",
            )
        ]
    )