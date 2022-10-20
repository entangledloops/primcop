from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd

from . import ids
from src.sample_data import samples_dict

import src.analysis as analysis

def render(app: Dash) -> html.Div:
    @app.callback(
        Output(ids.TABLE, "children"),
        Input(ids.SEQUENCE_DROPDOWN, "value")
    )
    def update_table(value, max_rows=10) -> html.Table:
        table_df = analysis.update_df(value).round(4).drop(columns=["Sequence_ID"])

        return html.Table([
            html.Thead(
                html.Tr([html.Th(col) for col in table_df.columns])
            ),
            html.Tbody([
                html.Tr([
                    html.Td(table_df.iloc[i][col]) for col in table_df.columns
                ]) for i in range(min(len(table_df), max_rows))
            ])
        ])
    return html.Div(id=ids.TABLE)