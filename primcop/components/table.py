import dash
import dash_core_components as dcc
import dash_html_components as html
from dash import Dash
from dash.dependencies import Input, Output
import pandas as pd

import primcop.analysis as analysis
from primcop.components.sequence_dropdown import SEQUENCE_DROPDOWN


TABLE = "table"


def render(app: Dash) -> html.Div:
    @app.callback(Output(TABLE, "children"), Input(SEQUENCE_DROPDOWN, "value"))
    def update_table(value, max_rows=10) -> html.Table:
        if not value:
            raise dash.exceptions.PreventUpdate()
        table_df = analysis.get_df(value).round(4)

        return html.Table(
            [
                html.Thead(html.Tr([html.Th(col) for col in table_df.columns])),
                html.Tbody(
                    [
                        html.Tr(
                            [html.Td(table_df.iloc[i][col]) for col in table_df.columns]
                        )
                        for i in range(min(len(table_df), max_rows))
                    ]
                ),
            ]
        )

    return html.Div(id=TABLE)
