import dash
import dash_core_components as dcc
import dash_html_components as html
from dash import Dash
from dash.dependencies import Input, Output
import plotly.graph_objs as go

import src.analysis as analysis
from src.sample_data import SAMPLES
from . import ids


def get_scatter(**kwargs):
    return go.Scatter(
        **kwargs,
        mode="lines+markers",
        fill="none",
        opacity=0.7,
        marker={"size": 10, "line": {"width": 0.5, "color": "white"}},
    )


def get_bar(**kwargs):
    return go.Bar(
        **kwargs,
        opacity=0.7,
        marker_color="red",
    )


SCORE_METHODS_UI = {
    analysis.PAPA: get_scatter,
    analysis.PRIMA: get_scatter,
    analysis.FOLD_INDEX: get_bar,
}


def render(app: Dash) -> html.Div:
    @app.callback(
        Output(ids.SCATTER_PLOT, "children"),
        Input(ids.SEQUENCE_DROPDOWN, "value"),
    )
    def update_scatter_plot(value) -> dcc.Graph:
        if not value:
            raise dash.exceptions.PreventUpdate()
        df = analysis.get_df(value)
        ranges = analysis.get_ranges(df)
        plots = []
        for i, k in enumerate(ranges.keys(), 1):
            plot = SCORE_METHODS_UI[k](
                x=df[analysis.SEQUENCE_POSITION],
                y=df[k],
                text=df[analysis.AMINO_ACID],
                name=k,
                yaxis=f"y{i}",
            )
            plots.append(plot)
        fig = {
            "data": plots,
            "layout": go.Layout(
                title="Measures of Potential Prion Activity",
                xaxis={"title": "Amino Acid Position in Sequence"},
                yaxis=dict(
                    title=analysis.PAPA,
                    range=ranges[analysis.PAPA],
                    titlefont=dict(color="#1f77b4"),
                    tickfont=dict(color="#1f77b4"),
                ),
                yaxis2=dict(
                    title=analysis.PRIMA,
                    range=ranges[analysis.PRIMA],
                    titlefont=dict(color="#ff7f0e"),
                    tickfont=dict(color="#ff7f0e"),
                    anchor="free",
                    overlaying="y",
                    side="left",
                    position=0.20,
                ),
                yaxis3=dict(
                    title=analysis.FOLD_INDEX,
                    range=ranges[analysis.FOLD_INDEX],
                    titlefont=dict(color="#FF0000"),
                    tickfont=dict(color="#FF0000"),
                    anchor="x",
                    overlaying="y",
                    side="right",
                ),
                legend={"x": -0.15, "y": 1},
                hovermode="closest",
            ),
        }
        return dcc.Graph(figure=fig)

    return html.Div(id=ids.SCATTER_PLOT)
