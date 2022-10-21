import dash
import dash_core_components as dcc
import dash_html_components as html
from dash import Dash
from dash.dependencies import Input, Output
import plotly.graph_objs as go

import primcop.analysis as analysis
from primcop.components.sequence_dropdown import SEQUENCE_DROPDOWN


SCATTER_PLOT = "scatter-plot"


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

def get_papa_axis(ranges):
    return dict(
        title=analysis.PAPA,
        range=ranges,
        titlefont=dict(color="#1f77b4"),
        tickfont=dict(color="#1f77b4"),
    )

def get_prima_axis(ranges):
    return dict(
        title=analysis.PRIMA,
        range=ranges,
        titlefont=dict(color="#ff7f0e"),
        tickfont=dict(color="#ff7f0e"),
        anchor="free",
        overlaying="y",
        side="left",
        position=0.20,
    )

def get_fold_index(ranges):
    return dict(
        title=analysis.FOLD_INDEX,
        range=ranges,
        titlefont=dict(color="#FF0000"),
        tickfont=dict(color="#FF0000"),
        anchor="x",
        overlaying="y",
        side="right",
    )

SCORE_METHODS_YAXIS = {
    analysis.PAPA: get_papa_axis,
    analysis.PRIMA: get_prima_axis,
    analysis.FOLD_INDEX: get_fold_index,
}

def render(app: Dash) -> html.Div:
    @app.callback(
        Output(SCATTER_PLOT, "children"),
        Input(SEQUENCE_DROPDOWN, "value"),
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
        method_names = [name for name in df if name in SCORE_METHODS_YAXIS]
        yaxes = {
            f"yaxis{i+1}": SCORE_METHODS_YAXIS[name](ranges[name])
            for i, name in enumerate(method_names)
        }
        fig = {
            "data": plots,
            "layout": go.Layout(
                title="Measures of Potential Prion Activity",
                xaxis={"title": "Amino Acid Position in Sequence"},
                **yaxes,
                legend={"x": -0.15, "y": 1},
                hovermode="closest",
            ),
        }
        return dcc.Graph(figure=fig)

    return html.Div(id=SCATTER_PLOT)
