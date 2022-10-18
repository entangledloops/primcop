from dash import Dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.graph_objs as go

import src.analysis as a
from src.sample_data import samples_dict

from . import ids

def render(app: Dash) -> html.Div:
    @app.callback(
        Output(ids.SCATTER_PLOT, "children"),
        Input(ids.SEQUENCE_DROPDOWN, "value"),
    )
    def update_scatter_plot(value) -> dcc.Graph:
        sequence = value
        sequence_id = [seq_id for seq_id, seq_val in samples_dict.items() if seq_val == value]
        df = a.analyze_sequence(sequence, sequence_id)
        papa_y_range, prima_y_range, fold_index_y_range = a.get_ranges(df)
        main_plot = [
            go.Scatter(
                x=df[df["Sequence_ID"] == i]["Sequence Position"],
                y=df[df["Sequence_ID"] == i]["PAPA score"],
                text=df[df["Sequence_ID"] == i]["Amino Acid"],
                mode="lines+markers",
                fill="none",
                opacity=0.7,
                marker={"size": 10, "line": {"width": 0.5, "color": "white"}},
                name="PAPA",
                yaxis="y1",
            )
            for i in df.Sequence_ID.unique()
        ]
        sub_plot1 = [
            go.Scatter(
                x=df[df["Sequence_ID"] == i]["Sequence Position"],
                y=df[df["Sequence_ID"] == i]["PRIMA score"],
                text=df[df["Sequence_ID"] == i]["Amino Acid"],
                mode="lines+markers",
                fill="none",
                opacity=0.7,
                marker={"size": 10, "line": {"width": 0.5, "color": "white"}},
                name="PRIMA",
                yaxis="y2",
            )
            for i in df.Sequence_ID.unique()
        ]
        sub_plot2 = [
            go.Bar(
                x=df[df["Sequence_ID"] == i]["Sequence Position"],
                y=df[df["Sequence_ID"] == i]["FoldIndex score"],
                text=df[df["Sequence_ID"] == i]["Amino Acid"],
                opacity=0.7,
                marker_color="red",
                name="FoldIndex",
                yaxis="y3",
            )
            for i in df.Sequence_ID.unique()
        ]
        main_plot.extend(sub_plot1)
        main_plot.extend(sub_plot2)

        fig = {
            "data": main_plot,
            "layout": go.Layout(
                title="Measures of Potential Prion Activity",
                xaxis={"title": "Amino Acid Position in Sequence"},
                yaxis=dict(
                    title="PAPA",
                    range=papa_y_range,
                    titlefont=dict(color="#1f77b4"),
                    tickfont=dict(color="#1f77b4"),
                ),
                yaxis2=dict(
                    title="PRIMA",
                    range=prima_y_range,
                    titlefont=dict(color="#ff7f0e"),
                    tickfont=dict(color="#ff7f0e"),
                    anchor="free",
                    overlaying="y",
                    side="left",
                    position=0.20,
                ),
                yaxis3=dict(
                    title="FoldIndex",
                    range=fold_index_y_range,
                    titlefont=dict(color="#FF0000"),
                    tickfont=dict(color="#FF0000"),
                    anchor="x",
                    overlaying="y",
                    side="right",
                ),
                legend={"x": -0.25, "y": 1},
                hovermode="closest",
            ),
        }
        return dcc.Graph(figure=fig, id=ids.SCATTER_PLOT)

    return html.Div(id=ids.SCATTER_PLOT)