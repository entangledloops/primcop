#!/usr/bin/env python

"""
Analyze amino acid sequences for prion activity using three different, sliding window-based algorithms.

This webapp is an outgrowth of the Prion Maintenance Collaborative Project (PRIMCOP) which aims to further our
understanding of prion proteins, using wet lab and computational approaches. The algorithms implemented are:
    Prion Aggregation Propensity Algorithm (PAPA) - measures the propensity of an amino acid sequence to assume a prion-like configuration
    FoldIndex - measures the intrinsic foldedness of an amino acid sequence
    Prion Maintenance Algorithm (PRIMA) - measures likelihood that a prion will 'maintain' or transfect another prion.
"""

__author__ = "William Rivera"
__credits__ = ["William Rivera", "Kyle Maclea", "Raghava Adusimilli", "Stephen Dunn"]
__license__ = "MIT"
__version__ = "0.1.1"
__maintainer__ = "William Rivera"
__email__ = "wto3@wildcats.unh.edu"
__status__ = "Development"


import dash
import dash_core_components as dcc
import dash_html_components as html
import hydra
import numpy as np
import pandas as pd
import plotly.graph_objs as go

# from prions import *  # import dictionaries with prion characteristic data for individual amino acids
from windows import *  # import functions for calculating window scores

# Create dictionary of sample data for generating sample plots
samples = dict(zip(sample_labels, sample_sequences))


def get_fold_index_scores(sequence: str) -> list:
    """
    Calculates FoldIndex scores for all windows in a given sequence.
    Positive(+) FoldIndex scores represent proteins/domains that are likely folded.
    Negative(-) FoldIndex scores represent proteins/domains that are likely intrinsically unfolded.

    Parameters:
        sequence: An amino acid sequence
    Returns: list of FoldIndex scores for all windows in the given sequence.
    """
    charges = window_scores(sequence, charge)
    hydrophobicities = window_scores(sequence, hydrophobicity)
    fold_index_list = [
        2.785 * hydrophobicities[i] - abs(charges[i]) - 1.151
        for i in range(len(sequence))
    ]
    fold_index_list = super_window_scores(sequence, fold_index_list)
    min_fold_index_score = min(fold_index_list)
    min_fold_index_position = fold_index_list.index(min_fold_index_score)
    return (min_fold_index_score, min_fold_index_position, fold_index_list)


def get_papa_scores(sequence: str, ignore_fold_index: bool = True) -> list:
    """
    Calculates PAPA scores for all windows in a given sequence.

    Parameters:
        sequence: An amino acid sequence
    Returns: list of PAPA scores for all windows in the given sequence.
    """
    fold_index_scores = get_fold_index_scores(sequence)[2]
    papa_scores = window_scores(sequence, propensity, True)
    if ignore_fold_index:
        papa_score_list = super_window_scores(sequence, papa_scores)
        print(papa_score_list)
    else:
        papa_score_list = super_window_scores(
            sequence, papa_scores, fold_index_scores=fold_index_scores
        )
    max_papa_score = max(papa_score_list)
    max_papa_position = papa_score_list.index(max_papa_score)
    # the case when no window had a negative foldIndex
    if max_papa_score is None:
        max_papa_score = -1.0
        max_papa_position = -1

    return (max_papa_score, max_papa_position, papa_score_list)


def get_prima_scores(sequence: str) -> list:
    """
    Calculates PRIMA scores for all windows in a given sequence.
    Parameters:
        sequence: An amino acid sequence
    Returns: list of PRIMA scores for all windows in the given sequence.
    """
    maintenance_scores = window_scores(sequence, maintenance)
    prima_score_list = [maintenance_scores[i] for i in range(len(sequence))]
    max_prima_score = max(prima_score_list)
    max_prima_position = prima_score_list.index(max_prima_score)
    # ? the case when no window had a negative foldIndex
    # TODO: Confirm if this is necessary with Maclea
    if max_prima_score is None:
        max_prima_score = -1.0
        max_prima_position = -1
    return (
        max_prima_score,
        max_prima_position,
        super_window_scores(sequence, prima_score_list),
    )


def assess_sequence(sequence: str) -> tuple:
    """
    Scan sequence for potential prion activity, assigning window- and sequence-level scores 
    """

    fold_index_list = get_fold_index_scores(sequence)
    prima_score_list = get_prima_scores(sequence)
    papa_score_list = get_papa_scores(sequence)

    return (
        get_fold_index_scores(sequence),
        get_papa_scores(sequence),
        get_prima_scores(sequence),
    )


def run_analysis(sequence: str, sequence_id: str) -> pd.DataFrame:
    """
    Create dataframe report of prionic activity to inform plots
    """
    sequence_id = [sequence_id]
    # ? Why are score, pos, prima_score and prima_position unused?

    elements = assess_sequence(sequence)
    seq_papa_score, papa_pos, papa_scores = elements[1]
    seq_fold_index_score, fold_index_pos, fold_index_scores = elements[0]
    seq_prima_score, prima_pos, prima_scores = elements[2]

    array = np.empty(((len(sequence) - window_size), 4))
    for i, s in enumerate(papa_scores):
        array[i, :] = [i + 1, s, prima_scores[i], fold_index_scores[i]]
    df = pd.DataFrame(
        data=array,
        columns=["Sequence Position", "PAPA score", "PRIMA score", "FoldIndex score"],
    )

    temp_list = list(sequence)
    aa_arr = temp_list[: len(sequence) - window_size]
    df2 = pd.DataFrame(data=aa_arr, columns=["Amino Acid"])
    seq_id_column = sequence_id * (len(sequence) - window_size)
    seq_df = pd.DataFrame(data=seq_id_column, columns=["Sequence_ID"])
    df3 = pd.concat([seq_df, df2, df], axis=1)
    print(df3)

    return df3


# TODO Attempt to refactor into lambda function
def setrange(lows: list, highs: list) -> tuple:
    """
    Establishes range of value scale from -|max dev| to + |max dev|. Scales graphs correctly.
    """
    ranges = []
    for a, b in zip(lows, highs):
        if abs(a) > abs(b):
            temp = [-abs(a) - 0.05, abs(a) + 0.05]
            ranges.append(temp)
        else:
            temp = [-abs(b) - 0.05, abs(b) + 0.05]
            ranges.append(temp)
    return ranges


app = dash.Dash()  # creates Dash app

df = run_analysis(default_seq, default_seq_id)

lower_bounds = df[["PAPA score", "PRIMA score", "FoldIndex score"]].min().tolist()
upper_bounds = df[["PAPA score", "PRIMA score", "FoldIndex score"]].max().tolist()


score_ranges = setrange(lower_bounds, upper_bounds)
papa_y_range = score_ranges[0]
prima_y_range = score_ranges[1]
fold_index_y_range = score_ranges[2]

sc1 = [
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
sc2 = [
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
sc3 = [
    go.Bar(
        x=df[df["Sequence_ID"] == i]["Sequence Position"],
        y=df[df["Sequence_ID"] == i]["FoldIndex score"],
        text=df[df["Sequence_ID"] == i]["Amino Acid"],
        # mode='lines+markers',
        # fill='none',
        opacity=0.7,
        # marker={
        #    'size': 10,
        #    'line': {'width': 0.5, 'color': 'white'}
        # },
        name="FoldIndex",
        yaxis="y3",
    )
    for i in df.Sequence_ID.unique()
]
sc1.extend(sc2)
sc1.extend(sc3)

app.layout = html.Div(
    children=[
        html.H1(children="Prion Maintenance Collaborative Project (PRIMCOP)"),
        html.Div(
            children="""
        A Dash-powered WebApp allowing for the prediction of the intrinsic foldedness (FoldIndex),
        prion aggregation propensity (PAPA), and prion maintenance propensity (PRIMA) of an amino acid sequence.
    """
        ),
        html.H4(children="Sample Sequences"),
        dcc.Dropdown(
            id="my-dropdown",
            placeholder="Please select sample protein sequence from the Alberti et al. dataset below. Ure2 is currently selected.",
            options=[{"label": k, "value": v} for k, v in samples.items()],
            value="placeholder",
        ),
        html.Div(id="output-container"),
        dcc.Graph(
            id="prion-analysis",
            figure={
                "data": sc1,
                "layout": go.Layout(
                    title="Prion Analysis Curves",
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
                        position=0.50,
                    ),
                    yaxis3=dict(
                        title="FoldIndex",
                        range=fold_index_y_range,
                        titlefont=dict(color="#6abc6a"),
                        tickfont=dict(color="#6abc6a"),
                        anchor="x",
                        overlaying="y",
                        side="right",
                    ),
                    # margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                    legend={"x": -0.1, "y": 1},
                    hovermode="closest",
                ),
            },
        ),
    ]
)


# TODO Find way to eliminate redundancy in trace scatter plot generation
@app.callback(
    dash.dependencies.Output("prion-analysis", "figure"),
    [dash.dependencies.Input("my-dropdown", "value")],
)
def update_figure(value):
    if value == "placeholder":
        value = default_seq
    seq_id = [seq_id for seq_id, seq_val in samples.items() if seq_val == value]
    upd_sequence_id = seq_id
    upd_sequence = value
    print(f"Update Sequence: {upd_sequence}")
    df = run_analysis(upd_sequence, upd_sequence_id)

    min_scores = df[["PAPA score", "PRIMA score", "FoldIndex score"]].min().tolist()
    max_scores = df[["PAPA score", "PRIMA score", "FoldIndex score"]].max().tolist()
    score_ranges = setrange(min_scores, max_scores)
    papa_y_range = score_ranges[0]
    prima_y_range = score_ranges[1]
    fold_index_y_range = score_ranges[2]

    traces = [
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
    traces2 = [
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
    traces3 = [
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
    traces.extend(traces2)
    traces.extend(traces3)

    return {
        "data": traces,
        "layout": go.Layout(
            title="Prion Analysis Curves",
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
            # margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={"x": -0.1, "y": 1},
            hovermode="closest",
        ),
    }


if __name__ == "__main__":
    app.run_server()
