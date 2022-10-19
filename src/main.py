#!/usr/bin/env python

"""
Analyze amino acid sequences for prion activity using three different, sliding window-based algorithms.

This webapp is an outgrowth of the Prion Maintenance Collaborative Project (PRIMCOP) which aims to further our
understanding of prion proteins, using wet lab and computational approaches. The algorithms implemented are:
    Prion Aggregation Propensity Algorithm (PAPA) - measures the propensity of an amino acid sequence to assume a prion-like configuration
    FoldIndex - measures the intrinsic foldedness of an amino acid sequence
    Prion Maintenance Algorithm (PRIMA) - measures likelihood that a prion will 'maintain' or transfect another prion.
"""

import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd
import plotly.graph_objs as go

# from prions import *  # import dictionaries with prion characteristic data for individual amino acids
from windows import *  # import functions for calculating window scores

import components

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
        fold_index_list,
        prima_score_list,
        papa_score_list,
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


app = dash.Dash()  # creates Dash app
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
            figure=components.get_analysis_figure(run_analysis(default_seq, default_seq_id))
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
    return components.get_analysis_figure(df)


if __name__ == "__main__":
    app.run_server()
