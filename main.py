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
__credits__ = ["William Rivera", "Kyle Maclea",
               "Raghava Adusimilli", "Stephen Dunn"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "William Rivera"
__email__ = "wto3@wildcats.unh.edu"
__status__ = "Development"

from prions import *  # import dictionaries with prion characteristic data for individual amino acids
import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd
import plotly.graph_objs as go

# Default window size
window_size = 41
half_window_size = int(window_size / 2)

aminoacidDict = dict(zip(sample_labels, sample_sequences))

default_seq = sample_sequences[-3]
default_seq_id = sample_labels[-3]

# Function Declarations


def get_window(sequence, position):
    """
    Determines sliding window position and contents.

    Parameters:
        sequence: An amino acid sequence
        position: Starting or left-most index of window along sequence

    Returns: A window (i.e., sequence slice) defined by its starting and ending position in the sequence.
    """
    start = max(position - half_window_size, 0)
    end = min(len(sequence), position + half_window_size + 1)
    return start, end


def window_score(sequence, position, aa_dict, ignore_consecutive_prolines=False):
    """
    Performs 1 of 3 algorithmic analyses (FoldIndex, PAPA, PRIMA) on a sequence slice, calculating a single window score.

    Parameters:
        sequence: An amino acid sequence
        position: ith or leftmost position of sliding window
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines:

    Returns: calculated window score for a window in a given amino acid sequence.
    """

    start, end = get_window(sequence, position)
    score = 0.0
    for i in range(start, end):
        if sequence[i] not in amino_acids:
            continue
        if not ignore_consecutive_prolines:
            score += aa_dict[sequence[i]]
        else:
            if sequence[i] != 'P':
                score += aa_dict[sequence[i]]
            elif (i > 0 and sequence[i - 1] == 'P') or (i > 1 and sequence[i - 2] == 'P'):
                pass
            else:
                score += aa_dict[sequence[i]]
    return score / (end - start)


def window_scores(sequence, aa_dict, ignore_consecutive_prolines=False):
    """
    Performs 1 of 3 algorithmic analyses (FoldIndex, PAPA, PRIMA) on a sequence and generates list of all window scores.

    Parameters:
        sequence: An amino acid sequence
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines: If consecutive residues in a sequence are proline (P), do not apply algorithm.

    Returns: a list containing calculated window scores for all windows in a given amino acid sequence.
    """

    return [window_score(sequence, i, aa_dict, ignore_consecutive_prolines) for i in range(len(sequence))]


def super_window_scores(sequence, window_scores, fold_index_scores=None):
    """
    Takes weighted average of window scores for a given sequence.

    Parameters:
        sequence: An amino acid sequence
        window_scores: The dictionary containing the log-odds value for the algorithm being performed on sequence.
        fold_index_scores: Ignore FoldIndex scores when performing analysis.

    Returns: list of weighted window scores of a particular algorithm in the given sequence.
    """
    scores = []
    for i in range(len(sequence) - window_size):
        if fold_index_scores is not None and fold_index_scores[i] > 0:
            scores.append(None)
        else:
            score = 0.0
            weights = 0.0
            for j in range(i, i + window_size):
                start, end = get_window(sequence, j)
                score += (end - start) * window_scores[j]
                weights += (end - start)
            scores.append(score / weights)
    return scores


def fold_index(sequence):
    """
    Calculates FoldIndex scores for all windows in a given sequence.

    Parameters:
        sequence: An amino acid sequence
    Returns: list of FoldIndex scores for all windows in the given sequence.
    """
    charges = window_scores(sequence, charge)
    hydrophobicities = window_scores(sequence, hydrophobicity)
    fold_index_list = [2.785 * hydrophobicities[i] -
                       abs(charges[i]) - 1.151 for i in range(len(sequence))]
    return super_window_scores(sequence, fold_index_list)


def prima_score(sequence):
    """
    Calculates PRIMA scores for all windows in a given sequence.

    Parameters:
        sequence: An amino acid sequence
    Returns: list of PRIMA scores for all windows in the given sequence.
    """
    maintenance_scores = window_scores(sequence, maintenance)
    prima_score_list = [maintenance_scores[i] for i in range(len(sequence))]
    return super_window_scores(sequence, prima_score_list)


def classify(sequence, ignore_fold_index=True):
    """
    Doc String Placeholder
    """

    fold_index_list = fold_index(sequence)
    prima_score_list = prima_score(sequence)

    window_propensities = window_scores(sequence, propensities, True)
    if ignore_fold_index:
        scores = super_window_scores(sequence, window_propensities)
    else:
        scores = super_window_scores(
            sequence, window_propensities, fold_index_list)
    max_score = max(scores)
    max_position = scores.index(max_score)
    # the case when no window had a negative foldIndex
    if max_score is None:
        max_score = -1.0
        max_position = -1

    max_prima_score = max(prima_score_list)
    max_prima_position = prima_score_list.index(max_prima_score)
    # the case when no window had a negative foldIndex
    if max_prima_score is None:
        max_prima_score = -1.0
        max_prima_position = -1

    return max_score, max_position, scores, fold_index_list, prima_score_list, max_prima_score, max_prima_position


# TODO Clean up documentation and rename this function


def analyze(sequence, sequence_id):
    """
    Doc String Placeholder
    """
    sequence_id = [sequence_id]
    # ? Why are score, pos, prima_score and prima_position unused?
    score, pos, scores, fold_indexes, prima_scores, prima_score, prima_position = classify(
        sequence)
    n_columns = 4
    arr = np.empty(((len(sequence) - window_size), n_columns))
    for i, s in enumerate(scores):
        arr[i, :] = [i + 1, s, prima_scores[i], fold_indexes[i]]
    df = pd.DataFrame(data=arr,
                      columns=['Sequence Position', 'PAPA score',
                               'PRIMA score', 'FoldIndex score']
                      )

    temp_list = list(sequence)
    aa_arr = temp_list[:len(sequence) - window_size]
    df2 = pd.DataFrame(data=aa_arr,
                       columns=['Amino Acid'])
    seq_id_column = sequence_id * (len(sequence) - window_size)
    seq_df = pd.DataFrame(data=seq_id_column,
                          columns=['Sequence_ID'])
    df3 = pd.concat([seq_df, df2, df], axis=1)

    return df3


# TODO Attempt to refactor into lambda function

def setrange(lows, highs):
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

df = analyze(default_seq, default_seq_id)

min_scores = df[['PAPA score', 'PRIMA score',
                 'FoldIndex score']].min().tolist()
max_scores = df[['PAPA score', 'PRIMA score',
                 'FoldIndex score']].max().tolist()


score_ranges = setrange(min_scores, max_scores)
papa_y_range = score_ranges[0]
prima_y_range = score_ranges[1]
fold_index_y_range = score_ranges[2]

sc1 = [go.Scatter(
    x=df[df['Sequence_ID'] == i]['Sequence Position'],
    y=df[df['Sequence_ID'] == i]['PAPA score'],
    text=df[df['Sequence_ID'] == i]['Amino Acid'],
    mode='lines+markers',
    fill='none',
    opacity=0.7,
    marker={
        'size': 10,
        'line': {'width': 0.5, 'color': 'white'}
    },
    name='PAPA',
    yaxis='y1'
) for i in df.Sequence_ID.unique()]
sc2 = [go.Scatter(
    x=df[df['Sequence_ID'] == i]['Sequence Position'],
    y=df[df['Sequence_ID'] == i]['PRIMA score'],
    text=df[df['Sequence_ID'] == i]['Amino Acid'],
    mode='lines+markers',
    fill='none',
    opacity=0.7,
    marker={
        'size': 10,
        'line': {'width': 0.5, 'color': 'white'}
    },
    name='PRIMA',
    yaxis='y2'
) for i in df.Sequence_ID.unique()]
sc3 = [go.Bar(
    x=df[df['Sequence_ID'] == i]['Sequence Position'],
    y=df[df['Sequence_ID'] == i]['FoldIndex score'],
    text=df[df['Sequence_ID'] == i]['Amino Acid'],
    # mode='lines+markers',
    # fill='none',
    opacity=0.7,
    # marker={
    #    'size': 10,
    #    'line': {'width': 0.5, 'color': 'white'}
    # },
    name='FoldIndex',
    yaxis='y3'
) for i in df.Sequence_ID.unique()]
sc1.extend(sc2)
sc1.extend(sc3)

app.layout = html.Div(children=[
    html.H1(children='Prion Maintenance Collaborative Project (PRIMCOP)'),

    html.Div(children='''
        A Dash-powered WebApp allowing for the prediction of the intrinsic foldedness (FoldIndex),
        prion aggregation propensity (PAPA), and prion maintenance propensity (PRIMA) of an amino acid sequence.
    '''),

    html.H4(children='Sample Sequences'),
    dcc.Dropdown(
        id='my-dropdown',
        placeholder='Please select sample protein sequence from the Alberti et al. dataset below. Ure2 is currently selected.',
        options=[{'label': k, 'value': v} for k, v in aminoacidDict.items()],
        value='placeholder'
    ),
    html.Div(id='output-container'),
    dcc.Graph(
        id='prion-analysis',
        figure={
            'data': sc1,
            'layout': go.Layout(
                title='Prion Analysis Curves',
                xaxis={'title': 'Amino Acid Position in Sequence'},
                yaxis=dict(
                    title='PAPA',
                    range=papa_y_range,
                    titlefont=dict(
                        color='#1f77b4'
                    ),
                    tickfont=dict(
                        color='#1f77b4'
                    )
                ),
                yaxis2=dict(
                    title='PRIMA',
                    range=prima_y_range,
                    titlefont=dict(
                        color='#ff7f0e'
                    ),
                    tickfont=dict(
                        color='#ff7f0e'
                    ),
                    anchor='free',
                    overlaying='y',
                    side='left',
                    position=0.50
                ),
                yaxis3=dict(
                    title='FoldIndex',
                    range=fold_index_y_range,
                    titlefont=dict(
                        color='#6abc6a'
                    ),
                    tickfont=dict(
                        color='#6abc6a'
                    ),
                    anchor='x',
                    overlaying='y',
                    side='right'
                ),
                # margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
                legend={'x': -0.1, 'y': 1},
                hovermode='closest'
            )
        }
    )
])


# TODO Find way to eliminate redundancy in trace scatter plot generation
@app.callback(
    dash.dependencies.Output('prion-analysis', 'figure'),
    [dash.dependencies.Input('my-dropdown', 'value')])
def update_figure(value):
    seq_id = [seq_id for seq_id, seq_val in aminoacidDict.items()
              if seq_val == value]
    upd_sequence_id = seq_id
    upd_sequence = value
    df = analyze(upd_sequence, upd_sequence_id)

    min_scores = df[['PAPA score', 'PRIMA score',
                     'FoldIndex score']].min().tolist()
    max_scores = df[['PAPA score', 'PRIMA score',
                     'FoldIndex score']].max().tolist()
    score_ranges = setrange(min_scores, max_scores)
    papa_y_range = score_ranges[0]
    prima_y_range = score_ranges[1]
    fold_index_y_range = score_ranges[2]

    traces = [go.Scatter(
        x=df[df['Sequence_ID'] == i]['Sequence Position'],
        y=df[df['Sequence_ID'] == i]['PAPA score'],
        text=df[df['Sequence_ID'] == i]['Amino Acid'],
        mode='lines+markers',
        fill='none',
        opacity=0.7,
        marker={
            'size': 10,
            'line': {'width': 0.5, 'color': 'white'}
        },
        name='PAPA',
        yaxis='y1'
    ) for i in df.Sequence_ID.unique()]
    traces2 = [go.Scatter(
        x=df[df['Sequence_ID'] == i]['Sequence Position'],
        y=df[df['Sequence_ID'] == i]['PRIMA score'],
        text=df[df['Sequence_ID'] == i]['Amino Acid'],
        mode='lines+markers',
        fill='none',
        opacity=0.7,
        marker={
            'size': 10,
            'line': {'width': 0.5, 'color': 'white'}
        },
        name='PRIMA',
        yaxis='y2'
    ) for i in df.Sequence_ID.unique()]
    traces3 = [go.Bar(
        # TODO change bar color to red
        x=df[df['Sequence_ID'] == i]['Sequence Position'],
        y=df[df['Sequence_ID'] == i]['FoldIndex score'],
        text=df[df['Sequence_ID'] == i]['Amino Acid'],
        # mode='lines+markers',
        # fill='none',
        opacity=0.7,
        # marker={
        #    'size': 10,
        #    'line': {'width': 0.5, 'color': 'white'}
        # },
        name='FoldIndex',
        yaxis='y3'
    ) for i in df.Sequence_ID.unique()]
    traces.extend(traces2)
    traces.extend(traces3)

    return {
        'data': traces,
        'layout': go.Layout(
            title='Prion Analysis Curves',
            xaxis={'title': 'Amino Acid Position in Sequence'},
            yaxis=dict(
                title='PAPA',
                range=papa_y_range,
                titlefont=dict(
                    color='#1f77b4'
                ),
                tickfont=dict(
                    color='#1f77b4'
                )
            ),
            yaxis2=dict(
                title='PRIMA',
                range=prima_y_range,
                titlefont=dict(
                    color='#ff7f0e'
                ),
                tickfont=dict(
                    color='#ff7f0e'
                ),
                anchor='free',
                overlaying='y',
                side='left',
                position=0.20
            ),
            yaxis3=dict(
                title='FoldIndex',
                range=fold_index_y_range,
                titlefont=dict(
                    color='#6abc6a'
                ),
                tickfont=dict(
                    color='#6abc6a'
                ),
                anchor='x',
                overlaying='y',
                side='right'
            ),
            # margin={'l': 40, 'b': 40, 't': 10, 'r': 10},
            legend={'x': -0.1, 'y': 1},
            hovermode='closest'
        )
    }


if __name__ == '__main__':
    app.run_server()
