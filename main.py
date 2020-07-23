"""
Author: William Rivera
Project: Prion Maintenance Collaborative Project

This is a webapp that analyzes amino acid sequences using three different algorithms.
Prion Aggregation Propensity Algorithm (PAPA) - measures the propensity of an amino acid sequence to assume a prion-like configuration
FoldIndex - measures the intrinsic foldedness of an amino acid sequence
Prion Maintenance Algorithm (PRIMA) - measures likelihood that a prion will 'maintain' or transfect another prion.

dashappserver.py - This file is responsible for launching and rendering our webapp.
"""

import dash  # import dash web app platform, allows for Dash app creation
import dash_core_components as dcc  # extends Dash capability for data visualizatino
import dash_html_components as html  # allows for rendering of HTML using Python code
import plotly.graph_objs as go  # allows for use of Plotly library plotting functionality
import pandas as pd  # used for data access and manipulation
import numpy as np  # ditto


# Amino Acid Attributes
# **********************************************************************************************************************
# **********************************************************************************************************************

# All possible amino acid residues
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# PAPA scores for all individual amino acids
propensities = {'A': -0.396490246, 'C': 0.415164505, 'D': -1.276997939, 'E': -0.605023827, 'F': 0.838732498,
                'G': -0.039220713, 'H': -0.278573356, 'I': 0.813697862, 'K': -1.576748587, 'L': -0.040005335,
                'M': 0.673729095, 'N': 0.080295334, 'P': -1.197447496, 'Q': 0.069168387, 'R': -0.405858577,
                'S': 0.133912418, 'T': -0.11457038, 'V': 0.813697862, 'W': 0.666735081, 'Y': 0.77865336}

# Hydrophobicity scores for all individual amino acids
hydrophobicity = {'A': 0.7, 'C': 0.777777778, 'D': 0.111111111, 'E': 0.111111111, 'F': 0.811111111,
                  'G': 0.455555556, 'H': 0.144444444, 'I': 1, 'K': 0.066666667, 'L': 0.922222222,
                  'M': 0.711111111, 'N': 0.111111111, 'P': 0.322222222, 'Q': 0.111111111, 'R': 0,
                  'S': 0.411111111, 'T': 0.422222222, 'V': 0.966666667, 'W': 0.4, 'Y': 0.355555556}

# Amino acid charge states
charge = {'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0,
          'G': 0, 'H': 0, 'I': 0, 'K': 1, 'L': 0,
          'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 1,
          'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

# Maintenance scores for all individual amino acids
maintenance = {'A': 0.28, 'C': 0.45, 'D': -0.19, 'E': 0.99, 'F': 0.33,
               'G': 0.047, 'H': -0.064, 'I': -0.57, 'K': 0, 'L': -0.48,
               'M': -1.80, 'N': 0.18, 'P': 0.065, 'Q': 0.11, 'R': -0.88,
               'S': 0.18, 'T': -0.059, 'V': -0.35, 'W': 1.40, 'Y': 1.03}

# default window size
window_size = 41
half_window_size = int(window_size / 2)

sequence_labels = ['Ure2', 'Sup35', 'Rnq1', 'New1', 'Nsp1', 'Puf2', 'Pub1', 'Swi1', 'Ybr016w', 'Cbk1', 'Lsm4',
                   'Ybl081w', 'Ksp1', 'Asm4', 'Gln3', 'Ypr022c', 'Rlm1', 'Nrp1']

sample_sequences = ['MMNNNGNQVSNLSNALRQVNIGNRNSNTTTDQSNINFEFSTGVNNNNNNNSSSNNNNVQNNNSGRNGSQNNDNENNIKNTLEQHRQQQQ',
                    'MSDSNQGNNQQNYQQYSQNGNQQQGNNRYQGYQAYNAQAQPAGGYYQNYQGYSGYQQGGYQQYNPDAGYQQQYNPQGGYQQYNPQGGYQQQFNPQGGRGNYKNFNYNNNLQGYQAGFQPQSQGMSLNDFQKQQKQ',
                    'SGSGGGSQSMGASGLAALASQFFKSGNNSQGQGQGQGQGQGQGQGQGQGSFTALASLASSFMNSNNNNQQGQNQSSGGSSFGALASMASSFMHSNNNQNSNNSQQGYNQSYQNGNQNSQGYNNQQYQGGNGGYQQQQGQSGGAFSSLASMAQSYLGGGQTQSNQQQYNQQGQNNQQQYQQQGQNYQHQQQGQQQQQGHSSSFSALASMASSYLGNNSNSNSSYGGQQQANEYGRPQQNGQQQSNEYGRPQYGGNQNSNGQHESFNFSGNFSQQNNNGNQNRY',
                    'GSNNASKKSSYQQQRNWKQGGNYQQGGYQSYNSNYNNYNNYNNYNNYNNYNNYNKYNGQGYQKSTYKQSAVTPNQSG',
                    'MNFNTPQQNKTPFSFGTANNNSNTTNQNSSTGAGAFGTGQSTFGFNNSAPNNTNNANSSITPAFGSNNTGNTAFGNSNPTSNVFGSNNSTTNTFGSNSAGTSLFGSSSAQQTKSNGTAGGNTFGSSSLFNNSTNSNTTKPAFGGLNFGGGNNTTPSSTGNANTSNNLFGATANAN',
                    'NSYFNNQQVVYSGNQNQNQNGNSNGLDELNSQFDSFRIANGTNLSLPIVNLPNVSNNNNNYNNSGYSSQMNPLSRSVSHNNNNNTNNYNNNDNDNNNNNNNNNNNNNNNNNNNNNSNNSNNNNNNDTSLYRYRSYGY',
                    'NNNNNNYQQRRNYGNNNRGGFRQYNSNNNNNMNMGMNMNMNMNMNNSRGMPPSSMGMPIGAMPLPSQGQPQQSQTIGLPPQVNPQ',
                    'MDFFNLNNNNNNNNTTTTTTTTNNNNTNNNNTNNNNNPANNTNNNNSTGHSSNTNNNTNNNNTNTGASGVDDFQNFFDPKPFDQNLDSNNNNSNSNNNDNNNSNTVASSTNFTSPTAVVNNAAPANVTGGKAANFIQNQSPQFNSPYDSNNSNTNLNSLSPQAILAKNSIIDSSNLPLQAQQQLYGGNNNNNSTGIANDNVITPHFITNVQSISQNSSSSTPNTNSNSTPNANQQFLPFNNSASNNGNLTSNQLISNYAASNSMDRSSSASNEFVPNTSDNNNNSNNHNMRNNSNNKTSNNNNVTAVPAATPANTNNSTSNANTVFSERAAMFAALQQKQQQRFQALQQQQQQQQNQQQQNQQPQQQQQQQQNPKFLQSQRQQQQRSILQSLNPA',
                    'MSANDYYGGTAGEKSQYSRPSNPPPSSAHQNKTQERGYPPQQQQQYYQQQQQHPGYYNQQGYNQQGYNQQGYNQQGYNQQGYNQQGYNQQGHQQPVYVQQQPPQR',
                    'MYNSSTNHHEGAPTSGHGYYMSQQQDQQHQQQQQYANEMNPYQQIPRPPAAGFSSNYMKEQGSHQSLQEHLQRETGNLGSGFTDVPALNYPATPPPHNNYAASNQMINTPPPSMGGLYRHNNNSQSMVQNGNGSGNAQLPQLSPGQYSIESEYNQNLNGSSSSSPFHQPQTLRSNGSYSSGLRSVKSFQRLQQEQENVQVQQQLSQAQQQNSRQQQQQLQYQQQQQQQQQQQHMQIQQQQQQQQQQQQSQSPVQSGFNNGTISNYMYFER',
                    'SNNNSNSNGPGHKRYYNNRDSNNNRGNYNRRNNNNGNSNRRPYSQNRQYNNSNSSNINNSINSINSNNQNMNNGLGGSVQHHFNSSSPQKVEF',
                    'NAPSHQSNYHPHYNHMKYNNTGSYYYYNNNNNSSVNPHNQAGLQSINRSIPSAPYGAYNQNRANDVPYMNTQKKHHRFSANNNLNQQKYKQYPQYTSNPMVTAHLKQTYPQLYYNSNVNAHNNNNNSNNNNNNNNNSNNNNNLYNQTQFSTRYFNSNSSPSLTSSTSNSSSPYNQSTFEYIL',
                    'PSVQHRYMEGFSNNNNKQYRQNRNYNNNNNNSNNNHGSNYNNFNNGNSYIKGWNKNFNKYRRPSSSSYTGKSPLSRYNMSYNHNNNSSINGY',
                    'MFGIRSGNNNGGFTNLTSQAPQTTQMEQSQSQLQPQPQPQPQQQQQHLQFNGSSDASSLRFGNSLSNTVNANNYSSNIGNNSINNNNIKNGTNNISQHGQGNNPSW',
                    'ILPKPSPNSANSQQFNMNMNLMNTTNNVSAGNSVASSPRIISSANFNSNSPLQQNLLSNSFQRQGMNIPRRKMSRNASYSSSFMAASLQQLHEQQQVDVNSNTNTNSNRQNWNSSNSVSTNSRSSNFVSQ',
                    'PPPSITQESKFRPFLQQAQQPQQVQQSQQPQQIQQLQQLQFPQQLRAPLQQPMLQQQMHPQQASPTFPSYDPRIRNNGQNGNQFFNLIFDNRTGVNGFEVDAANNNGNGNDQNMNINPAVQQQRYQDRNFASSSYQQPLQPLTQDQQQEQYFQQQKLAQQQQQQQQQQQQQQQLPPQNPFGDPLTSSSSGANLSV',
                    'QKGNNGRMVIKLPNANAPNGSNNGNGSNNNNHPYPFGSGSSPLFSATQPYIATPLQPSNIPGGPFQQNTSFLAQRQTQQYQQMSFKKQSQTVPL',
                    'PIIIADHFSGNNNIAPNYRYNNNINNNNNNINNMTNNRYNINNNINGNGNGNGNNSNNNNNHNNNHNNNHHNGSINSNSNTNNNNNNNNGNNSNNCNSNIGMGGCGSNM'
                    ]

aminoacidDict = dict(zip(sequence_labels, sample_sequences))

# sequence_id = ['']
# sequence = 'MSDSNQGNNQQNYQQYSQNGNQQQGNNRYQGYQAYNAQAQPAGGYYQNYQGYSGYQQGGYQQYNPDAGYQQQYNPQGGYQQYNPQGGYQQQFNPQGGRGNYKNFNYNNNLQGYQAGFQPQSQGMSLNDFQKQQKQ'
default_seq = sample_sequences[-3]
default_seq_id = sequence_labels[-3]
# **********************************************************************************************************************
# **********************************************************************************************************************


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
        window_scores: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        fold_index_scores: Ignore FoldIndex scores when performing analysis.

    Returns: list of weighted window scores of a particular algorithm in the given sequence.
    """
    scores = []
    for i in range(len(sequence) - window_size):
        if (fold_index_scores is not None and fold_index_scores[i] > 0):
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


"""
FoldIndex formula:
2.785(H) - |R| - 1.151
H = sum of the hydrophobicities across the window
R = net charge (where D/E=-1; K/R=+1; all others = neutral (including H))
"""

def fold_index(sequence):
    """
    Calculates FoldIndex scores for all windows in a given sequence.

    Parameters:
        sequence: An amino acid sequence
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines:

    Returns: list of FoldIndex scores for all windows in the given sequence.
    """
    charges = window_scores(sequence, charge)
    hydrophobicities = window_scores(sequence, hydrophobicity)
    fold_index_list = [2.785 * hydrophobicities[i] - abs(charges[i]) - 1.151 for i in range(len(sequence))]
    return super_window_scores(sequence, fold_index_list)


def prima_score(sequence):
    """
    Calculates PRIMA scores for all windows in a given sequence.

    Parameters:
        sequence: An amino acid sequence
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines:

    Returns: list of PRIMA scores for all windows in the given sequence.
    """
    maintenance_scores = window_scores(sequence, maintenance)
    prima_score_list = [maintenance_scores[i] for i in range(len(sequence))]
    return super_window_scores(sequence, prima_score_list)


def classify(sequence, ignore_fold_index=True):
    fold_index_list = fold_index(sequence)
    # print 'number of positions with negative fold index', sum([1 for f in fold_index_list if f < 0]), 'out of', len(sequence)

    prima_score_list = prima_score(sequence)

    window_propensities = window_scores(sequence, propensities, True)
    if ignore_fold_index:
        scores = super_window_scores(sequence, window_propensities)
    else:
        scores = super_window_scores(sequence, window_propensities, fold_index_list)
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


def analyze(sequence, sequence_id):
    sequence_id = [sequence_id]
    score, pos, scores, fold_indexes, prima_scores, prima_score, prima_position = classify(sequence)
    n_columns = 4
    arr = np.empty(((len(sequence) - window_size), n_columns))
    for i, s in enumerate(scores):
        arr[i, :] = [i + 1, s, prima_scores[i], fold_indexes[i]]
    df = pd.DataFrame(data=arr,
                      columns=['Sequence Position', 'PAPA score', 'PRIMA score', 'FoldIndex score']
                      )

    temp_list = list(sequence)
    aa_arr = temp_list[:len(sequence) - window_size]
    df2 = pd.DataFrame(data=aa_arr,
                       columns=['Amino Acid'])
    seq_id_column = sequence_id * (len(sequence) - window_size)
    seq_df = pd.DataFrame(data=seq_id_column,
                          columns=['Sequence_ID'])
    df3 = pd.concat([seq_df, df2, df], axis=1)

    return (df3)

def setrange(lows, highs):
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

min_scores = df[['PAPA score', 'PRIMA score', 'FoldIndex score']].min().tolist()
max_scores = df[['PAPA score', 'PRIMA score', 'FoldIndex score']].max().tolist()





score_ranges = setrange(min_scores, max_scores)
papa_y_range = score_ranges[0]
prima_y_range = score_ranges[1]
fold_index_y_range = score_ranges[2]
print(score_ranges)

sc1 = [go.Scatter(
    x=df[df['Sequence_ID'] == i]['Sequence Position'],
    y=df[df['Sequence_ID'] == i]['PAPA score'],
    text=df[df['Sequence_ID'] == i]['Amino Acid'],
    mode='lines+markers',
    fill='tozeroy',
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
sc3 = [go.Scatter(
    x=df[df['Sequence_ID'] == i]['Sequence Position'],
    y=df[df['Sequence_ID'] == i]['FoldIndex score'],
    text=df[df['Sequence_ID'] == i]['Amino Acid'],
    mode='lines+markers',
    fill='tozeroy',
    opacity=0.7,
    marker={
        'size': 10,
        'line': {'width': 0.5, 'color': 'white'}
    },
    name='FoldIndex',
    yaxis='y3'
) for i in df.Sequence_ID.unique()]
sc1.extend(sc2)
sc1.extend(sc3)

app.layout = html.Div(children=[
    html.H1(children='Prion Maintenance Collaborative Project (PRIMCOP)'),

    html.Div(children='''
        A Dash-powered WebApp allowing for the prediction of the intrinsic foldedness (FoldIndex), prion aggregation propensity (PAPA), and prion maintenance propensity(PRIMA) of an amino acid sequence.
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


@app.callback(
    dash.dependencies.Output('prion-analysis', 'figure'),
    [dash.dependencies.Input('my-dropdown', 'value')])
def update_figure(value):
    seq_id = [seq_id for seq_id, seq_val in aminoacidDict.items() if seq_val == value]
    upd_sequence_id = seq_id
    upd_sequence = value
    df = analyze(upd_sequence, upd_sequence_id)

    min_scores = df[['PAPA score', 'PRIMA score', 'FoldIndex score']].min().tolist()
    max_scores = df[['PAPA score', 'PRIMA score', 'FoldIndex score']].max().tolist()
    score_ranges = setrange(min_scores, max_scores)
    papa_y_range = score_ranges[0]
    prima_y_range = score_ranges[1]
    fold_index_y_range = score_ranges[2]
    print(score_ranges)

    #    print(dff)

    traces = [go.Scatter(
        x=df[df['Sequence_ID'] == i]['Sequence Position'],
        y=df[df['Sequence_ID'] == i]['PAPA score'],
        text=df[df['Sequence_ID'] == i]['Amino Acid'],
        mode='lines+markers',
        fill='tozeroy',
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
    traces3 = [go.Scatter(
        x=df[df['Sequence_ID'] == i]['Sequence Position'],
        y=df[df['Sequence_ID'] == i]['FoldIndex score'],
        text=df[df['Sequence_ID'] == i]['Amino Acid'],
        mode='lines+markers',
        fill='tozeroy',
        opacity=0.7,
        marker={
            'size': 10,
            'line': {'width': 0.5, 'color': 'white'}
        },
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
