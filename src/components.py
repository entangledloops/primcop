import plotly.graph_objs as go


AMINO_ACID = "Amino Acid"
FOLD_INDEX = "FoldIndex"
FOLD_INDEX_SCORE = "FoldIndex score"
PAPA = "PAPA"
PAPA_SCORE = "PAPA score"
PRIMA = "PRIMA"
PRIMA_SCORE = "PRIMA score"
SEQUENCE_ID = "Sequence_ID"
SEQUENCE_POSITION = "Sequence Position"


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
        # mode='lines+markers',
        # fill='none',
        # marker={
        #    'size': 10,
        #    'line': {'width': 0.5, 'color': 'white'}
        # },
    )

def get_analysis_components(df) -> list:
    papas = []
    primas = []
    folds = []
    for seq_id in df[SEQUENCE_ID].unique():
        # select the relevant subset of the dataframe
        seq = df[df[SEQUENCE_ID] == seq_id]
        pos = seq[SEQUENCE_POSITION]
        amino_acid = seq[AMINO_ACID]

        # PAPA
        papa_scatter = get_scatter(
            x=pos,
            y=seq[PAPA_SCORE],
            text=amino_acid,
            name=PAPA,
            yaxis="y1",
        )
        papas.append(papa_scatter)

        # PRIMA
        prima_scatter = get_scatter(
            x=pos,
            y=seq[PRIMA_SCORE],
            text=amino_acid,
            name=PRIMA,
            yaxis="y2",
        )
        primas.append(prima_scatter)

        # Fold Index
        fold_index_bar = get_bar(
            x=pos,
            y=seq[FOLD_INDEX_SCORE],
            text=amino_acid,
            name=FOLD_INDEX,
            yaxis="y3",
        )
        folds.append(fold_index_bar)
    return [*papas, *primas, *folds]


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


def get_analysis_figure(df):
    scores = df[[PAPA_SCORE, PRIMA_SCORE, FOLD_INDEX_SCORE]]
    lower_bounds = scores.min().tolist()
    upper_bounds = scores.max().tolist()
    score_ranges = setrange(lower_bounds, upper_bounds)
    papa_y_range = score_ranges[0]
    prima_y_range = score_ranges[1]
    fold_index_y_range = score_ranges[2]
    traces = get_analysis_components(df)
    return {
        "data": traces,
        "layout": go.Layout(
            title="Prion Analysis Curves",
            xaxis={"title": "Amino Acid Position in Sequence"},
            yaxis=dict(
                title=PAPA,
                range=papa_y_range,
                titlefont=dict(color="#1f77b4"),
                tickfont=dict(color="#1f77b4"),
            ),
            yaxis2=dict(
                title=PRIMA,
                range=prima_y_range,
                titlefont=dict(color="#ff7f0e"),
                tickfont=dict(color="#ff7f0e"),
                anchor="free",
                overlaying="y",
                side="left",
                position=0.20,
            ),
            yaxis3=dict(
                title=FOLD_INDEX,
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