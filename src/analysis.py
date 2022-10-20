from src.characteristics import AMINO_ACIDS, PROPENSITY, HYDROPHOBICITY, CHARGE, MAINTENANCE
from src.sample_data import SAMPLE_LABELS, SAMPLE_SEQUENCES, SAMPLES
import numpy as np
import pandas as pd

# Constants
WINDOW_SIZE = 41  # empirically identified default window size
AMINO_ACID = "Amino Acid"
FOLD_INDEX = "FoldIndex"
PAPA = "PAPA"
PRIMA = "PRIMA"
SEQUENCE_ID = "Sequence_ID"
SEQUENCE_POSITION = "Sequence Position"


# Function Declarations
def get_window_bounds(sequence: str, position: int) -> tuple:
    """
    Determines sliding window position and contents.

    Parameters:
        sequence: An amino acid sequence
        position: Starting or left-most index of window along sequence

    Returns: A window (i.e., sequence slice) defined by its starting and ending position in the sequence.
    """
    start = max(position - (WINDOW_SIZE // 2), 0)
    end = min(len(sequence), position + (WINDOW_SIZE // 2) + 1)
    return start, end


def calculate_window_score(
    sequence: str,
    position: int,
    aa_dict: dict,
    ignore_consecutive_prolines: bool = False,
) -> float:
    """
    Applies analysis algorithm (FoldIndex, PAPA or PRIMA) on a sequence slice, calculating a single window score.

    Parameters:
        sequence: An amino acid sequence
        position: ith or leftmost position of sliding window
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines: Flag to instruct sliding window to either include or omit scores for windows containing prolines.

    Returns: calculated window score for a window in a given amino acid sequence.
    """

    start, end = get_window_bounds(sequence, position)
    score = 0.0
    for i in range(start, end):
        if sequence[i] not in AMINO_ACIDS:
            continue
        if not ignore_consecutive_prolines:
            score += aa_dict[sequence[i]]
        else:
            if sequence[i] != "P":
                score += aa_dict[sequence[i]]
            elif (i > 0 and sequence[i - 1] == "P") or (
                i > 1 and sequence[i - 2] == "P"
            ):
                pass
            else:
                score += aa_dict[sequence[i]]
    return score / (end - start)


def calculate_window_scores(
    sequence: str, aa_dict: dict, ignore_consecutive_prolines: bool = False
) -> list:
    """
    Performs 1 of 3 algorithmic analyses (FoldIndex, PAPA, PRIMA) on a sequence and generates list of all window scores.

    Parameters:
        sequence: An amino acid sequence
        aa_dict: The dictionary containing the appropriate log-odds value for the algorithm being performed on sequence.
        ignore_consecutive_prolines: If consecutive residues in a sequence are proline (P), do not apply algorithm.

    Returns: a list containing calculated window scores for all windows in a given amino acid sequence.
    """

    return [
        calculate_window_score(sequence, i, aa_dict, ignore_consecutive_prolines)
        for i in range(len(sequence))
    ]


def calculate_super_window_scores(
    sequence: str, window_scores: list, fold_index_scores: bool = None
) -> list:
    """
    Takes weighted average of window scores for a given sequence.

    Parameters:
        sequence: An amino acid sequence
        window_scores: The dictionary containing the log-odds value for the algorithm being performed on sequence.
        fold_index_scores: Ignore FoldIndex scores when performing analysis.

    Returns: list of weighted window scores of a particular algorithm in the given sequence.
    """
    scores = []
    for i in range(len(sequence) - WINDOW_SIZE):
        if fold_index_scores is not None and fold_index_scores[i] > 0:
            scores.append(None)
        else:
            score = 0.0
            weights = 0.0
            for j in range(i, i + WINDOW_SIZE):
                start, end = get_window_bounds(sequence, j)
                score += (end - start) * window_scores[j]
                weights += end - start
            scores.append(score / weights)
    return scores


def get_fold_index_scores(sequence: str) -> list:
    """
    Calculates FoldIndex scores for all windows in a given sequence.
    Positive(+) FoldIndex scores represent proteins/domains that are likely folded.
    Negative(-) FoldIndex scores represent proteins/domains that are likely intrinsically unfolded.

    Parameters:
        sequence: An amino acid sequence
    Returns: list of FoldIndex scores for all windows in the given sequence.
    """
    charges = calculate_window_scores(sequence, CHARGE)
    hydrophobicities = calculate_window_scores(sequence, HYDROPHOBICITY)
    fold_index_list = [
        2.785 * hydrophobicities[i] - abs(charges[i]) - 1.151
        for i in range(len(sequence))
    ]
    fold_index_list = calculate_super_window_scores(sequence, fold_index_list)
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
    papa_scores = calculate_window_scores(sequence, PROPENSITY, True)
    if ignore_fold_index:
        papa_score_list = calculate_super_window_scores(sequence, papa_scores)
    else:
        papa_score_list = calculate_super_window_scores(
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
    maintenance_scores = calculate_window_scores(sequence, MAINTENANCE)
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
        calculate_super_window_scores(sequence, prima_score_list),
    )

SCORE_METHODS = {
    PAPA: get_papa_scores,
    PRIMA: get_prima_scores,
    FOLD_INDEX: get_fold_index_scores,
}

def analyze_sequence(sequence: str) -> pd.DataFrame:
    """
    Create dataframe with prionic activity scores to inform plots
    """
    elements = {k: v(sequence) for k, v in SCORE_METHODS.items()}
    _, _, papa_scores = elements[PAPA]
    _, _, fold_index_scores = elements[FOLD_INDEX]
    _, _, prima_scores = elements[PRIMA]

    score_array = np.empty(((len(sequence) - WINDOW_SIZE), 4))
    for idx, papa_score in enumerate(papa_scores):
        score_array[idx, :] = [idx + 1, papa_score, prima_scores[idx], fold_index_scores[idx]]
    score_df = pd.DataFrame(
        data=score_array,
        columns=[SEQUENCE_POSITION, *elements.keys()],
    )
    return score_df


def get_df(sequence) -> pd.DataFrame:
    score_df = analyze_sequence(sequence)
    seq_id = SAMPLE_LABELS[SAMPLE_SEQUENCES.index(sequence)]

    # build a new dataframe
    seq_id_col = pd.DataFrame(data=[seq_id] * (len(sequence) - WINDOW_SIZE), columns=[SEQUENCE_ID])
    seq_list = list(sequence)
    aa_df = pd.DataFrame(data=seq_list[: len(sequence) - WINDOW_SIZE], columns=[AMINO_ACID])
    df = pd.concat([seq_id_col, aa_df, score_df], axis=1)

    return df


def get_ranges(df: pd.DataFrame) -> dict:
    """
    Establishes range of value scale from -|max dev| to + |max dev|. Scales graphs correctly.
    """
    methods = list(SCORE_METHODS.keys())
    lower_bounds = df[methods].min().tolist()
    upper_bounds = df[methods].max().tolist()
    ranges = {}
    for method, lower, upper in zip(methods, lower_bounds, upper_bounds):
        largest = max(abs(lower), abs(upper))
        ranges[method] = [-abs(largest) - 0.05, abs(largest) + 0.05]
    return ranges
