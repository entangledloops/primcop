import collections

import numpy as np
import pandas as pd

from src.characteristics import (
    AMINO_ACIDS,
    PROPENSITY,
    HYDROPHOBICITY,
    CHARGE,
    MAINTENANCE,
)
from src.sample_data import SAMPLE_LABELS, SAMPLE_SEQUENCES


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

    Returns:
        The (start, end) of a window s.t. ``(end-start) <= WINDOW_SIZE``.
    """
    start = max(position - (WINDOW_SIZE // 2), 0)
    end = min(len(sequence), position + (WINDOW_SIZE // 2) + 1)
    return start, end


def calculate_window_scores(
    sequence: str,
    aa_dict: dict,
    ignore_consecutive_prolines: bool = False,
) -> float:
    """
    Applies a sliding window over sequence, computing the average score of the window
    as per the score mapping in ``aa_dict``. Returns a list of scores for all windows.

    Args:
        sequence: An amino acid sequence
        aa_dict: The dictionary containing the appropriate log-odds value for the
            algorithm being performed on sequence.
        ignore_consecutive_prolines: If consecutive residues in a sequence are
            proline (P), do not apply algorithm.

    Returns:
        A list of window scores for all windows in a given amino acid sequence.
    """
    scores = []
    for pos in range(len(sequence)):
        start, end = get_window_bounds(sequence, pos)
        recent_acids = collections.deque(maxlen=3)
        score = 0.0
        for acid in sequence[start:end]:
            recent_acids.append(acid)
            if acid not in AMINO_ACIDS:
                continue
            if (
                ignore_consecutive_prolines
                and "P" == acid
                and recent_acids.count("P") > 1
            ):
                continue
            score += aa_dict[acid]
        scores.append(score / (end - start))
    return scores


def calculate_super_window_scores(
    sequence: str, window_scores: list[float], fold_index_scores: bool = None
) -> list:
    """
    Takes weighted average of window scores for a given sequence.

    Args:
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
    Positive scores represent proteins/domains that are likely folded.
    Negative scores represent proteins/domains that are likely intrinsically unfolded.

    Args:
        sequence: An amino acid sequence

    Returns:
        List of FoldIndex scores for all windows in the given sequence.
    """
    charges = calculate_window_scores(sequence, CHARGE)
    hydrophobicities = calculate_window_scores(sequence, HYDROPHOBICITY)
    fold_index_scores = [
        2.785 * h - abs(c) - 1.151 for h, c, in zip(hydrophobicities, charges)
    ]
    fold_index_scores = calculate_super_window_scores(sequence, fold_index_scores)
    return fold_index_scores


def get_papa_scores(sequence: str, ignore_fold_index: bool = True) -> list:
    """
    Calculates PAPA scores for all windows in a given sequence.

    Args:
        sequence: An amino acid sequence

    Returns:
        List of PAPA scores for all windows in the given sequence.
    """
    papa_scores = calculate_window_scores(sequence, PROPENSITY, True)
    if ignore_fold_index:
        papa_scores = calculate_super_window_scores(sequence, papa_scores)
    else:
        fold_index_scores = get_fold_index_scores(sequence)
        papa_scores = calculate_super_window_scores(
            sequence, papa_scores, fold_index_scores=fold_index_scores
        )
    return papa_scores


def get_prima_scores(sequence: str) -> list:
    """
    Calculates PRIMA scores for all windows in a given sequence.
    Parameters:
        sequence: An amino acid sequence
    Returns: list of PRIMA scores for all windows in the given sequence.
    """
    maintenance_scores = calculate_window_scores(sequence, MAINTENANCE)
    return calculate_super_window_scores(sequence, maintenance_scores)


SCORE_METHODS = {
    PAPA: get_papa_scores,
    PRIMA: get_prima_scores,
    FOLD_INDEX: get_fold_index_scores,
}


def analyze_sequence(sequence: str) -> pd.DataFrame:
    """
    Create dataframe with prionic activity scores to inform plots
    """
    scores = {k: v(sequence) for k, v in SCORE_METHODS.items()}
    papa_scores = scores[PAPA]
    prima_scores = scores[PRIMA]
    fold_index_scores = scores[FOLD_INDEX]

    score_array = np.empty((len(papa_scores), 4))
    for idx, papa_score in enumerate(papa_scores):
        score_array[idx, :] = [
            idx + 1,
            papa_score,
            prima_scores[idx],
            fold_index_scores[idx],
        ]
    score_df = pd.DataFrame(
        data=score_array,
        columns=[SEQUENCE_POSITION, *scores.keys()],
    )
    return score_df


def get_df(sequence) -> pd.DataFrame:
    score_df = analyze_sequence(sequence)
    seq_id = SAMPLE_LABELS[SAMPLE_SEQUENCES.index(sequence)]

    # build a new dataframe
    seq_id_col = pd.DataFrame(
        data=[seq_id] * (len(sequence) - WINDOW_SIZE), columns=[SEQUENCE_ID]
    )
    seq_list = list(sequence)
    aa_df = pd.DataFrame(
        data=seq_list[: len(sequence) - WINDOW_SIZE], columns=[AMINO_ACID]
    )
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
