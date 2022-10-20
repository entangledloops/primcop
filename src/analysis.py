from src.characteristics import amino_acids, propensity, hydrophobicity, charge, maintenance
from src.sample_data import samples_dict
import numpy as np
import pandas as pd

# Variables
window_size = 41  # empirically identified default window size
half_window_size = int(window_size / 2)


# Function Declarations
def get_window(sequence: str, position: int) -> tuple:
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

    start, end = get_window(sequence, position)
    score = 0.0
    for i in range(start, end):
        if sequence[i] not in amino_acids:
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
    for i in range(len(sequence) - window_size):
        if fold_index_scores is not None and fold_index_scores[i] > 0:
            scores.append(None)
        else:
            score = 0.0
            weights = 0.0
            for j in range(i, i + window_size):
                start, end = get_window(sequence, j)
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
    charges = calculate_window_scores(sequence, charge)
    hydrophobicities = calculate_window_scores(sequence, hydrophobicity)
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
    papa_scores = calculate_window_scores(sequence, propensity, True)
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
    maintenance_scores = calculate_window_scores(sequence, maintenance)
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


def get_sequence_scores(sequence: str) -> tuple:
    """
    Scan sequence for potential prion activity, assigning window- and sequence-level scores 
    """
    return (
        get_fold_index_scores(sequence),
        get_papa_scores(sequence),
        get_prima_scores(sequence),
    )


def analyze_sequence(sequence: str, sequence_id: str) -> pd.DataFrame:
    """
    Create dataframe with prionic activity scores to inform plots
    """
    elements = get_sequence_scores(sequence)
    seq_papa_score, papa_pos, papa_scores = elements[1]
    seq_fold_index_score, fold_index_pos, fold_index_scores = elements[0]
    seq_prima_score, prima_pos, prima_scores = elements[2]

    score_array = np.empty(((len(sequence) - window_size), 4))
    for idx, papa_score in enumerate(papa_scores):
        score_array[idx, :] = [idx + 1, papa_score, prima_scores[idx], fold_index_scores[idx]]
    score_df = pd.DataFrame(
        data=score_array,
        columns=["Sequence Position", "PAPA score", "PRIMA score", "FoldIndex score"],
    )

    seq_id_col = pd.DataFrame(data=[sequence_id] * (len(sequence) - window_size), columns=["Sequence_ID"])
    seq_list = list(sequence)
    aa_df = pd.DataFrame(data=seq_list[: len(sequence) - window_size], columns=["Amino Acid"])
    df = pd.concat([seq_id_col, aa_df, score_df], axis=1)
    return df

def update_df(value) -> pd.DataFrame:
        sequence = value
        sequence_id = [seq_id for seq_id, seq_val in samples_dict.items() if seq_val == value]
        df = analyze_sequence(sequence, sequence_id)
        return df


def set_range(lows: list, highs: list) -> tuple:
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


def get_ranges(df: pd.DataFrame) -> tuple:
    lower_bounds = df[["PAPA score", "PRIMA score", "FoldIndex score"]].min().tolist()
    upper_bounds = df[["PAPA score", "PRIMA score", "FoldIndex score"]].max().tolist()
    return set_range(lower_bounds, upper_bounds)
