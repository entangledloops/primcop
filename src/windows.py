from prions import *

# Default window size
window_size = 41
half_window_size = int(window_size / 2)

aminoacidDict = dict(zip(sample_labels, sample_sequences))

default_seq = sample_sequences[-3]
default_seq_id = sample_labels[-3]

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


def window_score(
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


def window_scores(
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
        window_score(sequence, i, aa_dict, ignore_consecutive_prolines)
        for i in range(len(sequence))
    ]


def super_window_scores(
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
