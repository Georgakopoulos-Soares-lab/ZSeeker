from Bio.Seq import Seq
import numpy as np

from typing import Optional
from .constants import Params

class ZDNACalculatorSeq(Seq):

    def __init__(self, data,  params: Optional[Params] = None) -> None:
        super().__init__(data.upper())
        self.scoring_array = []  # Calculate scoring array in subarrays method do do only when necessary
        if params is None:
            self.params = Params()
        else:
            self.params = params


    def zdna_calculator_layered(self):
        """
        Calculates the scoring array for a given sequence, using seperate scorring for transitions, cadence and TA penalty
        :return: np.ndarray: Scoring array
        """
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0  # k is the number of mismatches
        cadence = 1  # indicating if the alternative purine-pyrimidine pattern is unbroken
        consecutive_AT_counter = 0  # Number of consecutive TAs or ATs
        for i in range(len(self) - 1):
            match self[i] + self[i + 1]:
                case "GC" | "CG":
                    scoring_array[i] = self.params.GC_weight
                    mismatches_counter = 0  # Reset the number of mismatches
                case "GT" | "TG":
                    scoring_array[i] = self.params.GT_weight
                    mismatches_counter = 0
                case "AC" | "CA":
                    scoring_array[i] = self.params.AC_weight
                    mismatches_counter = 0
                case "AT" | "TA":
                    adjusted_weight = self.params.AT_weight
                    # If we have less consecutive ATs than the provided AT scoring list, use the corresponding element
                    # Else use the last element of the list as penalty
                    if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                        adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                    else:
                        adjusted_weight += self.params.consecutive_AT_scoring[-1]
                    scoring_array[i] = adjusted_weight
                    consecutive_AT_counter += 1
                    mismatches_counter = 0
                case _:  # If there was a mismatch
                    mismatches_counter += 1
                    if self.params.mismatch_penalty_type == "exponential":
                        # Limit the maximum value of the mismatch penalty to avoid exploding values
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value ** mismatches_counter if mismatches_counter < 15 else -32000
                    elif self.params.mismatch_penalty_type == "linear":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value - self.params.mismatch_penalty_linear_delta * (mismatches_counter-1)
                    else:
                        raise ValueError(f"Mismatch penalty type not recognized. Valid options are {Params.mismatch_penalty_choices}")
            if self[i] + self[i + 1] in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
                scoring_array[i] += self.params.cadence_reward
        return scoring_array

    # TODO: Fix - This implementation results in the scoring array being generated every time we splice the sequence, which is not optimal.
    def zdna_calculator_transitions(self) -> np.ndarray:
        """
        Calculates the scoring array for a given sequence, treating each transition individually
        :return: np.ndarray: Scoring array
        """
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0  # mismatches_counter is the number of mismatches
        consecutive_AT_counter = 0
        for i in range(len(self) - 1):
            match self[i] + self[i + 1]:
                case "GC" | "CG":
                    scoring_array[i] = self.params.GC_weight
                    mismatches_counter = 0  # Reset the number of mismatches
                case "GT" | "TG":
                    scoring_array[i] = self.params.GT_weight
                    mismatches_counter = 0
                case "AC" | "CA":
                    scoring_array[i] = self.params.AC_weight
                    mismatches_counter = 0
                case "AT" | "TA":
                    adjusted_weight = self.params.AT_weight
                    # If we have less consecutive ATs than the provided AT scoring list, use the corresponding element
                    # Else use the last element of the list as penalty
                    if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                        adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                    else:
                        adjusted_weight += self.params.consecutive_AT_scoring[-1]
                    scoring_array[i] = adjusted_weight
                    consecutive_AT_counter += 1
                    mismatches_counter = 0
                case _:  # If there was a mismatch
                    mismatches_counter += 1
                    if self.params.mismatch_penalty_type == "exponential":
                        # Limit the maximum value of the mismatch penalty to avoid exploding values
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value ** mismatches_counter if mismatches_counter < 15 else -32000
                    elif self.params.mismatch_penalty_type == "linear":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value - self.params.mismatch_penalty_linear_delta * (mismatches_counter-1)
                    else:
                        raise ValueError(f"Mismatch penalty type not recognized. Valid options are {Params.mismatch_penalty_choices}")
        return scoring_array

    def zdna_calculator_coverage(self):
        """
        Adithya's proposal
        :return:
        """
        s = str(self)
        weights = {"gc": self.params.GC_weight, "gt": self.params.GT_weight, "ac": self.params.AC_weight, "at": self.params.AT_weight}
        scores = np.full(len(self), fill_value=-3, dtype=float)

        state = 0

        for i in range(len(s)):
            if i == 0:
                if s[i] in ['G', 'A']:
                    state = -1
                elif s[i] in ['C', 'T']:
                    state = 1

            elif state == -1:
                if s[i] in ['C', 'T']:
                    state = 2
                    sub_str = s[i - 1:i + 1]
                    match sub_str:
                        case "GC":
                            scores[i - 1:i + 1] = [weights["gc"]] * 2
                            last_weight = weights["gc"]
                        case "GT":
                            scores[i - 1:i + 1] = [weights["gt"]] * 2
                            last_weight = weights["gt"]
                        case "AC":
                            scores[i - 1:i + 1] = [weights["ac"]] * 2
                            last_weight = weights["ac"]
                        case "AT":
                            scores[i - 1:i + 1] = [weights["at"]] * 2
                            last_weight = weights["at"]

            elif state == 1:
                if s[i] in ['G', 'A']:
                    state = -2
                    sub_str = s[i - 1:i + 1]
                    match sub_str:
                        case "CG":
                            scores[i - 1:i + 1] = [weights["gc"]] * 2
                            last_weight = weights["gc"]
                        case "TG":
                            scores[i - 1:i + 1] = [weights["gt"]] * 2
                            last_weight = weights["gt"]
                        case "CA":
                            scores[i - 1:i + 1] = [weights["ac"]] * 2
                            last_weight = weights["ac"]
                        case "TA":
                            scores[i - 1:i + 1] = [weights["at"]] * 2
                            last_weight = weights["at"]

            elif state == -2:
                if s[i] in ['C', 'T']:
                    state = 2
                    sub_str = s[i - 1:i + 1]
                    match sub_str:
                        case "GC":
                            scores[i - 1] = max(scores[i - 1], weights["gc"])
                            scores[i] = weights["gc"]
                            last_weight = weights["gc"]
                        case "GT":
                            scores[i - 1] = max(scores[i - 1], weights["gt"])
                            scores[i] = weights["gt"]
                            last_weight = weights["gt"]
                        case "AC":
                            scores[i - 1] = max(scores[i - 1], weights["ac"])
                            scores[i] = weights["ac"]
                            last_weight = weights["ac"]
                        case "AT":
                            scores[i - 1] = max(scores[i - 1], weights["at"])
                            scores[i] = weights["at"]
                            last_weight = weights["at"]
                elif s[i] in ['G', 'A']:
                    state = -3

            elif state == 2:
                if s[i] in ['G', 'A']:
                    state = -2
                    sub_str = s[i - 1:i + 1]
                    match sub_str:
                        case "CG":
                            scores[i - 1] = max(scores[i - 1], weights["gc"])
                            scores[i] = weights["gc"]
                            last_weight = weights["gc"]
                        case "TG":
                            scores[i - 1] = max(scores[i - 1], weights["gt"])
                            scores[i] = weights["gt"]
                            last_weight = weights["gt"]
                        case "CA":
                            scores[i - 1] = max(scores[i - 1], weights["ac"])
                            scores[i] = weights["ac"]
                            last_weight = weights["ac"]
                        case "TA":
                            scores[i - 1] = max(scores[i - 1], weights["at"])
                            scores[i] = weights["at"]
                            last_weight = weights["at"]
                elif s[i] in ['C', 'T']:
                    state = 3

            elif state == -3:
                if s[i] in ['G', 'A']:
                    state = -1
                elif s[i] in ['C', 'T']:
                    sub_str = s[i - 1:i + 1]
                    state = 2
                    match sub_str:
                        case "GC":
                            if weights["gc"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["gc"]] * 2
                            else:
                                scores[i] = weights["gc"]
                            last_weight = weights["gc"]
                        case "GT":
                            if weights["gt"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["gt"]] * 2
                            else:
                                scores[i] = weights["gt"]
                            last_weight = weights["gt"]
                        case "AC":
                            if weights["ac"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["ac"]] * 2
                            else:
                                scores[i] = weights["ac"]
                            last_weight = weights["ac"]
                        case "AT":
                            if weights["at"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["at"]] * 2
                            else:
                                scores[i] = weights['at']
                            last_weight = weights["at"]

            elif state == 3:
                if s[i] in ['C', 'T']:
                    state = 1
                elif s[i] in ['G', 'A']:
                    sub_str = s[i - 1:i + 1]
                    state = -2
                    match sub_str:
                        case "CG":
                            if weights["gc"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["gc"]] * 2
                            else:
                                scores[i] = weights["gc"]
                            last_weight = weights["gc"]
                        case "TG":
                            if weights["gt"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["gt"]] * 2
                            else:
                                scores[i] = weights["gt"]
                            last_weight = weights["gt"]
                        case "CA":
                            if weights["ac"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["ac"]] * 2
                            else:
                                scores[i] = weights["ac"]
                            last_weight = weights["ac"]
                        case "TA":
                            if weights["at"] > last_weight:
                                scores[i - 2] = -3
                                scores[i - 1:i + 1] = [weights["at"]] * 2
                            else:
                                scores[i] = weights["at"]
                            last_weight = weights["at"]

        for i in range(1, len(scores)):
            if scores[i] == -3 and scores[i - 1] < 0:
                scores[i] = scores[i - 1] - 2

        return scores

    # TODO: Convert to generator
    def subarrays_above_threshold(self) -> list[tuple[int, int, int, str]]:
        # Calculate scoring array here instead of at initialization to avoid calculating it unecessarily
        match self.params.method:
            case "transitions":
                self.scoring_array = self.zdna_calculator_transitions()
            case "coverage":
                self.scoring_array = self.zdna_calculator_coverage()
            case "layered":
                self.scoring_array = self.zdna_calculator_layered()
            case _:  # If the method is not recognized
                raise ValueError(f"Method {self.params.method} not recognized. Valid options are: {Params.method_choices}")


        subarrays_above_threshold = []
        max_ending_here = self.scoring_array[0]
        start_idx = end_idx = 0
        current_max = 0  # Stores the maximum value of the currently considered candidate array
        candidate_array = None
        for i in range(1, len(self.scoring_array)):
            num = self.scoring_array[i]
            if num >= max_ending_here + num:
                start_idx = i
                end_idx = i + 1  # +1 because we are looking at the next element of the sequence to determine the weight
                max_ending_here = num
            else:
                max_ending_here += num
                end_idx = i + 1  # +1 because we are looking at the next element of the sequence to determine the weight
            if max_ending_here >= self.params.threshold and current_max < max_ending_here:
                candidate_array = (start_idx, end_idx, max_ending_here, str(self[start_idx: end_idx + 1]))
                current_max = max_ending_here

            # If we have a candidate array and the max_ending_here value dropped below 0, there is no way this element
            # will be party of a subarray above the threshold, so we add it to the list and reset the candidate array
            # Additionally, if the drop between the highest seen value and the current value is above the threshold
            # we also want to stop
            if candidate_array and (max_ending_here < 0 or current_max - max_ending_here >= self.params.threshold):
                subarrays_above_threshold.append(candidate_array)
                candidate_array = None
                max_ending_here = current_max = 0

        # If we have finished iterating without finding a negative value towards the end, we add the last candidate array
        if candidate_array:
            subarrays_above_threshold.append(candidate_array)

        return subarrays_above_threshold


# if __name__ == '__main__':
#     #test_string = "TTCGGCGAACTTCGGCGAGC"
#     #test_string = "TCGAGCGCCAGCGCCAAGCGCAGAGCGCAGGAGCGCAGCGCAGAGCGCCAG"
#     #test_string = "TCGCGCGATCGCGCGCGCGCGCGCCGCGATATATCG"
    
#     test_string = "AAAAACTAGCGGTCCG"

#     # test_string = "CGCGCGTATATATAGCCCCCCCCCCCCAGAGAGGGGAAGAGAGGGCCCCCAGAGTATATAGCGGAGAGC"
#     params = Params(method="layered", mismatch_penalty_type="linear", mismatch_penalty_starting_value=3, threshold=5)
#     seq = ZDNACalculatorSeq(test_string, params=params)

#     print([x for x in seq])
#     print(seq.scoring_array)
#     print(sum(seq.scoring_array))
#     print(seq.subarrays_above_threshold())

#     tests = False
#     if tests:
#         print("-------------------------")
#         params = Params(mismatch_penalty_type="linear", mismatch_penalty_starting_value=3, threshold=0)
#         seq = ZDNACalculatorSeq(test_string, params=params)
#         seq.scoring_array = [2, 10, 2, 2, -3, -2, 2, 2, 3]
#         for elem in seq.subarrays_above_threshold():
#             print(elem[0], elem[1], elem[2], seq.scoring_array[elem[0]: elem[1]])


