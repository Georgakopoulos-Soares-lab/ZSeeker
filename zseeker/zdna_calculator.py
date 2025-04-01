from Bio.Seq import Seq
import numpy as np
from typing import Optional
from .constants import Params

class ZDNACalculatorSeq(Seq):

    def __init__(self, data, params: Optional[Params] = None) -> None:
        super().__init__(data.upper())
        self.scoring_array = []
        if params is None:
            self.params = Params()
        else:
            self.params = params

    def zdna_calculator_layered(self):
        """
        Calculates the scoring array for a given sequence, using separate scoring 
        for transitions, cadence, and TA penalty.
        """
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0
        consecutive_AT_counter = 0  # Counter for consecutive AT/TA transitions

        for i in range(len(self) - 1):
            transition = self[i] + self[i + 1]
            match transition:
                case "GC" | "CG":
                    scoring_array[i] = self.params.GC_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0  # Reset because different transition
                case "GT" | "TG":
                    scoring_array[i] = self.params.GT_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0
                case "AC" | "CA":
                    scoring_array[i] = self.params.AC_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0
                case "AT" | "TA":
                    adjusted_weight = self.params.AT_weight
                    if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                        adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                    else:
                        adjusted_weight += self.params.consecutive_AT_scoring[-1]
                    scoring_array[i] = adjusted_weight
                    consecutive_AT_counter += 1
                    mismatches_counter = 0
                case _:
                    mismatches_counter += 1
                    consecutive_AT_counter = 0  # Reset because transition is not AT/TA
                    if self.params.mismatch_penalty_type == "exponential":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value ** mismatches_counter \
                                           if mismatches_counter < 15 else -32000
                    elif self.params.mismatch_penalty_type == "linear":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value \
                                           - self.params.mismatch_penalty_linear_delta * (mismatches_counter - 1)
                    else:
                        raise ValueError(f"Mismatch penalty type not recognized. "
                                         f"Valid options are {Params.mismatch_penalty_choices}")

            # cadence reward
            if transition in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
                scoring_array[i] += self.params.cadence_reward

        return scoring_array

    def zdna_calculator_transitions(self) -> np.ndarray:
        """
        Calculates the scoring array for a given sequence, treating each transition individually.
        """
        scoring_array = np.empty(len(self) - 1, dtype=float)
        mismatches_counter = 0
        consecutive_AT_counter = 0  # Counter for consecutive AT/TA transitions

        for i in range(len(self) - 1):
            transition = self[i] + self[i + 1]
            match transition:
                case "GC" | "CG":
                    scoring_array[i] = self.params.GC_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0  # Reset counter
                case "GT" | "TG":
                    scoring_array[i] = self.params.GT_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0
                case "AC" | "CA":
                    scoring_array[i] = self.params.AC_weight
                    mismatches_counter = 0
                    consecutive_AT_counter = 0
                case "AT" | "TA":
                    adjusted_weight = self.params.AT_weight
                    if consecutive_AT_counter < len(self.params.consecutive_AT_scoring):
                        adjusted_weight += self.params.consecutive_AT_scoring[consecutive_AT_counter]
                    else:
                        adjusted_weight += self.params.consecutive_AT_scoring[-1]
                    scoring_array[i] = adjusted_weight
                    consecutive_AT_counter += 1
                    mismatches_counter = 0
                case _:
                    mismatches_counter += 1
                    consecutive_AT_counter = 0  # Reset because not an AT/TA transition
                    if self.params.mismatch_penalty_type == "exponential":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value ** mismatches_counter \
                                           if mismatches_counter < 15 else -32000
                    elif self.params.mismatch_penalty_type == "linear":
                        scoring_array[i] = -self.params.mismatch_penalty_starting_value \
                                           - self.params.mismatch_penalty_linear_delta * (mismatches_counter - 1)
                    else:
                        raise ValueError(f"Mismatch penalty type not recognized. "
                                         f"Valid options are {Params.mismatch_penalty_choices}")

        return scoring_array

    def zdna_calculator_coverage(self):
        """
        Adithya's coverage-based approach (placeholder).
        """
        s = str(self)
        weights = {
            "gc": self.params.GC_weight,
            "gt": self.params.GT_weight,
            "ac": self.params.AC_weight,
            "at": self.params.AT_weight
        }
        scores = np.full(len(self), fill_value=-3, dtype=float)
        state = 0
        last_weight = 0

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

        # A quick fix for leftover -3 values
        for i in range(1, len(scores)):
            if scores[i] == -3 and scores[i - 1] < 0:
                scores[i] = scores[i - 1] - 2

        return scores

    def subarrays_above_threshold(self) -> list[tuple[int, int, int, str]]:
        """
        Finds subarrays above the threshold using a variation of Kadane's approach,
        but stopping or "dropping" the subarray if the drop from its max is >= drop_threshold.
        """
        match self.params.method:
            case "transitions":
                self.scoring_array = self.zdna_calculator_transitions()
            case "coverage":
                self.scoring_array = self.zdna_calculator_coverage()
            case "layered":
                self.scoring_array = self.zdna_calculator_layered()
            case _:
                raise ValueError(f"Method {self.params.method} not recognized. "
                                 f"Valid options are: {Params.method_choices}")

        subarrays_above_threshold = []
        max_ending_here = self.scoring_array[0]
        start_idx = end_idx = 0
        current_max = 0
        candidate_array = None

        for i in range(1, len(self.scoring_array)):
            num = self.scoring_array[i]
            if num >= max_ending_here + num:
                start_idx = i
                end_idx = i + 1
                max_ending_here = num
            else:
                max_ending_here += num
                end_idx = i + 1

            if max_ending_here >= self.params.threshold and current_max < max_ending_here:
                candidate_array = (start_idx, end_idx, max_ending_here, str(self[start_idx: end_idx + 1]))
                current_max = max_ending_here

            # Stop or drop the subarray if the drop from its max exceeds drop_threshold.
            if candidate_array and (max_ending_here < 0 or
                                    current_max - max_ending_here >= self.params.drop_threshold):
                subarrays_above_threshold.append(candidate_array)
                candidate_array = None
                max_ending_here = current_max = 0

        if candidate_array:
            subarrays_above_threshold.append(candidate_array)

        return subarrays_above_threshold
