from dataclasses import dataclass, field


@dataclass(frozen=True, kw_only=True)
class Params:

    GC_weight: float = 7.0
    GT_weight: float = 1.0
    AC_weight: float = 1.0
    AT_weight: float = 0.5
    display_sequence_score: int = 0
    mismatch_penalty_starting_value: int = 1.0 
    mismatch_penalty_linear_delta: int = 2.0
    mismatch_penalty_type: str = "linear"
    mismatch_penalty_choices: tuple[str] = ("linear", "exponential")
    method_choices: tuple[str] = ("transitions", "coverage", "layered")
    consecutive_AT_scoring: tuple[int] = (0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0)
    n_jobs: int = 8
    cadence_reward: float = 1
    method: str = "transitions"
    threshold: int = 50
    headers: list[str] = field(repr=False,
                               default_factory=lambda : ["Chromosome",
                                                         "Start",
                                                         "End",
                                                         "Z-DNA Score",
                                                         "Sequence",
                                                         "Total Sequence Score"])

    @property
    def __new_dict__(self) -> dict[str, int]:
        representation = self.__dict__.copy()
        representation.pop("headers")
        return representation
