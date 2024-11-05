from dataclasses import dataclass, field

@dataclass(frozen=True, kw_only=True)
class Params:
    GC_weight: float = 7.0
    GT_weight: float = 1.25
    AC_weight: float = 1.25
    AT_weight: float = 0.5
    display_sequence_score: int = 0
    mismatch_penalty_starting_value: int = 3
    mismatch_penalty_linear_delta: int = 3
    mismatch_penalty_type: str = "linear"
    mismatch_penalty_choices: tuple[str] = ("linear", "exponential")
    method_choices: tuple[str] = ("transitions", "coverage", "layered")
    consecutive_AT_scoring: tuple[float] = (0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0)
    n_jobs: int = 8
    output_dir: str = "zdna_extractions"
    cadence_reward: float = 1
    method: str = "transitions"
    threshold: int = 50
    headers: list[str] = field(repr=False,
                               default_factory=lambda : ["Chromosome",
                                                         "Start",
                                                         "End",
                                                         "Z-DNA Score",
                                                         "Sequence",
                                                         "totalSequenceScore"]
                               )

    @property
    def __new_dict__(self) -> dict[str, int]:
        representation = self.__dict__.copy()
        representation.pop("headers")
        return representation
