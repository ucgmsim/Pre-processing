import dataclasses
import numpy as np


@dataclasses.dataclass
class Asperity:
    """
    Represents an asperity within a fault.
    x0, x1 should be in the range [-L/2, L/2] when L is the total fault length
    y0, y1 should be in the range [0, W] where W is the fault width
    v is the base slip value to be applied in the asperity region. This will be scaled relative to the background value.
    """

    x0: int
    x1: int
    y0: int
    y1: int
    v: float

    def __post_init__(self):
        if self.x0 > self.x1:
            self.x0, self.x1 = self.x1, self.x0
        if self.y0 > self.y1:
            self.y0, self.y1 = self.y1, self.y0
        assert not np.isnan(self.v)

    def to_asperity_file_format(self):
        return f"{self.v} {self.x0} {self.y0} {self.x1} {self.y1}"

    def __str__(self):
        return f"An asperity covering the region from {self.x0} to {self.x1} along strike, and from {self.y0} to {self.y1} down dip. The asperity has a value of {self.v}."
