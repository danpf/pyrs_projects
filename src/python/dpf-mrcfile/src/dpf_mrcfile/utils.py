from copy import deepcopy
from typing import Sequence

import numpy as np

from . import MRC


def add_MRCs(maps: Sequence[MRC]) -> MRC:
    if not maps:
        raise RuntimeError("Cannot add_MRCs that don't exist!")
    map_0 = deepcopy(maps[0])
    for map_ in maps[1:]:
        map_0.volume_data += map_.volume_data
    return map_0


def zero_volume_data(map_in: MRC) -> MRC:
    map_in.volume_data = np.zeros(map_in.volume_data.shape)
    return map_in
