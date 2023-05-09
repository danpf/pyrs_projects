import numpy as np

from . import MRC


def match_histogram(base: np.ndarray, to_change: np.ndarray) -> np.ndarray:
    X, Y, Z = np.meshgrid(*map(np.arange, base.shape[::-1]))
    base_table = np.vstack((base.ravel(), X.ravel(), Y.ravel(), Z.ravel())).T
    base_table = base_table[base_table[:, 0].argsort()]

    X, Y, Z = np.meshgrid(*map(np.arange, to_change.shape[::-1]))
    to_change_table = np.vstack((to_change.ravel(), Z.ravel(), Y.ravel(), X.ravel())).T
    to_change_table = to_change_table[to_change_table[:, 0].argsort()]

    assert base_table.shape == to_change_table.shape

    to_change_table[:, 0] = base_table[:, 0]
    tct = to_change_table
    tct = tct[np.lexsort((tct[:, 1], tct[:, 2], tct[:, 3]))]
    tct = tct[:, 0].reshape(base.shape).swapaxes(0, 1)
    return tct


def homogenize_and_combine_volumes(mrcs_in: list[MRC]) -> MRC:
    SMALL_DIFF = 0.0000001

    maps = []
    for mrc in mrcs_in:
        map_ = mrc.volume_data.copy()
        map_[map_ == 0] = SMALL_DIFF
        map_ /= np.max(map_)
        maps.append(map_)

    for i in range(len(maps)):
        if i == 0:
            continue
        maps[i] = match_histogram(maps[0], maps[i])

    sum_map = np.sum(maps, axis=0)

    sum_map[sum_map == 0] = SMALL_DIFF
    sum_map = np.stack([sum_map] * len(maps), axis=3)

    weights_map = np.stack(maps, axis=3)
    weights_map /= sum_map

    avgs = np.average(sum_map, weights=weights_map, axis=3)
    avgs -= SMALL_DIFF

    dens_map = MRC()
    dens_map.volume_data = avgs
    dens_map.header_from_data()
    return dens_map
