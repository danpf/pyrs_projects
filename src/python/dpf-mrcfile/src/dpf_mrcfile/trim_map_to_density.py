import math

import numpy as np

from . import MRC


def new_max_and_min(mi, ma, cubesize):
    dist = abs(ma - mi)
    remaineder = cubesize - dist
    new_min = mi - math.ceil(remaineder / 2)
    new_max = ma + math.floor(remaineder / 2)
    return new_min, new_max


def trim_map_to_density(map_in: MRC, lower_limit: float, padding: int, force_cube: bool = False) -> MRC:
    xs, ys, zs = np.where(map_in.volume_data >= lower_limit)

    x_min = min(xs) - padding
    x_max = max(xs) + padding

    y_min = min(ys) - padding
    y_max = max(ys) + padding

    z_min = min(zs) - padding
    z_max = max(zs) + padding

    if force_cube:
        x_dist = abs(x_max - x_min)
        y_dist = abs(y_max - y_min)
        z_dist = abs(z_max - z_min)
        box_size = max(x_dist, y_dist, z_dist)
        new_map = np.zeros((box_size, box_size, box_size))
        xyz_shape = map_in.volume_data[min(xs) : max(xs), min(ys) : max(ys), min(zs) : max(zs)].shape
        xslice, yslice, zslice = [slice(int((box_size - x) / 2), int((box_size - x) / 2) + x) for x in xyz_shape]
        new_map[xslice, yslice, zslice] = map_in.volume_data[min(xs) : max(xs), min(ys) : max(ys), min(zs) : max(zs)]
        originx = map_in.originx - (xslice.start - min(xs)) * map_in.A_per_pixel
        originy = map_in.originy - (yslice.start - min(ys)) * map_in.A_per_pixel
        originz = map_in.originz - (zslice.start - min(zs)) * map_in.A_per_pixel
    else:
        new_map = np.zeros((x_max - x_min, y_max - y_min, z_max - z_min))
        print("left sides", x_min, y_min, z_min, new_map.shape)
        new_map = map_in.volume_data[
            x_min - padding : x_max + padding, y_min - padding : y_max + padding, z_min - padding : z_max + padding
        ]
        originx = map_in.originx - padding * map_in.A_per_pixel
        originy = map_in.originy - padding * map_in.A_per_pixel
        originz = map_in.originz - padding * map_in.A_per_pixel

    map_in.volume_data = new_map

    a_per_pix = map_in.A_per_pixel

    map_in.header_from_data()
    if not force_cube:
        map_in.originx = a_per_pix * x_min + originx
        map_in.originy = a_per_pix * y_min + originy
        map_in.originz = a_per_pix * z_min + originz
    else:
        map_in.originx = originx
        map_in.originy = originy
        map_in.originz = originz

    map_in.A_per_pixel = a_per_pix
    return map_in
