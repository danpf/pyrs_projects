"""standard map conversions."""

import numpy as np

from . import MRC


def mrc_to_situs(my_map: MRC) -> str:
    """convert mrc object to situs string."""
    situs_str = (
        f"{my_map.A_per_pixel:.6f} {my_map.originx:.6f} {my_map.originy:.6f}"
        f" {my_map.originz:.6f} {my_map.nx} {my_map.ny} {my_map.nz}\n\n"
    )
    volume_data_to_write = my_map.volume_data.flatten("F")
    f_string = [f"{val:11.6f}" for val in volume_data_to_write]
    count = 0
    for f in f_string:
        situs_str += f
        count += 1
        if count == 10:
            situs_str += " \n"
            count = 0
        else:
            situs_str += " "
    return situs_str


def situs_to_mrc(my_situs: str) -> MRC:
    data = []
    a_per_pix = None
    originx = None
    originy = None
    originz = None
    nx = None
    ny = None
    nz = None
    for i, line in enumerate(my_situs.split("\n")):
        line = line.strip()
        split_line = line.split()
        if not line:
            continue
        if i == 0:
            a_per_pix = float(split_line[0])
            originx = float(split_line[1])
            originy = float(split_line[2])
            originz = float(split_line[3])
            nx = int(split_line[4])
            ny = int(split_line[5])
            nz = int(split_line[6])
        else:
            data.extend([float(x) for x in split_line])
    error_text = f"Unable to obtain the following fields from situs file: {my_situs} --"
    og_len = len(error_text)
    for x in (a_per_pix, originx, originy, originz, nx, ny, nz):
        if x is None:
            error_text += f" {x=}"
    if len(error_text) != og_len:
        raise RuntimeError(error_text)

    data = np.reshape(np.array(data), (nx, ny, nz), order="F")
    out_mrc = MRC()
    out_mrc.volume_data = data
    out_mrc.header_from_data()
    out_mrc.A_per_pixel = a_per_pix
    out_mrc.originx = originx
    out_mrc.originy = originy
    out_mrc.originz = originz
    return out_mrc


def convert_map_to_mrc(map_in):
    """MRC files use orgin(x,y,z) to set map origin coordinates in space MAP
    files use voxel indexes n(x,y,z)start to set map origin coordinates in
    space."""
    if map_in.nxstart or map_in.nystart or map_in.nzstart:
        map_in.originx = map_in.xlen / map_in.mx * map_in.nxstart
        map_in.originy = map_in.ylen / map_in.my * map_in.nystart
        map_in.originz = map_in.zlen / map_in.mz * map_in.nzstart
        map_in.nxstart = 0
        map_in.nystart = 0
        map_in.nzstart = 0
    return map_in


def convert_mrc_to_map(mrc_in):
    """MRC files use orgin(x,y,z) to set map origin coordinates in space MAP
    files use voxel indexes n(x,y,z)start to set map origin coordinates in
    space."""
    if mrc_in.originx or mrc_in.originy or mrc_in.originz:
        mrc_in.nxstart = round(mrc_in.originx * mrc_in.mx / mrc_in.xlen)
        mrc_in.nystart = round(mrc_in.originy * mrc_in.my / mrc_in.ylen)
        mrc_in.nzstart = round(mrc_in.originz * mrc_in.mz / mrc_in.zlen)
        mrc_in.originx = 0
        mrc_in.originy = 0
        mrc_in.originz = 0
    return mrc_in
