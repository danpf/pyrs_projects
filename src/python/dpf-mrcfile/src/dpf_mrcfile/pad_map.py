import numpy as np

from . import MRC


def _check_z_slice(dens_map, z, xlen, ylen):
    for x in range(xlen):
        for y in range(ylen):
            if dens_map.volume_data[x, y, z] > 0:
                return True
    return False


def _check_x_slice(dens_map, x, ylen, zlen):
    for y in range(ylen):
        for z in range(zlen):
            if dens_map.volume_data[x, y, z] > 0:
                return True
    return False


def _check_y_slice(dens_map, y, xlen, zlen):
    for x in range(xlen):
        for z in range(zlen):
            if dens_map.volume_data[x, y, z] > 0:
                return True
    return False


def pad_map(dens_map: MRC, padding: float) -> MRC:
    """Smartly add padding to em maps.

    modifies the input MRC
    """
    xfwdpad = 0
    xbackpad = 0
    yfwdpad = 0
    ybackpad = 0
    zfwdpad = 0
    zbackpad = 0
    for z in range(dens_map.volume_data.shape[2]):
        if _check_z_slice(dens_map, z, dens_map.volume_data.shape[0], dens_map.volume_data.shape[1]):
            zfwdpad = z
            break
    for z in range(dens_map.volume_data.shape[2] - 1, 0, -1):
        if _check_z_slice(dens_map, z, dens_map.volume_data.shape[0], dens_map.volume_data.shape[1]):
            zbackpad = dens_map.volume_data.shape[2] - 1 - z
            break
    for x in range(dens_map.volume_data.shape[0]):
        if _check_x_slice(dens_map, x, dens_map.volume_data.shape[1], dens_map.volume_data.shape[2]):
            xfwdpad = x
            break
    for x in range(dens_map.volume_data.shape[0] - 1, 0 - 1, -1):
        if _check_x_slice(dens_map, x, dens_map.volume_data.shape[1], dens_map.volume_data.shape[2]):
            xbackpad = dens_map.volume_data.shape[0] - 1 - x
            break
    for y in range(dens_map.volume_data.shape[1]):
        if _check_y_slice(dens_map, y, dens_map.volume_data.shape[0], dens_map.volume_data.shape[2]):
            yfwdpad = y
            break
    for y in range(dens_map.volume_data.shape[1] - 1, 0 - 1, -1):
        if _check_y_slice(dens_map, y, dens_map.volume_data.shape[0], dens_map.volume_data.shape[2]):
            ybackpad = dens_map.volume_data.shape[1] - 1 - y
            break
    padding = int(round(1 / dens_map.A_per_pixel * padding))
    xfwdpad = max(0, padding - xfwdpad)
    xbackpad = max(0, padding - xbackpad)
    yfwdpad = max(0, padding - yfwdpad)
    ybackpad = max(0, padding - ybackpad)
    zfwdpad = max(0, padding - zfwdpad)
    zbackpad = max(0, padding - zbackpad)

    xlen = dens_map.xlen / dens_map.nx
    new_o = (
        dens_map.originx - ((xfwdpad + xbackpad) * dens_map.A_per_pixel / 2),
        dens_map.originy - ((yfwdpad + ybackpad) * dens_map.A_per_pixel / 2),
        dens_map.originz - ((zfwdpad + zbackpad) * dens_map.A_per_pixel / 2),
    )

    dens_map.volume_data = np.pad(
        dens_map.volume_data, ((xfwdpad, xbackpad), (yfwdpad, ybackpad), (zfwdpad, zbackpad)), "constant"
    )

    dens_map.header_from_data()
    dens_map.xlen = xlen * dens_map.nx
    dens_map.ylen = xlen * dens_map.ny
    dens_map.zlen = xlen * dens_map.nz
    dens_map.originx = new_o[0]
    dens_map.originy = new_o[1]
    dens_map.originz = new_o[2]

    return dens_map
