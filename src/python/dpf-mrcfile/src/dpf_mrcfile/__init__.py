import math
import struct
import time
from pathlib import Path
from typing import Any

import numpy as np


def read_single_entity(type_: str, file_handle, n_to_read: int, endian: str) -> Any:
    return struct.unpack(f"{endian}{type_}", file_handle.read(n_to_read))[0]


class MRC:
    nx: int
    "number of columns in 3d data array"
    ny: int
    "number of rows in 3d data array"
    nz: int
    "number of sections in 3d data array"

    imod: int
    """Mode:
    0 8-bit signed integer (range -128 to 127)
    1 16-bit signed integer
    2 32-bit signed real
    3 transform : complex 16-bit integers
    4 transform : complex 32-bit reals
    6 16-bit unsigned integer
    """

    nxstart: int
    "location of first column in unit cell"
    nystart: int
    "location of first row in unit cell"
    nzstart: int
    "location of first section in unit cell"

    mx: int
    "sampling along X axis of unit cell"
    my: int
    "sampling along Y axis of unit cell"
    mz: int
    "sampling along Z axis of unit cell"

    xlen: float
    "cell x length in angstroms"
    ylen: float
    "cell y length in angstroms"
    zlen: float
    "cell z length in angstroms"

    alpha: float
    "alpha cell angle in degrees"
    beta: float
    "beta cell angle in degrees"
    gamma: float
    "gamma cell angle in degrees"

    mapc: int
    "axis that correspons to cols (1,2,3 for X,Y,Z)"
    mapr: int
    "axis that correspons to rows (1,2,3 for X,Y,Z)"
    mapr: int
    "axis that correspons to sections (1,2,3 for X,Y,Z)"

    amin: float
    "minimum density value"
    amax: float
    "maximum density value"
    amean: float
    "mean density value"

    ispg: int
    """Space Group Number
    Spacegroup 0 implies a 2D image or image stack. For crystallography,
    ISPG represents the actual spacegroup. For single volumes from EM/ET,
    the spacegroup should be 1. For volume stacks, we adopt the convention
    that ISPG is the spacegroup number + 400, which in EM/ET will typically
    be 401.
    """

    NSYMBT: int
    """
    NSYMBT specifies the size of the extended header in bytes, whether it
    contains symmetry records (as in the original format definition) or any
    other kind of additional metadata.
    """

    extra_space: Any
    "extra space that could be used for anything. 0 by default"

    """

    For transforms (Mode 3 or 4), ORIGIN is the phase origin of the
    transformed image in pixels, e.g. as used in helical processing of the
    MRC package. For a transform of a padded image, this value corresponds
    to the pixel position in the padded image of the center of the unpadded
    image.

    For other modes, ORIGIN specifies the real space location of a
    subvolume taken from a larger volume. In the (2-dimensional) example
    shown above, the header of the map containing the subvolume (red
    rectangle) would contain ORIGIN = 100, 120 to specify its position
    with respect to the original volume (assuming the original volume has
    its own ORIGIN set to 0, 0).
    """

    originx: float
    "phase origin (pixels) or subvolume (A) for X"
    originy: float
    "phase origin (pixels) or subvolume (A) for Y"
    originz: float
    "phase origin (pixels) or subvolume (A) for Z"

    mapstring: str
    "character string 'MAP ' to identify file type"

    """
    Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating the
    representation of float, complex, integer and character datatypes.
    Bytes 215 and 216 are unused. The CCP4 library contains a general
    representation of datatypes, but in practice it is safe to use
    0x44 0x44 0x00 0x00 for little endian machines, and
    0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library uses
    this information to automatically byte-swap data if appropriate, when
    tranferring data files between machines.
    """

    machine_stamp: Any
    "machine stamp encoding byte order of data"

    rms_dev: float
    "rms deviation of map from mean density"

    num_labels: int
    "number of labels being used"

    text_header: list[str]
    "10x80 character text labels"

    all_images: np.ndarray
    "Holds all image data in np array"

    volume_data: np.ndarray
    "Holds all volume data in np array"

    # seperate data
    A_per_pixel: float

    def __init__(self):
        """Init all variables to zero.

        mainly used for reference/record keeping
        all information taken from:
        http://www.ccpem.ac.uk/mrc_format/mrc2014.php

        We could have used the CCP-EM library. But this is a much easier to
        understand object, which will be necessary when manipulating these
        files as much as we are going to.
        """

        """
        The data block of an MRC format file holds a 3D array of data (of type
        specified by MODE). NX, NY, NZ specify the dimensions (in grid points)
        of this array. In EM, this will correspond to the dimensions of a
        volume/map, or the combined size of an image/volume stack. In
        crystallography, this will correspond to the dimensions of a map,
        which may cover a crystallographic unit cell or may cover some fraction
        or multiple of a unit cell.
        """
        # number of columns in 3d data array
        self.nx = 0
        # number of rows in 3d data array
        self.ny = 0
        # number of sections in 3d data array
        self.nz = 0

        """
        0 8-bit signed integer (range -128 to 127)
        1 16-bit signed integer
        2 32-bit signed real
        3 transform : complex 16-bit integers
        4 transform : complex 32-bit reals
        6 16-bit unsigned integer
        """
        # Mode,
        self.imod = 0

        # location of first column in unit cell
        self.nxstart = 0
        # location of first row in unit cell
        self.nystart = 0
        # location of first section in unit cell
        self.nzstart = 0

        """
        In crystallographic usage, MZ represents the number of intervals,
        or sampling grid, along Z in a crystallographic unit cell. This need
        not be the same as NZ, if the map doesn't cover exactly a single unit
        cell. For microscopy, where there is no unit cell, MZ represents the
        number of sections in a single volume. For a volume stack, NZ/MZ will
        be the number of volumes in the stack. For images, MZ = 1.
        """
        # sampling along X axis of unit cell
        self.mx = 0
        # sampling along Y axis of unit cell
        self.my = 0
        # sampling along Z axis of unit cell
        self.mz = 0

        """
        cell dimensions in angstroms
        """
        # cell x length in angstroms
        self.xlen = 0
        # cell y length in angstroms
        self.ylen = 0
        # cell z length in angstroms
        self.zlen = 0

        """
        cell angles in degrees
        """
        # alpha angle
        self.alpha = 0
        # beta angle
        self.beta = 0
        # gamma angle
        self.gamma = 0

        """
        In EM MAPC,MAPR,MAPS = 1,2,3 so that sections and images are
        perpendicular to the Z axis. In crystallography, other orderings
        are possible. For example, in some spacegroups it is convenient to
        section along the Y axis (i.e. where this is the polar axis).
        """
        # axis that correspons to cols (1,2,3 for X,Y,Z)
        self.mapc = 0
        # axis that correspons to rows (1,2,3 for X,Y,Z)
        self.mapr = 0
        # axis that correspons to sections (1,2,3 for X,Y,Z)
        self.maps = 0

        """
        Density statistics may not be kept up-to-date for image/volume stacks,
        since it is expensive to recalculate these every time a new
        image/volume is added/deleted. We have proposed the following
        convention: DMAX < DMIN, DMEAN < (smaller of DMIN and DMAX), RMS < 0
        each indicate that the quantity in question is not well determined
        """
        # minimum density value
        self.amin = 0
        # maximum density value
        self.amax = 0
        # mean density value
        self.amean = 0

        """
        Spacegroup 0 implies a 2D image or image stack. For crystallography,
        ISPG represents the actual spacegroup. For single volumes from EM/ET,
        the spacegroup should be 1. For volume stacks, we adopt the convention
        that ISPG is the spacegroup number + 400, which in EM/ET will typically
        be 401.
        """
        # space group number
        self.ispg = 0

        """
        NSYMBT specifies the size of the extended header in bytes, whether it
        contains symmetry records (as in the original format definition) or any
        other kind of additional metadata.
        """
        # size of extended header in bytes
        self.NSYMBT = 0

        # extra space that could be used for anything. 0 by default
        self.extra_space = 0

        """
        For transforms (Mode 3 or 4), ORIGIN is the phase origin of the
        transformed image in pixels, e.g. as used in helical processing of the
        MRC package. For a transform of a padded image, this value corresponds
        to the pixel position in the padded image of the center of the unpadded
        image.

        For other modes, ORIGIN specifies the real space location of a
        subvolume taken from a larger volume. In the (2-dimensional) example
        shown above, the header of the map containing the subvolume (red
        rectangle) would contain ORIGIN = 100, 120 to specify its position
        with respect to the original volume (assuming the original volume has
        its own ORIGIN set to 0, 0).
        """
        # phase origin (pixels) or subvolume (A) for X
        self.originx = 0
        # phase origin (pixels) or subvolume (A) for Y
        self.originy = 0
        # phase origin (pixels) or subvolume (A) for Z
        self.originz = 0

        # character string 'MAP ' to identify file type
        self.mapstring = ""

        """
        Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating the
        representation of float, complex, integer and character datatypes.
        Bytes 215 and 216 are unused. The CCP4 library contains a general
        representation of datatypes, but in practice it is safe to use
        0x44 0x44 0x00 0x00 for little endian machines, and
        0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library uses
        this information to automatically byte-swap data if appropriate, when
        tranferring data files between machines.
        """
        # machine stamp encoding byte order of data
        self.machine_stamp = 0

        # rms deviation of map from mean density
        self.rms_dev = 0

        # number of labels being used
        self.num_labels = 0

        # 10x80 character text labels
        self.text_header = []

        self.all_images = np.array([])
        self.volume_data = np.array([])

        # seperate data
        # This is temporary. I'm not sure how i should set this
        self.A_per_pixel = 0

    def header_from_data(self):
        """Attempt to build a header from given data.

        It is impossible to calculate all of the specifics from just the
        image data, however, we can at least build something that is
        default/will load in chimera with this function.
        """
        assert len(self.volume_data) != 0 or len(self.all_images) != 0
        base_data = None
        if len(self.volume_data):
            base_data = self.volume_data
        else:
            base_data = np.dstack(self.all_images)

        self.nx = base_data.shape[0]
        self.ny = base_data.shape[1]
        self.nz = base_data.shape[2]

        self.imod = 2  # TODO should this be an option?

        self.nxstart = 0
        self.nystart = 0
        self.nzstart = 0

        self.mx = base_data.shape[0]
        self.my = base_data.shape[1]
        self.mz = base_data.shape[2]

        self.xlen = self.A_per_pixel * self.nx
        self.ylen = self.A_per_pixel * self.ny
        self.zlen = self.A_per_pixel * self.nz

        self.alpha = 90
        self.beta = 90
        self.gamma = 90

        self.mapc = 1
        self.mapr = 2
        self.maps = 3

        self.amin = np.amin(base_data)
        self.amax = np.amax(base_data)
        self.amean = np.mean(base_data)

        self.ispg = 1

        self.NSYMBT = 0

        empty_space_list = 25 * [0]
        self.extra_space = struct.pack("25i", *empty_space_list)

        self.originx = 0
        self.originy = 0
        self.originz = 0

        self.mapstring = "MAP "
        self.machine_stamp = (68, 65, 0, 0)

        self.rms_dev = math.sqrt(np.sum(np.square(base_data - self.amean)) / base_data.size)

        self.num_labels = 1

        my_80 = f"MRC file written by dimaio-daimyo (author: Danny Farrell) {time.strftime('%H:%M:%S-%d/%m/%Y')}   ".encode(  # noqa: E501
            "utf-8"
        )
        blank_80 = 80 * " ".encode("utf-8")
        self.text_header = [my_80] + 9 * [blank_80]

    def header_as_string(self):
        s = (
            f"nx,ny,nz: {self.nx, self.ny, self.nz}\n"
            f"imod: {self.imod}\n"
            f"nxstart, nystart, nzstart: {self.nxstart, self.nystart, self.nzstart}\n"
            f"mx, my, mz: {self.mx, self.my, self.mz}\n"
            f"xlen, ylen, zlen: {self.xlen, self.ylen, self.zlen}\n"
            f"alpha, beta, gamma: {self.alpha, self.beta, self.gamma}\n"
            f"mapc, mapr, maps: {self.mapc, self.mapr, self.maps}\n"
            f"amin, amax, amean: {self.amin, self.amax, self.amean}\n"
            f"ispg: {self.ispg}\n"
            f"NSYMBT: {self.NSYMBT}\n"
            f"extra_space: not shown\n"
            f"originx, originy, originz: {self.originx, self.originy, self.originz}\n"
            f"mapstring: >{self.mapstring}<\n"
            f"machine stamp = {self.machine_stamp}\n"
            f"rms dev = {self.rms_dev}\n"
            f"num labels = {self.num_labels}\n"
            f"text header = {self.text_header}\n"
        )
        return s

    def _import_file(self, filename: str | Path) -> None:
        """import a mrc/mrcs file and store the contents in this class.

        The data for a mrc and mrcs file are exactly the same, but the image
        stacks are just each Z-direction. so we import files exactly the same,
        but we store the actual data with a different naming scheme.

        NOTE:
            Since the data is just read in one at a time, other functions will
        bear the responsibility of figuring that shit out.
        """
        endian = "<"
        with open(filename, "rb") as f:
            in_x = f.read(4)
            self.nx = struct.unpack(f"{endian}i", in_x)[0]
            if not (0 < self.nx < 65536):
                endian = ">"
                self.nx = struct.unpack(f"{endian}i", in_x)[0]
            self.ny = read_single_entity("i", f, 4, endian)
            self.nz = read_single_entity("i", f, 4, endian)

            # def read_single_entity(type_: str, file_handle, n_to_read: int, endian: str) -> Any:

            self.imod = read_single_entity("i", f, 4, endian)

            self.nxstart = read_single_entity("i", f, 4, endian)
            self.nystart = read_single_entity("i", f, 4, endian)
            self.nzstart = read_single_entity("i", f, 4, endian)

            self.mx = read_single_entity("i", f, 4, endian)
            self.my = read_single_entity("i", f, 4, endian)
            self.mz = read_single_entity("i", f, 4, endian)

            self.xlen = read_single_entity("f", f, 4, endian)
            self.ylen = read_single_entity("f", f, 4, endian)
            self.zlen = read_single_entity("f", f, 4, endian)

            self.alpha = read_single_entity("f", f, 4, endian)
            self.beta = read_single_entity("f", f, 4, endian)
            self.gamma = read_single_entity("f", f, 4, endian)

            self.mapc = read_single_entity("i", f, 4, endian)
            self.mapr = read_single_entity("i", f, 4, endian)
            self.maps = read_single_entity("i", f, 4, endian)

            self.amin = read_single_entity("f", f, 4, endian)
            self.amax = read_single_entity("f", f, 4, endian)
            self.amean = read_single_entity("f", f, 4, endian)

            self.ispg = read_single_entity("i", f, 4, endian)

            self.NSYMBT = read_single_entity("i", f, 4, endian)

            self.extra_space = f.read(100)

            self.originx = read_single_entity("f", f, 4, endian)
            self.originy = read_single_entity("f", f, 4, endian)
            self.originz = read_single_entity("f", f, 4, endian)

            mapstring = struct.unpack("4c", f.read(4))
            self.mapstring = "".join([x.decode("utf-8") for x in mapstring])

            self.machine_stamp = struct.unpack("4b", f.read(4))

            self.rms_dev = read_single_entity("f", f, 4, endian)

            self.num_labels = read_single_entity("i", f, 4, endian)

            self.text_header = [struct.unpack("80s", f.read(80))[0] for _ in range(10)]

            if math.isclose(self.xlen / self.nx, self.ylen / self.ny, rel_tol=0.01, abs_tol=0) and math.isclose(
                self.ylen / self.ny, self.zlen / self.nz, rel_tol=0.01, abs_tol=0
            ):
                self.A_per_pixel = self.xlen / self.nx

            self.alldata = np.frombuffer(f.read(4 * self.nx * self.ny * self.nz), dtype="<f4")

    def read_mrc(self, filename: str | Path):
        """Read an mrc file and store the contents.

        Inputs:
            filename: The mrc filename to read in.

        NOTE:
            We don't read this in as fortran.. perhaps we should?
        Not entirely sure that would be useful at this point.
        """
        try:
            self._import_file(filename)
        except Exception:
            print(f"failure to read {filename}")
            raise
        self.volume_data = np.reshape(self.alldata, (self.nx, self.ny, self.nz), order="F")
        self.alldata = 0

    def read_mrcs(self, filename: str | Path) -> None:
        """Read an mrcs file and store the contents.

        Inputs:
            filename: The mrcs filename to read in.
        """
        self._import_file(filename)
        self.alldata = np.reshape(self.alldata, (self.nx, self.ny, self.nz), order="F")
        self.all_images = np.dsplit(self.alldata, self.nz)
        self.alldata = 0

    def set_mrcs(self, all_images):
        """Sets the all_images variable and sets header data.

        Inputs:
            List of images in numpy format
        """
        self.all_images = all_images
        self.header_from_data()

    def set_mrc(self, volume_data):
        """Sets the volume_data variable and sets header data.

        Inputs:
            Volume image in numpy format
        """
        self.volume_data = volume_data
        self.header_from_data()

    def export_header(self, f):
        """Export just the header of the mrc file.

        Inputs:
            f: file already opened in binary

        TODO:
            Recalculate things like min, max, and mean on rewrite
        """
        f.write(struct.pack("i", self.nx))
        f.write(struct.pack("i", self.ny))
        f.write(struct.pack("i", self.nz))

        f.write(struct.pack("i", self.imod))

        f.write(struct.pack("i", self.nxstart))
        f.write(struct.pack("i", self.nystart))
        f.write(struct.pack("i", self.nzstart))

        f.write(struct.pack("i", self.mx))
        f.write(struct.pack("i", self.my))
        f.write(struct.pack("i", self.mz))

        f.write(struct.pack("f", self.xlen))
        f.write(struct.pack("f", self.ylen))
        f.write(struct.pack("f", self.zlen))

        f.write(struct.pack("f", self.alpha))
        f.write(struct.pack("f", self.beta))
        f.write(struct.pack("f", self.gamma))

        f.write(struct.pack("i", self.mapc))
        f.write(struct.pack("i", self.mapr))
        f.write(struct.pack("i", self.maps))

        f.write(struct.pack("f", self.amin))
        f.write(struct.pack("f", self.amax))
        f.write(struct.pack("f", self.amean))

        f.write(struct.pack("i", self.ispg))

        f.write(struct.pack("i", self.NSYMBT))

        f.write(self.extra_space)

        f.write(struct.pack("f", self.originx))
        f.write(struct.pack("f", self.originy))
        f.write(struct.pack("f", self.originz))

        for x in self.mapstring:
            f.write(struct.pack("c", x.encode("utf-8")))

        for x in self.machine_stamp:
            f.write(struct.pack("b", x))

        f.write(struct.pack("f", self.rms_dev))

        f.write(struct.pack("i", self.num_labels))

        for x in self.text_header:
            f.write(struct.pack("80s", x))

    def write_mrc_file(self, filename):
        """Write this class to an mrc file.

        Inputs:
            filename: filename to write to
        """
        with open(filename, "wb") as f:
            self.export_header(f)
            volume_data_to_write = self.volume_data.flatten("F")
            f.write(struct.pack(f"={volume_data_to_write.size}f", *volume_data_to_write))

    def write_mrcs_file(self, filename):
        """Write this class to an mrcs file.

        Inputs:
            filename: filename to write to
        """
        with open(filename, "wb") as f:
            self.export_header(f)
            image_data_to_write = np.dstack(self.all_images).flatten("F")
            f.write(struct.pack(f"={image_data_to_write.size}f", *image_data_to_write))

    def get_origin(self):
        return np.array([self.originx, self.originy, self.originz])

    def set_origin(self, new_origin: np.ndarray):
        self.originx, self.originy, self.originz = new_origin


def read_mrc(filename: str) -> MRC:
    mrc = MRC()
    mrc.read_mrc(filename)
    return mrc
