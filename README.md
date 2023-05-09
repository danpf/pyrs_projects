# pyrs_projects
A collection of personal projects

# Striped Smith-Waterman implementation in Rust with python bindings -- dpf-ssw-aligner-rs
This project was my first look into rust, how rust interacts with pyo3, and with SIMD programming.

[Rust Library](src/rust/dpf-ssw-aligner-rs)
[Python Library](src/python/dpf-ssw-aligner-rs)

I was especially impressed with the seamless integration between rust and python via pyo3, especially after many
years of fighting with pybind11.


# MRC/MRCS file handling in a simple and easy package - dpf-mrcfile
MRC is the file format that store the voxel-based data of electron microscopy.
MRCS is the same file format as for MRC files, but instead of holding voxel based data, it holds 2d images
stacked along a matrix's Z-axis.

This project is a very simple library built to make mrc file modifications, creation, and interpretation,
of MRC/MRCS files as easy and simple as possible.

[Python Library](src/python/dpf-mrcfile)
