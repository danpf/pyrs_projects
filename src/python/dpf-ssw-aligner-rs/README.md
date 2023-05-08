
# dpf_ssw_aligner_rspy

This is python + rust bindings of a Striped Smith Waterman algorithm.

It is a reimplementation of this library: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/ using rust.

It takes about 2x the amount of time as the C implementation takes to run, but it is not built for speed, it is built to mirror
the C implementation as closely as possible.  I'm interested in optimizing the code, but I'm not actively working on that.

## Description
SSW is a fast implementation of the Smith-Waterman algorithm, which uses the Single-Instruction Multiple-Data (SIMD)
instructions to parallelize the algorithm at the instruction level.  It can return the Smith-Waterman score,
alignment location and traceback path (cigar) of the optimal alignment accurately; and return the sub-optimal
alignment score and location heuristically.


## Contributing
Please note that changes to the underlying rust library should be done against `src/rust/dpf-ssw-aligner-rs/src/lib.rs` as 
this library only implements the python wrapping part. I just copy that file into this folder and rename it.

## Usage:
`TODO`


## Thanks to Mengyao Zhao for the C implementation and Yongan Zhao for the python wrapper


## boston college license
License: MIT

Copyright (c) 2012-2015 Boston College

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
