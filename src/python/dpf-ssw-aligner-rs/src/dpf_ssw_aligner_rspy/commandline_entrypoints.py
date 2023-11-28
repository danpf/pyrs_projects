from __future__ import annotations

import argparse
import sys
import time
from collections.abc import Iterator

from .aligning import Aligner
from .tracing import TraceResult


def parse_args(cmdline_args: list[str]):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--match-score",
        type=int,
        default=2,
        help="a positive integer as the score for a match in genome sequence alignment. [default: 2]",
    )
    parser.add_argument(
        "-x",
        "--mismatch-score",
        type=int,
        default=2,
        help="a positive integer as the score for a mismatch in genome sequence alignment. (is multiplied by -1)",
    )
    parser.add_argument(
        "-o",
        "--gap-open-penalty",
        type=int,
        default=3,
        help="a positive integer as the penalty for the gap opening in genome sequence alignment. [default: 3], (is multiplied by -1)",
    )
    parser.add_argument(
        "-e",
        "--gap-extension-penalty",
        type=int,
        default=1,
        help="a positive integer as the penalty for the gap extension in genome sequence alignment. [default: 1], (is multiplied by -1)",
    )
    parser.add_argument(
        "-p",
        "--protein",
        action="store_true",
        help="Do protein sequence alignment. default is genome alignment",
    )
    parser.add_argument(
        "--matrix-file",
        default="",
        help="a file for either Blosum or Pam weight matrix.",
    )
    parser.add_argument(
        "--matrix",
        choices=("BLOSUM62", "BLOSUM50"),
        help="The built in matrix to use - overridden by --matrix-file",
        default="BLOSUM50",
    )
    parser.add_argument(
        "-f",
        "--nThr",
        default=0,
        help="a positive integer. Only output the alignments with the Smith-Waterman score >= N.",
    )
    parser.add_argument(
        "-r",
        "--try-rc-and-use-best",
        action="store_true",
        help="The best alignment will be picked between the original read alignment and the reverse complement read alignment. [default: False]",
    )
    parser.add_argument("-t", "--target", help="target file", required=True)
    parser.add_argument("-q", "--query", help="query file", required=True)
    return parser.parse_args(cmdline_args)


def cmdline_main(cmdline_args: list[str]) -> Iterator[TraceResult]:
    args = parse_args(cmdline_args)
    aligner = Aligner(
        is_protein=args.protein,
        matrix=args.matrix,
        matrix_file=args.matrix_file,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_open_penalty=args.gap_open_penalty,
        gap_extension_penalty=args.gap_extension_penalty,
        try_rc_and_use_best=args.try_rc_and_use_best,
        flag=2,
        mat=[],
    )
    t1 = time.time()
    for r in aligner.run_from_files(args.query, args.target):
        print(r.target_aln)
        print(r.visual_aln)
        print(r.query_aln)
        print(r.cigar_aln)
        yield r
    t2 = time.time()
    print(f"CPU time: {t2-t1} seconds")


def cmdline_wrapper():
    for _ in cmdline_main(sys.argv[1:]):
        pass
