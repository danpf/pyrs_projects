from pathlib import Path
import gzip
from typing import Iterator
from dataclasses import dataclass, field
import os
import sys
import math
import argparse
import time

# import the contents of the Rust library into the Python extension
# optional: include the documentation from the Rust module
from .dpf_ssw_aligner_rspy import PyCigar as PyCigar, py_ssw_align as py_ssw_align, PyAlign as PyAlign, PyProfile as PyProfile


# fmt: off
BLOSUM50 = [
     #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
        5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,    # A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,    # R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,    # N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,    # D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,    # C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,    # Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,    # E
        0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,    # G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,    # H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,    # I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,    # L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,    # K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,    # M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,    # F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,    # P
        1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,    # S
        0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5,    # T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5,    # W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5,    # Y
        0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5,    # V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5,    # B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5,    # Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,    # X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1,     # *
       ]

BLOSUM62 = [
    #    A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
         4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,   # A
        -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,   # R
        -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,   # N
        -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,   # D
         0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,   # C
        -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,   # Q
        -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,   # E
         0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,   # G
        -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,   # H
        -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,   # I
        -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,   # L
        -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,   # K
        -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,   # M
        -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,   # F
        -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,   # P
         1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,   # S
         0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,   # T
        -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,   # W
        -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,   # Y
         0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,   # V
        -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,   # B
        -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,   # Z
         0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,   # X
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1,   # *
]
# fmt: on



def read_matrix(filename: Path) -> tuple[list[str], dict[str, int], dict[int, str], list[int]]:
    """
    read a score matrix for either DNA or protein
    assume the format of the input score matrix is the same as that of http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

    """
    with open(filename, "r") as f:
        lines = f.read().split("\n")
        dEle2Int = {}
        dInt2Ele = {}
        lScore = []
        lEle = []
        for i, x in enumerate(lines):
            if x.startswith("#"):
                continue
            else:
                if not dEle2Int:
                    lEle = x.strip().split()
                    for j, ele in enumerate(lEle):
                        dEle2Int[ele] = j
                        dEle2Int[ele.lower()] = j
                        dInt2Ele[i] = ele
                else:
                    lScore.extend([int(y) for y in x.strip().split()[1:]])

        assert lEle
        return lEle, dEle2Int, dInt2Ele, lScore


def read_fasta_and_fastq_files(filename: Path) -> Iterator[tuple[str, str, str]]:
    """
    read a sequence file
    @param  sFile   sequence file
    """

    def read_one_fasta_from_lines(lines: list[str]) -> Iterator[tuple[str, str, str]]:
        """
        read a fasta file
        """
        sId = ""
        sSeq = ""
        for l in lines:
            if l.startswith(">"):
                if sSeq:
                    yield sId, sSeq, ""
                sId = l.strip()[1:].split()[0]
                sSeq = ""
            else:
                sSeq += l.strip()

        yield sId, sSeq, ""

    def read_one_fastaq_from_lines(lines: list[str]) -> Iterator[tuple[str, str, str]]:
        """
        read a fastq file
        @param  f   file handler
        """
        sId = ""
        sSeq = ""
        # s3 = ""
        sQual = ""
        for i, l in enumerate(lines):
            if l.startswith("@"):
                sId = l.strip()[1:].split()[0]
                sSeq = lines[i + 1].strip()
                # s3 = lines[i + 2].strip()
                sQual = lines[i + 3].strip()

            yield sId, sSeq, sQual

    ext = os.path.splitext(filename)[1][1:].strip().lower()
    if ext == "gz" or ext == "gzip":
        with gzip.open(filename, "r") as f:
            # l = f.readline().decode()
            text = f.read().decode()
    else:
        with open(filename, "r") as f:
            text = f.read()

    # read
    if text.startswith(">"):
        for sId, sSeq, sQual in read_one_fasta_from_lines(text.split("\n")):
            yield sId, sSeq, sQual
    elif text.startswith("@"):
        for sId, sSeq, sQual in read_one_fastaq_from_lines(text.split("\n")):
            yield sId, sSeq, sQual
    else:
        raise RuntimeError(f"File format cannot be recognized {filename=} -- {ext=} -- {text[:20]=}")


def seq_to_int_representation(seq: str, lEle: list[str], dEle2Int: dict[str, int]) -> list[int]:
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    print("TODO should convert to c_int8")
    # TODO
    # num_decl = len(seq) * ct.c_int8
    # num = num_decl()
    num = []
    for ele in seq:
        n = dEle2Int[lEle[-1]]
        try:
            n = dEle2Int[ele]
        except KeyError:
            pass
        finally:
            num.append(n)
    return num


def align_one(qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = py_ssw_align(qProfile, rNum, nRLen, nOpen, nExt, nFlag, 0, 0, nMaskLen)
    if res is None:
        raise RuntimeError("Problem in running ssw_align - bindings returned None")

    nScore = res.get_score1()
    nScore2 = res.get_score2()
    nRefBeg = res.get_ref_begin1()
    nRefEnd = res.get_ref_end1()
    nQryBeg = res.get_read_begin1()
    nQryEnd = res.get_read_end1()
    nRefEnd2 = res.get_ref_end2()
    cigar = res.get_cigar()
    lCigar = cigar.get_seq()
    nCigarLen = cigar.get_length()

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)


def buildPath(q, r, nQryBeg, nRefBeg, lCigar):
    """
    build cigar string and align path based on cigar array returned by ssw_align
    @param  q   query sequence
    @param  r   reference sequence
    @param  nQryBeg   begin position of query sequence
    @param  nRefBeg   begin position of reference sequence
    @param  lCigar   cigar array
    """
    sCigarInfo = "MIDNSHP=X"
    sCigar = ""
    sQ = ""
    sA = ""
    sR = ""
    nQOff = nQryBeg
    nROff = nRefBeg
    for _, x in enumerate(lCigar):
        n = x >> 4
        m = x & 15
        if m > 8:
            c = "M"
        else:
            c = sCigarInfo[m]
        sCigar += str(n) + c

        if c == "M":
            sQ += q[nQOff : nQOff + n]
            sA += "".join(["|" if q[nQOff + j] == r[nROff + j] else "*" for j in range(n)])
            sR += r[nROff : nROff + n]
            nQOff += n
            nROff += n
        elif c == "I":
            sQ += q[nQOff : nQOff + n]
            sA += " " * n
            sR += "-" * n
            nQOff += n
        elif c == "D":
            sQ += "-" * n
            sA += " " * n
            sR += r[nROff : nROff + n]
            nROff += n
    return sCigar, sQ, sA, sR


def build_default_matrices(matrix_name: str):
    dEle2Int = {}
    dInt2Ele = {}
    if matrix_name == "BLOSUM50":
        matrix = BLOSUM50
    elif matrix_name == "BLOSUM62":
        matrix = BLOSUM62
    else:
        raise RuntimeError(f"Unrecognized built-in matrix name {matrix_name}")
    lEle = (
        "A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *".split()
    )
    for i, ele in enumerate(lEle):
        dEle2Int[ele] = i
        dEle2Int[ele.lower()] = i
        dInt2Ele[i] = ele
    # nEleNum = len(lEle)
    return lEle, dEle2Int, dInt2Ele, matrix

@dataclass
class Aligner:
    protein: bool
    matrix: str
    matrix_file: None | str
    match_score: int
    mismatch_score: int
    try_rc_and_use_best: bool
    mat: list[int]
    dRc: dict[str, str] = field(default_factory=dict)


    def set_params_from_matrices(self):
        self.lEle = []
        self.dRc = {}
        self.dEle2Int = {}
        self.dInt2Ele = {}
        if not self.protein:
            # init DNA score matrix
            if not self.matrix_file:
                self.lEle = ["A", "C", "G", "T", "N"]
                dRc = {"A": "T", "C": "G", "G": "C", "T": "A", "a": "T", "c": "G", "g": "C", "t": "A"}
                for i, ele in enumerate(self.lEle):
                    self.dEle2Int[ele] = i
                    self.dEle2Int[ele.lower()] = i
                    self.dInt2Ele[i] = ele
                self.mat = [0] * (len(self.lEle)**2)
                for i in range(len(self.lEle) - 1):
                    for j in range(len(self.lEle) - 1):
                        if self.lEle[i] == self.lEle[j]:
                            self.mat[i * len(self.lEle) + j] = self.match_score
                        else:
                            self.mat[i * len(self.lEle) + j] -= self.mismatch_score
            else:
                self.lEle, self.dEle2Int, self.dInt2Ele, self.mat = read_matrix(Path(self.matrix_file))
        else:
            # load AA score matrix
            if not self.matrix_file:
                self.lEle, self.dEle2Int, self.dInt2Ele, self.mat = build_default_matrices(self.matrix)
            else:
                # assume the format of the input score matrix is the same as that of http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
                self.lEle, self.dEle2Int, self.dInt2Ele, self.mat = read_matrix(Path(self.matrix_file))

    def __init__(self):
        self.set_params_from_matrices()
        assert self.match_score >= 0
        assert self.mismatch_score >= 0
        if self.try_rc_and_use_best and self.protein:
            raise RuntimeError("Reverse complement alignment is not available for protein sequences.")



def main(args):
    lEle = []
    dRc = {}
    dEle2Int = {}
    dInt2Ele = {}
    if not args.bProtein:
        # init DNA score matrix
        if not args.matrix_file:
            lEle = ["A", "C", "G", "T", "N"]
            dRc = {"A": "T", "C": "G", "G": "C", "T": "A", "a": "T", "c": "G", "g": "C", "t": "A"}
            for i, ele in enumerate(lEle):
                dEle2Int[ele] = i
                dEle2Int[ele.lower()] = i
                dInt2Ele[i] = ele
            nEleNum = len(lEle)
            lScore = [0] * (nEleNum**2)
            for i in range(nEleNum - 1):
                for j in range(nEleNum - 1):
                    if lEle[i] == lEle[j]:
                        lScore[i * nEleNum + j] = args.match_score.
                    else:
                        lScore[i * nEleNum + j] = -args.mismatch_score.
        else:
            lEle, dEle2Int, dInt2Ele, lScore = read_matrix(args.matrix_file)
    else:
        # load AA score matrix
        if not args.matrix_file:
            lEle, dEle2Int, dInt2Ele, lScore = read_matrix(args.matrix)
        else:
            # assume the format of the input score matrix is the same as that of http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
            lEle, dEle2Int, dInt2Ele, lScore = read_matrix(args.matrix_file)

    if args.bBest and args.bProtien:
        sys.stderr.write("Reverse complement alignment is not available for protein sequences.\n")

    mat = list(lScore)
    # set flag
    nFlag = 0
    if args.bPath:
        nFlag = 2
    # print sam head
    if args.bSam and args.bHeader and args.bPath:
        print("@HD\tVN:1.4\tSO:queryname")
        for sRId, sRSeq, _ in read_fasta_and_fastq_files(args.target):
            print("@SQ\tSN:{}\tLN:{}".format(sRId, len(sRSeq)))
    elif args.bSam and not args.bPath:
        sys.stderr.write("SAM format output is only available together with option -c.\n")
        args.bSam = False

    # ssw = ssw_lib.CSsw(args.sLibPath)
    # iterate query sequence
    for sQId, sQSeq, sQQual in read_fasta_and_fastq_files(args.query):
        # build query profile
        qNum = seq_to_int_representation(sQSeq, lEle, dEle2Int)
        qProfile = PyProfile(qNum, len(sQSeq), mat, len(lEle), 2)
        # build rc query profile
        if args.bBest and not args.bProtein:
            sQRcSeq = "".join([dRc[x] for x in sQSeq[::-1]])
            qRcNum = seq_to_int_representation(sQRcSeq, lEle, dEle2Int)
            qRcProfile = PyProfile(qRcNum, len(sQSeq), mat, len(lEle), 2)
        # set mask len
        nMaskLen = len(sQSeq) // 2

        # iter target sequence
        for sRId, sRSeq, _ in read_fasta_and_fastq_files(args.target):
            rNum = seq_to_int_representation(sRSeq, lEle, dEle2Int)
            # format of res: (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
            res = align_one(qProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)
            # align rc query
            resRc = None
            if args.bBest and not args.bProtein:
                resRc = align_one(qRcProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)

            # build cigar and trace back path
            strand = 0
            if resRc == None or res[0] > resRc[0]:
                resPrint = res
                strand = 0
                sCigar, sQ, sA, sR = buildPath(sQSeq, sRSeq, res[4], res[2], res[8])
            else:
                resPrint = resRc
                strand = 1
                sCigar, sQ, sA, sR = buildPath(sQRcSeq, sRSeq, resRc[4], resRc[2], resRc[8])
            print("query/target")
            print(sQ)
            print(sR)

            # print results
            if not args.bSam:
                print("target_name: {}\nquery_name: {}\noptimal_alignment_score: {}\t".format(sRId, sQId, resPrint[0])),
                if resPrint[1] > 0:
                    print("suboptimal_alignment_score: {}\t".format(resPrint[1])),
                if strand == 0:
                    print("strand: +\t"),
                else:
                    print("strand: -\t"),
                if resPrint[2] + 1:
                    print("target_begin: {}\t".format(resPrint[2] + 1)),
                print("target_end: {}\t".format(resPrint[3] + 1)),
                if resPrint[4] + 1:
                    print("query_begin: {}\t".format(resPrint[4] + 1)),
                print("query_end: {}\n".format(resPrint[5] + 1))
                if resPrint[-2] > 0:
                    n1 = 1 + resPrint[2]
                    n2 = min(60, len(sR)) + resPrint[2] - sR.count("-", 0, 60)
                    n3 = 1 + resPrint[4]
                    n4 = min(60, len(sQ)) + resPrint[4] - sQ.count("-", 0, 60)
                    for i in range(0, len(sQ), 60):
                        print("Target:{:>8}\t{}\t{}".format(n1, sR[i : i + 60], n2))
                        n1 = n2 + 1
                        n2 = n2 + min(60, len(sR) - i - 60) - sR.count("-", i + 60, i + 120)

                        print("{: ^15}\t{}".format("", sA[i : i + 60]))

                        print("Query:{:>9}\t{}\t{}\n".format(n3, sQ[i : i + 60], n4))
                        n3 = n4 + 1
                        n4 = n4 + min(60, len(sQ) - i - 60) - sQ.count("-", i + 60, i + 120)
            else:
                print("{}\t".format(sQId)),
                if resPrint[0] == 0:
                    print("4\t*\t0\t255\t*\t*\t0\t0\t*\t*"),
                else:
                    mapq = int(-4.343 * math.log(1 - abs(resPrint[0] - resPrint[1]) / float(resPrint[0])))
                    mapq = int(mapq + 4.99)
                    if mapq >= 254:
                        mapq = 254
                    if strand == 1:
                        print("16\t"),
                    else:
                        print("0\t"),
                    print("{}\t{}\t{}\t".format(sRId, resPrint[2] + 1, mapq)),
                    print(sCigar),
                    print("\t*\t0\t0\t"),
                    print(sQSeq[resPrint[4] : resPrint[5] + 1]) if strand == 0 else print(
                        sQRcSeq[resPrint[4] : resPrint[5] + 1]
                    ),
                    print("\t"),
                    if sQQual:
                        if strand == 0:
                            print(sQQual[resPrint[4] : resPrint[5] + 1]),
                        else:
                            print(sQQual[-resPrint[4] - 1 : -resPrint[5] - 1 : -1])
                    else:
                        print("*"),

                    print("\tAS:i:{}".format(resPrint[0])),
                    print("\tNM:i:{}\t".format(len(sA) - sA.count("|"))),
                    if resPrint[1] > 0:
                        print("ZS:i:{}".format(resPrint[1]))

        # ssw.init_destroy(qProfile)
        # if args.bBest and not args.bProtein:
        #     ssw.init_destroy(qRcProfile)


def parse_args(cmdline_args: list[str]):
    parser = argparse.ArgumentParser()
    # parser.add_argument('-l', '--sLibPath', default='./', help='path of libssw.so')
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
        "--nOpen",
        type=int,
        default=3,
        help="a positive integer as the penalty for the gap opening in genome sequence alignment. [default: 3]",
    )
    parser.add_argument(
        "-e",
        "--nExt",
        type=int,
        default=1,
        help="a positive integer as the penalty for the gap extension in genome sequence alignment. [default: 1]",
    )
    parser.add_argument(
        "-p",
        "--bProtein",
        action="store_true",
        help="Do protein sequence alignment. Without this option, pyssw will do genome sequence alignment. [default: False]",
    )
    parser.add_argument(
        "--matrix-file", default="", help="a file for either Blosum or Pam weight matrix."
    )
    parser.add_argument(
        "--matrix", choices=("BLOSUM62", "BLOSUM50"), help="The built in matrix to use - overridden by --matrix-file", default="BLOSUM50"
    )
    parser.add_argument("-c", "--bPath", action="store_true", help="Return the alignment path. [default: False]")
    parser.add_argument(
        "-f",
        "--nThr",
        default=0,
        help="a positive integer. Only output the alignments with the Smith-Waterman score >= N.",
    )
    parser.add_argument(
        "-r",
        "--bBest",
        action="store_true",
        help="The best alignment will be picked between the original read alignment and the reverse complement read alignment. [default: False]",
    )
    parser.add_argument("-s", "--bSam", action="store_true", help="Output in SAM format. [default: no header]")
    parser.add_argument(
        "-header", "--bHeader", action="store_true", help="If -s is used, include header in SAM output."
    )
    parser.add_argument("target", help="target file")
    parser.add_argument("query", help="query file")
    args = parser.parse_args(cmdline_args)
    return args


def cmdline_main(cmdline_args: list[str]):
    t1 = time.time()
    main(parse_args(cmdline_args))
    t2 = time.time()
    print(f"CPU time: {t2-t1} seconds")


if __name__ == "__main__":
    cmdline_main(sys.argv[1:])
