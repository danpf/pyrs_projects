from typing import Iterator
from pathlib import Path
import os
import gzip


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
