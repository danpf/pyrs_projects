from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ._rs_bind import (
    PyProfile,
    py_ssw_align,
)
from .builtin_matrices import build_default_matrices
from .file_io import read_fasta_and_fastq_files, read_matrix
from .tracing import TraceResult


@dataclass
class AlignResult:
    score1: int
    score2: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    target_end2: int
    cigar_seq: list[int]
    is_rc: bool


@dataclass
class SSWSeq:
    id_: str
    seq: str
    int_seq: list[int]
    quality: str
    profile: PyProfile
    is_protein: bool
    rc_seq: str
    rc_int_seq: list[int]
    rc_quality: str
    rc_profile: PyProfile
    is_rc: bool

    @classmethod
    def build_seq(
        cls,
        is_protein: bool,
        id_: str,
        seq: str,
        quality: str,
        elements: list[str],
        element_to_int: dict[str, int],
        reverse_complement_map: dict[str, str],
        mat: list[int],
    ) -> SSWSeq:
        int_seq = seq_to_int_representation(seq, elements, element_to_int)
        profile = PyProfile(int_seq, len(int_seq), mat, len(elements), 2)

        if is_protein:
            rc_seq = ""
            rc_quality = ""
            rc_int_seq = []
            rc_profile = None
        else:
            rc_seq = "".join([reverse_complement_map[x] for x in seq[::-1]])
            rc_quality = quality[::-1]
            rc_int_seq = seq_to_int_representation(rc_seq, elements, element_to_int)
            rc_profile = PyProfile(rc_int_seq, len(int_seq), mat, len(elements), 2)

        ret = cls(
            id_=id_,
            seq=seq,
            int_seq=int_seq,
            quality=quality,
            profile=profile,
            is_protein=is_protein,
            rc_seq=rc_seq,
            rc_int_seq=rc_int_seq,
            rc_quality=rc_quality,
            rc_profile=rc_profile,  # type: ignore
            is_rc=False,
        )
        return ret

    def reverse_complement(self) -> None:
        if self.is_protein:
            raise RuntimeError(
                "Reverse complement alignment is not available for protein sequences."
            )
        self.seq, self.rc_seq = self.rc_seq, self.seq
        self.int_seq, self.rc_int_seq = self.rc_int_seq, self.int_seq
        self.quality, self.rc_quality = self.rc_quality, self.quality
        self.is_rc = not self.is_rc


def seq_to_int_representation(
    seq: str, lEle: list[str], dEle2Int: dict[str, int]
) -> list[int]:
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    # TODO does this need to be a c_int8?
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


def align_one(
    query: SSWSeq,
    target: SSWSeq,
    gap_open_penalty: int,
    gap_extension_penalty: int,
    flag: int,
    mask_len: int,
    rc: bool,
) -> AlignResult:
    """
    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)
    """
    if not rc:
        res = py_ssw_align(
            query.profile,
            target.int_seq,
            len(target.int_seq),
            gap_open_penalty,
            gap_extension_penalty,
            flag,
            0,
            0,
            mask_len,
        )
    else:
        res = py_ssw_align(
            query.rc_profile,
            target.rc_int_seq,
            len(target.rc_int_seq),
            gap_open_penalty,
            gap_extension_penalty,
            flag,
            0,
            0,
            mask_len,
        )
    if res is None:
        raise RuntimeError("Problem in running ssw_align - bindings returned None")

    ret = AlignResult(
        score1=res.get_score1(),
        score2=res.get_score2(),
        query_start=res.get_read_begin1(),
        query_end=res.get_read_end1(),
        target_start=res.get_ref_begin1(),
        target_end=res.get_ref_end1(),
        target_end2=res.get_ref_end2(),
        cigar_seq=res.get_cigar().get_seq(),
        is_rc=rc,
    )
    return ret


@dataclass
class Aligner:
    is_protein: bool
    matrix: str
    matrix_file: None | str
    match_score: int
    mismatch_score: int
    gap_open_penalty: int
    gap_extension_penalty: int
    try_rc_and_use_best: bool
    flag: int
    mat: list[int]
    reverse_complement_map: dict[str, str] = field(default_factory=dict)
    elements: list[str] = field(default_factory=list)
    element_to_int: dict[str, int] = field(default_factory=dict)
    int_to_element: dict[int, str] = field(default_factory=dict)

    def _set_dna_params(self):
        self.elements = ["A", "C", "G", "T", "N"]
        self.reverse_complement_map = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "a": "T",
            "c": "G",
            "g": "C",
            "t": "A",
        }
        if not self.matrix_file:
            for i, ele in enumerate(self.elements):
                self.element_to_int[ele] = i
                self.element_to_int[ele.lower()] = i
                self.int_to_element[i] = ele
            self.mat = [0] * (len(self.elements) ** 2)
            for i in range(len(self.elements) - 1):
                for j in range(len(self.elements) - 1):
                    if self.elements[i] == self.elements[j]:
                        self.mat[i * len(self.elements) + j] = self.match_score
                    else:
                        self.mat[i * len(self.elements) + j] -= self.mismatch_score
        else:
            (
                self.elements,
                self.element_to_int,
                self.int_to_element,
                self.mat,
            ) = read_matrix(Path(self.matrix_file))

    def _set_aa_params(self):
        # load AA score matrix
        if not self.matrix_file:
            (
                self.elements,
                self.element_to_int,
                self.int_to_element,
                self.mat,
            ) = build_default_matrices(self.matrix)
        else:
            # assume the format of the input score matrix is the same as that of http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
            (
                self.elements,
                self.element_to_int,
                self.int_to_element,
                self.mat,
            ) = read_matrix(Path(self.matrix_file))

    def _set_params_from_matrices(self):
        self.reverse_complement_map = {}
        if not self.is_protein:
            self._set_dna_params()
        else:
            self._set_aa_params()

    def __post_init__(self):
        self._set_params_from_matrices()
        assert self.match_score >= 0
        assert self.mismatch_score >= 0
        if self.try_rc_and_use_best and self.is_protein:
            raise RuntimeError(
                "Reverse complement alignment is not available for protein sequences."
            )

    def _align_and_build_traceback(
        self, target: SSWSeq, query: SSWSeq, mask_len: int
    ) -> TraceResult:
        res = align_one(
            query,
            target,
            self.gap_open_penalty,
            self.gap_extension_penalty,
            self.flag,
            mask_len,
            False,
        )
        rc_res = None
        if self.try_rc_and_use_best:
            rc_res = align_one(
                query,
                target,
                self.gap_open_penalty,
                self.gap_extension_penalty,
                self.flag,
                mask_len,
                True,
            )

        if rc_res is None or res.score1 > rc_res.score2:
            best_res = res
        else:
            best_res = rc_res

        trace_result = TraceResult.from_align_result(
            query.rc_seq if best_res.is_rc else query.seq,
            target.seq,
            best_res.query_start,
            best_res.target_start,
            best_res.cigar_seq,
        )
        return trace_result

    def run_from_sequences(
        self, query_seqs: Sequence[tuple[str, str]], target_seqs: Sequence[tuple[str, str]]
    ) -> Iterator[TraceResult]:
        for _query_id, _query_seq in query_seqs:
            query_sswseq = SSWSeq.build_seq(
                self.is_protein,
                id_=_query_id,
                seq=_query_seq,
                quality="",
                elements=self.elements,
                element_to_int=self.element_to_int,
                reverse_complement_map=self.reverse_complement_map,
                mat=self.mat,
            )
            mask_len = len(query_sswseq.seq) // 2
            for _target_id, _target_seq in target_seqs:
                target_sswseq = SSWSeq.build_seq(
                    self.is_protein,
                    id_=_target_id,
                    seq=_target_seq,
                    quality="",
                    elements=self.elements,
                    element_to_int=self.element_to_int,
                    reverse_complement_map=self.reverse_complement_map,
                    mat=self.mat,
                )
                yield self._align_and_build_traceback(
                    target_sswseq, query_sswseq, mask_len
                )

    def run_from_files(
        self, query_file: str, target_file: str
    ) -> Iterator[TraceResult]:
        for _query_id, _query_seq, _query_quality in read_fasta_and_fastq_files(
            Path(query_file)
        ):
            query_sswseq = SSWSeq.build_seq(
                self.is_protein,
                id_=_query_id,
                seq=_query_seq,
                quality=_query_quality,
                elements=self.elements,
                element_to_int=self.element_to_int,
                reverse_complement_map=self.reverse_complement_map,
                mat=self.mat,
            )
            mask_len = len(query_sswseq.seq) // 2
            for _target_id, _target_seq, _target_quality in read_fasta_and_fastq_files(
                Path(target_file)
            ):
                target_sswseq = SSWSeq.build_seq(
                    self.is_protein,
                    id_=_target_id,
                    seq=_target_seq,
                    quality=_target_quality,
                    elements=self.elements,
                    element_to_int=self.element_to_int,
                    reverse_complement_map=self.reverse_complement_map,
                    mat=self.mat,
                )
                yield self._align_and_build_traceback(
                    target_sswseq, query_sswseq, mask_len
                )
