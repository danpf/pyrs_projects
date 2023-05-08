from dataclasses import dataclass


@dataclass
class TraceResult:
    target_aln: str
    query_aln: str
    cigar_aln: str
    visual_aln: str

    @classmethod
    def from_align_result(
        cls, query_seq: str, target_seq: str, query_start: int, target_start: int, cigar_seq: list[int]
    ) -> "TraceResult":
        sCigarInfo = "MIDNSHP=X"
        cigar_aln = ""
        query_aln = ""
        visual_aln = ""
        target_aln = ""
        current_query_offset = query_start
        current_target_offset = target_start
        for _, x in enumerate(cigar_seq):
            n = x >> 4
            m = x & 15
            if m > 8:
                c = "M"
            else:
                c = sCigarInfo[m]
            cigar_aln += str(n) + c

            if c == "M":
                query_aln += query_seq[current_query_offset : current_query_offset + n]
                visual_aln += "".join(
                    [
                        "|" if query_seq[current_query_offset + j] == target_seq[current_target_offset + j] else "*"
                        for j in range(n)
                    ]
                )
                target_aln += target_seq[current_target_offset : current_target_offset + n]
                current_query_offset += n
                current_target_offset += n
            elif c == "I":
                query_aln += query_seq[current_query_offset : current_query_offset + n]
                visual_aln += " " * n
                target_aln += "-" * n
                current_query_offset += n
            elif c == "D":
                query_aln += "-" * n
                visual_aln += " " * n
                target_aln += target_seq[current_target_offset : current_target_offset + n]
                current_target_offset += n
        return cls(
            target_aln=target_aln,
            query_aln=query_aln,
            cigar_aln=cigar_aln,
            visual_aln=visual_aln,
        )
