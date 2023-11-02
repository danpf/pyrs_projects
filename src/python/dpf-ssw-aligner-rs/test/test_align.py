
from dpf_ssw_aligner_rspy import cmdline_main

def test_align_01():
    cmdline_main(["--matrix-file", "/Users/dan.farrell/git/Complete-Striped-Smith-Waterman-Library/demo/blosum62.txt", "-c", "/Users/dan.farrell/git/pyaln/t.fa", "/Users/dan.farrell/git/pyaln/rq.fa"])
    # cmdline_main(["-c", "t.fa", "q.fa"])

if __name__ == "__main__":
    test_align_01()

