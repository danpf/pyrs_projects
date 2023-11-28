#!/usr/bin/env python
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from shutil import copytree

from dpf_ssw_aligner_rspy.aligning import Aligner
from dpf_ssw_aligner_rspy.commandline_entrypoints import cmdline_main

SSW_TEST_DIR = Path(__file__).parent


class TestPyPulchra(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_data_dir = Path(self.temp_dir.name) / "test_data"
        copytree(SSW_TEST_DIR / "test_data", self.test_data_dir)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_align_via_cmdline(self):
        query_seq_file = self.test_data_dir / "r1_query.fq"
        target_seq_file = self.test_data_dir / "r1.fa"
        ret = list(
            cmdline_main(["-t", str(target_seq_file), "-q", str(query_seq_file)])
        )
        assert len(ret) == 1
        aln = ret[0]
        target_expected = "ggttcgacagcgacgccgcgagtc-cgagagaggagccgcgggcg-ccgtgg---atagagcaggaggggccg-gagtat------tgggaccgg---aacacac"
        visual_expected = "***** ***  *** *** ***** ** * ********* ***** ******   ********** ******* ** ***      **** ****   *******"
        query_expected = "GGATC-CCA--GAC-CCG-GAGACTCG-G-GAGTACCCG-GGGCGTCCGTGGGGCATGGAGGTGG-GGGGTCGTGA-TCTGCGCCCTGGG-CCGGGGTTACTCAC"

        assert aln.query_aln == query_expected
        assert aln.target_aln == target_expected
        assert aln.visual_aln == visual_expected

    def test_align_via_api(self):
        query = [
            (
                "sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens GN=HBA1 PE=1 SV=2",
                "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
            )
        ]
        ref = [
            (
                ">sp|P01942|HBA_MOUSE Hemoglobin subunit alpha OS=Mus musculus GN=Hba PE=1 SV=2",
                "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR",
            )
        ]
        aligner = Aligner(
            is_protein=True,
            matrix="BLOSUM62",
            matrix_file="",
            match_score=1,
            mismatch_score=1,
            gap_open_penalty=1,
            gap_extension_penalty=1,
            try_rc_and_use_best=False,
            flag=2,
            mat=[],
        )
        ret_l = list(aligner.run_from_sequences(query, ref))
        assert len(ret_l) == 1
        ret = ret_l[0]
        assert (
            ret.target_aln
            == "MVLS-GEDKSNIKAAWGKIGGH-GAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASA-AGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASH-HPADFTPAVHASLDKFLASVSTVLTSKYR"
        )
        assert (
            ret.query_aln
            == "MVLSPA-DKTNVKAAWGKVGAHAG-EYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVA-HVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHL-PAEFTPAVHASLDKFLASVSTVLTSKYR"
        )


if __name__ == "__main__":
    unittest.main()
