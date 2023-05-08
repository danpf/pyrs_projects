from pathlib import Path
from dpf_ssw_aligner_rspy.commandline_entrypoints import cmdline_main
from dpf_ssw_aligner_rspy.aligning import Aligner


def test_align_via_cmdline(test_data_dir):
    query_seq_file = Path(test_data_dir) / "r1_query.fq"
    target_seq_file = Path(test_data_dir) / "r1.fa"
    ret = list(cmdline_main(["-t", str(target_seq_file), "-q", str(query_seq_file)]))
    assert len(ret) == 1
    aln = ret[0]
    target_expected = (
        "ggttcgacagcgacgccgcgagtc-cgagagaggagccgcgggcg-ccgtgg---atagagcaggaggggccg-gagtat------tgggaccgg---aacacac"
    )
    visual_expected = (
        "***** ***  *** *** ***** ** * ********* ***** ******   ********** ******* ** ***      **** ****   *******"
    )
    query_expected = (
        "GGATC-CCA--GAC-CCG-GAGACTCG-G-GAGTACCCG-GGGCGTCCGTGGGGCATGGAGGTGG-GGGGTCGTGA-TCTGCGCCCTGGG-CCGGGGTTACTCAC"
    )

    assert aln.query_aln == query_expected
    assert aln.target_aln == target_expected
    assert aln.visual_aln == visual_expected


def test_align_via_api():
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
