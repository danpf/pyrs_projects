#![feature(portable_simd)]

use std::cmp::{max, min};
use std::mem;

use std::simd::Which::{First, Second};
use std::simd::{
    i16x8, simd_swizzle, u8x16, SimdInt, SimdOrd, SimdPartialEq, SimdPartialOrd, SimdUint, Which,
};

pub const MAPSTR: &str = "MIDNSHP=X";

#[derive(Clone)]
pub struct Profile {
    profile_byte: Option<Vec<u8x16>>,
    profile_word: Option<Vec<i16x8>>,
    read: Vec<i8>,
    mat: Vec<i8>,
    read_len: i32,
    n: i32,
    bias: u8,
}

#[derive(Clone, Debug)]
pub struct Align {
    pub score1: u16,
    pub score2: u16,
    pub ref_begin1: i32,
    pub ref_end1: i32,
    pub read_begin1: i32,
    pub read_end1: i32,
    pub ref_end2: i32,
    pub cigar: Cigar,
    pub flag: u16,
}

#[derive(Clone, Debug)]
pub struct AlignmentEnd {
    score: u16,
    ref_position: i32,  // 0-based position
    read_position: i32, // alignment ending position on read, 0-based
}

#[derive(Clone, Debug)]
pub struct Cigar {
    pub seq: Vec<u32>,
    pub length: usize,
}

impl Profile {
    pub fn ssw_init(read: Vec<i8>, read_len: i32, mat: Vec<i8>, n: i32, score_size: i8) -> Self {
        let mut p = Profile {
            profile_byte: None,
            profile_word: None,
            bias: 0,
            read: read.clone(),
            mat: mat.clone(),
            read_len,
            n,
        };

        if score_size == 0 || score_size == 2 {
            let mut bias = 0;
            for i in 0..n * n {
                if mat[i as usize] < bias {
                    bias = mat[i as usize];
                }
            }

            p.bias = bias.abs() as u8;
            p.profile_byte = Some(query_profile_byte(&read, &mat, read_len, n, p.bias));
        }
        if score_size == 1 || score_size == 2 {
            p.profile_word = Some(query_profile_word(
                &read,
                &mat,
                read_len as usize,
                n as usize,
            ));
        }

        p
    }
}

pub fn ssw_align(
    prof: Profile,
    ref_seq: Vec<i8>,
    ref_len: i32,
    weight_gap_o: u8,
    weight_gap_e: u8,
    flag: u8,
    filters: u16,
    filterd: i32,
    mask_len: i32,
) -> Option<Align> {
    let mut bests: [AlignmentEnd; 2] = [
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
    ];
    let mut bests_reverse: [AlignmentEnd; 2] = [
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
    ];
    let mut word = 0;
    let read_len = prof.read_len;
    let mut read_reverse: Vec<i8> = Vec::new();
    let mut path: Option<Cigar> = None;

    let mut r = Align {
        score1: 0,
        score2: 0,
        ref_begin1: -1,
        ref_end1: -1,
        read_begin1: -1,
        read_end1: -1,
        ref_end2: -1,
        cigar: Cigar {
            seq: Vec::new(),
            length: 0,
        },
        flag: 0,
    };

    if mask_len < 15 {
        eprintln!(
            "When maskLen < 15, the function ssw_align doesn't return 2nd best alignment information."
        );
    }

    // Find the alignment scores and ending positions
    if let Some(profile_byte) = &prof.profile_byte {
        bests = sw_sse2_byte(
            &ref_seq,
            0,
            ref_len,
            read_len,
            weight_gap_o,
            weight_gap_e,
            &profile_byte,
            u8::MAX,
            prof.bias,
            mask_len,
        );
        if let Some(profile_word) = &prof.profile_word {
            if bests[0].score == 255 {
                bests = sw_sse2_word(
                    &ref_seq,
                    0,
                    ref_len,
                    read_len,
                    weight_gap_o,
                    weight_gap_e,
                    &profile_word,
                    u16::MAX,
                    mask_len,
                );
                word = 1;
            } else if bests[0].score == 255 {
                eprintln!("Please set 2 to the score_size parameter of the function ssw_init, otherwise the alignment results will be incorrect.");
                // TODO maybe return None?
                return None
            }
        }
    } else if let Some(profile_word) = &prof.profile_word {
        bests = sw_sse2_word(
            &ref_seq,
            0,
            ref_len,
            read_len,
            weight_gap_o,
            weight_gap_e,
            &profile_word,
            u16::MAX,
            mask_len,
        );
        word = 1;
    } else {
        eprintln!("Please call the function ssw_init before ssw_align.");
        return None
    }

    r.score1 = bests[0].score;
    r.ref_end1 = bests[0].ref_position;
    r.read_end1 = bests[0].read_position;

    if mask_len >= 15 {
        r.score2 = bests[1].score;
        r.ref_end2 = bests[1].ref_position;
    } else {
        r.score2 = 0;
        r.ref_end2 = -1;
    }

    if flag == 0 || (flag == 2 && r.score1 < filters as u16) {
        eprintln!("notsure 9n.");
        return Some(r)
    }

    // Find the beginning position of the best alignment.
    read_reverse = seq_reverse(&prof.read, r.read_end1);
    if word == 0 {
        if prof.profile_byte.is_some() {
            let v_p8 =
                query_profile_byte(&read_reverse, &prof.mat, r.read_end1 + 1, prof.n, prof.bias);
            bests_reverse = sw_sse2_byte(
                &ref_seq,
                1,
                r.ref_end1 + 1,
                r.read_end1 + 1,
                weight_gap_o,
                weight_gap_e,
                &v_p8,
                r.score1 as u8,
                prof.bias,
                mask_len,
            );
        } else {
            eprintln!("Error: ssw_align given prof.profile_byte that was None.");
            return None
        }
    } else {
        if prof.profile_word.is_some() {
            let v_p16 = query_profile_word(
                &read_reverse,
                &prof.mat,
                (r.read_end1 + 1) as usize,
                prof.n as usize,
            );
            bests_reverse = sw_sse2_word(
                &ref_seq,
                1,
                r.ref_end1 + 1,
                r.read_end1 + 1,
                weight_gap_o,
                weight_gap_e,
                &v_p16,
                r.score1,
                mask_len,
            );
        } else {
            eprintln!("Error: ssw_align given prof.profile_word that was None.");
            return None
        }
    }
    r.ref_begin1 = bests_reverse[0].ref_position;
    r.read_begin1 = r.read_end1 - bests_reverse[0].read_position;

    if r.score1 > bests_reverse[0].score {
        eprintln!("Warning: The alignment path of one pair of sequences may miss a small part. ssw_align");
        r.flag = 2;
    }

    if (flag & 7) == 0
        || ((flag & 2) != 0 && r.score1 < filters as u16)
        || ((flag & 4) != 0
            && (r.ref_end1 - r.ref_begin1 > filterd || r.read_end1 - r.read_begin1 > filterd))
    {
        eprintln!("notsure 7n.");
        return Some(r);
    }

    // Generate cigar.
    let ref_len = r.ref_end1 - r.ref_begin1 + 1;
    let read_len = r.read_end1 - r.read_begin1 + 1;
    let band_width = (ref_len - read_len).abs() + 1;
    path = banded_sw(
        &ref_seq[(r.ref_begin1 as usize)..=(r.ref_end1 as usize)],
        &prof.read[(r.read_begin1 as usize)..=(r.read_end1 as usize)],
        ref_len,
        read_len,
        r.score1 as i32,
        weight_gap_o as u32,
        weight_gap_e as u32,
        band_width,
        &prof.mat,
        prof.n,
    );

    if let Some(p) = path {
        r.cigar = p;
        // r.cigar_len = p.length;
    } else {
        r.flag = 1; // banded_sw failed
    }
    Some(r)
}

impl Align {
    pub fn mark_mismatch(
        &mut self,
        ref_begin1: i32,
        read_begin1: i32,
        read_end1: i32,
        ref_seq: &[i8],
        read: &[i8],
        read_len: i32,
        cigar: &mut Cigar,
    ) -> u32 {
        let mut mismatch_length: u32 = 0;
        let mut p: i32 = 0;
        let mut length: u32;
        let mut s: i32 = (cigar.length as i32) + 2;
        let mut new_cigar = Cigar {
            seq: vec![0; s as usize],
            length: 0,
        };
        let mut length_m = 0;
        let mut length_x = 0;
        let mut op;

        let mut ref_seq = &ref_seq[ref_begin1 as usize..];
        let mut read = &read[read_begin1 as usize..];
        if read_begin1 > 0 {
            new_cigar.seq[p as usize] = to_cigar_int(read_begin1 as u32, 'S');
            p += 1;
        }
        for i in 0..cigar.length {
            op = cigar_int_to_op(cigar.seq[i as usize]);
            length = cigar_int_to_len(cigar.seq[i as usize]);
            if op == 'M' {
                for j in 0..length {
                    if ref_seq[0] != read[0] {
                        mismatch_length += 1;
                        new_cigar.seq = store_previous_m(
                            2,
                            &mut length_m,
                            &mut length_x,
                            &mut p,
                            &mut s,
                            new_cigar.seq,
                        );
                        length_x += 1;
                    } else {
                        new_cigar.seq = store_previous_m(
                            1,
                            &mut length_m,
                            &mut length_x,
                            &mut p,
                            &mut s,
                            new_cigar.seq,
                        );
                        length_m += 1;
                    }
                    ref_seq = &ref_seq[1..];
                    read = &read[1..];
                }
            } else if op == 'I' {
                read = &read[length as usize..];
                mismatch_length += length;
                new_cigar.seq = store_previous_m(
                    0,
                    &mut length_m,
                    &mut length_x,
                    &mut p,
                    &mut s,
                    new_cigar.seq,
                );
                new_cigar.seq = add_cigar(new_cigar.seq, &mut p, &mut s, length as u32, 'I');
            } else if op == 'D' {
                ref_seq = &ref_seq[length as usize..];
                mismatch_length += length;
                new_cigar.seq = store_previous_m(
                    0,
                    &mut length_m,
                    &mut length_x,
                    &mut p,
                    &mut s,
                    new_cigar.seq,
                );
                new_cigar.seq = add_cigar(new_cigar.seq, &mut p, &mut s, length as u32, 'D');
            }
        }
        new_cigar.seq = store_previous_m(
            0,
            &mut length_m,
            &mut length_x,
            &mut p,
            &mut s,
            new_cigar.seq,
        );

        length = (read_len - read_end1 - 1) as u32;
        if length > 0 {
            new_cigar.seq = add_cigar(new_cigar.seq, &mut p, &mut s, length as u32, 'S');
        }

        new_cigar.length = p as usize;
        *cigar = new_cigar;
        mismatch_length
    }
}

pub fn to_cigar_int(length: u32, op_letter: char) -> u32 {
    (length << 4) | ENCODED_OPS[op_letter as usize] as u32
}

pub fn cigar_int_to_op(cigar_int: u32) -> char {
    let index = (cigar_int & 0xf) as usize;
    if index > 8 {
        'M'
    } else {
        MAPSTR.chars().nth(index).unwrap()
    }
}

pub fn cigar_int_to_len(cigar_int: u32) -> u32 {
    cigar_int >> 4
}

pub fn seq_reverse(seq: &[i8], end: i32) -> Vec<i8> {
    let mut reverse = vec![0; (end + 1) as usize];
    let mut start = 0;
    let mut end = end;

    while start <= end {
        reverse[start as usize] = seq[end as usize];
        reverse[end as usize] = seq[start as usize];
        start += 1;
        end -= 1;
    }

    reverse
}
//
pub fn add_cigar(
    mut new_cigar: Vec<u32>,
    p: &mut i32,
    s: &mut i32,
    length: u32,
    op: char,
) -> Vec<u32> {
    if p >= s {
        *s += 1;
        new_cigar.resize(*s as usize, 0);
    }
    new_cigar[*p as usize] = to_cigar_int(length, op);
    *p += 1;
    new_cigar
}
//
pub fn store_previous_m(
    choice: i8,
    length_m: &mut u32,
    length_x: &mut u32,
    p: &mut i32,
    s: &mut i32,
    mut new_cigar: Vec<u32>,
) -> Vec<u32> {
    if *length_m > 0 && (choice == 2 || choice == 0) {
        new_cigar = add_cigar(new_cigar, p, s, *length_m, '=');
        *length_m = 0;
    } else if *length_x > 0 && (choice == 1 || choice == 0) {
        new_cigar = add_cigar(new_cigar, p, s, *length_x, 'X');
        *length_x = 0;
    }
    new_cigar
}

// fn kroundup32(x: i32) -> i32 {
pub fn kroundup32<
    T: std::ops::Shr<i32, Output = T>
        + std::ops::SubAssign
        + std::ops::AddAssign
        + std::ops::BitOrAssign
        + Copy
        + PartialEq
        + From<u8>,
>(
    mut x: T,
) -> T {
    // let mut x = x - 1;
    x -= T::from(1);
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x += T::from(1);
    x
}
// Convert the coordinate in the scoring matrix into the coordinate in one line of the band.
pub fn set_u(w: i32, i: i32, j: i32) -> i32 {
    let x = i - w;
    let x = if x > 0 { x } else { 0 };
    j - x + 1
}

// Convert the coordinate in the direction matrix into the coordinate in one line of the band.
pub fn set_d(w: i32, i: i32, j: i32, p: i32) -> i32 {
    let x = i - w;
    let x = if x > 0 { x } else { 0 };
    let x = j - x;
    x * 3 + p
}

/* array index is an ASCII character value from a CIGAR,
element value is the corresponding integer opcode between 0 and 8 */
const ENCODED_OPS: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, /*   */
    0, /* ! */
    0, /* " */
    0, /* # */
    0, /* $ */
    0, /* % */
    0, /* & */
    0, /* ' */
    0, /* ( */
    0, /* ) */
    0, /* * */
    0, /* + */
    0, /* , */
    0, /* - */
    0, /* . */
    0, /* / */
    0, /* 0 */
    0, /* 1 */
    0, /* 2 */
    0, /* 3 */
    0, /* 4 */
    0, /* 5 */
    0, /* 6 */
    0, /* 7 */
    0, /* 8 */
    0, /* 9 */
    0, /* : */
    0, /* ; */
    0, /* < */
    7, /* = */
    0, /* > */
    0, /* ? */
    0, /* @ */
    0, /* A */
    0, /* B */
    0, /* C */
    2, /* D */
    0, /* E */
    0, /* F */
    0, /* G */
    5, /* H */
    1, /* I */
    0, /* J */
    0, /* K */
    0, /* L */
    0, /* M */
    3, /* N */
    0, /* O */
    6, /* P */
    0, /* Q */
    0, /* R */
    4, /* S */
    0, /* T */
    0, /* U */
    0, /* V */
    0, /* W */
    8, /* X */
    0, /* Y */
    0, /* Z */
    0, /* [ */
    0, /* \ */
    0, /* ] */
    0, /* ^ */
    0, /* _ */
    0, /* ` */
    0, /* a */
    0, /* b */
    0, /* c */
    0, /* d */
    0, /* e */
    0, /* f */
    0, /* g */
    0, /* h */
    0, /* i */
    0, /* j */
    0, /* k */
    0, /* l */
    0, /* m */
    0, /* n */
    0, /* o */
    0, /* p */
    0, /* q */
    0, /* r */
    0, /* s */
    0, /* t */
    0, /* u */
    0, /* v */
    0, /* w */
    0, /* x */
    0, /* y */
    0, /* z */
    0, /* { */
    0, /* | */
    0, /* } */
    0, /* ~ */
    0, /*  */
];

pub fn query_profile_byte(
    read_num: &[i8],
    mat: &[i8],
    read_len: i32,
    n: i32, // the edge length of the square matrix mat
    bias: u8,
) -> Vec<u8x16> {
    let seg_len = (read_len + 15) / 16; // Split the 128 bit register into 16 pieces.
                                        // Each piece is 8 bit. Split the read into 16 segments.
                                        // Calculate 16 segments in parallel.
    let mut v_profile = vec![u8x16::splat(0); (n * seg_len) as usize];
    let mut t = 0;
    let mut tmpv: [u8; 16] = [0; 16];

    // Generate query profile, rearrange query sequence & calculate the weight of match/mismatch
    for nt in 0..n {
        for i in 0..seg_len {
            let mut j = i;
            for _seg_num in 0..16 {
                tmpv[_seg_num] = if j >= read_len {
                    bias
                } else {
                    mat[(nt as usize) * (n as usize) + (read_num[j as usize] as usize)] as u8 + bias
                };
                j += seg_len;
            }
            v_profile[t] = u8x16::from_array(tmpv);
            t += 1;
        }
    }

    v_profile
}
//
// // static zero_u8x16: u8x16 = u8x16::splat(0);
// // const ZERO_I16X8: i16x8 = i16x8::from_array([0,0,0,0,0,0,0,0]);
// // const _0XFF_U8X16: u8x16 = u8x16::from_array([0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF]);
// // const MSB_MASK: u8x16 = u8x16::from_array([0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80]);
//
// // const RSHIFT_16_8: [Which; 16] = [Second(0), Second(1), Second(2), Second(3), Second(4), Second(5), Second(6), Second(7), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7)];
// // const RSHIFT_16_4: [Which; 16] = [Second(0), Second(1), Second(2), Second(3), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11)];
// // const RSHIFT_16_2: [Which; 16] = [Second(0), Second(1), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11), First(12), First(13)];
// // const RSHIFT_16_1: [Which; 16] = [Second(0), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11), First(12), First(13), First(14)];
//
// // const LSHIFT_16_8: [Which; 16] = [First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1), Second(2), Second(3), Second(4), Second(5), Second(6), Second(7)];
// // const LSHIFT_16_4: [Which; 16] = [First(4), First(5), First(6), First(7),  First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1), Second(2), Second(3)];
// // const LSHIFT_16_2: [Which; 16] = [First(2), First(3), First(4), First(5), First(6), First(7),  First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1)];
//
// // const RSHIFT_i16X8_8: [Which; 16] = [Second(0), Second(1), Second(2), Second(3), Second(4), Second(5), Second(6), Second(7), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7)];
// // const RSHIFT_i16X8_4: [Which; 16] = [Second(0), Second(1), Second(2), Second(3), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11)];
// // const RSHIFT_i16X8_2: [Which; 16] = [Second(0), Second(1), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11), First(12), First(13)];
// // const RSHIFT_i16X8_1: [Which; 16] = [Second(0), First(0), First(1), First(2), First(3), First(4), First(5), First(6), First(7), First(8), First(9), First(10), First(11), First(12), First(13), First(14)];
//
// // const LSHIFT_i16X8_8: [Which; 16] = [First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1), Second(2), Second(3), Second(4), Second(5), Second(6), Second(7)];
// // const LSHIFT_i16X8_4: [Which; 16] = [First(4), First(5), First(6), First(7),  First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1), Second(2), Second(3)];
// // const LSHIFT_i16X8_2: [Which; 16] = [First(2), First(3), First(4), First(5), First(6), First(7),  First(8), First(9), First(10), First(11), First(12), First(13), First(14), First(15), Second(0), Second(1)];
//
pub fn sw_sse2_byte(
    ref_seq: &[i8],
    ref_dir: i8, // 0: forward ref; 1: reverse ref
    ref_len: i32,
    read_len: i32,
    weight_gap_o: u8, // will be used as -
    weight_gap_e: u8, // will be used as -
    v_profile: &Vec<u8x16>,
    terminate: u8, // the best alignment score: used to terminate
    // the matrix calculation when locating the
    // alignment beginning point. If this score
    // is set to 0, it will not be used
    bias: u8, // Shift 0 point to a positive value.
    mask_len: i32,
) -> [AlignmentEnd; 2] {
    // Some helper macros to be used later

    let mut this_max = 0; // the max alignment score
    let mut end_ref = -1; // 0_based best alignment ending point; Initialized as isn't aligned -1.
    let seg_len = (read_len + 15) / 16; // number of segment
                                        //
    const lshift_16_1: [Which; 16] = [
        Second(0),
        First(0),
        First(1),
        First(2),
        First(3),
        First(4),
        First(5),
        First(6),
        First(7),
        First(8),
        First(9),
        First(10),
        First(11),
        First(12),
        First(13),
        First(14),
    ];

    // array to record the largest score of each reference position
    let mut max_column = vec![0u8; ref_len as usize];

    let mut pv_h_store = vec![u8x16::splat(0); seg_len as usize];
    let mut pv_h_load = vec![u8x16::splat(0); seg_len as usize];
    let mut pv_e = vec![u8x16::splat(0); seg_len as usize];
    let mut pv_hmax = vec![u8x16::splat(0); seg_len as usize];

    // 16 byte insertion begin vector
    let v_gap_o = u8x16::splat(weight_gap_o);

    // 16 byte insertion extension vector
    let v_gap_e = u8x16::splat(weight_gap_e);

    // 16 byte bias vector
    let v_bias = u8x16::splat(bias);
    let zero_u8x16: u8x16 = u8x16::splat(0);

    let mut v_max_score = zero_u8x16; // Trace the highest score of the whole SW matrix.
    let mut v_max_mark = zero_u8x16; // Trace the highest score till the previous column.
    let mut begin: i32 = 0;
    let mut end: i32 = ref_len;
    let mut step: i32 = 1;

    // outer loop to process the reference sequence
    if ref_dir == 1 {
        begin = ref_len - 1;
        end = -1;
        step = -1;
    }
    for i in (begin..=end).step_by(step as usize) {
        let (mut e, mut v_f, mut v_max_column) = (zero_u8x16, zero_u8x16, zero_u8x16);
        let mut v_h = pv_h_store[(seg_len - 1) as usize];
        v_h = simd_swizzle!(v_h, zero_u8x16, lshift_16_1);
        let v_p = &v_profile[ref_seq[i as usize] as usize * seg_len as usize..]; // Right part of the vProfile
                                                                                 //
                                                                                 // Swap the 2 H buffers.
        mem::swap(&mut pv_h_load, &mut pv_h_store);

        // inner loop to process the query sequence
        for j in 0..seg_len as usize {
            v_h = v_h.saturating_add(v_p[j]);
            v_h = v_h.saturating_sub(v_bias); // vH will be always > 0

            // Get max from vH, vE and vF.
            e = pv_e[j];
            v_h = v_h.simd_max(e);
            v_h = v_h.simd_max(v_f);
            v_max_column = v_max_column.simd_max(v_h);

            // Save vH values.
            pv_h_store[j] = v_h;

            // Update vE value.
            v_h = v_h.saturating_sub(v_gap_o); // saturation arithmetic, result >= 0
            e = e.saturating_sub(v_gap_e);
            e = e.simd_max(v_h);
            pv_e[j] = e;

            // Update vF value.
            v_f = v_f.saturating_sub(v_gap_e);
            v_f = v_f.simd_max(v_h);

            // Load the next vH.
            v_h = pv_h_load[j];
        }

        // Lazy_F loop: has been revised to disallow adjecent insertion and then deletion, so don't update E(i, j), learn from SWPS3
        let mut shouldbreak = false;
        for k in 0..16 {
            v_f = simd_swizzle!(v_f, zero_u8x16, lshift_16_1);
            for j in 0..seg_len as usize {
                v_h = pv_h_store[j];
                v_h = v_h.simd_max(v_f);
                v_max_column = v_max_column.simd_max(v_h); // newly added line
                pv_h_store[j] = v_h;
                v_h = v_h.saturating_sub(v_gap_o);
                v_f = v_f.saturating_sub(v_gap_e);
                if (v_f.saturating_sub(v_h)).simd_eq(zero_u8x16).all() {
                    shouldbreak = true;
                    break;
                }
            }
            if shouldbreak {
                break;
            }
        }

        v_max_score = v_max_score.simd_max(v_max_column);
        let v_temp = v_max_mark.simd_eq(v_max_score);
        if !v_temp.all() {
            let temp = v_max_score.reduce_max();
            v_max_mark = v_max_score;
            if temp > this_max {
                this_max = temp;
                if this_max + bias >= 255 {
                    break;
                } //overflow
                end_ref = i;
                // Store the column with the highest alignment score in order to trace the alignment ending position on read.
                for j in 0..seg_len as usize {
                    pv_hmax[j] = pv_h_store[j];
                }
            }
        } else {
        }

        // Record the max score of the current column.
        max_column[i as usize] = v_max_column.reduce_max();
        if max_column[i as usize] == terminate {
            break;
        }
    }
    // Trace the alignment ending position on read.
    let t: Vec<u8> = pv_hmax
        .iter()
        .flat_map(|x| x.as_array().iter().cloned())
        .collect();
    let column_len = seg_len * 16;
    let mut end_read = read_len - 1;
    for (i, t_item) in (0..column_len).zip(t.iter()) {
        if *t_item == this_max {
            let temp = i / 16 + i % 16 * seg_len;
            if temp < end_read {
                end_read = temp;
            }
        }
    }

    // Find the most possible 2nd best alignment.
    let mut bests: [AlignmentEnd; 2] = [
        AlignmentEnd {
            score: if this_max + bias >= 255 {
                255
            } else {
                this_max as u16
            },
            ref_position: end_ref,
            read_position: end_read,
        },
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
    ];

    let mut edge: i32 = if end_ref - mask_len > 0 {
        end_ref - mask_len
    } else {
        0
    };
    for i in 0..edge {
        if max_column[i as usize] as u16 > bests[1].score {
            bests[1].score = max_column[i as usize] as u16;
            bests[1].ref_position = i;
        }
    }
    edge = if end_ref + mask_len > ref_len {
        ref_len
    } else {
        end_ref + mask_len
    };
    for i in edge + 1..ref_len {
        if max_column[i as usize] as u16 > bests[1].score {
            bests[1].score = max_column[i as usize] as u16;
            bests[1].ref_position = i;
        }
    }
    bests
}

pub fn query_profile_word(read_num: &[i8], mat: &[i8], read_len: usize, n: usize) -> Vec<i16x8> {
    let seg_len = (read_len + 7) / 8;
    let mut v_profile = vec![i16x8::splat(0); n * seg_len];
    let mut tmpv: [i16; 8] = [0; 8];
    let mut t = 0;
    for nt in 0..n {
        for i in 0..seg_len {
            let mut j = i;
            for seg_num in 0..8 {
                if j < read_len {
                    tmpv[seg_num] = mat[nt * n + read_num[j] as usize] as i16;
                } else {
                    tmpv[seg_num] = 0;
                }
                j += seg_len;
            }
            v_profile[t] = i16x8::from_array(tmpv);
            t += 1;
        }
    }
    v_profile
}

pub fn sw_sse2_word(
    ref_seq: &[i8],
    ref_dir: i8, // 0: forward ref; 1: reverse ref
    ref_len: i32,
    read_len: i32,
    weight_gap_o: u8, // will be used as -
    weight_gap_e: u8, // will be used as -
    v_profile: &[i16x8],
    terminate: u16,
    mask_len: i32,
) -> [AlignmentEnd; 2] {
    let mut _max: u16 = 0;
    let mut end_ref = 0;
    let seg_len = (read_len + 7) / 8;
    let zero_i16x8: i16x8 = i16x8::splat(0);

    let mut max_column = vec![0u16; ref_len as usize];

    let mut pv_h_store = vec![zero_i16x8; seg_len as usize];
    let mut pv_h_load = vec![zero_i16x8; seg_len as usize];
    let mut pv_e = vec![zero_i16x8; seg_len as usize];
    let mut pv_h_max = vec![zero_i16x8; seg_len as usize];

    // let v_gap_o = i16x8::splat(weight_gapo as i16);
    // let v_gap_e = i16x8::splat(weight_gape as i16);
    let v_gap_o = i16x8::splat(weight_gap_o as i16);
    let v_gap_e = i16x8::splat(weight_gap_e as i16);

    const lshift_8_2: [Which; 8] = [
        Second(0),
        Second(1),
        First(0),
        First(1),
        First(2),
        First(3),
        First(4),
        First(5),
    ];

    let mut v_max_score = zero_i16x8;
    let mut v_max_mark = zero_i16x8;

    let (edge, begin, end, step) = if ref_dir == 1 {
        (0, ref_len - 1, -1, -1)
    } else {
        (ref_len, 0, ref_len, 1)
    };
    const lshift_8_1: [Which; 8] = [
        Second(0),
        First(0),
        First(1),
        First(2),
        First(3),
        First(4),
        First(5),
        First(6),
    ];

    // Outer loop to process the reference sequence
    let mut i = begin;
    while i != end {
        let mut e: i16x8;
        let mut v_f = zero_i16x8;
        let mut v_h = pv_h_store[seg_len as usize - 1];
        v_h = simd_swizzle!(v_h, zero_i16x8, lshift_8_1);

        let mut v_max_column = zero_i16x8;
        let v_p = &v_profile[ref_seq[i as usize] as usize * seg_len as usize..]; // Right part of the vProfile

        mem::swap(&mut pv_h_load, &mut pv_h_store);

        for j in 0..seg_len {
            v_h = v_h.saturating_add(v_p[j as usize]);

            e = pv_e[j as usize];
            v_h = v_h.simd_max(e);
            v_h = v_h.simd_max(v_f);
            v_max_column = v_max_column.simd_max(v_h);

            pv_h_store[j as usize] = v_h;

            v_h = v_h.saturating_sub(v_gap_o).simd_max(zero_i16x8);
            e = e.saturating_sub(v_gap_e).simd_max(zero_i16x8);
            e = e.simd_max(v_h);
            pv_e[j as usize] = e;

            v_f = v_f.saturating_sub(v_gap_e).simd_max(zero_i16x8);
            v_f = v_f.simd_max(v_h);

            v_h = pv_h_load[j as usize];
        }
        // Lazy_F loop
        let mut shouldbreak = false;
        for _ in 0..8 {
            v_f = simd_swizzle!(v_f, zero_i16x8, lshift_8_1);
            for j in 0..seg_len {
                v_h = pv_h_store[j as usize];
                v_h = v_h.simd_max(v_f);
                v_max_column = v_max_column.simd_max(v_h);
                pv_h_store[j as usize] = v_h;
                v_h = v_h.saturating_sub(v_gap_o).simd_max(zero_i16x8);
                v_f = v_f.saturating_sub(v_gap_e).simd_max(zero_i16x8);
                if v_f.simd_le(v_h).all() {
                    shouldbreak = true;
                    break;
                }
            }
            if shouldbreak {
                break;
            }
        }

        v_max_score = v_max_score.simd_max(v_max_column);
        let v_temp = v_max_mark.simd_eq(v_max_score);
        if !v_temp.all() {
            let temp = v_max_score.reduce_max() as u16;
            v_max_mark = v_max_score;

            if temp > _max {
                _max = temp;
                end_ref = i;
                for j in 0..seg_len {
                    pv_h_max[j as usize] = pv_h_store[j as usize];
                }
            }
        }

        max_column[i as usize] = max(v_max_column.reduce_max() as i16, 0) as u16;
        if max_column[i as usize] as u16 == terminate {
            break;
        }
        i += step;
    }

    let t: Vec<i16> = pv_h_max
        .iter()
        .flat_map(|x| x.as_array().iter().cloned())
        .collect();
    let column_len = seg_len * 8;
    let mut end_read = read_len - 1;
    for (i, t_item) in (0..column_len).zip(t.iter()) {
        if *t_item as u16 == _max {
            let temp = i / 8 + i % 8 * seg_len;
            if temp < end_read {
                end_read = temp;
            }
        }
    }

    let mut bests: [AlignmentEnd; 2] = [
        AlignmentEnd {
            score: _max,
            ref_position: end_ref,
            read_position: end_read,
        },
        AlignmentEnd {
            score: 0,
            ref_position: 0,
            read_position: 0,
        },
    ];

    let edge = if (end_ref - mask_len) > 0 {
        end_ref - mask_len
    } else {
        0
    };
    for i in 0..edge {
        if max_column[i as usize] as u16 > bests[1].score {
            bests[1].score = max_column[i as usize] as u16;
            bests[1].ref_position = i;
        }
    }
    let edge = if (end_ref + mask_len) > ref_len {
        ref_len
    } else {
        end_ref + mask_len
    };
    for i in edge..ref_len {
        if max_column[i as usize] as u16 > bests[1].score {
            bests[1].score = max_column[i as usize] as u16;
            bests[1].ref_position = i;
        }
    }
    bests
}

pub fn banded_sw(
    ref_: &[i8],
    read: &[i8],
    ref_len: i32,
    read_len: i32,
    score: i32,
    weight_gap_o: u32,
    weight_gap_e: u32,
    mut band_width: i32,
    mat: &[i8],
    n: i32,
) -> Option<Cigar> {
    let mut c: Vec<u32> = vec![0; 16];
    let mut i;
    let mut j;
    let mut f;
    let mut e = 0;
    let mut temp1 = 0;
    let mut temp2;
    let mut s = 16;
    let mut s1: i32 = 8;
    let mut s2: u32 = 1024;
    // let mut l;
    let mut max_val = 0;
    let mut op;
    let mut prev_op;
    let mut h_b: Vec<i32>;
    let mut e_b: Vec<i32>;
    let mut h_c: Vec<i32>;
    let mut direction: Vec<i8>;
    let mut direction_line_start: usize = 0;

    // let mut direction_line: &[i8];
    // let mut direction: Vec<i8>;
    // let &mut direction_line: Vec<i8>;

    let len = max(ref_len, read_len);
    h_b = vec![0; s1 as usize];
    e_b = vec![0; s1 as usize];
    h_c = vec![0; s1 as usize];

    let mut width_d;

    loop {
        let width = band_width * 2 + 3;
        width_d = band_width * 2 + 1;

        while width >= s1 {
            s1 += 1;
            s1 = kroundup32(s1);
            h_b.resize(s1 as usize, 0);
            e_b.resize(s1 as usize, 0);
            h_c.resize(s1 as usize, 0);
        }

        // if width_d * read_len * 3 >= direction.len() as i32 {
        //     direction.resize((width_d * read_len * 3) as usize, 0);
        // }
        while (width_d * read_len * 3) as u32 >= s2 {
            // width_d*readLen* overflow before s2
            s2 += 1;
            s2 = kroundup32(s2);
            // direction = (int8_t*)realloc(direction, s2 * sizeof(int8_t));
        }

        direction = vec![0; s2 as usize];

        for j in 1..(width - 1) {
            h_b[j as usize] = 0;
        }

        for i in 0..read_len {
            let (mut beg, mut end) = (0, ref_len - 1);
            let mut u: i32 = 0;

            j = i - band_width;
            beg = max(beg, j);
            j = i + band_width;
            end = min(end, j);
            let edge = min(end + 1, width - 1);

            f = 0;
            h_b[0] = 0;
            e_b[0] = 0;
            h_b[edge as usize] = 0;
            e_b[edge as usize] = 0;
            h_c[0] = 0;

            direction_line_start = (width_d * i * 3) as usize;

            // TODO Might not be right?
            for j in beg..=end {
                u = set_u(band_width, i, j);
                e = set_u(band_width, i - 1, j);
                let b = set_u(band_width, i, j - 1);
                let d = set_u(band_width, i - 1, j - 1);
                let de = set_d(band_width, i, j, 0);
                let df = set_d(band_width, i, j, 1);
                let dh = set_d(band_width, i, j, 2);

                temp1 = if i == 0 {
                    -(weight_gap_o as i32)
                } else {
                    h_b[e as usize] - weight_gap_o as i32
                };
                temp2 = if i == 0 {
                    -(weight_gap_e as i32)
                } else {
                    e_b[e as usize] - weight_gap_e as i32
                };
                e_b[u as usize] = max(temp1, temp2);
                direction[direction_line_start + de as usize] = if temp1 > temp2 { 3 } else { 2 };

                temp1 = h_c[b as usize] - weight_gap_o as i32;
                temp2 = f - weight_gap_e as i32;
                f = max(temp1, temp2);
                direction[direction_line_start + df as usize] = if temp1 > temp2 { 5 } else { 4 };

                let e1 = max(e_b[u as usize], 0);
                let f1 = max(f, 0);
                temp1 = max(e1, f1);
                temp2 = h_b[d as usize]
                    + mat[ref_[j as usize] as usize * n as usize + (read[i as usize]) as usize]
                        as i32;
                h_c[u as usize] = max(temp1, temp2);

                if h_c[u as usize] > max_val {
                    max_val = h_c[u as usize];
                }

                if temp1 <= temp2 {
                    direction[direction_line_start + dh as usize] = 1;
                } else {
                    direction[direction_line_start + dh as usize] = if e1 > f1 {
                        direction[direction_line_start + de as usize]
                    } else {
                        direction[direction_line_start + df as usize]
                    };
                }
            }

            for j in 1..=u {
                h_b[j as usize] = h_c[j as usize];
            }
        }

        band_width *= 2;

        if max_val >= score || band_width > len {
            break;
        }
    }

    band_width /= 2;

    // Trace back
    i = read_len - 1;
    j = ref_len - 1;
    e = 0;
    let mut l = 0;
    op = 'M';
    prev_op = 'M';
    temp2 = 2;

    while i >= 0 && j > 0 {
        temp1 = set_d(band_width, i, j, temp2);
        match direction[direction_line_start + temp1 as usize] {
            1 => {
                i -= 1;
                j -= 1;
                temp2 = 2;
                direction_line_start -= (width_d * 3) as usize;
                op = 'M';
            }
            2 => {
                i -= 1;
                temp2 = 0;
                direction_line_start -= (width_d * 3) as usize;
                op = 'I';
            }
            3 => {
                i -= 1;
                temp2 = 2;
                direction_line_start -= (width_d * 3) as usize;
                op = 'I';
            }
            4 => {
                j -= 1;
                temp2 = 1;
                op = 'D';
            }
            5 => {
                j -= 1;
                temp2 = 2;
                op = 'D';
            }
            _ => {
                eprintln!("Trace back error: {}.", direction[direction_line_start]);
                return None;
            }
        }

        if op == prev_op {
            e += 1;
        } else {
            l += 1;
            if l >= s {
                s += 1;
                s = kroundup32(s);
            }
            c.resize(s as usize, 0);
            c[l as usize - 1] = to_cigar_int(e as u32, prev_op);
            prev_op = op;
            e = 1;
        }
    }

    if op == 'M' {
        l += 1;
        if l >= s {
            s += 1;
            s = kroundup32(s);
        }
        c.resize(s as usize, 0);
        c[l as usize - 1] = to_cigar_int((e + 1) as u32, op);
    } else {
        l += 2;
        if l >= s {
            s += 1;
            s = kroundup32(s);
        }
        c.resize(s as usize, 0);
        c[l as usize - 2] = to_cigar_int(e as u32, op);
        c[l as usize - 1] = to_cigar_int(1 as u32, 'M');
    }

    // Reverse cigar
    let mut c1 = vec![0; l as usize];
    let mut s = 0;
    let mut e = l - 1;
    while s <= e {
        c1[s as usize] = c[e as usize];
        c1[e as usize] = c[s as usize];
        s += 1;
        e -= 1;
    }

    let result = Cigar {
        seq: c1,
        length: l as usize,
    };

    // let result = Cigar {
    //     seq: Vec::new(),
    //     length: 0,
    // };
    Some(result)
}

#[cfg(test)]
mod tests {
    use std::simd::SimdPartialOrd;

    use super::*;
    // use portable_simd::u8x16;
    // use indexmap::IndexMap;
    // static USAGE_EXAMPLE_CIF_TEXT: &[u8] = include_bytes!("../tests/data/usage-example.cif.json");
    // static USAGE_EXAMPLE_CIF_TEXT: &str  = include_str!("../tests/data/usage-example.cif");
    // static USAGE_EXAMPLE_CIF_JSON_TEXT: &[u8] = include_bytes!("../tests/data/usage-example.cif.json");

    #[test]
    fn simple_01() {
        let vec_u8x16: Vec<u8x16> = vec![
            u8x16::from_array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]),
            u8x16::from_array([
                17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
            ]),
        ];

        let vec_u8: Vec<u8> = vec_u8x16
            .iter()
            .flat_map(|x| x.as_array().iter().cloned())
            .collect();

        println!("{:?}", vec_u8);
        let a = u8x16::from_array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]);
        let b = u8x16::from_array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        // let shift = 4; // Shift by 4 bytes (32 bits)
        let result = simd_swizzle!(
            a,
            b,
            [
                Second(0),
                Second(1),
                Second(2),
                Second(3),
                Second(4),
                Second(5),
                Second(6),
                Second(7),
                First(0),
                First(1),
                First(2),
                First(3),
                First(4),
                First(5),
                First(6),
                First(7)
            ]
        );
        let v_0x80 = u8x16::splat(0x80);
        // let result = a.rotate_lanes_right::<8>().as_array().clone(); // Shift right by 4 bytes
        println!("{:?}", result);
        println!("{:?}", a.reduce_max());
        let aa = u8x16::from_array([
            0x80, 0x7F, 0x80, 0x7F, 0x80, 0x7F, 0x80, 0x7F, 0x80, 0x7F, 0x80, 0x7F, 0x80, 0x7F,
            0x80, 0x7F,
        ]);
    }
}
