#![feature(portable_simd)]

use pyo3::prelude::*;
mod dpf_ssw_aligner;

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyCigar {
    inner: dpf_ssw_aligner::Cigar,
}

#[pymethods]
impl PyCigar {
    pub fn get_seq(&self) -> PyResult<Vec<u32>> {
        Ok(self.inner.seq.clone())
    }
    pub fn get_length(&self) -> PyResult<usize> {
        Ok(self.inner.length)
    }
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct PyAlign {
    inner: dpf_ssw_aligner::Align,
}

#[pymethods]
impl PyAlign {
    pub fn get_score1(&self) -> PyResult<u16> {
        Ok(self.inner.score1)
    }
    pub fn get_score2(&self) -> PyResult<u16> {
        Ok(self.inner.score2)
    }
    pub fn get_ref_begin1(&self) -> PyResult<i32> {
        Ok(self.inner.ref_begin1)
    }
    pub fn get_ref_end1(&self) -> PyResult<i32> {
        Ok(self.inner.ref_end1)
    }
    pub fn get_read_begin1(&self) -> PyResult<i32> {
        Ok(self.inner.read_begin1)
    }
    pub fn get_read_end1(&self) -> PyResult<i32> {
        Ok(self.inner.read_end1)
    }
    pub fn get_ref_end2(&self) -> PyResult<i32> {
        Ok(self.inner.ref_end2)
    }
    pub fn get_cigar(&self) -> PyResult<PyCigar> {
        Ok(PyCigar {
            inner: self.inner.cigar.clone(),
        })
    }
    pub fn get_flag(&self) -> PyResult<u16> {
        Ok(self.inner.flag)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyProfile {
    pub inner: dpf_ssw_aligner::Profile,
}

#[pymethods]
impl PyProfile {
    #[new]
    pub fn ssw_init(read: Vec<i8>, read_len: i32, mat: Vec<i8>, n: i32, score_size: i8) -> Self {
        PyProfile {
            inner: dpf_ssw_aligner::Profile::ssw_init(read, read_len, mat, n, score_size),
        }
    }
}

#[pyfunction]
pub fn py_ssw_align(
    prof: PyProfile,
    ref_seq: Vec<i8>,
    ref_len: i32,
    weight_gap_o: u8,
    weight_gap_e: u8,
    flag: u8,
    filters: u16,
    filterd: i32,
    mask_len: i32,
) -> PyResult<Option<PyAlign>> {
    // does this need to be pyresult?
    let ret = dpf_ssw_aligner::ssw_align(
            prof.inner,
            ref_seq,
            ref_len,
            weight_gap_o,
            weight_gap_e,
            flag,
            filters,
            filterd,
            mask_len,
    );
    if ret.is_some() {
        return Ok(Some( PyAlign { inner: ret.unwrap() } ))
    } else {
        return Ok(None)
    }
}

/// This module is implemented in Rust.
#[pymodule]
#[pyo3(name="_rs_bind")]
fn dpf_ssw_aligner_rspy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCigar>()?;
    m.add_class::<PyAlign>()?;
    m.add_class::<PyProfile>()?;
    m.add_function(wrap_pyfunction!(py_ssw_align, m)?)?;
    Ok(())
}
