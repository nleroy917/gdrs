use pyo3::prelude::*;

use gdrs::{calc_gc_content, models::RegionSet};

use std::path::Path;

use crate::models::PyGenomeAssembly;

#[pyfunction(name="calc_gc_content")]
pub fn py_calc_gc_content(file: String, genome: &PyGenomeAssembly) -> f64 {    
    let path = Path::new(&file);
    let rs = RegionSet::from_bed(path).unwrap();

    let gc_content = calc_gc_content(&rs, &genome.genome_assembly).unwrap();

    gc_content
}