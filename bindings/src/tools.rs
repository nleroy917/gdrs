use pyo3::prelude::*;

use std::path::Path;

use crate::models::PyGenomeAssembly;

#[pyfunction(name = "calc_gc_content")]
pub fn py_calc_gc_content(file: String, genome: &PyGenomeAssembly) -> anyhow::Result<f64> {
    let path = Path::new(&file);
    let rs = gdrs::models::RegionSet::from_bed(path)?;

    gdrs::calc_gc_content(&rs, &genome.genome_assembly)
}

#[pyfunction(name = "calc_neighbor_distances")]
pub fn py_calc_neighbor_distances(file: String) -> anyhow::Result<Vec<u32>> {
    let path = Path::new(&file);
    let rs = gdrs::models::RegionSet::from_bed(path)?;

    gdrs::calc_neighbor_distances(&rs)
}
