use pyo3::prelude::*;

use std::collections::HashMap;
use std::path::Path;

use crate::models::PyGenomeAssembly;

#[pyfunction(name = "calc_gc_content")]
pub fn py_calc_gc_content(file: String, genome: &PyGenomeAssembly) -> anyhow::Result<Vec<f64>> {
    let path = Path::new(&file);
    let rs = gdrs::models::RegionSet::from_bed(path)?;

    gdrs::calc_gc_content(&rs, &genome.genome_assembly)
}

#[pyfunction(name = "calc_neighbor_distances")]
pub fn py_calc_neighbor_distances(file: String) -> anyhow::Result<Vec<u32>> {
    let path = Path::new(&file);
    let rs = gdrs::models::RegionSet::from_bed(path)?;

    let rs = rs.into_sorted();

    gdrs::calc_neighbor_distances(&rs)
}

#[pyfunction(name = "calc_dincleotide_frequency")]
pub fn py_calc_dinucleotide_frequency(
    file: String,
    genome: &PyGenomeAssembly,
) -> anyhow::Result<HashMap<String, f64>> {
    let path = Path::new(&file);
    let rs = gdrs::models::RegionSet::from_bed(path)?;

    let frequencies = gdrs::calc_dinucl_freq(&rs, &genome.genome_assembly)?;

    let mut freq_map: HashMap<String, f64> = HashMap::new();

    // Convert Dinucleotide to String and push to HashMap
    for (di, freq) in frequencies {
        freq_map.insert(di.to_string()?, freq);
    }

    Ok(freq_map)
}
