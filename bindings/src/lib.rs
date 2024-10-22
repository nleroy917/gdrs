use pyo3::prelude::*;

mod models;
mod tools;

#[pymodule]
fn gdrs(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<models::PyGenomeAssembly>()?;
    m.add_function(wrap_pyfunction!(tools::py_calc_widths, m)?)?;
    m.add_function(wrap_pyfunction!(tools::py_calc_gc_content, m)?)?;
    m.add_function(wrap_pyfunction!(tools::py_calc_neighbor_distances, m)?)?;
    m.add_function(wrap_pyfunction!(tools::py_calc_dinucleotide_frequency, m)?)?;
    Ok(())
}
