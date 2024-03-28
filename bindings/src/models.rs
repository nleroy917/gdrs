use pyo3::prelude::*;
use gdrs::prelude::*;

use std::path::Path;

#[pyclass(name="GenomeAssembly")]
pub struct PyGenomeAssembly {
    pub genome_assembly: GenomeAssembly
}

#[pymethods]
impl PyGenomeAssembly {
    #[new]
    pub fn new(path: String) -> Self {
        let path = Path::new(&path);
        let genome_assembly = GenomeAssembly::from_fasta(path).unwrap();
        PyGenomeAssembly {
            genome_assembly
        }
    }
}