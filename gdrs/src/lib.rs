use std::collections::HashMap;

use anyhow::{ensure, Result};

pub mod models;

use models::{Dinucleotide, GenomeAssembly, RegionSet};

pub fn calc_neighbor_distances(region_set: &RegionSet) -> Result<Vec<u32>> {
    // make sure that the regions are sorted
    ensure!(
        region_set.is_sorted(),
        "RegionSet must be sorted to compute neighbor distances!"
    );

    let mut distances = vec![];

    // iterate over all chromosomes
    for chr in region_set.iter_chroms() {
        // if there is only one region on the chromosome, skip it, can't calculate distance between one region
        if region_set.iter_regions(chr).count() < 2 {
            continue;
        }

        let regions: Vec<_> = region_set.iter_regions(chr).collect();
        for window in regions.windows(2) {
            let distance = window[1].start - window[0].end;
            distances.push(distance);
        }
    }

    Ok(distances)
}

pub fn calc_gc_content(region_set: &RegionSet, genome: &GenomeAssembly) -> Result<f64> {
    let mut gc_count: u32 = 0;
    let mut total_count: u32 = 0;

    // for region in region_set
    for chr in region_set.iter_chroms() {
        for region in region_set.iter_regions(chr) {
            let seq = genome.seq_from_region(region);
            match seq {
                Ok(seq) => {
                    for base in seq {
                        match base.to_ascii_lowercase() {
                            b'g' | b'c' => {
                                gc_count += 1;
                            }
                            _ => {}
                        }
                        total_count += 1;
                    }
                }
                Err(e) => {
                    return Err(anyhow::anyhow!(
                        "Error getting sequence for region {}:{}-{}: {}",
                        region.chr.to_string(),
                        region.start,
                        region.end,
                        e
                    ));
                }
            }
        }
    }

    Ok(gc_count as f64 / total_count as f64)
}

pub fn calc_widths(region_set: &RegionSet) -> Result<Vec<u32>> {
    let mut widths: Vec<u32> = Vec::new();
    for chr in region_set.iter_chroms() {
        for region in region_set.iter_regions(chr) {
            widths.push(region.end - region.start)
        }
    }
    Ok(widths)
}

pub fn calc_dinucl_freq(
    region_set: &RegionSet,
    genome: &GenomeAssembly,
) -> Result<HashMap<Dinucleotide, f64>> {

    let mut dinucl_freqs: HashMap<Dinucleotide, f64> = HashMap::new();

    for chr in region_set.iter_chroms() {
        for region in region_set.iter_regions(chr) {
            let seq = genome.seq_from_region(region)?;
            for aas in seq.windows(2) {
                let diucl = Dinucleotide::from_bytes(aas);
                match diucl {
                    Some(dinucl) => {
                        let current_freq = dinucl_freqs.entry(dinucl).or_insert(0.0);
                        *current_freq += 1.0;
                    }
                    None => continue,
                }
            }
        }
    }

    Ok(dinucl_freqs)
}

pub mod prelude {
    pub use super::calc_gc_content;
    pub use super::calc_neighbor_distances;
    pub use super::models::{GenomeAssembly, Region, RegionSet};
}
