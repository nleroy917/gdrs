// use std::fs::File;
// use std::path::Path;

use anyhow::{ensure, Result};
// use bio::io::fasta;

pub mod models;

use models::RegionSet;

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

// pub fn gc_content(region_set: &RegionSet, genome: &Path) -> Result<f64> {
//     let mut gc_count: u32 = 0;
//     let mut total_count: u32 = 0;

//     // open the genome file
//     let file = File::open(genome)?;
//     let genome = fasta::Reader::new(file);

//     let records = genome.records();

//     // for region in region_set

//     // get sequence from chr:start-end

//     // count num Gs and num Cs, increment gc_count

//     // increment total count

//     Ok(gc_count as f64 / total_count as f64)
// }

pub mod prelude {
    pub use super::calc_neighbor_distances;
    pub use super::models::{Region, RegionSet};
}
