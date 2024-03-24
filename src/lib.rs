use anyhow::Result;

pub mod models;

use models::RegionSet;

pub fn calc_neighnor_distances(region_set: &mut RegionSet) -> Result<Vec<u32>> {
    
    // make sure that the regions are sorted
    if !region_set.is_sorted(){
        region_set.sort();
    }

    let mut distances = vec![];

    // iterate over all chromosomes
    for chr in region_set.iter_choms() {
        // if there is only one region on the chromosome, skip it, can't calculate distance between one region
        if region_set.iter_regions(chr).count() < 2 {
            continue;
        }
        for i in 0..region_set.iter_regions(chr).count() - 1 {
            let region = region_set.iter_regions(chr).nth(i).unwrap();
            let next_region = region_set.iter_regions(chr).nth(i + 1).unwrap();

            let distance = next_region.start - region.end;

            distances.push(distance);
            
        } 
    }

    Ok(distances)
}

pub mod prelude {
    pub use super::models::{Region, RegionSet};
    pub use super::calc_neighnor_distances;
}