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
        for i in 0..region_set.iter_regions(chr).count() - 1 {
            let region = region_set.iter_regions(chr).nth(i).unwrap();
            let next_region = region_set.iter_regions(chr).nth(i + 1).unwrap();

            let distance = next_region.start - region.end;
            distances.push(distance);
            
        } 
    }

    Ok(distances)
}

#[cfg(test)]
mod tests {
    use super::*;

    
}
