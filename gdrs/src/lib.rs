use std::collections::HashMap;

use anyhow::{ensure, Result};

pub mod models;

use models::{Dinucleotide, GenomeAssembly, RegionSet, TSSIndex};

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

pub fn calc_gc_content(
    region_set: &RegionSet,
    genome: &GenomeAssembly,
    ignore_unk_chroms: bool,
) -> Result<Vec<f64>> {
    // for region in region_set
    let mut gc_contents: Vec<f64> = vec![];
    for chr in region_set.iter_chroms() {

        // check if the chrom is even in genome
        if ignore_unk_chroms && !genome.contains_chr(chr) {
            continue;
        }

        for region in region_set.iter_regions(chr) {

            let mut gc_count: u32 = 0;
            let mut total_count: u32 = 0;
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
                    gc_contents.push(gc_count as f64 / total_count as f64);
                }
                Err(e) => {
                    if ignore_unk_chroms {
                        continue;
                    } else {
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
    }

    Ok(gc_contents)
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

pub fn calc_tss_dist(region_set: &RegionSet, tss_index: &TSSIndex) -> Result<Vec<u32>> {
    let mut tss_dists: Vec<u32> = Vec::with_capacity(region_set.len());

    for chr in region_set.iter_chroms() {
        if !tss_index.has_chr(chr) {
            continue;
        }
        for region in region_set.iter_regions(chr) {
            let tsses = tss_index.query(region);
            ensure!(
                tsses.is_some(),
                format!(
                    "No TSS's found for region: {}:{}-{}. Double-check your index!",
                    region.chr, region.start, region.end
                )
            );

            let midpoint = region.end - region.start;

            let dists = tsses.unwrap().into_iter().map(|tss| {
                let tss_midpoint = tss.end - tss.start;
                midpoint - tss_midpoint
            });

            tss_dists.push(dists.min().unwrap());
        }
    }

    Ok(tss_dists)
}

pub mod prelude {
    pub use super::calc_dinucl_freq;
    pub use super::calc_gc_content;
    pub use super::calc_neighbor_distances;
    pub use super::calc_tss_dist;
    pub use super::calc_widths;
    pub use super::models::{GenomeAssembly, Region, RegionSet, TSSIndex};
}
