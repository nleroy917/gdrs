use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;

pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32,
}

impl PartialEq for Region {
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.start == other.start && self.end == other.end
    }
}

pub struct RegionSet {
    regions: HashMap<String, Vec<Region>>,
    sorted: bool,
}

impl RegionSet {
    pub fn from_bed_file(value: &Path) -> Result<RegionSet> {
        let file = File::open(value)?;
        let reader = BufReader::new(file);

        let mut regions = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let fields = line.split('\t').collect::<Vec<&str>>();
            let chr = fields[0];
            let start = fields[1].parse::<u32>()?;
            let end = fields[2].parse::<u32>()?;

            let region = Region {
                chr: chr.to_string(),
                start,
                end,
            };

            regions
                .entry(chr.to_string())
                .or_insert(Vec::new())
                .push(region);
        }

        Ok(RegionSet {
            regions,
            sorted: false,
        })
    }

    pub fn iter_choms(&self) -> impl Iterator<Item = &String> {
        self.regions.keys()
    }

    pub fn iter_regions(&self, chr: &str) -> impl Iterator<Item = &Region> {
        self.regions.get(chr).unwrap().iter()
    }

    pub fn sort(&mut self) {
        for regions in self.regions.values_mut() {
            regions.sort_by(|a, b| a.start.cmp(&b.start));
        }

        // set the sorted flag to true
        self.sorted = true;
    }

    pub fn is_sorted(&self) -> bool {
        self.sorted
    }
}
