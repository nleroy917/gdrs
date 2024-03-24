use std::collections::HashMap;
use std::path::Path;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;

use anyhow::Result;

pub struct Region {
    pub chr: String,
    pub start: u32,
    pub end: u32
}

impl PartialEq for Region {
    fn eq(&self, other: &Self) -> bool {
        self.chr == other.chr && self.start == other.start && self.end == other.end
    }
}

pub struct RegionSet {
    regions: HashMap<String, Vec<Region>>
}

impl RegionSet {
    fn from_bed_file(value: &Path) -> Result<RegionSet> {

        let file = File::open(value)?;
        let mut reader = BufReader::new(file);

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
                end
            };

            regions.entry(chr.to_string()).or_insert(Vec::new()).push(region);
        }

        Ok(RegionSet {
            regions
        })
        
    }
}