use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use anyhow::ensure;
use flate2::read::GzDecoder;

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
    pub fn from_bed(value: &Path) -> Result<RegionSet> {
        let is_gzipped = value.extension() == Some(OsStr::new("gz"));
        let file = File::open(value)?;

        let file: Box<dyn Read> = match is_gzipped {
            true => Box::new(GzDecoder::new(file)),
            false => Box::new(file),
        };

        let reader = BufReader::new(file);

        let mut regions = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let fields = line.split('\t').collect::<Vec<&str>>();

            ensure!(fields.len() >= 3, "Invalid BED file format!");

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

    pub fn iter_chroms(&self) -> impl Iterator<Item = &String> {
        self.regions.keys()
    }

    pub fn iter_regions(&self, chr: &str) -> impl Iterator<Item = &Region> {
        self.regions.get(chr).unwrap().iter()
    }

    pub fn into_sorted(self) -> RegionSet {
        let mut regions = self.regions;
        for region_vec in regions.values_mut() {
            region_vec.sort_by(|a, b| a.start.cmp(&b.start));
        }

        RegionSet {
            regions,
            sorted: true,
        }
    }

    pub fn is_sorted(&self) -> bool {
        self.sorted
    }
}
