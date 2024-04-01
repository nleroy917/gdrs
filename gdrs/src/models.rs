use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;

use anyhow::ensure;
use anyhow::Result;
use bio::io::fasta;
use flate2::read::GzDecoder;

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

            ensure!(
                fields.len() >= 3,
                "Invalid BED file format found. File lacks the three necessary columns."
            );

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

pub struct GenomeAssembly {
    seq_map: HashMap<String, Vec<u8>>,
}

impl GenomeAssembly {
    pub fn from_fasta(path: &Path) -> Result<GenomeAssembly> {
        let file = File::open(path)?;
        let genome = fasta::Reader::new(file);

        let records = genome.records();

        // store the genome in a hashmap
        let mut seq_map: HashMap<String, Vec<u8>> = HashMap::new();
        for record in records {
            match record {
                Ok(record) => {
                    seq_map.insert(record.id().to_string(), record.seq().to_owned());
                }
                Err(e) => {
                    return Err(anyhow::anyhow!("Error reading genome file: {}", e));
                }
            }
        }

        Ok(GenomeAssembly { seq_map })
    }

    pub fn seq_from_region<'a>(&'a self, coords: &Region) -> Result<&'a [u8]> {
        let chr = &coords.chr;
        let start = coords.start;
        let end = coords.end;

        let seq = self.seq_map.get(chr);

        match seq {
            Some(seq) => {
                let seq = &seq[start as usize..end as usize];
                Ok(seq)
            }
            None => Err(anyhow::anyhow!(
                "Unknown chromosome found in region set: {}",
                chr.to_string()
            )),
        }
    }
}

#[derive(PartialEq, Eq, Hash)]
pub enum Dinucleotide {
    Aa,
    Ac,
    Ag,
    At,
    Ca,
    Cc,
    Cg,
    Ct,
    Ga,
    Gc,
    Gg,
    Gt,
    Ta,
    Tc,
    Tg,
    Tt,
}

impl Dinucleotide {
    pub fn from_bytes(bytes: &[u8]) -> Option<Dinucleotide> {
        match bytes {
            b"Aa" => Some(Dinucleotide::Aa),
            b"Ac" => Some(Dinucleotide::Ac),
            b"Ag" => Some(Dinucleotide::Ag),
            b"At" => Some(Dinucleotide::At),
            b"Ca" => Some(Dinucleotide::Ca),
            b"Cc" => Some(Dinucleotide::Cc),
            b"Cg" => Some(Dinucleotide::Cg),
            b"Ct" => Some(Dinucleotide::Ct),
            b"Ga" => Some(Dinucleotide::Ga),
            b"Gc" => Some(Dinucleotide::Gc),
            b"Gg" => Some(Dinucleotide::Gg),
            b"Gt" => Some(Dinucleotide::Gt),
            b"Ta" => Some(Dinucleotide::Ta),
            b"Tc" => Some(Dinucleotide::Tc),
            b"Tg" => Some(Dinucleotide::Tg),
            b"Tt" => Some(Dinucleotide::Tt),
            _ => None,
        }
    }
}