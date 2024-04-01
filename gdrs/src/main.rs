use anyhow::{ensure, Context, Result};
use gdrs::models::GenomeAssembly;
use gdrs::prelude::*;
use gdrs::{calc_gc_content, calc_neighbor_distances};
use std::io::stdout;
use std::io::Write;
use std::path::Path;

use clap::{arg, Command};

pub mod consts {
    pub const VERSION: &str = env!("CARGO_PKG_VERSION");
    pub const PKG_NAME: &str = env!("CARGO_PKG_NAME");
    pub const BIN_NAME: &str = env!("CARGO_PKG_NAME");
    pub const ND_CMD: &str = "nd";
    pub const GC_CMD: &str = "gc";
}

fn build_neighbor_distances_cli() -> Command {
    Command::new(consts::ND_CMD)
        .author("Nathan LeRoy")
        .about("Calculate distances between consecutive regions in a query region set")
        .arg(arg!(<path> "Path to bed file OR folder of bed files").required(true))
}

fn build_gc_content_cli() -> Command {
    Command::new(consts::GC_CMD)
        .author("Nathan LeRoy")
        .about("Compute the gc content of a query region set")
        .arg(arg!(<path> "Path to bed file OR folder of bed files").required(true))
        .arg(arg!(-g --genome <GENOME> "genome assembly file").required(true))
}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Nathan LeRoy")
        .about("A command line tool that provides functions for 1) calculating and 2) visualizing a variety of statistics for a collection of genomic ranges.")
        .subcommand_required(true)
        .subcommand(build_neighbor_distances_cli())
        .subcommand(build_gc_content_cli())
}

fn main() -> Result<()> {
    // parse the cli
    let app = build_parser();
    let matches = app.get_matches();

    // build handler for stdout
    let stdout = stdout();
    let mut handle = stdout.lock();

    match matches.subcommand() {
        Some((consts::ND_CMD, matches)) => {
            let path_to_data = matches
                .get_one::<String>("path")
                .expect("Path to data is required.");

            let path_to_data = Path::new(path_to_data);

            ensure!(
                !path_to_data.is_dir(),
                "Please provide a path to a file, not a directory"
            );

            let rs = RegionSet::from_bed(path_to_data)
                .with_context(|| {
                    format!(
                        "Failed to parse bedfile: '{}'",
                        path_to_data.to_string_lossy()
                    )
                })?
                .into_sorted();

            let distances = calc_neighbor_distances(&rs)
                .with_context(|| "Error calculating neighbor distances")?;

            for dist in distances {
                // write to stdout
                handle.write_all(format!("{}\n", dist).as_bytes()).unwrap();
            }

            Ok(())
        }

        Some((consts::GC_CMD, matches)) => {
            // parse cli matches
            let path_to_data = matches
                .get_one::<String>("path")
                .expect("Path to data is required.");
            let genome = matches
                .get_one::<String>("genome")
                .expect("Please specify a genome assembly file");

            // parse given region set
            let region_set = RegionSet::from_bed(Path::new(path_to_data))
                .with_context(|| format!("Failed to parse bedfile: '{}'", path_to_data))?;

            // read in the genome file
            let genome = Path::new(genome);
            let genome = GenomeAssembly::from_fasta(genome).with_context(|| {
                format!("Error reading genome file: '{}'", genome.to_string_lossy())
            })?;

            // compute gc content
            let gc_content = calc_gc_content(&region_set, &genome).unwrap();

            // dump to std-out
            handle.write_all(format!("{}\n", gc_content).as_bytes())?;

            Ok(())
        }
        _ => unreachable!("Subcommand not found"),
    }
}
