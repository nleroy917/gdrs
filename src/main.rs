use gdrs::calc_gc_content;
use glob::glob;
use std::io::stdout;
use std::io::Write;
use std::path::Path;

use clap::{arg, Command};
use indicatif::{ProgressBar, ProgressStyle};

use gdrs::calc_neighbor_distances;
use gdrs::prelude::*;

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

fn main() {
    let app = build_parser();
    let matches = app.get_matches();

    match matches.subcommand() {
        Some((consts::ND_CMD, matches)) => {
            let path_to_data = matches
                .get_one::<String>("path")
                .expect("Path to data is required.");

            let path_to_data = Path::new(path_to_data);

            // if dir, assume folder of bed files
            if path_to_data.is_dir() {
                let n_files = glob(&format!("{}/*.bed*", path_to_data.to_str().unwrap()))
                    .expect("Failed to read glob pattern")
                    .count() as u64;

                let files = glob(&format!("{}/*.bed*", path_to_data.to_str().unwrap()))
                    .expect("Failed to read glob pattern");

                let progress_bar = ProgressBar::new(n_files);
                progress_bar.set_style(
                    ProgressStyle::with_template(
                        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
                    )
                    .unwrap()
                    .progress_chars("##-"),
                );

                for entry in files {
                    let entry = entry.unwrap();

                    let path = entry.as_path();
                    let rs = RegionSet::from_bed(path);

                    match rs {
                        Ok(rs) => {
                            let distances = calc_neighbor_distances(&rs.into_sorted()).unwrap();

                            let stdout = stdout();
                            let mut handle = stdout.lock();

                            for dist in distances {
                                // write to stdout
                                handle.write_all(format!("{}\n", dist).as_bytes()).unwrap();
                            }
                        }
                        Err(_e) => {
                            eprintln!("Error reading file: {:?}... skipping", path);
                        }
                    }

                    progress_bar.inc(1);
                }

            // else assume its a bed file
            } else {
                let rs = RegionSet::from_bed(path_to_data).unwrap().into_sorted();

                let stdout = stdout();
                let mut handle = stdout.lock();

                let distances = calc_neighbor_distances(&rs).unwrap();

                for dist in distances {
                    // write to stdout
                    handle.write_all(format!("{}\n", dist).as_bytes()).unwrap();
                }
            }
        }

        Some((consts::GC_CMD, matches)) => {
            let path_to_data = matches
                .get_one::<String>("path")
                .expect("Path to data is required.");
            let genome = matches
                .get_one::<String>("genome")
                .expect("Please specify a genome assembly file");

            let region_set = RegionSet::from_bed(Path::new(path_to_data)).unwrap();
            let genome = Path::new(genome);

            let gc_content = calc_gc_content(&region_set, genome).unwrap();

            println!("{}", gc_content)
        }
        _ => unreachable!("Subcommand not found"),
    }
}
