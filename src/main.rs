use glob::glob;
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
}

fn build_neighbor_distances_cli() -> Command {
    Command::new(consts::ND_CMD)
        .author("Nathan LeRoy")
        .about("Calculate distances between consecutive regions in a query region set")
        .arg(arg!(<path> "Path to bed file OR folder of bed files").required(true))
}

fn build_parser() -> Command {
    Command::new(consts::BIN_NAME)
        .bin_name(consts::BIN_NAME)
        .version(consts::VERSION)
        .author("Databio")
        .about("Performance critical tools for working with genomic interval data with an emphasis on preprocessing for machine learning pipelines.")
        .subcommand_required(true)
        .subcommand(build_neighbor_distances_cli())
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
                            let _distances = calc_neighbor_distances(&rs.into_sorted()).unwrap();
                        }
                        Err(e) => {
                            eprintln!("Error reading file: {:?}... skipping", path);
                        }
                    }
                    
                    progress_bar.inc(1);
                }

            // else assume its a bed file
            } else {
                let rs = RegionSet::from_bed(path_to_data).unwrap().into_sorted();
                let _distances = calc_neighbor_distances(&rs).unwrap();
            }
            println!("Done.")
        }
        _ => unreachable!("Subcommand not found"),
    }
}
