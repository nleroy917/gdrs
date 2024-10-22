use gdrs::prelude::*;
use rstest::*;

use std::path::Path;

mod tests {

    use gdrs::{calc_dinucl_freq, models::GenomeAssembly};

    use super::*;

    #[rstest]
    fn test_calc_neighbor_distances() {
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed"))
            .unwrap()
            .into_sorted();
        let actual = [
            202_684_741,
            32_120_636,
            3_824_685,
            43_450_165,
            714_607,
            105_011_571,
        ];

        let distances = calc_neighbor_distances(&region_set).unwrap();

        // make sure that the distances are the same, doesnt matter in which order
        for distance in distances {
            assert!(actual.contains(&distance));
        }
    }

    #[rstest]
    fn test_calc_neighbor_distances_zipped() {
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed.gz"))
            .unwrap()
            .into_sorted();
        let actual = [
            202_684_741,
            32_120_636,
            3_824_685,
            43_450_165,
            714_607,
            105_011_571,
        ];

        let distances = calc_neighbor_distances(&region_set).unwrap();

        // make sure that the distances are the same, doesnt matter in which order
        for distance in distances {
            assert!(actual.contains(&distance));
        }
    }

    #[rstest]
    fn test_gc_content() {
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed.gz")).unwrap();
        let genome = Path::new("/Users/nathanleroy/genomes/hg38/hg38.fa");
        let genome = GenomeAssembly::from_fasta(genome).unwrap();

        let gc_content = calc_gc_content(&region_set, &genome, true).unwrap();
        println!("{:?}", gc_content);
        assert!(gc_content[0] >= 0.0);
    }

    #[rstest]
    fn test_dinucleotide_freq() {
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed.gz")).unwrap();
        let genome = Path::new("/Users/nathanleroy/genomes/hg38/hg38.fa");
        let genome = GenomeAssembly::from_fasta(genome).unwrap();

        let freqs = calc_dinucl_freq(&region_set, &genome);
        assert!(freqs.is_ok());
        let freqs = freqs.unwrap();
        for (_di, freq) in freqs {
            assert!(freq >= 0.0);
        }
    }
}
