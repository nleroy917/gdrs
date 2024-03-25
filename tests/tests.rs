use gdrs::prelude::*;
use rstest::*;

use std::path::Path;

mod tests {
    use super::*;

    #[rstest]
    fn test_calc_neighbor_distances() {
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed")).unwrap().into_sorted();
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
        let region_set = RegionSet::from_bed(Path::new("tests/data/test.bed.gz")).unwrap().into_sorted();
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
}
