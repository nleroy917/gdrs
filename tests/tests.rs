use gdrs::prelude::*;
use rstest::*;

use std::path::Path;

mod tests {
    use super::*;

    #[rstest]
    fn test_calc_neighnor_distances() {
        let mut region_set = RegionSet::from_bed_file(Path::new("tests/data/test.bed")).unwrap();
        let actual = [
            202_684_741,
            32_120_636,
            3_824_685,
            43_450_165,
            714_607,
            105_011_571
        ];

        let distances = calc_neighnor_distances(&mut region_set).unwrap();
        
        // make sure that the distances are the same, doesnt matter in which order
        for distance in distances {
            assert!(actual.contains(&distance));
        }
        
    }
    
}