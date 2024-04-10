# gdrs
`gdrs` is a rust implementation of the [GenomicDistributions](https://github.com/databio/GenomicDistributions) package.

## Roadmap
Ideally, we would like to implement all the features of the original package. However, we are starting with the basic features and will add more features as we go along. There currently are no plans to implement the plotting features of the original package, just the data generation features:

- [ ] Chromosome distribution
- [x] Neighbor distances
- [x] GC content
- [ ] Partition calculations
- [ ] Cumulative partition calculations
- [ ] Distance to TSS
- [x] Dinucleotide frequency

## Installation
```bash
cargo add gdrs
```
