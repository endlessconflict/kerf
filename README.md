# kerf

CRISPR off-target search over haplotypes, in Zig. Reference-only tools miss sites that a person's variants create or destroy, so kerf scans each haplotype as its own sequence.

## How it works

A haplotype is the reference with a set of variants applied (`haplotype`). Variants are SNVs, insertions and deletions, sorted by position, with a map from every haplotype base back to its reference position. Guides of 1 to 62 bases are scanned along the haplotype and along its reverse complement in search mode with unit edit costs, so mismatches and bulges (extra or missing bases) both cost 1. A hit at cost at most K counts as a site when an NGG PAM follows the protospacer.

The scan is the packed kernel of [bitdp](https://github.com/endlessconflict/bitdp), which puts several guides into each 64-bit word with a spacer row between them and derives its carry-chain program at compile time. On one core of a Ryzen 7 8840U it scans 96 guides over chromosome 22 (50.8 Mbp, both strands) in 11 s, about 18 GCUPS.

The test compares every reported site with a plain dynamic program over the same haplotype. The inputs are random genomes, random SNVs, insertions and deletions, guides with planted PAMs, more guides than fit one word, 40 rounds on the reference and on the variant haplotype. The two agree exactly.

## Example on real data

Chromosome 22 (GRCh38) and one 1000 Genomes high-coverage sample (HG00096), one haplotype at a time, with 40 626 and 42 200 biallelic variants applied. Guides are random reference 20-mers followed by NGG, K = 4.

| Haplotype | Guides | Sites on the reference | Sites on the haplotype | Gained | Lost |
|---|---|---|---|---|---|
| 1 | 96 | 262 771 | 262 723 | 974 | 1 027 |
| 2 | 96 | 262 771 | 262 710 | 988 | 1 047 |

A site counts as gained when the haplotype has one and the reference has none for that guide within 8 bases, and as lost in the opposite case. About one site in 260 changes, in each direction. Of the gained sites, 118 and 148 score below K, which means a variant made them strictly better matches. Costs at exactly K also flip in and out of the limit. The random guides here are not a curated panel, so read the counts as a scale, not as a result about any guide.

## Limits

- The PAM is NGG. Other PAMs and the guide-start rules of other nucleases are not there yet.
- Costs are unit edit costs. Position-weighted scoring (such as CFD) is not implemented.
- Sites are found on each haplotype separately. There is no population index, so a cohort costs one scan per distinct haplotype.
- Phased variants must be sorted and spaced. A variant that overlaps the previous one is skipped.
- Ends are reported per hit. Neighbouring end positions of one site appear as separate hits, and the demo merges hits within 3 bases.

## Usage

Requires Zig 0.16.0.

```zig
const kerf = @import("kerf");

var hap = try kerf.haplotype(gpa, reference, variants);
defer hap.deinit(gpa);
var sites: std.ArrayList(kerf.Site) = .empty;
try kerf.scan(gpa, hap, guides, 4, &sites); // guides: []const []const u8
```

```sh
zig build test
zig build demo -Doptimize=ReleaseFast
./zig-out/bin/kerf-demo chr22.fa hap1.vcf 96 4 1   # guides, K, seed
```

## References

See [REFERENCES.md](REFERENCES.md).
