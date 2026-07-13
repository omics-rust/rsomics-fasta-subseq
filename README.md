# rsomics-fasta-subseq

Extract a subsequence from every record in a FASTA by 1-based region — a
value-exact Rust port of `seqkit subseq`.

## Install

```
cargo install rsomics-fasta-subseq
```

## Usage

```
# bases 1..100 of every record
rsomics-fasta-subseq -r 1:100 genome.fa

# the last 50 bases (negative, from the end)
rsomics-fasta-subseq -r -50:-1 genome.fa

# append the coordinate to each id, write to a file
rsomics-fasta-subseq -r 1:100 -R genome.fa -o subset.fa
```

- `-r, --region` — 1-based region, e.g. `1:100`, `-50:-1`, `5:-5` (required).
- `-R, --region-coord` — append `:start-end` to each output sequence id.
- `-w, --line-width` — output wrap width (default `60`; `0` = no wrap).
- `-o, --output` — output path (`-` = stdout).

## Origin

Independent Rust reimplementation of `seqkit subseq` (region mode), based on
documented behaviour and black-box byte-for-byte comparison against `seqkit`
v2.13: forward and from-end coordinates, clamping past the sequence end, the
`-R` coordinate suffix, and the line-wrap width are matched against the upstream
binary. Frozen goldens live in `tests/golden/`.

License: MIT OR Apache-2.0.
Upstream credit: [seqkit](https://github.com/shenwei356/seqkit) (MIT).
