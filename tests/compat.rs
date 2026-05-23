//! Compatibility tests: byte-identical output vs seqkit subseq v2.13.0.
//!
//! Golden fixtures generated with:
//!   seqkit subseq -r 1:5   tests/golden/test.fa              → first5.fa
//!   seqkit subseq -r -5:-1 tests/golden/test.fa              → last5.fa
//!   seqkit subseq -r 2:-2  tests/golden/test.fa              → inner.fa
//!   seqkit subseq -r 1:100 tests/golden/test.fa              → clamp.fa
//!   seqkit subseq -r 1:5   -R tests/golden/test.fa           → `first5_coord.fa`
//!   seqkit subseq -r -5:-1 -R tests/golden/test.fa           → `last5_coord.fa`
//!   seqkit subseq -r 1:5   -w 0 tests/golden/test.fa         → nowrap.fa

use std::path::Path;

use rsomics_fasta_subseq::{Region, SubseqOptions, subseq};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn run(region: &str, append_coord: bool, line_width: usize) -> String {
    let region = Region::parse(region).expect("valid region");
    let opts = SubseqOptions {
        region,
        append_coord,
        line_width,
    };
    let mut buf = Vec::new();
    subseq(&golden("test.fa"), &opts, &mut buf).expect("subseq failed");
    String::from_utf8(buf).expect("valid utf8")
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(golden(name)).expect("golden file missing")
}

#[test]
fn first5() {
    assert_eq!(
        run("1:5", false, 60),
        expected("first5.fa"),
        "first5 mismatch"
    );
}

#[test]
fn last5() {
    assert_eq!(
        run("-5:-1", false, 60),
        expected("last5.fa"),
        "last5 mismatch"
    );
}

#[test]
fn inner() {
    assert_eq!(
        run("2:-2", false, 60),
        expected("inner.fa"),
        "inner mismatch"
    );
}

#[test]
fn clamp_beyond_end() {
    assert_eq!(
        run("1:100", false, 60),
        expected("clamp.fa"),
        "clamp mismatch"
    );
}

#[test]
fn first5_with_coord() {
    assert_eq!(
        run("1:5", true, 60),
        expected("first5_coord.fa"),
        "first5 coord mismatch"
    );
}

#[test]
fn last5_with_coord() {
    assert_eq!(
        run("-5:-1", true, 60),
        expected("last5_coord.fa"),
        "last5 coord mismatch"
    );
}

#[test]
fn no_line_wrap() {
    assert_eq!(
        run("1:5", false, 0),
        expected("nowrap.fa"),
        "nowrap mismatch"
    );
}
