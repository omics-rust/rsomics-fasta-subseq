//! Compatibility tests: byte-identical output vs seqkit subseq v2.13.0.
//!
//! Golden fixtures generated with:
//!   seqkit subseq -r 1:5    tests/golden/test.fa             → first5.fa
//!   seqkit subseq -r -5:-1  tests/golden/test.fa             → last5.fa
//!   seqkit subseq -r 2:-2   tests/golden/test.fa             → inner.fa
//!   seqkit subseq -r 1:100  tests/golden/test.fa             → clamp.fa
//!   seqkit subseq -r 1:5   -R tests/golden/test.fa           → `first5_coord.fa`
//!   seqkit subseq -r -5:-1 -R tests/golden/test.fa           → `last5_coord.fa`
//!   seqkit subseq -r 1:5   -w 0 tests/golden/test.fa         → nowrap.fa
//!   seqkit subseq -r -21:-1 tests/golden/test.fa             → negclamp.fa
//!
//! `negclamp.fa` covers the "last N bases" idiom where N exceeds a short
//! sequence's length: seqkit clamps the negative start to position 1 and emits
//! the whole sequence (chr2/chr3), while a longer sequence (chr1) yields a true
//! last-21 window. A region with `start < 0 && end > 0` is rejected loudly, so
//! it has no output golden — see `mixed_sign_rejected`.

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

/// `-r -21:-1` on `test.fa`: chr1 (len 28) yields a genuine last-21 window;
/// chr2 (len 20) and chr3 (len 12) are shorter than 21, so the negative start
/// clamps to position 1 and the whole sequence is emitted — matching seqkit.
#[test]
fn neg_clamp_short_seq() {
    assert_eq!(
        run("-21:-1", false, 60),
        expected("negclamp.fa"),
        "neg-clamp mismatch"
    );
}

/// seqkit rejects `start < 0 && end > 0` with a non-zero exit; we mirror that at
/// parse time so no ambiguous window is ever emitted.
#[test]
fn mixed_sign_rejected() {
    for region in ["-1:20", "-5:10", "-20:1"] {
        let err = match Region::parse(region) {
            Err(e) => e,
            Ok(_) => panic!("mixed-sign region {region} must be rejected"),
        };
        assert!(
            err.to_string()
                .contains("when start < 0, end should not > 0"),
            "unexpected error for {region}: {err}"
        );
    }
}

/// Fail-loud at the binary boundary: a mixed-sign region exits non-zero and
/// names the reason on stderr. Under `--json` the framework additionally emits
/// a single machine-readable error envelope (also on stderr; stdout stays
/// empty because no record was produced).
#[test]
fn mixed_sign_binary_fails_loud() {
    let bin = env!("CARGO_BIN_EXE_rsomics-fasta-subseq");
    let input = golden("test.fa");

    let out = std::process::Command::new(bin)
        .args(["-r", "-1:20"])
        .arg(&input)
        .output()
        .expect("spawn binary");
    assert!(!out.status.success(), "expected non-zero exit");
    let stderr = String::from_utf8(out.stderr).expect("utf8 stderr");
    assert!(
        stderr.contains("when start < 0, end should not > 0"),
        "stderr missing reason: {stderr}"
    );

    let json = std::process::Command::new(bin)
        .args(["-r", "-1:20", "--json"])
        .arg(&input)
        .output()
        .expect("spawn binary");
    assert!(!json.status.success(), "expected non-zero exit with --json");
    assert!(
        json.stdout.is_empty(),
        "no record was produced, so stdout must be empty: {:?}",
        String::from_utf8_lossy(&json.stdout)
    );
    let stderr = String::from_utf8(json.stderr).expect("utf8 stderr");
    let envelope = stderr
        .lines()
        .find(|l| l.starts_with('{'))
        .expect("a JSON error envelope on stderr");
    assert!(
        envelope.contains("\"status\":\"error\"")
            && envelope.contains("when start < 0, end should not > 0"),
        "malformed error envelope: {envelope}"
    );
}
