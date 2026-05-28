use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_fasta_subseq(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fasta-subseq");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fa = manifest.join("tests/golden/clamp.fa");
    c.bench_function("rsomics-fasta-subseq golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([fa.to_str().unwrap(), "-r", "1:5"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fasta_subseq);
criterion_main!(benches);
