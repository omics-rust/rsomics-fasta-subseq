use std::io::Write;
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

/// 1-based, inclusive region specification (negatives count from end).
#[derive(Clone, Copy)]
pub struct Region {
    pub start: i64,
    pub end: i64,
}

impl Region {
    /// Parse `"start:end"` where each side is a non-zero integer.
    pub fn parse(s: &str) -> Result<Self> {
        let (left, right) = s.split_once(':').ok_or_else(|| {
            RsomicsError::InvalidInput(format!("invalid region `{s}` — expect start:end"))
        })?;

        let start: i64 = left
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid region start `{left}`")))?;
        let end: i64 = right
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("invalid region end `{right}`")))?;

        if start == 0 || end == 0 {
            return Err(RsomicsError::InvalidInput(
                "both start and end should not be 0".to_string(),
            ));
        }

        // A start counted from the end paired with an end counted from the
        // front spans an ambiguous window; seqkit rejects it rather than guess.
        if start < 0 && end > 0 {
            return Err(RsomicsError::InvalidInput(
                "when start < 0, end should not > 0".to_string(),
            ));
        }

        Ok(Self { start, end })
    }

    /// Resolve region against a sequence of `seq_len` bases.
    ///
    /// Returns `(start0, end0_excl)` — 0-based half-open — clamped to
    /// `[0, seq_len]` so slicing `seq[start0..end0_excl]` is always in-bounds.
    /// Returns `(0, 0)` when the region misses the sequence entirely.
    #[must_use]
    pub fn resolve(&self, seq_len: usize) -> (usize, usize) {
        // Sequence lengths in bioinformatics are well within i64 range (max ~3 GB
        // for a human chromosome); the cast cannot wrap in practice.
        #[allow(clippy::cast_possible_wrap)]
        let len = seq_len as i64;

        // Convert 1-based signed to 0-based inclusive. A negative start that
        // counts past the front of a short sequence clamps to position 1
        // (index 0) — seqkit's "last N bases" idiom yields the whole sequence
        // rather than nothing. A positive end past the sequence clamps to its
        // last base.
        let start0 = if self.start > 0 {
            self.start - 1
        } else {
            (len + self.start).max(0)
        };
        let end0_incl = if self.end > 0 {
            self.end.min(len) - 1
        } else {
            len + self.end
        };

        if end0_incl < 0 || start0 > end0_incl {
            return (0, 0);
        }

        // Both values are now non-negative and ≤ seq_len, so the cast to usize
        // is safe on any platform where usize ≥ 32 bits.
        #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
        let s = start0.clamp(0, len) as usize;
        #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
        let e = (end0_incl + 1).clamp(0, len) as usize;
        (s, e)
    }
}

pub struct SubseqOptions {
    pub region: Region,
    /// Append `:{start}-{end}` to the sequence ID in the output header.
    pub append_coord: bool,
    /// 0 = no wrap.
    pub line_width: usize,
}

pub fn subseq(input: &Path, opts: &SubseqOptions, out: &mut dyn Write) -> Result<u64> {
    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let mut count = 0u64;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| RsomicsError::InvalidInput(format!("{e}")))?;

        let full_header = std::str::from_utf8(rec.id())
            .map_err(|e| RsomicsError::InvalidInput(format!("non-UTF-8 header: {e}")))?;

        let (id, desc) = match full_header.find(' ') {
            Some(pos) => (&full_header[..pos], Some(&full_header[pos..])),
            None => (full_header, None),
        };

        let seq = rec.seq();
        let (s, e) = opts.region.resolve(seq.len());
        let window = &seq[s..e];

        write_record(
            id,
            desc,
            window,
            &opts.region,
            opts.append_coord,
            opts.line_width,
            out,
        )?;
        count += 1;
    }
    Ok(count)
}

fn write_record(
    id: &str,
    desc: Option<&str>,
    seq: &[u8],
    region: &Region,
    append_coord: bool,
    line_width: usize,
    out: &mut dyn Write,
) -> Result<()> {
    // Header reconstruction mirrors seqkit's output:
    //   with -R and no description:  ">id:s-e \n"   (trailing space)
    //   with -R and description:     ">id:s-e desc\n"
    //   without -R:                  ">id desc\n"    (original header intact)
    if append_coord {
        match desc {
            Some(d) => {
                write!(out, ">{id}:{}-{}{d}", region.start, region.end)
                    .map_err(RsomicsError::Io)?;
            }
            None => {
                write!(out, ">{id}:{}-{} ", region.start, region.end).map_err(RsomicsError::Io)?;
            }
        }
    } else {
        match desc {
            Some(d) => write!(out, ">{id}{d}").map_err(RsomicsError::Io)?,
            None => write!(out, ">{id}").map_err(RsomicsError::Io)?,
        }
    }
    out.write_all(b"\n").map_err(RsomicsError::Io)?;

    if seq.is_empty() {
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    } else if line_width == 0 || seq.len() <= line_width {
        out.write_all(seq).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    } else {
        for chunk in seq.chunks(line_width) {
            out.write_all(chunk).map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
        }
    }
    Ok(())
}
