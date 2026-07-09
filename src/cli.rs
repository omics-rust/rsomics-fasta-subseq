use std::io::BufWriter;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, HelpSpec, Origin};

use rsomics_fasta_subseq::{Region, SubseqOptions, subseq};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser)]
#[command(
    name = "rsomics-fasta-subseq",
    version,
    about = "Extract FASTA subsequences by 1-based region",
    disable_help_flag = true
)]
pub struct Cli {
    /// Input FASTA file
    input: PathBuf,

    /// Region to extract, e.g. 1:100, -50:-1, 5:-5
    #[arg(short = 'r', long, allow_hyphen_values = true)]
    region: String,

    /// Append :{start}-{end} coordinate to sequence IDs in the output
    #[arg(short = 'R', long = "region-coord")]
    region_coord: bool,

    /// Output FASTA line width (0 = no wrapping)
    #[arg(short = 'w', long = "line-width", default_value_t = 60)]
    line_width: usize,

    /// Output file (default: stdout)
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let region = Region::parse(&self.region)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                std::fs::File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };

        let opts = SubseqOptions {
            region,
            append_coord: self.region_coord,
            line_width: self.line_width,
        };

        subseq(&self.input, &opts, &mut out)?;
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Extract FASTA subsequences by 1-based region.",
    origin: Some(Origin {
        upstream: "seqkit subseq",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &["[OPTIONS] -r <region> <input>"],
    sections: &[],
    examples: &[
        Example {
            description: "First 100 bases of each sequence",
            command: "rsomics-fasta-subseq -r 1:100 genome.fa",
        },
        Example {
            description: "Last 50 bases",
            command: "rsomics-fasta-subseq -r -50:-1 genome.fa",
        },
        Example {
            description: "Bases 101 onwards (trim first 100)",
            command: "rsomics-fasta-subseq -r 101:-1 genome.fa",
        },
        Example {
            description: "Append coordinates to IDs",
            command: "rsomics-fasta-subseq -r 1:100 -R genome.fa",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
