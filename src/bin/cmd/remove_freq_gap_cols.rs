use std::io::{self};

use famlib::fastaio::format_from_string;

use crate::data::{DataSink, DataSource};

use super::{datasink, datasource, Command, ToError};

pub struct RemoveFreqGapColumns{}

impl RemoveFreqGapColumns {
    pub fn rm_columns(
        fs: DataSource,
        fo: DataSink,
        min_freq: f64
    ) -> io::Result<()> {
        fs
            .get_sequence_collection()
            .ok_or_else(
                || "Cannot get sequence collection from input.\n"
                    .to_io_error()
            )?
            .to_msa()
            .map_err(|_| "Input is not an alignment.\n".to_io_error())
            .map( |mut msa| {msa.remove_frq_gap_columns(min_freq); msa})
            .map(|msa| fo.write_fasta(&msa))?
    }
}

impl Command for RemoveFreqGapColumns {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("rm-gap-cols-frq") {
            let format = match m.value_of("format") {
                Some(format) => format_from_string(format)?,
                None => {
                    eprintln!("[WARN] No format provided, assuming fasta");
                    format_from_string("fasta")?
                }
            };
            let default_min_freq = 0.5f64;
            let min_freq = match m.value_of("min_freq") {
                Some(freq) => freq.parse::<f64>()
                    .unwrap_or(default_min_freq)
                    .max(0.0f64)
                    .min(1.0f64),
                None => default_min_freq
            };
            let input = datasource(m, format);
            let output = datasink(m);
            Self::rm_columns(input, output, min_freq)?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("rm-gap-cols-frq")
            .is_some()
    }
}