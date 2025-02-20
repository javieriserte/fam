use std::io::{self};

use famlib::fastaio::format_from_string;

use crate::data::{DataSink, DataSource};

use super::{datasink, datasource, Command, ToError};

pub struct RemoveAllGapColumns{}

impl RemoveAllGapColumns {
    pub fn remove_all_gap_columns_command(
        fs: DataSource,
        fo: DataSink,
    ) -> io::Result<()> {
        fs
            .get_sequence_collection()
            .ok_or_else(
                || "Cannot get sequence collection from input.\n".to_io_error()
            )?
            .to_msa()
            .map_err(|_| "Input is not an alignment.\n".to_io_error())
            .map(|mut msa| {msa.remove_all_gap_columns(); msa})
            .map(|msa| fo.write_fasta(&msa))?
    }
}

impl Command for RemoveAllGapColumns {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("rm-gap-cols") {
            let format = match m.value_of("format") {
                Some(format) => {
                    format_from_string(format)?
                },
                None => {
                    eprintln!("[WARN] No format provided, assuming fasta");
                    format_from_string("fasta")?
                }
            };
            let input = datasource(m, format);
            let output = datasink(m);
            Self::remove_all_gap_columns_command(input, output)?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("rm-gap-cols")
            .is_some()
    }
}