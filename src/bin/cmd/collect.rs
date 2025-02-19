
use std::io::{self};
use crate::data::DataSink;
use super::Command;
use clap::ArgMatches;
use famlib::fastaio::{format_from_string, sequence_collection_from_stdin, InputFormats};

pub struct Collect {}
impl Collect {
    pub fn collect_command(ds: DataSink, format: InputFormats) -> io::Result<()> {
        let msa = sequence_collection_from_stdin(format)?;
        ds.write_fasta(&msa)
    }
}
impl Command for Collect {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("collect") {
            let input = m.value_of("output").unwrap();
            let ds = DataSink::FilePath(input.to_string());
            let format = match m.value_of("format") {
                Some(format) => {
                    format_from_string(format)?
                },
                None => {
                    eprintln!("[WARN] No format provided, assuming fasta");
                    format_from_string("fasta")?
                }
            };
            return Collect::collect_command(ds, format);
        };
        Ok(())
    }

    fn works_with(&self, matches: &ArgMatches) -> bool {
        matches
            .subcommand_matches("collect")
            .is_some()
    }
}

