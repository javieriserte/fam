
use std::io::{self};
use crate::data::{DataSink};
use super::Command;
use clap::{ArgMatches};
use famlib::fastaio::sequence_collection_from_stdin;

pub struct Collect {}
impl Collect {
    pub fn collect_command(ds: DataSink) -> io::Result<()> {
        let msa = sequence_collection_from_stdin()?;
        ds.write_fasta(&msa)
    }
}
impl Command for Collect {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("collect") {
            let input = m.value_of("output").unwrap();
            let ds = DataSink::FilePath(input.to_string());
            return Collect::collect_command(ds);
        };
        Ok(())
    }
}

