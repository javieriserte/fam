
use std::io::{self, ErrorKind};
use crate::data::{DataSink, DataSource};
use super::Command;
use clap::{ArgMatches};

pub struct Gapstrip {}

impl Gapstrip {
    pub fn gapstrip_command(fs: DataSource, fo: DataSink) -> io::Result<()> {
        let input = fs.get_sequence_collection().unwrap();
        let msa = match input.to_msa() {
            Ok(x) => x,
            Err(_) => {
                return Err(std::io::Error::new(
                    ErrorKind::Other,
                    "Input needs to be an alignment. \
                        All sequences should be the same length.\n",
                ));
            }
        };
        let msa = msa.gapstrip();
        fo.write_fasta(&msa)
    }
}
impl Command for Gapstrip {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(gsmatches) = matches.subcommand_matches("gapstrip") {
            let input = match gsmatches.value_of("input") {
                None => DataSource::stdin(),
                Some(x) => DataSource::from(&x),
            };
            let output = match gsmatches.value_of("output") {
                None => DataSink::StdOut,
                Some(x) => DataSink::FilePath(String::from(x)),
            };
            Self::gapstrip_command(input, output)?
        };
        Ok(())
    }
}