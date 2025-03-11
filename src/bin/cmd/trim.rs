use std::io::{Result, Error, ErrorKind::InvalidData};
use famlib::trim::TrimBufferedSequenceCollection;

use crate::data::{DataSink, DataSource};

use super::{datasink, datasource, Command};

pub struct Trim{ }

impl Trim {
    pub fn trim_fixed(
        fs: DataSource,
        fo: DataSink,
        right: usize,
        left: usize,
    ) -> Result<()> {
        let input = fs
            .get_buffered_sequence_collection()
            .unwrap();
        let result = TrimBufferedSequenceCollection::trim_fixed(
            Box::new(input),
            right,
            left
        );
        fo
            .write_buffered_to_fasta(&result)
            .map_err(Into::into)

    }
    pub fn trim_by_gaps(
        fs: DataSource,
        fo: DataSink,
        right: bool,
        left: bool,
    ) -> Result<()> {
        let input = fs
            .get_sequence_collection()
            .unwrap();
        let result = famlib::trim::Trim::trim_by_gaps(&input, right, left);
        fo.write_fasta(&result)
            .map_err(Into::into)
    }
    pub fn trim_by_terminal_gaps(
        fs: DataSource,
        fo: DataSink,
        right: bool,
        left: bool,
    ) -> Result<()> {
        let input = fs
            .get_sequence_collection()
            .unwrap();
        let result = famlib::trim::Trim::trim_by_terminal_gaps(
            &input,
            right,
            left
        );
        fo.write_fasta(&result)
            .map_err(Into::into)
    }
}

impl Command for Trim {
    fn run(&self, matches: &clap::ArgMatches) ->  std::io::Result<()> {
        if let Some(m) = matches.subcommand_matches("trim") {
            if let Some(m1) = m.subcommand_matches("fixed") {
                let input = datasource(m1);
                let sink = datasink(m1);
                let right = m1
                    .value_of("right")
                    .map(|x| x.parse::<usize>())
                    .unwrap_or(Ok(0))
                    .map_err(|e| Error::new(InvalidData, e))?;
                let left = m1
                    .value_of("left")
                    .map(|x| x.parse::<usize>())
                    .unwrap_or(Ok(0))
                    .map_err(|e| Error::new(InvalidData, e))?;
                Self::trim_fixed(input, sink, left, right)?;
            }
            if let Some(m1) = m.subcommand_matches("by-gaps") {
                let input = datasource(m1);
                let sink = datasink(m1);
                let mut right = m1.is_present("right");
                let mut left = m1.is_present("left");
                if !(right || left) {
                    right = true;
                    left = true
                }
                Self::trim_by_gaps(input, sink, left, right)?;
            }
            if let Some(m1) = m.subcommand_matches("by-terminal-gaps") {
                let input = datasource(m1);
                let sink = datasink(m1);
                let mut right = m1.is_present("right");
                let mut left = m1.is_present("left");
                if !(right || left) {
                    right = true;
                    left = true
                }
                Self::trim_by_terminal_gaps(input, sink, left, right)?;
            }
        }
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("trim")
            .is_some()
    }
}