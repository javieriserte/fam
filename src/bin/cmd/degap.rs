
use std::io;

use famlib::degap::DegapBufferedSequenceCollection;
use crate::data::{DataSink, DataSource};
use super::{datasink, datasource, Command};

pub struct Degap {}

impl Degap {
    pub fn degap(
        input: DataSource,
        output: DataSink,
    ) -> io::Result<()> {
        let bsq = input
            .get_buffered_sequence_collection()
            .unwrap();
        let result = DegapBufferedSequenceCollection::degap(
            Box::new(bsq)
        );
        output
            .write_buffered_to_fasta(&result)
            .map_err(Into::into)
    }
}

impl Command for Degap {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("degap") {
            let input = datasource(m);
            let output = datasink(m);
            Self::degap(
                input,
                output
            )?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("degap")
            .is_some()
    }
}