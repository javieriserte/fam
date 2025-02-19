
use std::io;

use famlib::{degap::DegapBufferedSequenceCollection, fastaio::format_from_string};
use crate::data::{DataSink, DataSource};
use super::{datasink, datasource, Command};

pub struct Degap {}

impl Degap {
    pub fn degap(
        input: DataSource,
        output: DataSink,
        accept_dots: bool
    ) -> io::Result<()> {
        let bsq = input
            .get_buffered_sequence_collection()
            .unwrap();
        let result = DegapBufferedSequenceCollection::degap(
            Box::new(bsq),
            accept_dots
        );
        output
            .write_buffered_to_fasta(&result)
            .map_err(Into::into)
    }
}

impl Command for Degap {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("degap") {
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
            let accetps_dots = m.is_present("accept-dots");
            Self::degap(
                input,
                output,
                accetps_dots
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