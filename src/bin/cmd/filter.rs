use std::io;

use famlib::{fastaio::format_from_string, filter::FilterBufferedSequenceCollection};
use crate::data::{DataSink, DataSource};
use super::{datasink, datasource, Command};

pub struct Filter{}

impl Filter {
    pub fn filter(
        input: DataSource,
        output: DataSink,
        ignore_case:bool,
        keep: bool,
        pattern: &str
    ) -> io::Result<()> {
        let bsq = input
            .get_buffered_sequence_collection()
            .unwrap();
        let result = FilterBufferedSequenceCollection::filter_regex_id(
            Box::new(bsq),
            ignore_case,
            keep,
            pattern
        );
        output
            .write_buffered_to_fasta(&result)
            .map_err(Into::into)
    }
}

impl Command for Filter {
    fn run (&self, matches: &clap::ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("filter") {
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
            let case_insentitive = m.is_present("ignore_case");
            let keep = !m.is_present("exclude");
            let pattern = m.value_of("pattern").unwrap();
            Self::filter(
                input,
                output,
                case_insentitive,
                keep,
                pattern
            )?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("filter")
            .is_some()
    }
}
