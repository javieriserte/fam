use std::io;

use famlib::filter::Filter as SqFilter;
use crate::data::{DataSink, DataSource};
use super::{datasink, datasource, Command};

pub struct Filter{}

impl Filter {
    pub fn filter(
        input: DataSource,
        output: DataSink,
        ignore_case:bool,
        pattern: &str
    ) -> io::Result<()> {
        let input = input.get_sequence_collection().unwrap();
        let result = match ignore_case {
            true => {
                input.filter_regex_id_ignore_case(pattern)
            },
            false => {
                input.filter_regex_id(pattern)
            }
        };
        output
            .write_fasta(&result)
            .map_err(Into::into)
    }
}

impl Command for Filter {
    fn run (&self, matches: &clap::ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("filter") {
            let input = datasource(m);
            let output = datasink(m);
            let case_insentitive = m.is_present("ignore-case");
            let pattern = m.value_of("pattern").unwrap();
            Self::filter(
                input,
                output,
                case_insentitive,
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
