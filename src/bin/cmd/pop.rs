use std::io::{self, ErrorKind};
use famlib::seqs::SequenceAccesors;
use crate::data::{DataSink, DataSource};
use super::{Command, datasink, datasource};

pub struct Pop{}

impl Pop {
    pub fn pop_command(fs: DataSource, fo: DataSink, id: String) -> io::Result<()> {
        let mut input = fs.get_sequence_collection().unwrap();
        let seqs = input.move_up(&id);
        match seqs {
            Ok(_) => fo.write_fasta(&input),
            Err(x) => {
                return Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("Could not pop reference {}.\n", x),
                ))
            }
        }
    }
}

impl Command for Pop {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("pop") {
            let input = datasource(m);
            let output = datasink(m);
            let id = m.value_of("id").unwrap();
            Self::pop_command(input, output, String::from(id))?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("pop")
            .is_some()
    }
}

