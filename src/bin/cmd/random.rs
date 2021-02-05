use std::io::{self, ErrorKind};
use clap::ArgMatches;
use famlib::{random::RandomGen};
use crate::data::{DataSink, DataSource};
use super::Command;

pub struct Random{}

impl Random {
    fn _get_input(matches: &ArgMatches) -> DataSource {
        match matches.value_of("input") {
            None => DataSource::stdin(),
            Some(x) => DataSource::from(&x),
        }
    }
    fn _get_output(matches: &ArgMatches) -> DataSink {
        match matches.value_of("output") {
            None => DataSink::StdOut,
            Some(x) => DataSink::FilePath(String::from(x)),
        }
    }
    pub fn shuffle_command(fs: DataSource, fo: DataSink, fixed: bool) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle(fixed);
                fo.write_fasta(msa)
            }
            Err(_) => {
                Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("A MSA is required.\n")))
                }
            }
    }
    pub fn shuffle_rows_command(fs: DataSource, fo: DataSink) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle_rows();
                fo.write_fasta(msa)
            }
            Err(_) => {
                Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("A MSA is required.\n")))
                }
        }

    }
    pub fn shuffle_cols_command(fs: DataSource, fo: DataSink) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle_rows();
                fo.write_fasta(msa)
            }
            Err(_) => {
                Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("A MSA is required.\n")))
                }
        }
    }
}

impl Command for Random {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("shuffle") {
            match m.subcommand() {
                ("all", sm) => {
                    let fixed = sm.map_or_else(
                        || false, |x| x.is_present("fixed"));
                    Self::shuffle_command(
                        Self::_get_input(m),
                        Self::_get_output(m),
                        fixed)?
                },
                ("rows", _) => {
                    Self::shuffle_rows_command(
                        Self::_get_input(m),
                        Self::_get_output(m))?
                },
                ("cols", _) => {
                    Self::shuffle_cols_command(
                        Self::_get_input(m),
                        Self::_get_output(m))?
                },
                (_, _) => {}
            };
        }
        Ok(())
    }
}

