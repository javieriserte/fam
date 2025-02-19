use std::io::{self, ErrorKind};
use famlib::{fastaio::format_from_string, random::RandomGen};
use crate::data::{DataSink, DataSource};
use super::{Command, datasink, datasource};

pub struct Random{}

impl Random {
    pub fn shuffle_command(
        fs: DataSource,
        fo: DataSink,
        fixed: bool,
    ) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle(fixed);
                fo.write_fasta(&msa)
            }
            Err(_) => {
                Err(
                    std::io::Error::new(
                        ErrorKind::Other,
                        "A MSA is required.\n".to_string()
                    )
                )
            }
        }
    }
    pub fn shuffle_rows_command(
        fs: DataSource,
        fo: DataSink
    ) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle_rows();
                fo.write_fasta(&msa)
            }
            Err(_) => {
                Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("A MSA is required.\n")))
                }
        }

    }
    pub fn shuffle_cols_command(
        fs: DataSource,
        fo: DataSink,
    ) -> io::Result<()> {
        let msa = fs.get_sequence_collection().unwrap().to_msa();
        match msa {
            Ok(mut msa) => {
                msa.shuffle_cols();
                fo.write_fasta(&msa)
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
            let format = match m.value_of("format") {
                Some(format) => {
                    format_from_string(format)?
                },
                None => {
                    eprintln!("[WARN] No format provided, assuming fasta");
                    format_from_string("fasta")?
                }
            };
            match m.subcommand() {
                ("all", sm) => {
                    let fixed = sm.map_or_else(
                        || false, |x| x.is_present("fixed"));
                    Self::shuffle_command(
                        datasource(m, format),
                        datasink(m),
                        fixed,
                    )?
                },
                ("rows", sm) => {
                    let format = match m.value_of("format") {
                        Some(format) => {
                            format_from_string(format)?
                        },
                        None => {
                            eprintln!("[WARN] No format provided, assuming fasta");
                            format_from_string("fasta")?
                        }
                    };
                    Self::shuffle_rows_command(
                        datasource(sm.unwrap(), format),
                        datasink(sm.unwrap()),
                    )?
                },
                ("cols", sm) => {
                    Self::shuffle_cols_command(
                        datasource(sm.unwrap(), format),
                        datasink(sm.unwrap())
                    )?
                },
                (_, _) => {}
            };
        }
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("shuffle")
            .is_some()
    }
}

