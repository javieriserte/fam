
use std::io::{self, ErrorKind};
use crate::data::{DataSink, DataSource};
use super::{datasink, datasource, Command, ToError};
use clap::ArgMatches;
use famlib::degap::DegapBufferedSequenceCollection;

pub struct Gap {}

impl Gap {
    pub fn gapstrip_command(
        fs: DataSource,
        fo: DataSink
    ) -> io::Result<()> {
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

    pub fn remove_all_gap_columns_command(
        fs: DataSource,
        fo: DataSink,
    ) -> io::Result<()> {
        fs
            .get_sequence_collection()
            .ok_or_else(
                || "Cannot get sequence collection from input.\n".to_io_error()
            )?
            .to_msa()
            .map_err(|_| "Input is not an alignment.\n".to_io_error())
            .map(|mut msa| {msa.remove_all_gap_columns(); msa})
            .map(|msa| fo.write_fasta(&msa))?
    }

    pub fn rm_columns(
        fs: DataSource,
        fo: DataSink,
        min_freq: f64
    ) -> io::Result<()> {
        fs
            .get_sequence_collection()
            .ok_or_else(
                || "Cannot get sequence collection from input.\n"
                    .to_io_error()
            )?
            .to_msa()
            .map_err(|_| "Input is not an alignment.\n".to_io_error())
            .map( |mut msa| {msa.remove_frq_gap_columns(min_freq); msa})
            .map(|msa| fo.write_fasta(&msa))?
    }
}
impl Command for Gap {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("gap") {
            match m.subcommand() {
                ("strip", Some(m1)) => {
                    let input = datasource(m1);
                    let output = datasink(m1);
                    Self::gapstrip_command(input, output)?
                }
                ("degap", Some(m1)) => {
                    let input = datasource(m1);
                    let output = datasink(m1);
                    let accetps_dots = m1.is_present("accept-dots");
                    Self::degap(input, output, accetps_dots)?
                },
                ("remove-columns", Some(m1)) => {
                    let input = datasource(m1);
                    let output = datasink(m1);
                    let freq = m1
                        .value_of("by-freq")
                        .map(
                            |freq| freq.parse::<f64>()
                                .unwrap_or(0.5f64)
                                .max(0.0f64)
                                .min(1.0f64)
                        );
                    match freq {
                        Some(freq) =>
                            Self::rm_columns(input, output, freq)?,
                        None =>
                            Self::remove_all_gap_columns_command(input, output)?
                    }
                },
                _ => eprintln!("No subcommand provided"),
            }
        };
        Ok(())
    }

    fn works_with(&self, matches: &ArgMatches) -> bool {
        matches
            .subcommand_matches("gap")
            .is_some()
    }
}
