use std::io::{self};
use crate::data::{DataSink, DataSource};
use super::{datasink, inputformat, Command, ToError};
use clap::ArgMatches;
use famlib::{
    combine::concat,
    seqs::SequenceAccesors,
    combine::join,
    combine::merge
};

pub struct Combine{}

impl Combine {
    pub fn concat_command(
        dss: Vec<DataSource>,
        sink: DataSink
    ) -> io::Result<()> {
        let seqcols = dss.iter().fold(
            Result::Ok(vec![]),
            |acc, x| {
                match x.get_sequence_collection() {
                    Some(y) => acc.map(|mut z| { z.push(y); z }),
                    None => Err(
                        "All inputs must be sequence collections".to_io_error()
                    )
                }
            }
        );
        seqcols
            .and_then( |x| concat(x).map_err(|x| x.into()))
            .and_then( |x| sink.write_fasta(&x))
    }
    pub fn join_command(
        dss:Vec<DataSource>,
        sink: DataSink
    ) -> io::Result<()> {
        dss
        .iter()
        .fold(
            Result::Ok(vec![]),
            |acc, x| {
                match x.get_sequence_collection() {
                    Some(y) => acc.map(|mut z| { z.push(y); z }),
                    None => Err(
                        "All inputs must be sequence collections".to_io_error()
                    )
                }
            }
        )
        .and_then(
            |x| {
                match x.iter().all(|y| y.size() == x[0].size()) {
                    true => Ok(x),
                    false => Err(
                        "All sequence collections should have the same size."
                            .to_io_error()
                    )
                }
            }
        )
        .and_then(|x| join(x).map_err(|x| x.into()))
        .and_then(|x| sink.write_fasta(&x))
    }
    pub fn merge_command(
        dss: Vec<DataSource>,
        sink: DataSink,
        outer: bool
    ) -> io::Result<()> {
        dss
            .iter()
            .fold(
                Result::Ok(vec![]),
                |acc, x| {
                    match x.get_sequence_collection() {
                        Some(y) => acc.map(|mut z| { z.push(y); z }),
                        None => Err(
                            "All inputs must be sequence collections"
                                .to_io_error()
                        )
                    }
                }
            )
            .and_then(|x| merge(x, outer).map_err(|x| x.into()))
            .and_then(|x| sink.write_fasta(&x))
    }
}

impl Command for Combine {
    fn run(&self, matches: &ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("combine") {
            match m.subcommand() {
                ("concat", Some(m)) => {
                    let format = inputformat(&m);
                    let inputs = m.values_of("infiles").unwrap();
                    let files: Vec<DataSource> = inputs
                        .map(
                            |x| DataSource::from(x, format)
                        )
                        .collect();
                    let sink = datasink(m);
                    return Self::concat_command(files, sink);
                },
                ("join", Some(m)) => {
                    let format = inputformat(&m);
                    let inputs = m.values_of("infiles").unwrap();
                    let files: Vec<DataSource> = inputs
                        .map(
                            |x| DataSource::from(x, format)
                        )
                        .collect();
                    let sink = datasink(m);
                    return Self::join_command(files, sink);
                },
                ("merge", Some(m)) => {
                    let format = inputformat(&m);
                    let inputs = m.values_of("infiles").unwrap();
                    let files: Vec<DataSource> = inputs
                        .map(
                            |x| DataSource::from(x, format)
                        )
                        .collect();
                    let outer = m.is_present("outer");
                    let sink = datasink(m);
                    return Self::merge_command(files, sink, outer);
                },
                _ => {
                    eprintln!("Invalid subcommand");
                }
            }
        }
        Ok(())
    }

    fn works_with(&self, matches: &ArgMatches) -> bool {
        matches
            .subcommand_matches("combine")
            .is_some()
    }
}