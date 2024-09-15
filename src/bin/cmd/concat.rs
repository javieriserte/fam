
use std::io::{self};
use crate::data::{DataSink, DataSource};
use super::{Command, datasink};
use clap::{ArgMatches};
use famlib::merge::{concat};

pub struct Concat {}

impl Concat {
    pub fn concat_command(dss: Vec<DataSource>, sink: DataSink) -> io::Result<()> {
        let seqcols = dss
            .iter()
            .map(|x| x.get_sequence_collection().unwrap())
            .collect::<Vec<_>>();
        match concat(seqcols) {
            Ok(x) => sink.write_fasta(&x)?,
            Err(x) => eprintln!("{}", x),
        }
        Ok(())
    }
}
impl Command for Concat {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("concat") {
            let inputs = m.values_of("input").unwrap();
            let files: Vec<DataSource> =
            inputs.map(DataSource::from).collect();
            let sink = datasink(m);
            return Self::concat_command(files, sink);
        }
        Ok(())
    }
}
