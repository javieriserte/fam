
use std::io::{self};
use crate::data::{DataSink, DataSource};
use super::{Command, datasink};
use clap::ArgMatches;
use famlib::merge::join;
use famlib::fastaio::format_from_string;
pub struct Join {}

impl Join {
    pub fn join_command(
        dss: Vec<DataSource>,
        sink: DataSink
    ) -> io::Result<()> {
        let seqcols = dss
            .iter()
            .map(|x| x.get_sequence_collection().unwrap())
            .collect::<Vec<_>>();
        match join(seqcols) {
            Ok(x) => sink.write_fasta(&x)?,
            Err(x) => eprintln!("{}", x),
        }
        Ok(())
    }
}
impl Command for Join {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(m) = matches.subcommand_matches("join") {
            let format = match m.value_of("format") {
                Some(format) => {
                    format_from_string(format)?
                },
                None => {
                    eprintln!("[WARN] No format provided, assuming fasta");
                    format_from_string("fasta")?
                }
            };
            let inputs = m.values_of("input").unwrap();
            let files: Vec<DataSource> = inputs
                .map(|x| DataSource::from(x, format))
                .collect();
            let sink = datasink(m);
            return Join::join_command(files, sink);
        }
        Ok(())
    }

    fn works_with(&self, matches: &ArgMatches) -> bool {
        matches
            .subcommand_matches("join")
            .is_some()
    }
}
