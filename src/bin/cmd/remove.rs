
use std::io::{self, ErrorKind};
use crate::data::{DataSink, DataSource};
use super::{Command, datasink, datasource};
use clap::{ArgMatches, Values};
use famlib::{edit_msa::EditMSA, seqs::SequenceAccesors};

pub struct Remove {}

/// Remove Sequence command
impl Remove {
    /// Remove Rows and columns from a MSA
    pub fn remove_command(
            fs: DataSource,
            fo: DataSink,
            rows: Vec<usize>,
            columns: Vec<usize>)
            -> io::Result<()> {
        let mut input = fs
            .get_sequence_collection()
            .unwrap();
        for i in rows {
            input.remove(i);
        }
        if !columns.is_empty() {
            let mut msa = match input.to_msa() {
                Ok(x) => x,
                Err(_) => {
                    return Err(std::io::Error::new(
                        ErrorKind::Other,
                        "Input needs to be an alignment. \
                            All sequences should be the same length.\n",
                    ));
                }
            };
            match msa.remove_columns(columns) {
                Ok(_) => {},
                Err(x) => {
                    return  Err(std::io::Error::new(
                        ErrorKind::Other,
                        format!("{}", x)
                    ))
                }
            }
            input = msa.seq_col_owned();
        }
        fo.write_fasta(&input)
    }
}
impl Command for Remove {
    fn run(&self, matches: &ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("remove") {
            let input = datasource(m);
            let sink = datasink(m);
            let val_to_vec = |x:Values| x.filter_map(|y|
                y.parse::<usize>().ok())
                .map(|x|x-1)
                .collect::<Vec<_>>();
            let rows = m
                .values_of("rows")
                .map_or_else( Vec::new, |x| val_to_vec(x));
            let cols = m
                .values_of("cols")
                .map_or_else(Vec::new, |x| val_to_vec(x));
            Self::remove_command(input, sink, rows, cols)?
        }
        Ok(())
    }
}

