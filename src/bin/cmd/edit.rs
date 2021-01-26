use std::io::{self, ErrorKind};
use famlib::seqs::SequenceAccesors;
use famlib::edit::EditSequence;
use crate::data::{DataSink, DataSource};
use super::Command;

pub struct Edit{}

impl Edit {
    pub fn edit_command(
            fs: DataSource,
            fo: DataSink,
            at: Vec<usize>,
            content: Vec<&str>)
            -> io::Result<()> {
        let mut input = fs.get_sequence_collection().unwrap();
        content.iter().enumerate().for_each(
            |(x, c)| {
                let seq = input.get_mut(at[0]-1+x).unwrap();
                let new=c.chars().collect::<Vec<_>>();
                let count = new.len();
                match seq.edit_replace(new, at[1]-1, count) {
                    Ok(_) => {}
                    Err(x) => {println!("{}", x)}
                }
                println!("hola");
            }
        );
        fo.write_fasta(input)
    }
}

impl Command for Edit {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("edit") {
            let input = match m.value_of("input") {
                None => DataSource::stdin(),
                Some(x) => DataSource::from(&x),
            };
            let output = match m.value_of("output") {
                None => DataSink::StdOut,
                Some(x) => DataSink::FilePath(String::from(x)),
            };
            let id = m.value_of("id").unwrap();
            let at = m.values_of("at")
                .unwrap()
                .map(|x| x.parse::<usize>().ok())
                .collect::<Vec<Option<usize>>>();
            if at.iter().any(|x| x.is_none()) || at.len() != 2 {
                return Err(std::io::Error::new(
                    ErrorKind::Other,
                    format!("Can not parse X, Y edit positions.\n",),
                ))
            }
            let at = at
                .iter()
                .take(2)
                .map(|x| x.unwrap())
                .collect::<Vec<usize>>();
            let content = m.values_of("content")
                .unwrap()
                .collect::<Vec<&str>>();
            Self::edit_command(input, output, at, content)?
        };
        Ok(())
    }
}

