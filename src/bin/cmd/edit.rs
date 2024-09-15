use std::io::{self, ErrorKind};
use famlib::seqs::SequenceAccesors;
use famlib::edit::EditSequence;
use crate::data::{DataSink, DataSource};
use super::{Command, datasink, datasource};

pub struct Edit{}

impl Edit {
    pub fn edit_replace(
            fs: DataSource,
            fo: DataSink,
            at: Vec<usize>,
            content: Vec<&str>)
            -> io::Result<()> {
        let mut input = fs.get_sequence_collection().unwrap();
        let col_idx = at[1]-1;
        for (x, c) in content.iter().enumerate() {
            let row_idx = at[0]-1+x;
            match input.get_mut(row_idx) {
                Some(seq) => {
                    let new=c.chars().collect::<Vec<_>>();
                    let count = new.len();
                    match seq.edit_replace(new, col_idx, count) {
                        Ok(_) => {}
                        Err(x) => {
                            return  Err(std::io::Error::new(
                                ErrorKind::Other,
                                format!("{}.\n", x))
                            )
                        }
                    }
                }
                None => {
                    return Err(std::io::Error::new(
                        ErrorKind::Other,
                        format!("Row index out of bounds: {}.", row_idx))
                    )
                }
            }
        };
        fo.write_fasta(&input)
    }
    pub fn edit_insert(
            fs: DataSource,
            fo: DataSink,
            at: Vec<usize>,
            content: Vec<&str>)
            -> io::Result<()> {
        let mut input = fs.get_sequence_collection().unwrap();
        let col_idx = at[1]-1;
        for (x, c) in content.iter().enumerate() {
            let row_idx = at[0]-1+x;
            match input.get_mut(row_idx) {
                Some(seq) => {
                    let new=c.chars().collect::<Vec<_>>();
                    match seq.edit_insert(new, col_idx) {
                        Ok(_) => {}
                        Err(x) => {
                            return  Err(std::io::Error::new(
                                ErrorKind::Other,
                                format!("{}.\n", x))
                            )
                        }
                    }
                }
                None => {
                    return Err(std::io::Error::new(
                        ErrorKind::Other,
                        format!("Row index out of bounds: {}.", row_idx))
                    )
                }
            }
        };
        fo.write_fasta(&input)
    }
    pub fn edit_delete(
            fs: DataSource,
            fo: DataSink,
            at: Vec<usize>,
            width: usize,
            height: usize)
            -> io::Result<()> {
        let mut input = fs.get_sequence_collection().unwrap();
        let col_idx = at[1]-1;
        println!("Col Index: {}", col_idx);
        for x in 0..height {
            let row_idx = at[0]-1+x;
            input
                .get_mut(row_idx)
                .ok_or_else(
                    || std::io::Error::new(
                        ErrorKind::Other,
                        format!("Row index out of bounds: {}.", row_idx)
                    )
                )?
                .edit_delete(col_idx, width)?;
        };
        fo.write_fasta(&input)
    }
}

impl Command for Edit {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("edit") {
            if let Some(m1) = m.subcommand_matches("replace") {
                let input = datasource(m1);
                let output = datasink(m1);
                let at = m1.values_of("at")
                    .unwrap()
                    .map(|x| x.parse::<usize>().ok())
                    .collect::<Vec<Option<usize>>>();
                if at.iter().any(|x| x.is_none()) || at.len() != 2 {
                    return Err(
                        std::io::Error::new(
                            ErrorKind::Other,
                            "Can not parse X, Y edit positions.".to_string(),
                        )
                    )
                }
                let at = at
                    .iter()
                    .take(2)
                    .map(|x| x.unwrap())
                    .collect::<Vec<usize>>();
                let content = m1.values_of("content")
                    .unwrap()
                    .collect::<Vec<&str>>();
                Self::edit_replace(input, output, at, content)?
            };
            if let Some(m1) = m.subcommand_matches("insert") {
                let input = datasource(m1);
                let output = datasink(m1);
                let at = m1.values_of("at")
                    .unwrap()
                    .map(|x| x.parse::<usize>().ok())
                    .collect::<Vec<Option<usize>>>();
                if at.iter().any(|x| x.is_none()) || at.len() != 2 {
                    return Err(std::io::Error::new(
                        ErrorKind::Other,
                        "Can not parse X, Y edit positions.".to_string(),
                    ))
                }
                let at = at
                    .iter()
                    .take(2)
                    .map(|x| x.unwrap())
                    .collect::<Vec<usize>>();
                let content = m1.values_of("content")
                    .unwrap()
                    .collect::<Vec<&str>>();
                Self::edit_insert(input, output, at, content)?
            };
            if let Some(m1) = m.subcommand_matches("delete") {
                let input = datasource(m1);
                let output = datasink(m1);
                let at = m1.values_of("at")
                    .unwrap()
                    .map(|x| x.parse::<usize>().ok())
                    .collect::<Vec<Option<usize>>>();
                if at.iter().any(|x| x.is_none()) || at.len() != 2 {
                    return Err(std::io::Error::new(
                        ErrorKind::Other,
                        "Can not parse X, Y edit positions.".to_string(),
                    ))
                }
                let at = at
                    .iter()
                    .take(2)
                    .map(|x| x.unwrap())
                    .collect::<Vec<usize>>();
                let err_gen = |_| Err(
                    std::io::Error::new(
                        ErrorKind::Other,
                        "Can not parse Width or Height.".to_string()
                    )
                );
                let width = m1.value_of("width").unwrap().parse::<usize>().or_else(err_gen)?;
                let height = m1.value_of("height").unwrap().parse::<usize>().or_else(err_gen)?;
                Self::edit_delete(input, output, at, width, height)?
            };
        }
        Ok(())
    }
}

