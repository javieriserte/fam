use famlib::{
    fastaio::{
        write_sequence_collection,
        sequence_collection_from_file,
        sequence_collection_from_stdin},
    seqs::{SequenceCollection,SequenceAccesors}
};
use std::fs::File;
use std::path::Path;
use std::io;
use std::io::stdout;

#[derive(Debug)]
pub struct DataSource {
    pub stdin: bool,
    pub filepath: Option<String>,
}

impl DataSource {
    pub fn get_sequence_collection(&self) -> Option<SequenceCollection> {
        match self.stdin {
            true => sequence_collection_from_stdin().ok(),
            false => match &self.filepath {
                Some(x) => sequence_collection_from_file(&Path::new(&x)).ok(),
                None => None,
            },
        }
    }
    pub fn source_name(&self) -> String {
        match self.stdin {
            true => String::from("StdIn"),
            false => format!("File: {}", self.filepath.clone().unwrap()),
        }
    }
    pub fn from(path: &str) -> Self {
        DataSource {
            stdin: false,
            filepath: Some(String::from(path)),
        }
    }
    pub fn stdin() -> Self {
        DataSource {
            stdin: true,
            filepath: None,
        }
    }
}

pub enum DataSink {
    StdOut,
    FilePath(String),
}

impl DataSink {
    pub fn write_fasta<T: SequenceAccesors>(&self, seqs: T) -> io::Result<()> {
        match self {
            DataSink::StdOut => {
                write_sequence_collection(seqs, stdout().lock())
            }
            DataSink::FilePath(x) => {
                let file = File::create(x).unwrap();
                write_sequence_collection(seqs, file)
            }
        }
    }
}