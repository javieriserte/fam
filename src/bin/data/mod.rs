use famlib::{
    fastaio::{
        buffered_sequence_collection_from_file,
        buffered_sequence_collection_from_stdin,
        sequence_collection_from_file,
        sequence_collection_from_stdin,
        write_buffered_sequence_collection,
        write_sequence_collection
    },
    seqs::{
        BufferedSeqCollection,
        BufferedSeqCollectionFromRead,
        SequenceAccesors,
        SequenceCollection
    }
};
use std::fs::File;
use std::path::Path;
use std::io;
use std::io::stdout;

#[derive(Debug)]
/// Representation of the reading input of a MSA or sequence collection.
pub enum DataSource {
    StdIn,
    FilePath(String),
}

impl DataSource {
    pub fn get_buffered_sequence_collection(&self)
        -> Option<BufferedSeqCollectionFromRead> {
        match self {
            DataSource::StdIn => buffered_sequence_collection_from_stdin().ok(),
            DataSource::FilePath(file) => {
                buffered_sequence_collection_from_file(&Path::new(&file)).ok()
            }
        }
    }
    pub fn get_sequence_collection(&self) -> Option<SequenceCollection> {
        match self {
            DataSource::StdIn => sequence_collection_from_stdin().ok(),
            DataSource::FilePath(file) => sequence_collection_from_file(
                &Path::new(&file)).ok()
        }
    }
    pub fn source_name(&self) -> String {
        match self {
            DataSource::StdIn => String::from("StdIn"),
            DataSource::FilePath(file) => format!(
                "File: {}", file
            ),
        }
    }
    pub fn from(path: &str) -> Self {
        DataSource::FilePath(String::from(path))
    }
}

/// Enum representation of the writing output for a MSA or sequence collection.
/// Possible values are:
/// - StdOut -> writes to the standard output.
/// - FilePath(String) -> writes to file on disk.
pub enum DataSink {
    StdOut,
    FilePath(String),
}

impl DataSink {
    /// Writes a SequenceAccessors to an output.
    pub fn write_fasta<T: SequenceAccesors>(
        &self,
        seqs: &T
    ) -> io::Result<()> {
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
    /// Writes a Buffered Sequence collection to fasta file.
    pub fn write_buffered_to_fasta<T: BufferedSeqCollection>(
        &self,
        seqs: &T
    ) -> io::Result<()> {
        match self {
            DataSink::StdOut => {
                write_buffered_sequence_collection(seqs, stdout().lock())
            }
            DataSink::FilePath(x) => {
                let file = File::create(x).unwrap();
                write_buffered_sequence_collection(seqs, file)
            }
        }
    }
}