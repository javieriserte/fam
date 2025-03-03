use std::io;
use clap::ArgMatches;
use famlib::fastaio::{format_from_string, InputFormats};
use std::result::Result::Err;

use crate::data::{DataSink, DataSource};
pub mod pop;
pub mod dimension;
pub mod collect;
pub mod remove;
pub mod edit;
pub mod random;
pub mod onepixel;
pub mod filter;
pub mod pad;
pub mod gap;
pub mod combine;
pub mod trim;

/// A trait to encapsulate command line execution code.
pub trait Command {
    fn run(&self, matches: &ArgMatches) ->  io::Result<()>;
    fn works_with(&self, matches: &ArgMatches) -> bool;
}

/// Creates a DataSink struct from the commandline arguments
pub fn datasink(matches: &ArgMatches) -> DataSink {
    match matches.value_of("output") {
        None => DataSink::StdOut,
        Some(x) => DataSink::FilePath(String::from(x)),
    }
}

/// Creates a DataSource struct from the commandline arguments
fn datasource(matches: &ArgMatches) -> DataSource {
    let format = inputformat(matches);
    match matches.value_of("input") {
        None => DataSource::StdIn(format),
        Some(x) => DataSource::from(&x, format),
    }
}
// Creates a format input from the commandline arguments
fn inputformat(matches: &ArgMatches) -> InputFormats {
    match matches.value_of("format") {
        Some(format) => {
            match format_from_string(format) {
                Ok(x) => x,
                Err(_) => {
                    eprintln!("[WARN] Invalid format provided, assuming fasta");
                    format_from_string("fasta").unwrap()
                }
            }
        },
        None => {
            eprintln!("[WARN] No format provided, assuming fasta");
            format_from_string("fasta").unwrap()
        }
    }
}

#[allow(dead_code)]
pub trait ToError {
    fn to_error(&self) -> Result<(), io::Error>;
    fn to_error_of_kind(&self, kind: io::ErrorKind) -> Result<(), io::Error>;
    fn to_io_error(&self) -> io::Error;
}

impl ToError for str {
    fn to_error(&self) -> Result<(), io::Error> {
        Err(io::Error::new(io::ErrorKind::Other, self))
    }
    fn to_error_of_kind(&self, kind: io::ErrorKind) -> Result<(), io::Error> {
        Err(io::Error::new(kind, self))
    }
    fn to_io_error(&self) -> io::Error {
        io::Error::new(io::ErrorKind::Other, self)
    }
}

impl ToError for String {
    fn to_error(&self) -> Result<(), io::Error> {
        Err(io::Error::new(io::ErrorKind::Other, self.clone()))
    }
    fn to_error_of_kind(&self, kind: io::ErrorKind) -> Result<(), io::Error> {
        Err(io::Error::new(kind, self.clone()))
    }
    fn to_io_error(&self) -> io::Error {
        io::Error::new(io::ErrorKind::Other, self.clone())
    }
}