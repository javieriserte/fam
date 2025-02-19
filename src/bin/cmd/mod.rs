use std::io;
use clap::ArgMatches;
use famlib::fastaio::InputFormats;

use crate::data::{DataSink, DataSource};
pub mod gs;
pub mod pop;
pub mod dimension;
pub mod collect;
pub mod join;
pub mod concat;
pub mod remove;
pub mod edit;
pub mod random;
pub mod onepixel;
pub mod filter;
pub mod degap;
pub mod pad;

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
fn datasource(matches: &ArgMatches, format: InputFormats) -> DataSource {
    match matches.value_of("input") {
        None => DataSource::StdIn(format),
        Some(x) => DataSource::from(&x, format),
    }
}
