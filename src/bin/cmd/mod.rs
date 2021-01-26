use std::io;
use clap::ArgMatches;
pub mod gs;
pub mod pop;
pub mod dimension;
pub mod collect;
pub mod join;
pub mod concat;
pub mod remove;
pub mod edit;

/// A trait to encapsulate command line execution code.
pub trait Command {
    fn run(&self, matches: &ArgMatches) ->  io::Result<()>;
}