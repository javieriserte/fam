extern crate clap;
mod cmd;
mod data;
use clap::{App, Arg, SubCommand};
use cmd::{
    Command, collect::Collect, concat::Concat, dimension::Dimension,
    gs::Gapstrip, join::Join, pop::Pop, remove::Remove
};
use std::io;

pub fn main() -> io::Result<()> {
    let matches = App::new("Fasta Alignment Manipulator")
        .version("0.0.2-dev")
        .author("Javier A. Iserte <jiserte@leloir.org.ar>")
        .about("Does many common manipulation of fasta files.")
        .subcommand(SubCommand::with_name("dimensions")
                    .about("Get the dimensions of the fasta file")
                    .arg(Arg::with_name("input")
                        .help("Input file")
                        .takes_value(true))
                    .arg(Arg::with_name("expanded")
                        .short("e")
                        .long("expanded")
                        .help("Show width of all sequences")))
        .subcommand(SubCommand::with_name("collect")
                    .about("Get sequences from stdin and writes to a file.")
                    .arg(Arg::with_name("output")
                        .help("Output file")
                        .takes_value(true)
                        .required(true)))
        .subcommand(SubCommand::with_name("gapstrip")
                    .about("Gapstrip an alignment")
                    .arg(Arg::with_name("input")
                        .short("i")
                        .long("in")
                        .takes_value(true)
                        .help("Input file to be gapstripped"))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("Output file to be gapstripped")))
        .subcommand(SubCommand::with_name("pop")
                    .about("Moves a sequence to the top")
                    .arg(Arg::with_name("input")
                        .short("i")
                        .long("in")
                        .takes_value(true)
                        .help("Input file"))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("Output file"))
                    .arg(Arg::with_name("id")
                        .long("id")
                        .takes_value(true)
                        .required(true)
                        .help("The Id of the sequence to be popped")))
        .subcommand(SubCommand::with_name("join")
                    .about("Joins (merge vertically) two sequence collections into one file")
                    .arg(Arg::with_name("input")
                        .min_values(2)
                        .takes_value(true)
                        .help("A list of comma separated input files"))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("Output file")))
        .subcommand(SubCommand::with_name("concat")
                    .about("Concatenates (merge horizontally) two sequence collections into one file")
                    .arg(Arg::with_name("input")
                        .min_values(2)
                        .takes_value(true)
                        .help("A list of comma separated input files"))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("Output file")))
        .subcommand(SubCommand::with_name("remove")
                    .about("Delete some rows or columns from the sequence collection")
                    .arg(Arg::with_name("input")
                        .short("i")
                        .long("in")
                        .takes_value(true)
                        .help("The input file"))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("The output file"))
                    .arg(Arg::with_name("rows")
                        .short("r")
                        .long("rows")
                        .takes_value(true)
                        .min_values(1)
                        .help("The order index of the row to be removed (first row has index 1)"))
                    .arg(Arg::with_name("cols")
                        .short("c")
                        .long("cols")
                        .takes_value(true)
                        .min_values(1)
                        .help("The id of the rows to be deleted")))
        .get_matches();

    let commands: Vec<Box<dyn Command>>= vec![
        Box::new(Gapstrip{}),
        Box::new(Pop{}),
        Box::new(Dimension{}),
        Box::new(Collect{}),
        Box::new(Join{}),
        Box::new(Concat{}),
        Box::new(Remove{})
    ];
    for cmd in commands {
        cmd.run(&matches)?;
    }
    Ok(())
}
