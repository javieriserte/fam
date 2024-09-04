extern crate clap;
mod cmd;
mod data;
use clap::{App, Arg, SubCommand};
use cmd::{
  Command,
  collect::Collect,
  concat::Concat,
  dimension::Dimension,
  edit::Edit,
  gs::Gapstrip,
  join::Join,
  onepixel::OnePixel,
  pop::Pop,
  random::Random,
  remove::Remove,
  filter::Filter
};
use std::io;

pub fn main() -> io::Result<()> {
    let matches = App::new("Fasta Alignment Manipulator")
        .version("0.0.4")
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
        .subcommand(SubCommand::with_name("edit")
            .about("Edit MSA content")
            .subcommand(SubCommand::with_name("replace")
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
                .arg(Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point"))
                .arg(Arg::with_name("content")
                    .short("c")
                    .long("content")
                    .required(true)
                    .takes_value(true)
                    .min_values(1)
                    .help("The replacement content of the edit cells")))
            .subcommand(SubCommand::with_name("insert")
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
                .arg(Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point"))
                .arg(Arg::with_name("content")
                    .short("c")
                    .long("content")
                    .required(true)
                    .takes_value(true)
                    .min_values(1)
                    .help("The content to be inserted.")))
            .subcommand(SubCommand::with_name("delete")
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
                .arg(Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point"))
                .arg(Arg::with_name("width")
                    .short("w")
                    .long("width")
                    .required(true)
                    .takes_value(true)
                    .help("."))
                .arg(Arg::with_name("height")
                    .short("h")
                    .long("height")
                    .required(true)
                    .takes_value(true)
                    .help(".")))
            )
        .subcommand(SubCommand::with_name("shuffle")
            .about("Randomize MSA")
            .subcommand(SubCommand::with_name("all")
                .about("Shuffle content of rows")
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
                .arg(Arg::with_name("fixed")
                    .short("f")
                    .long("fixed")
                    .takes_value(false)
                    .help("Keep gaps in a fixed position")))
            .subcommand(SubCommand::with_name("rows")
                .about("Shuffle row order")
                .arg(Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("The input file"))
                .arg(Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("The output file")))
            .subcommand(SubCommand::with_name("cols")
                .about("Shuffle column order")
                .arg(Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("The input file"))
                .arg(Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("The output file"))))
        .subcommand(SubCommand::with_name("plot")
            .about("Create a simple plot of the MSA")
            .arg(Arg::with_name("input")
                .short("i")
                .long("in")
                .takes_value(true)
                .help("The input file"))
            .arg(Arg::with_name("output")
                .short("o")
                .long("out")
                .required(true)
                .takes_value(true)
                .help("The output PNG file"))
            .arg(Arg::with_name("is_protein")
                .long("is_protein")
                .takes_value(false)
                .help("Use a protein color scheme [default]"))
            .arg(Arg::with_name("is_nucleic")
                .long("is_nucleic")
                .takes_value(false)
                .help("Use a nucleic acid color scheme")
                .conflicts_with("is_protein"))
            .arg(Arg::with_name("pixel_size")
                .long("--pixel_size")
                .takes_value(true)
                .default_value("1")
                .help("Keep gaps in a fixed position")))
        .subcommand(
            SubCommand::with_name("filter")
                .about("filter sequence matching a regex to the sequence id.")
                .arg(
                    Arg::with_name("input")
                        .short("i")
                        .long("in")
                        .takes_value(true)
                        .help("The input file")
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("out")
                        .takes_value(true)
                        .help("The output file")
                )
                .arg(
                    Arg::with_name("ignore_case")
                        .short("c")
                        .long("ignore-case")
                        .takes_value(false)
                        .help("The search is case insensitive")
                )
                .arg(
                    Arg::with_name("pattern")
                        .short("p")
                        .long("pattern")
                        .takes_value(true)
                        .required(true)
                        .help("The search pattern to match")
                    
                )
            
        )    
        .get_matches();

    let commands: Vec<Box<dyn Command>>= vec![
        Box::new(Gapstrip{}),
        Box::new(Pop{}),
        Box::new(Dimension{}),
        Box::new(Collect{}),
        Box::new(Join{}),
        Box::new(Concat{}),
        Box::new(Remove{}),
        Box::new(Edit{}),
        Box::new(Random{}),
        Box::new(OnePixel{}),
        Box::new(Filter{})
    ];
    let is_there_any_command = commands
        .iter()
        .any(|cmd| cmd.works_with(&matches));
    match is_there_any_command {
        true => {
            for cmd in commands {
                match cmd.run(&matches) {
                    Ok(_) => {}
                    Err(x) => {println!("Error: {}", x)}
                }
            }
        },
        false => {
            println!("{}", matches.usage())
        },
    }
    Ok(())
}
