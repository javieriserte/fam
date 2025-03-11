extern crate clap;

#[macro_use]
extern crate rust_i18n;

mod cmd;
mod data;
use clap::{App, Arg, SubCommand};
use cmd::{
    collect::Collect,
    dimension::Dimension,
    edit::Edit,
    filter::Filter,
    gap::Gap,
    onepixel::OnePixel,
    pad::PadWithGapsCommand,
    pop::Pop,
    random::Random,
    combine::Combine,
    remove::Remove,
    trim::Trim,
    Command,
    ToError
};
use std::{collections::HashMap, io};
mod config;
use config::config::get_config;

i18n!("locales");

fn add_dimensions_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let new_app = app.subcommand(
        SubCommand::with_name("dimensions")
            .about(_map.msg("example!"))
            .about("Get the dimensions of the fasta file")
            // This is how the mapping should be used.
            .arg(
                Arg::with_name("input")
                    .help("Input file")
                    .takes_value(true)
            )
            .arg(
                Arg::with_name("expanded")
                    .short("e")
                    .long("expanded")
                    .help("Show width of all sequences")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return new_app;
}

fn add_collect_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("collect")
            .about("Get sequences from stdin and writes to a file.")
            .arg(
                Arg::with_name("output")
                    .help("Output file")
                    .takes_value(true)
                    .required(true)
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app
}

fn add_pop_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("pop")
            .about("Moves a sequence to the top")
            .arg(
                Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("Input file")
                )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("Output file")
            )
            .arg(
                Arg::with_name("id")
                    .long("id")
                    .takes_value(true)
                    .required(true)
                    .help("The Id of the sequence to be popped")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app;
}


fn add_remove_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("remove")
            .about("Delete some rows or columns from the sequence collection")
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
                Arg::with_name("rows")
                    .short("r")
                    .long("rows")
                    .takes_value(true)
                    .min_values(1)
                    .help("The order index of the row to be removed (first row has index 1)")
            )
            .arg(
                Arg::with_name("cols")
                    .short("c")
                    .long("cols")
                    .takes_value(true)
                    .min_values(1)
                    .help("The id of the rows to be deleted")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app;
}

fn add_replace_subcommand<'a>(
    edit: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let edit = edit.subcommand(
        SubCommand::with_name("replace")
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
                Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point")
            )
            .arg(
                Arg::with_name("content")
                    .short("c")
                    .long("content")
                    .required(true)
                    .takes_value(true)
                    .min_values(1)
                    .help("The replacement content of the edit cells")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return edit;
}

fn add_insert_subcommand<'a>(
    edit: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let edit = edit.subcommand(
        SubCommand::with_name("insert")
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
                Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point")
            )
            .arg(
                Arg::with_name("content")
                    .short("c")
                    .long("content")
                    .required(true)
                    .takes_value(true)
                    .min_values(1)
                    .help("The content to be inserted.")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return edit;
}

fn add_delete_subcommand<'a>(
    edit: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let edit = edit.subcommand(
        SubCommand::with_name("delete")
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
                Arg::with_name("at")
                    .long("at")
                    .takes_value(true)
                    .required(true)
                    .min_values(2)
                    .help("X,Y positions of the edit point")
            )
            .arg(
                Arg::with_name("width")
                    .short("w")
                    .long("width")
                    .required(true)
                    .takes_value(true)
                    .help(".")
            )
            .arg(
                Arg::with_name("height")
                    .short("h")
                    .long("height")
                    .required(true)
                    .takes_value(true)
                    .help(".")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return edit;
}

fn add_edit_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &'a Messages
) -> App<'a, 'a> {
    let edit = SubCommand::with_name("edit")
        .about("Edit MSA content");
    let edit = add_replace_subcommand(edit, _map);
    let edit = add_insert_subcommand(edit, _map);
    let edit = add_delete_subcommand(edit, _map);
    let app = app.subcommand(edit);
    return app;
}


fn add_shuffle_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("shuffle")
            .about("Randomize MSA")
            .subcommand(
                SubCommand::with_name("all")
                    .about("Shuffle content of rows")
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
                        Arg::with_name("fixed")
                            .long("fixed")
                            .takes_value(false)
                            .help("Keep gaps in a fixed position")
                    )
                    .arg(
                        Arg::with_name("format")
                            .short("f")
                            .long("format")
                            .help("Specify the input format: [Fasta, Plain]")
                            .default_value("fasta")
                    )
            )
            .subcommand(
                SubCommand::with_name("rows")
                    .about("Shuffle row order")
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
                        Arg::with_name("format")
                            .short("f")
                            .long("format")
                            .help("Specify the input format: [Fasta, Plain]")
                            .default_value("fasta")
                    )
            )
            .subcommand(
                SubCommand::with_name("cols")
                    .about("Shuffle column order")
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
                        Arg::with_name("format")
                            .short("f")
                            .long("format")
                            .help("Specify the input format: [Fasta, Plain]")
                            .default_value("fasta")
                    )
            )
    );
    return app;
}

fn app_plot_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("plot")
            .about("Create a simple plot of the MSA")
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
                    .required(true)
                    .takes_value(true)
                    .help("The output PNG file")
            )
            .arg(
                Arg::with_name("is_protein")
                    .long("is_protein")
                    .takes_value(false)
                    .help("Use a protein color scheme [default]")
            )
            .arg(
                Arg::with_name("is_nucleic")
                    .long("is_nucleic")
                    .takes_value(false)
                    .help("Use a nucleic acid color scheme")
                    .conflicts_with("is_protein")
            )
            .arg(
                Arg::with_name("pixel_size")
                    .long("--pixel_size")
                    .takes_value(true)
                    .default_value("1")
                    .help("Keep gaps in a fixed position")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app;
}

fn add_filter_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
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
            .arg(
                Arg::with_name("exclude")
                    .short("e")
                    .long("exclude")
                    .takes_value(false)
                    .help("Exclude the matching sequences")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app;
}

fn add_pad_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("pad")
            .about("Pad sequence collection with gaps")
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
                Arg::with_name("width")
                    .short("w")
                    .long("width")
                    .takes_value(true)
                    .help("The final width of the sequences")
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
            )
    );
    return app;
}

fn add_combine_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("combine")
            .about("Combine two sequence collections")
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("Output file")
                    .global(true)
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
                    .global(true)
            )
            .arg(
                Arg::with_name("infiles")
                    .min_values(2)
                    .takes_value(true)
                    .help("A list of comma separated input files")
                    .global(true)
            )
            .subcommand(
                SubCommand::with_name("concat")
                    .about("Concatenate the sequences (Vertically)")
            )
            .subcommand(
                SubCommand::with_name("join")
                    .about("Join the sequences (Horizontally)")
            )
            .subcommand(
                SubCommand::with_name("merge")
                    .about("Merge the sequences (Join Horizontally by index)")
                    .arg(
                        Arg::with_name("outer")
                            .long("outer")
                            .takes_value(false)
                            .help(
                                "Include all sequences, even if they are not in all files"
                            )
                    )
            )
    );
    return app;
}

fn add_gap_subcommand<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("gap")
            .about("Gap operations")
            .arg(
                Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("The input file")
                    .global(true)
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("The output file")
                    .global(true)
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
                    .global(true)
            )
            .subcommand(
                SubCommand::with_name("strip")
                    .about("Gapstrip an alignment")
            )
            .subcommand(
                SubCommand::with_name("degap")
                    .about("Remove columns with gaps in all positions")
            )
            .subcommand(
                SubCommand::with_name("remove-columns")
                    .about("Remove columns with gaps")
                    .arg(
                        Arg::with_name("by-freq")
                            .long("by-freq")
                            .takes_value(true)
                            .help("The threshold to remove columns")
                    )
            )
    );
    return app;
}

fn add_trim_command<'a>(
    app: App<'a, 'a>,
    _map: &Messages
) -> App<'a, 'a> {
    app.subcommand(
        SubCommand::with_name("trim")
            .about("Gap operations")
            .arg(
                Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("The input file")
                    .global(true)
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("The output file")
                    .global(true)
            )
            .arg(
                Arg::with_name("format")
                    .short("f")
                    .long("format")
                    .help("Specify the input format: [Fasta, Plain]")
                    .default_value("fasta")
                    .global(true)
            )
            .subcommand(
                SubCommand::with_name("by-gaps")
                    .arg(
                        Arg::with_name("right")
                            .long("right")
                            .takes_value(false)
                            .help(
                                "Trim from the right until no gaps are present"
                            )
                    )
                    .arg(
                        Arg::with_name("left")
                            .long("left")
                            .takes_value(false)
                            .help(
                                "Trim from the left until no gaps are present"
                            )
                    )
            )
            .subcommand(
                SubCommand::with_name("fixed")
                    .arg(
                        Arg::with_name("right")
                            .long("right")
                            .takes_value(true)
                            .help(
                                "Trim from the right a fixed number of columns"
                            )
                    )
                    .arg(
                        Arg::with_name("left")
                            .long("left")
                            .takes_value(true)
                            .help(
                                "Trim from the left a fixed number of columns"
                            )
                    )
            )
            .subcommand(
                SubCommand::with_name("by-terminal-gaps")
                    .arg(
                        Arg::with_name("right")
                            .long("right")
                            .takes_value(true)
                            .help(
                                "Trim from the right by terminal gaps"
                            )
                    )
                    .arg(
                        Arg::with_name("left")
                            .long("left")
                            .takes_value(true)
                            .help(
                                "Trim from the left by terminal gaps"
                            )
                    )
            )

    )
}

fn create_app<'a>(map: &'a Messages)-> App<'a, 'a> {
    let mut app = App::new("Fasta Alignment Manipulator")
        .version("0.0.11")
        .author("Javier A. Iserte <javiserte@gmail.com>")
        .about("Does many common manipulation of fasta files.");
    app = add_dimensions_subcommand(app, &map);
    app = add_collect_subcommand(app, &map);
    app = app_plot_subcommand(app, &map);
    app = add_remove_subcommand(app, &map);
    app = add_shuffle_subcommand(app, &map);
    app = add_edit_subcommand(app, &map);
    app = add_pop_subcommand(app, &map);
    app = add_filter_subcommand(app, &map);
    app = add_pad_subcommand(app, &map);
    app = add_gap_subcommand(app, &map);
    app = add_combine_subcommand(app, &map);
    app = add_trim_command(app, &map);
    return app;
}

fn  create_translation_map<'a>() -> HashMap<&'a str, String> {
    let translation_keys = vec!["a"];
    let map = translation_keys
        .into_iter()
        .map(|x| (x, t!(x).into_owned()))
        .collect::<HashMap<_,_>>();
    map
}

struct Messages<'a, 'b> {
    mapping: &'b HashMap<&'a str, String>
}

impl <'a, 'b> Messages<'a, 'b> {
    pub fn msg(&self, key: &'a str) -> &'b str {
        self
            .mapping
            .get(key)
            .map(|x| x.as_str())
            .unwrap_or(key)
    }
}

pub fn main() -> io::Result<()> {
    let config = get_config();
    rust_i18n::set_locale(&config.locale.lang);
    let map = create_translation_map();
    let messages = Messages {mapping: &map};
    let app = create_app(&messages);
    let mut help_message = Vec::new();
    app
        .write_help(&mut help_message)
        .map_err(|_| t!("Could not write help message").to_io_error())?;
    let matches = app.get_matches();
    let commands: Vec<Box<dyn Command>> = vec![
        Box::new(Dimension{}),
        Box::new(Collect{}),
        Box::new(OnePixel{}),
        Box::new(Remove{}),
        Box::new(Random{}),
        Box::new(Edit{}),
        Box::new(Pop{}),
        Box::new(Filter{}),
        Box::new(PadWithGapsCommand{}),
        Box::new(Gap{}),
        Box::new(Combine{}),
        Box::new(Trim{}),
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
            println!("{}", String::from_utf8_lossy(&help_message));
        },
    }
    Ok(())
}
