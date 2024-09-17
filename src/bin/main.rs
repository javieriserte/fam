extern crate clap;
mod cmd;
mod data;
use clap::{App, Arg, SubCommand};
use cmd::{
    collect::Collect, concat::Concat, degap::Degap, dimension::Dimension, edit::Edit, filter::Filter, gs::Gapstrip, join::Join, onepixel::OnePixel, pad::PadWithGapsCommand, pop::Pop, random::Random, remove::Remove, Command
};
use std::io;


fn add_dimensions_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let new_app = app.subcommand(
        SubCommand::with_name("dimensions")
            .about("Get the dimensions of the fasta file")
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
    );
    return new_app;
}

fn add_collect_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("collect")
            .about("Get sequences from stdin and writes to a file.")
            .arg(
                Arg::with_name("output")
                    .help("Output file")
                    .takes_value(true)
                    .required(true)
            )
    );
    return app
}

fn add_gapstrip_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("gapstrip")
            .about("Gapstrip an alignment")
            .arg(
                Arg::with_name("input")
                    .short("i")
                    .long("in")
                    .takes_value(true)
                    .help("Input file to be gapstripped")
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("Output file to be gapstripped")
            )
    );
    return app;
}

fn add_pop_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
    );
    return app;
}

fn add_join_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("join")
            .about(
                "Joins (merge vertically) two sequence collections into one file"
            )
            .arg(
                Arg::with_name("input")
                    .min_values(2)
                    .takes_value(true)
                    .help("A list of comma separated input files")
            )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("Output file")
            )
    );
    return app;
}

fn add_concat_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("concat")
            .about(
                "Concatenates (merge horizontally) two sequence collections \
                into one file"
            )
            .arg(
                Arg::with_name("input")
                    .min_values(2)
                    .takes_value(true)
                    .help("A list of comma separated input files")
                )
            .arg(
                Arg::with_name("output")
                    .short("o")
                    .long("out")
                    .takes_value(true)
                    .help("Output file")
            )
    );
    return app;
}

fn add_remove_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
    );
    return app;
}

fn add_edit_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("edit")
            .about("Edit MSA content")
    );
    return app;
}

fn add_replace_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
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
    );
    return app;
}

fn add_insert_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
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
    );
    return app;
}

fn add_delete_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
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
    );
    return app;
}

fn add_shuffle_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
                            .short("f")
                            .long("fixed")
                            .takes_value(false)
                            .help("Keep gaps in a fixed position")
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
            )
    );
    return app;
}

fn app_plot_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
    );
    return app;
}

fn add_filter_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
    );
    return app;
}

fn add_degap_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
    let app = app.subcommand(
        SubCommand::with_name("degap")
            .about("Remove gaps from the alignment")
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
    );
    return app;
}

fn add_pad_subcommand<'a>(app: App<'a, 'a>) -> App<'a, 'a> {
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
    );
    return app;
}

fn create_app() -> App<'static, 'static> {
    let mut app = App::new("Fasta Alignment Manipulator")
        .version("0.0.5")
        .author("Javier A. Iserte <jiserte@leloir.org.ar>")
        .about("Does many common manipulation of fasta files.");
    app = add_dimensions_subcommand(app);
    app = add_collect_subcommand(app);
    app = add_gapstrip_subcommand(app);
    app = add_pop_subcommand(app);
    app = add_join_subcommand(app);
    app = add_concat_subcommand(app);
    app = add_remove_subcommand(app);
    app = add_edit_subcommand(app);
    app = add_replace_subcommand(app);
    app = add_insert_subcommand(app);
    app = add_delete_subcommand(app);
    app = add_shuffle_subcommand(app);
    app = app_plot_subcommand(app);
    app = add_filter_subcommand(app);
    app = add_degap_subcommand(app);
    app = add_pad_subcommand(app);
    return app;
}

pub fn main() -> io::Result<()> {
    let app = create_app();
    let matches = app.get_matches();
    let commands: Vec<Box<dyn Command>> = vec![
        Box::new(Dimension{}),
        Box::new(Gapstrip{}),
        Box::new(OnePixel{}),
        Box::new(Collect{}),
        Box::new(Concat{}),
        Box::new(Remove{}),
        Box::new(Random{}),
        Box::new(Join{}),
        Box::new(Edit{}),
        Box::new(Pop{}),
        Box::new(Filter{}),
        Box::new(Degap{}),
        Box::new(PadWithGapsCommand{}),
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
