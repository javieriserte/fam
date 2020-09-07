extern crate clap;
use std::io::stdout;
use std::io::Write;
use std::path::Path;
use clap::{Arg, App, SubCommand};
use std::io;
use std::fs::File;
use std::io::ErrorKind;
use std::io::BufWriter;
use famlib::seqs::{
    SequenceCollection,
    SequenceAccesors};
use famlib::fastaio::{
    sequence_collection_from_stdin,
    sequence_collection_from_file,
    write_sequence_collection};
use famlib::merge::{join, concat};

#[derive(Debug)]
pub struct DataSource{
    pub stdin: bool,
    pub filepath: Option<String>
}


impl DataSource {
    pub fn get_sequence_collection(&self) -> Option<SequenceCollection> {
        match self.stdin{
            true => sequence_collection_from_stdin().ok(),
            false => match &self.filepath {
                Some(x) => sequence_collection_from_file(&Path::new(&x)).ok(),
                None => None
            }
        }
    }
    pub fn source_name(&self) -> String {
        match self.stdin {
            true => String::from("StdIn"),
            false => format!("File: {}", self.filepath.clone().unwrap())
        }
    }
    pub fn from(path :&str) -> Self {
        DataSource{
            stdin:false,
            filepath: Some(String::from(path))
        }
    }
    pub fn stdin() -> Self {
        DataSource{
            stdin:true,
            filepath:None
        }
    }
}

pub enum DataSink{
    StdOut,
    FilePath(String)
}

impl DataSink {
    pub fn write_fasta<T: SequenceAccesors>(&self, seqs: T) -> io::Result<()> {
        match self {
            DataSink::StdOut => {
                write_sequence_collection(seqs, stdout().lock())
            },
            DataSink::FilePath(x) => {
                let file = File::create(x).unwrap();
                write_sequence_collection(seqs, file)
            }
        }
    }
    pub fn write_to_file<T: SequenceAccesors>(&self, seqs: T) -> io::Result<()>{
        match self {
            DataSink::StdOut => panic!("This should be not reachable"),
            DataSink::FilePath(x) => {
                let file = File::create(x).unwrap();
                write_sequence_collection(seqs, file)
            }
        }
    }
}

pub fn pop_command (
    fs: DataSource,
    fo: DataSink,
    id: String
) -> io::Result<()> {
    let mut input = fs.get_sequence_collection().unwrap();
    let seqs = input.move_up(&id);
    match seqs {
        Ok(_) => {fo.write_fasta(input)},
        Err(x) => {
            return Err(
                std::io::Error::new(
                    ErrorKind::Other,
                    format!("Could not pop reference {}.\n", x)
                )
            )
        }
    }
}

pub fn gapstrip_command(fs: DataSource, fo: DataSink) -> io::Result<()> {
    let input = fs.get_sequence_collection().unwrap();
    let msa = match input.to_msa() {
        Ok(x) => x,
        Err(_) => {
            return Err(
                std::io::Error::new(
                    ErrorKind::Other,
                    "Input needs to be an alignment. \
                    All sequences should be the same length.\n"
                )
            );
        }
    };
    let msa = msa.gapstrip();
    fo.write_fasta(msa)
}

pub fn dimension_command(fs: DataSource, expanded:bool) -> io::Result<()> {
    let seqcol = fs.get_sequence_collection().unwrap();
    let out = stdout();
    let mut writer = BufWriter::new(out.lock());
    writer.write_fmt(format_args!("Source: {}", fs.source_name()))?;
    writer.write_fmt(format_args!("Number of sequences: {}", seqcol.size()))?;
    let mut widths:Vec<usize> = seqcol.iter().map(|x| x.len()).collect();
    widths.sort();
    let widths_strings:Vec<String> = widths.into_iter().map(|x| x.to_string()).collect();
    if expanded {
        writer.write_fmt(format_args!("Sequences width: [{:?}]", widths_strings))?;
    } else {
        let max = widths_strings.last().unwrap();
        let min = widths_strings.first().unwrap();
        if min == max {
            writer.write_fmt(format_args!("Sequences width: [{}]", min))?;
        } else {
            writer.write_fmt(format_args!("Sequences width: [{} - {}]", min, max))?;
        }
    }
    Ok(())
}

pub fn collect_command(ds: DataSink) -> io::Result<()> {
    let msa = sequence_collection_from_stdin()?;
    ds.write_to_file(msa)
}

/// Join two or more sequence collections into one.
/// 
pub fn join_command(dss: Vec<DataSource>, sink: DataSink) -> io::Result<()> {
    let seqcols = dss
        .iter()
        .map(|x| 
            x.get_sequence_collection().unwrap())
        .collect::<Vec<_>>();
    match join(seqcols) {
        Ok(x) => sink.write_fasta(x)?,
        Err(x) => eprintln!("{}", x)
    }
    Ok(())
}

//
pub fn concat_command(dss: Vec<DataSource>, sink: DataSink) -> io::Result<()> {
    let seqcols = dss
        .iter()
        .map(|x|
            x.get_sequence_collection().unwrap())
        .collect::<Vec<_>>();
    match concat(seqcols) {
        Ok(x) => sink.write_fasta(x)?,
        Err(x) => eprintln!("{}", x)
    }
    Ok(())
}

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
            .get_matches();

    if let Some(gsmatches) = matches.subcommand_matches("gapstrip") {
        let input = match gsmatches.value_of("input") {
            None => DataSource::stdin(),
            Some(x) => DataSource::from(&x)
        };
        let output = match gsmatches.value_of("output") {
            None => {DataSink::StdOut},
            Some(x) => DataSink::FilePath(String::from(x))
        };
        return gapstrip_command(input, output);
    };
    if let Some(dimatches) = matches.subcommand_matches("dimensions") {
        let input = match dimatches.value_of("input") {
            None => DataSource::stdin(),
            Some(x) => DataSource::from(&x)
        };
        let expanded = dimatches.is_present("expanded");
        return dimension_command(input, expanded);
    };
    if let Some(dimatches) = matches.subcommand_matches("collect") {
        let input = dimatches.value_of("output").unwrap();
        let ds = DataSink::FilePath(input.to_string());
        return collect_command(ds);
    };

    if let Some(m) = matches.subcommand_matches("join") {
        let inputs = m.values_of("input").unwrap();
        let files: Vec<DataSource> = inputs
            .map(|x| DataSource::from(x))
            .collect();
        let sink = match m.value_of("output") {
            None => DataSink::StdOut,
            Some(x) => DataSink::FilePath(x.to_string())
        };
        return join_command(files, sink);
    }

    if let Some(m) = matches.subcommand_matches("concat") {
        let inputs = m.values_of("input").unwrap();
        let files: Vec<DataSource> = inputs
            .map(|x| DataSource::from(x))
            .collect();
        let sink = match m.value_of("output") {
            None => DataSink::StdOut,
            Some(x) => DataSink::FilePath(x.to_string())
        };
        return concat_command(files, sink);
    }
    Ok(())
}
