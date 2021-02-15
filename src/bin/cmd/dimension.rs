
use std::io::{self, BufWriter, Write, stdout};
use crate::data::{DataSource};
use super::{Command, datasource};
use famlib::seqs::SequenceAccesors;
use clap::{ArgMatches};

pub struct Dimension {}

impl Dimension {
    pub fn dimension_command(fs: DataSource, expanded: bool) -> io::Result<()> {
        let seqcol = fs.get_sequence_collection().unwrap();
        let out = stdout();
        let mut writer = BufWriter::new(out.lock());
        writer.write_fmt(format_args!("Source: {}", fs.source_name()))?;
        writer.write_fmt(format_args!("Number of sequences: {}", seqcol.size()))?;
        let mut widths: Vec<usize> = seqcol.iter().map(|x| x.len()).collect();
        widths.sort();
        let widths_strings: Vec<String> =
            widths.into_iter().map(|x| x.to_string()).collect();
        if expanded {
            writer.write_fmt(format_args!(
                "Sequences width: [{:?}]",
                widths_strings
            ))?;
        } else {
            let max = widths_strings.last().unwrap();
            let min = widths_strings.first().unwrap();
            if min == max {
                writer.write_fmt(format_args!("Sequences width: [{}]", min))?;
            } else {
                writer.write_fmt(format_args!(
                    "Sequences width: [{} - {}]",
                    min, max
                ))?;
            }
        }
        Ok(())
    }}
impl Command for Dimension {
    fn run(&self, matches: &ArgMatches) -> io::Result<()> {
        if let Some(dimatches) = matches.subcommand_matches("dimensions") {
            let input = datasource(dimatches);
            let expanded = dimatches.is_present("expanded");
            return Dimension::dimension_command(input, expanded);
        };
    Ok(())
    }
}