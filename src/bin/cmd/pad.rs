use super::Command;
use famlib::seqs::{BufferedSeqCollection, SequenceCollection};
use famlib::gapping::{PadWithGapsOnce, PadWithGaps};

pub struct PadWithGapsCommand { }

impl Command for PadWithGapsCommand {
  fn run(&self, matches: &clap::ArgMatches) ->  std::io::Result<()> {
    if let Some(m) = matches.subcommand_matches("pad") {
      let input = super::datasource(m);
      let output = super::datasink(m);
      match m.value_of("width") {
        Some(w) => {
          let bsc = Box::new(
            input
              .get_buffered_sequence_collection()
              .unwrap()
          ) as Box<dyn BufferedSeqCollection>;
          let width:usize = w.parse().unwrap();
          let x = bsc.pad_with_gaps(width);
          output
            .write_buffered_to_fasta(&x)
            .map_err(|e| -> std::io::Error { e.into() })?;
        }
        None => {
          let bsc = input
              .get_sequence_collection()
              .unwrap();
          let x: SequenceCollection = bsc.pad_with_gaps_to_max_length();
          output
            .write_fasta(&x)
            .map_err(|e| -> std::io::Error { e.into() })?;
        }
      }
    };
    Ok(())
  }

  fn works_with(&self, matches: &clap::ArgMatches) -> bool {
    matches
      .subcommand_matches("pad")
      .is_some()
  }
}

