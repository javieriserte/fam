use std::{io::{self, ErrorKind}, path::Path};
use famlib::plotting::OnePixelMsaPlotter;
use crate::data::DataSource;
use super::{Command, datasource};

pub struct OnePixel{}

impl OnePixel {
    pub fn plot_command(
        fs: DataSource,
        outfile: &str,
        is_protein: bool,
        pixel_size: usize
    )
        -> io::Result<()> {
        let input = fs.get_sequence_collection().unwrap();
        match input.to_msa() {
            Ok(msa) => {
                let mut plotter = OnePixelMsaPlotter::new(&msa)
                    .as_dna()
                    .with_pixel_size(pixel_size);
                if is_protein {
                    plotter = plotter.as_protein();
                };
                let outpath = Path::new(outfile);
                match plotter.save_png(outpath) {
                    Err(e) => Err(
                        std::io::Error::new(
                            ErrorKind::Other,
                            format!("There was a problem generating the \
                                    image file {}.\n", e)
                        )
                    ),
                    Ok(_) =>Ok(())
                }
            }
            Err(e) => {
                Err(
                    std::io::Error::new(
                        ErrorKind::Other,
                        format!("Input sequences are not a MSA {}.\n", e)
                    )
                )
            }
        }
    }
}

impl Command for OnePixel {
    fn run(&self, matches: &clap::ArgMatches) ->  io::Result<()> {
        if let Some(m) = matches.subcommand_matches("plot") {
            let input = datasource(m);
            let output = m.value_of("output").unwrap();
            let is_protein = ! m.is_present("is_dna");
            let pixel_size = m.value_of("pixel_size")
                .unwrap()
                .parse::<usize>()
                .unwrap_or(0);
            Self::plot_command(
                input,
                output,
                is_protein,
                pixel_size
            )?
        };
        Ok(())
    }

    fn works_with(&self, matches: &clap::ArgMatches) -> bool {
        matches
            .subcommand_matches("plot")
            .is_some()
    }
}

