extern crate graphics;
extern crate tempfile;

use std::{io::Result, path::Path, time};
use graphics_buffer::*;
use graphics::rectangle;
use time::SystemTime;

use crate::seqs::{Alignment, SequenceAccesors};

trait ColorScheme{
    fn color(&self, char: &char) -> [f32; 4];
}

struct ProteinColors {colors: [[f32; 4]; 14]}

impl ProteinColors {
    pub fn new() -> Self {
        ProteinColors{
            colors: [
                [0.901, 0.039, 0.039, 1.000],  // ASP, GLU
                [0.901, 0.901, 0.000, 1.000],  // CYS, MET
                [0.078, 0.352, 1.000, 1.000],  // LYS, ARG
                [0.980, 0.588, 0.000, 1.000],  // SER, THR
                [0.196, 0.196, 0.666, 1.000],  // PHE, TYR
                [0.000, 0.862, 0.862, 1.000],  // ASN, GLN
                [0.921, 0.921, 0.921, 1.000],  // GLY
                [0.058, 0.509, 0.058, 1.000],  // LEU, VAL, IL
                [0.784, 0.784, 0.784, 1.000],  // ALA
                [0.705, 0.352, 0.705, 1.000],  // TRP
                [0.509, 0.509, 0.823, 1.000],  // HIS
                [0.862, 0.588, 0.509, 1.000],  // PRO
                [0.500, 0.500, 0.500, 1.000],  // GAP
                [0.000, 0.000, 0.000, 1.000],  // UNKNOWN
                ]
            }
        }
    }

    impl ColorScheme for ProteinColors {
        fn color(&self, char: &char) -> [f32; 4] {
            match char {
                'D'| 'E' => self.colors[0],
                'C'| 'M' => self.colors[1],
                'K'| 'R' => self.colors[2],
                'S'| 'T' => self.colors[3],
                'F'| 'Y' => self.colors[4],
                'N'| 'Q' => self.colors[5],
                'G' => self.colors[6],
                'L'| 'V'| 'I' => self.colors[7],
                'A' => self.colors[8],
                'W' => self.colors[9],
                'H' => self.colors[10],
                'P' => self.colors[11],
                '-' => self.colors[12],
                _ => self.colors[13]
            }
        }
    }
    struct NucleicAcidColors {colors: [[f32; 4]; 6]}

    impl NucleicAcidColors {
        pub fn new() -> Self {
            NucleicAcidColors{
                colors: [
                    [0.901, 0.039, 0.039, 1.000],  // A
                    [0.901, 0.901, 0.000, 1.000],  // C
                    [0.078, 0.352, 1.000, 1.000],  // T
                    [0.980, 0.588, 0.000, 1.000],  // G
                    [0.500, 0.500, 0.500, 1.000],  // GAP
                    [0.000, 0.000, 0.000, 1.000],  // UNKNOWN
                ]
            }
        }
    }

    impl ColorScheme for NucleicAcidColors {
        fn color(&self, char: &char) -> [f32; 4] {
            match char {
                'A' => self.colors[0],
                'C' => self.colors[1],
                'T'| 'U' => self.colors[2],
                'G' => self.colors[3],
                '-' => self.colors[4],
                _ => self.colors[5]
            }
        }
    }
    pub struct OnePixelMsaPlotter<'a> {
        msa: &'a Alignment,
        is_protein: bool,
        pixel_size: usize
    }

    impl <'a> OnePixelMsaPlotter<'a> {
        pub fn new(msa: &'a Alignment) -> Self {
            OnePixelMsaPlotter{
                msa,
                is_protein: true,
                pixel_size: 1usize
            }
        }
        pub fn as_dna(mut self) -> Self {
            self.is_protein=false;
            self
        }
        pub fn as_protein(mut self) -> Self {
            self.is_protein=true;
            self
        }
        pub fn with_pixel_size(mut self, pixel_size: usize) -> Self {
            if pixel_size<50 {
                self.pixel_size=pixel_size;
            };
            self
        }
        pub fn save_png(&self, outfile: &Path) -> Result<()> {
            let margin:u32 = 20;
            let width:u32 = self.msa.length() as u32 * self.pixel_size as u32 + 2 * margin;
            let height:u32 = self.msa.size() as u32 * self.pixel_size as u32 + 2 * margin;
            let mut buffer = RenderBuffer::new(width, height);
            buffer.clear([1.0, 1.0, 1.0, 1.0]);
            let color_scheme: Box<dyn ColorScheme> = match self.is_protein {
                true => Box::new(ProteinColors::new()),
                false => Box::new(NucleicAcidColors::new())
            };
            let now = SystemTime::now();
            let mut last = 0u128;
            self.msa.iter().enumerate().for_each(
                |(i, x)| {
                    if i % 100 == 0 {
                        let current = now.elapsed().unwrap().as_millis();
                        println!(
                            "Row:{} = {:?} secs | Delta = {}",
                            i,
                            current,
                            current-last
                        );
                        last = current;
                    }
                    x.seq().unwrap().iter().enumerate().for_each(
                    |(j, c)| {
                        let color = color_scheme.color(c);
                        let x = margin as f64 + j as f64*self.pixel_size as f64;
                        let y = margin as f64 + i as f64*self.pixel_size as f64;
                        rectangle(
                            color,
                            [x, y, self.pixel_size as f64, self.pixel_size as f64],
                            IDENTITY,
                            &mut buffer,
                        );
                    }
                )}
            );
            // Save the buffer
            let x = outfile.display().to_string();
            buffer.save(&x[..]).unwrap();
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::Read};
    use tempfile::tempdir;
    use crate::seqs::{AnnotatedSequence, SequenceAccesors, SequenceCollection};
    use super::*;

    fn sample_display_msa() -> Alignment {
        let mut sq = SequenceCollection::new();
        let chars = "ATGAAAGTGCATTAA";
        for i in 0..10 {
            let ca = AnnotatedSequence::from_string(
                format!("Seq_{}", i), String::from(&chars[i..i+4]));
                sq.add(ca).unwrap();
            }
            sq.to_msa().ok().unwrap()
        }

    #[test]
    fn test_plotter_creation() {
        let msa = sample_display_msa();
        let p1 = OnePixelMsaPlotter::new(&msa);
        assert_eq!(p1.msa.size(), msa.size());
        assert_eq!(p1.msa.length(), msa.length());
        assert!(p1.is_protein);
        assert_eq!(p1.pixel_size, 1usize);
        let p1 = OnePixelMsaPlotter::new(&msa).as_dna();
        assert!(!p1.is_protein);
        let p1 = OnePixelMsaPlotter::new(&msa).as_dna().as_protein();
        assert!(p1.is_protein);
        let p1 = OnePixelMsaPlotter::new(&msa).as_dna().with_pixel_size(3);
        assert!(!p1.is_protein);
        assert!(p1.pixel_size==3);
    }
    #[test]
    fn test_plotter_images() {
        let msa = sample_display_msa();
        let mut plotter = OnePixelMsaPlotter::new(&msa).with_pixel_size(3);
        let tdir = tempdir().unwrap().into_path();
        let cpath = tdir.join("temp.png");
        plotter.save_png(cpath.as_path()).unwrap();
        let mut open_file = File::open(cpath).unwrap();
        let expected: Vec<u8> = vec![
            137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0,
            0, 52, 0, 0, 0, 70, 8, 6, 0, 0, 0, 131, 146, 183, 33, 0, 0, 1, 2,
            73, 68, 65, 84, 120, 156, 237, 216, 33, 178, 194, 64, 16, 69, 209,
            143, 249, 22, 21, 4, 54, 59, 64, 133, 133, 176, 9, 22, 196, 38, 88,
            8, 81, 236, 32, 6, 129, 32, 10, 11, 2, 186, 170, 5, 109, 166, 218,
            140, 120, 12, 247, 138, 39, 48, 169, 35, 134, 154, 100, 241, 178,
            254, 26, 10, 144, 122, 128, 212, 3, 164, 30, 32, 245, 0, 169, 7, 72,
            61, 64, 234, 1, 82, 15, 144, 122, 128, 178, 198, 113, 180, 245, 54,
            231, 173, 173, 119, 223, 221, 108, 189, 105, 154, 108, 189, 97, 24,
            108, 235, 5, 40, 171, 57, 208, 227, 176, 176, 245, 74, 136, 24, 160,
            36, 64, 89, 243, 60, 219, 122, 37, 68, 12, 80, 18, 160, 172, 248,
            183, 93, 170, 239, 123, 91, 175, 235, 58, 219, 122, 1, 202, 250, 25,
            80, 68, 44, 143, 43, 91, 239, 127, 95, 245, 241, 128, 210, 154, 6,
            149, 16, 241, 6, 241, 85, 103, 8, 80, 133, 170, 131, 226, 77, 161,
            132, 120, 62, 63, 191, 175, 215, 85, 31, 15, 40, 173, 57, 80, 233,
            125, 40, 34, 46, 151, 147, 173, 39, 127, 57, 5, 164, 14, 138, 103,
            168, 132, 136, 95, 131, 228, 111, 10, 128, 212, 65, 215, 235, 231,
            12, 149, 16, 49, 64, 73, 128, 212, 3, 164, 30, 32, 245, 0, 169, 7,
            72, 61, 64, 234, 1, 82, 15, 144, 122, 128, 212, 3, 164, 94, 115,
            160, 55, 102, 230, 236, 62, 98, 45, 95, 193, 0, 0, 0, 0, 73, 69, 78,
            68, 174, 66, 96, 130];
        let mut observed:Vec<u8> = vec![];
        open_file.read_to_end(&mut observed).ok();
        assert_eq!(expected, observed);

        plotter = plotter.as_dna().with_pixel_size(1);
        let cpath = tdir.join("temp2.png");
        plotter.save_png(cpath.as_path()).unwrap();
        let mut open_file = File::open(cpath).unwrap();
        let mut observed:Vec<u8> = vec![];
        open_file.read_to_end(&mut observed).ok();
        let expected = vec![137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73,
        72, 68, 82, 0, 0, 0, 44, 0, 0, 0, 50, 8, 6, 0, 0, 0, 39, 88, 57, 234,
        0, 0, 0, 201, 73, 68, 65, 84, 120, 156, 237, 212, 177, 13, 2, 49, 12,
        133, 97, 167, 65, 226, 234, 171, 82, 83, 51, 1, 67, 92, 205, 8, 212,
        140, 65, 205, 8, 76, 194, 44, 169, 82, 7, 137, 230, 176, 37, 114, 183,
        128, 159, 228, 72, 239, 239, 146, 34, 250, 20, 89, 78, 171, 38, 3, 69,
        48, 58, 130, 209, 17, 140, 142, 96, 116, 4, 163, 35, 24, 29, 193, 232,
        8, 70, 71, 48, 58, 130, 203, 52, 201, 121, 105, 82, 46, 73, 234, 253,
        40, 185, 53, 189, 245, 203, 29, 60, 95, 21, 253, 199, 90, 225, 193,
        223, 231, 142, 181, 194, 131, 109, 36, 122, 243, 227, 35, 135, 155, 235,
        243, 56, 176, 97, 243, 123, 149, 250, 210, 131, 99, 16, 112, 199, 218,
        44, 135, 255, 97, 155, 225, 142, 173, 139, 206, 112, 118, 125, 222, 31,
        188, 109, 9, 197, 202, 105, 128, 181, 102, 63, 220, 177, 182, 143, 227,
        207, 112, 73, 27, 214, 10, 15, 70, 71, 48, 58, 130, 209, 17, 140, 142,
        96, 116, 4, 163, 35, 24, 29, 193, 232, 8, 70, 71, 48, 186, 225, 192, 63,
        49, 69, 0, 122, 236, 10, 16, 179, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66,
        96, 130];
        assert_eq!(expected, observed);
    }
}