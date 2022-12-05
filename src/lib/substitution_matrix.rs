use std::{error::Error, io::ErrorKind};
use crate::seqs::Alignment;
use crate::conservation::amino_index;
use itertools::Itertools;

pub trait SubstitutionMatrix {
    fn build_from(&mut self, msa: &Alignment) -> Result<(), Box<dyn Error>>;
}
pub struct ProteinMatrix {
    pub data: [usize; Self::SIZE],
}
impl ProteinMatrix {
    const ALPHABET_SIZE: usize = 20;
    const SIZE:usize = Self::ALPHABET_SIZE * (Self::ALPHABET_SIZE +1)/2;
    pub fn new() -> Self {
        Self{
            data: [0; 210]
        }
    }
    pub fn set(
            &mut self,
            x: usize,
            y: usize,
            value: usize)
            -> Result<(), Box<dyn Error>> {
        let pos = ProteinMatrix::_from_xy(x, y)?;
        self.data[pos] = value;
        Ok(())
    }
    pub fn get(&self, x:usize, y:usize) -> Result<usize, Box<dyn Error>> {
        let pos = ProteinMatrix::_from_xy(x, y)?;
        Ok(self.data[pos])
    }
    pub fn _from_xy(x: usize, y: usize) -> Result<usize, Box<dyn Error>> {
        let pos = match y>=x {
            true => x*(x+1)/2+y,
            false => y*(y+1)/2+x
        };
        if pos < Self::SIZE &&
            x < Self::ALPHABET_SIZE &&
            y < Self::ALPHABET_SIZE {
            Ok(pos)
        } else {
            Err(
                Box::new(
                    std::io::Error::new(
                        ErrorKind::Other,
                        format!("Index out of range")
                    )
                )
            )
        }
    }
}
impl SubstitutionMatrix for ProteinMatrix {
    fn build_from(&mut self, msa: &Alignment) -> Result<(), Box<dyn Error>> {
        let mut data = [0; Self::SIZE];
        for col in msa.columns() {
            for pair in col.iter().combinations(2) {
                println!("{:?}", pair);
                let c1 = *pair[0];
                let c2 = *pair[1];
                let x = amino_index(*c1);
                let y = amino_index(*c2);
                match Self::_from_xy(x, y) {
                    Ok(pos) => data[pos] += 1,
                    Err(_) => {}
                }
            }
        }
        self.data = data;
        Ok(())
    }
}
#[cfg(test)]
mod test {
    use crate::seqs::{Alignment, AnnotatedSequence, SequenceAccesors};
    use super::{ProteinMatrix, SubstitutionMatrix};

    fn build_msa(seqs: Vec<&str>) -> Alignment {
        let mut msa = Alignment::new();
        for (i, s) in seqs.iter().enumerate(){
            let seqid = format!("Seq_{}", i);
            let ann = AnnotatedSequence::new(seqid, s.chars().collect());
            msa.add(ann).unwrap();
        }
        msa
    }
    #[test]
    fn test_protmatrix_with_one_column_all_equals() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["Q", "Q", "Q", "Q"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(6));
    }
    #[test]
    fn test_protmatrix_with_one_column_one_diff() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["Q", "Q", "Q", "A"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(3));
        assert_eq!(pm.get(13, 0).ok(), Some(3));
        assert_eq!(pm.get(0, 13).ok(), Some(3));
        assert_eq!(pm.get(0, 0).ok(), Some(0));
    }
    #[test]
    fn test_protmatrix_with_one_column_with_pairs() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["Q", "Q", "A", "A"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(1));
        assert_eq!(pm.get(13, 0).ok(), Some(4));
        assert_eq!(pm.get(0, 13).ok(), Some(4));
        assert_eq!(pm.get(0, 0).ok(), Some(1));
    }
    #[test]
    fn test_protmatrix_with_one_column_one_pair() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["Q", "Q", "A", "C"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(1));
        assert_eq!(pm.get(13, 0).ok(), Some(2));
        assert_eq!(pm.get(0, 13).ok(), Some(2));
        assert_eq!(pm.get(13, 1).ok(), Some(2));
        assert_eq!(pm.get(1, 13).ok(), Some(2));
        assert_eq!(pm.get(0, 1).ok(), Some(1));
        assert_eq!(pm.get(0, 0).ok(), Some(0));
        assert_eq!(pm.get(1, 1).ok(), Some(0));
    }
    #[test]
    fn test_protmatrix_with_one_column_all_different() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["E", "Q", "A", "C"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(3, 3).ok(), Some(0));
        assert_eq!(pm.get(3, 13).ok(), Some(1));
        assert_eq!(pm.get(3, 0).ok(), Some(1));
        assert_eq!(pm.get(3, 1).ok(), Some(1));
        assert_eq!(pm.get(13, 13).ok(), Some(0));
        assert_eq!(pm.get(13, 0).ok(), Some(1));
        assert_eq!(pm.get(13, 1).ok(), Some(1));
        assert_eq!(pm.get(0, 0).ok(), Some(0));
        assert_eq!(pm.get(0, 1).ok(), Some(1));
        assert_eq!(pm.get(1, 1).ok(), Some(0));
    }
    #[test]
    fn test_protmatrix_with_two_columns() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["QQ", "QA", "QQ", "QC"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(7));
        assert_eq!(pm.get(13, 0).ok(), Some(2));
        assert_eq!(pm.get(13, 1).ok(), Some(2));
        assert_eq!(pm.get(0, 1).ok(), Some(1));
    }
    #[test]
    fn test_protmatrix_with_two_columns_with_gaps() {
        let mut pm = ProteinMatrix::new();
        let seqs = vec!["QQ", "Q-", "Q-", "QA"];
        let msa = build_msa(seqs);
        pm.build_from(&msa).unwrap();
        assert_eq!(pm.get(13, 13).ok(), Some(6));
        assert_eq!(pm.get(13, 0).ok(), Some(1));
        assert_eq!(pm.get(0, 0).ok(), Some(0));
    }
}
