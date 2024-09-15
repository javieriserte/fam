
use std::error::Error;

use itertools::Itertools;

use crate::conservation::amino_index;
use crate::conservation::dna_index;
use crate::matrices::triangular_matrix::TriangularMatrix;
use crate::matrices::triangular_matrix::Num;
use crate::seqs::Alignment;

pub trait SubstitutionMatrix<T> where T: Num
{
    type Output;
    fn build_from(&mut self, msa: &Alignment) -> Result<(), Box<dyn Error>>;
    fn get(&self, x:usize, y:usize) -> Result<T, Box<dyn Error>>;
    fn set(&mut self, x:usize, y:usize, value: T) -> Result<(), Box<dyn Error>>;
    fn normalize(&self) -> Result<Self::Output, Box<dyn Error>>;
}

pub struct ProteinMatrix<T> where T: Num
{
    pub data: TriangularMatrix<T>
}

impl <T>ProteinMatrix<T> where T: Num
{
    const ALPHABET_SIZE: usize = 20;
    pub fn new() -> Self {
        Self{
            data: TriangularMatrix::new(Self::ALPHABET_SIZE)
        }
    }
}

impl <T>SubstitutionMatrix<T> for ProteinMatrix<T>
where T: Num
{
    type Output = ProteinMatrix<f64>;
    fn build_from(&mut self, msa: &Alignment) -> Result<(), Box<dyn Error>> {
        let mut data = TriangularMatrix::new(Self::ALPHABET_SIZE);
        for col in msa.columns() {
            for pair in col.iter().combinations(2) {
                let c1 = *pair[0];
                let c2 = *pair[1];
                let x = amino_index(*c1);
                let y = amino_index(*c2);
                data.increment(x, y).ok();
            }
        }
        self.data = data;
        Ok(())
    }
    fn set(
        &mut self,
        x: usize,
        y: usize,
        value: T)
    -> Result<(), Box<dyn Error>> {
        self.data.set(x, y, value)?;
        Ok(())
    }
    fn get(&self, x:usize, y:usize) -> Result<T, Box<dyn Error>> {
        self.data.get(x, y)
    }
    fn normalize(&self) -> Result<Self::Output, Box<dyn Error>> {
        let data = TriangularMatrix::new(Self::ALPHABET_SIZE);
        let mut pm = ProteinMatrix::new();
        pm.data = data;
        Ok(pm)
    }
}

pub struct DNAMatrix<T>
where T: Num {
    data: TriangularMatrix<T>
}

impl <T> DNAMatrix<T> where T:Num {
    const ALPHABET_SIZE: usize = 20;
    pub fn new() -> Self {
        Self{
            data: TriangularMatrix::new(Self::ALPHABET_SIZE)
        }
    }
}
impl <T> SubstitutionMatrix<T> for DNAMatrix<T> where T:Num {
    type Output = DNAMatrix<f64>;
    fn build_from(&mut self, msa: &Alignment) -> Result<(), Box<dyn Error>> {
        let mut data = TriangularMatrix::new(Self::ALPHABET_SIZE);
        for col in msa.columns() {
            for pair in col.iter().combinations(2) {
                let c1 = *pair[0];
                let c2 = *pair[1];
                let x = dna_index(*c1);
                let y = dna_index(*c2);
                data.increment(x, y)?;
            }
        }
        self.data = data;
        Ok(())
    }
    fn set(
        &mut self,
        x: usize,
        y: usize,
        value: T
    )
    -> Result<(), Box<dyn Error>> {
        self.data.set(x, y, value)?;
        Ok(())
    }
    fn get(&self, x:usize, y:usize) -> Result<T, Box<dyn Error>> {
        self.data.get(x, y)
    }

    fn normalize(&self) -> Result<Self::Output, Box<dyn Error>> {
        todo!()
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
        assert_eq!(pm.get(13, 13).ok(), Some(6.0));
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
