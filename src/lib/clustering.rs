use std::{cmp::max, error::Error};
use std::io::ErrorKind;

use crate::seqs::{Alignment, SequenceAccesors};

#[derive(Clone)]
pub struct Cluster {
    pub representative: usize,
    pub members: Vec<usize>
}

impl Cluster {
    pub fn new(representative: usize) -> Self {
        Cluster {
            representative,
            members:vec![]
        }
    }
    pub fn with_members(mut self, members: Vec<usize>) -> Self {
        self.members = members;
        self
    }
}

pub trait Clusterer {
    fn clusterize(
        &mut self)
        -> Result<(), Box<dyn Error>>;
    fn weights(&self) -> Result<Vec<f64>, Box<dyn Error>>;
    fn take_clusters(&mut self) -> Option<Vec<Cluster>>;
    fn clear(&mut self);
    fn clusters(&self) -> Option<Vec<Cluster>>;


}
/// Computes the similarity of two sequences as the fraction of
/// identical characters in both sequences.
///
/// Returns an error is sequences length is not the same.
///
/// If exclude_gaps is true, positions in which both sequences are
/// gaps are ignored to compute the total length.
///
/// Example:
/// ```
/// use famlib::clustering::identity_fraction;
/// let seq1 = vec!['A', 'B', 'C', '-', 'E'];
/// let seq2 = vec!['A', 'B', 'C', '-', '-'];
/// let id = identity_fraction(&seq1, &seq2, false).unwrap();
/// assert_eq!(id, 0.6);
/// let id = identity_fraction(&seq1, &seq2, true).unwrap();
/// assert_eq!(id, 0.75);
/// ```
pub fn identity_fraction(
        seq1: &Vec<char>,
        seq2: &Vec<char>,
        exclude_gaps: bool)
        -> Result<f64, Box<dyn Error>> {
    if seq1.len() == seq2.len() {
        let identical:usize = seq1
            .iter()
            .zip(seq2)
            .map(|(a, b)| {
                match b==a && *a != '-' {
                    true => 1,
                    false => 0
                }})
            .sum();
        let total = match exclude_gaps {
            true => seq1
                .iter()
                .zip(seq2)
                .map(|(a, b)| {
                    match *a!='-' || *b!='-' {
                        true => 1,
                        false => 0
                    }})
                .sum(),
            false => seq1.len()
        };
        Ok(identical as f64 / max(1, total) as f64)
    } else {
        Err(Box::new(std::io::Error::new(
            ErrorKind::Other,
            format!("Sequences have different length"),
        )))
    }
}

pub struct Hobohm1<'a> {
    msa: &'a Alignment,
    similarity: f64,
    exclude_gaps: bool,
    clusters: Option<Vec<Cluster>>
}

impl <'a> Hobohm1<'a> {
    pub fn new(msa: &'a Alignment) -> Self {
        Hobohm1{
            msa,
            similarity: 0.62,
            exclude_gaps: false,
            clusters: None
        }
    }
    pub fn with_similarity(mut self, similarity: f64) -> Self {
        self.similarity = similarity;
        self
    }
    pub fn exclude_gaps(mut self) -> Self {
        self.exclude_gaps = true;
        self
    }
    pub fn include_gaps(mut self) -> Self {
        self.exclude_gaps = false;
        self
    }
}

impl <'a> Clusterer for Hobohm1<'a> {
    fn clusterize(&mut self) -> Result<(), Box<dyn Error>> {
        let candidates = 0..self.msa.size();
        let mut members:Vec<Vec<usize>> = vec![];
        let mut repr_indexes:Vec<usize> = vec![];
        for index in candidates {
            let cseq = self.msa.get(index).unwrap().seq().unwrap();
            let mut add_cluster = true;
            for (i, repr) in repr_indexes.iter().enumerate() {
                match identity_fraction(
                        cseq,
                        self.msa.get(*repr).unwrap().seq().unwrap(),
                        false) {
                    Ok(sim) => {
                        if sim >= self.similarity {
                            // Add index to current cluster
                            members[i].push(index);
                            add_cluster = false;
                            break;
                        };
                    }
                    Err(_) => {}
                }
            }
            if add_cluster {
                repr_indexes.push(index);
                members.push(vec![index]);
            }
        }
        self.clusters = Some(repr_indexes.iter().zip(members).map(
            |(r, m)| Cluster::new(*r).with_members(m)
        ).collect());
        Ok(())
    }

    fn weights(&self) -> Result<Vec<f64>, Box<dyn Error>> {
        match &self.clusters {
            Some(clusters) => {
                let max_index = clusters
                    .iter()
                    .flat_map(|x| x.members.iter())
                    .max()
                    .unwrap();
                let mut ws = vec![0f64; *max_index+1];
                clusters
                    .iter()
                    .for_each(|x| {
                        x.members.iter().for_each(
                            |y| ws[*y] = 1f64/x.members.len() as f64
                        )
                    });
                Ok(ws)
            }
            None => {
                Err(Box::new(std::io::Error::new(
                    ErrorKind::Other,
                    format!("There are no clusters to compute weights"),
                )))
            }
        }
    }
    fn clear(&mut self) {
        self.clusters = None;
    }
    fn take_clusters(&mut self) -> Option<Vec<Cluster>> {
        self.clusters.take()
    }
    fn clusters(&self) -> Option<Vec<Cluster>> {
        self.clusters.clone()
    }

}
#[cfg(test)]
mod test {
    use crate::seqs::AnnotatedSequence;
    use crate::seqs::SequenceCollection;

    use super::*;
    #[test]
    fn test_identity_fraction() {
        let seq1 = vec!['A'; 100];
        let seq2 = vec!['A'; 100];
        let seq3 = vec!['A'; 98]
            .into_iter()
            .chain(vec!['B'; 2].into_iter())
            .collect::<Vec<char>>();
        let seq4 = vec!['A'; 62]
            .into_iter()
            .chain(vec!['B'; 38].into_iter())
            .collect::<Vec<char>>();
        let seq5 = vec!['A'; 50]
            .into_iter()
            .chain(vec!['-'; 50].into_iter())
            .collect::<Vec<char>>();
        let seq6 = vec!['A'; 75]
            .into_iter()
            .chain(vec!['-'; 25].into_iter())
            .collect::<Vec<char>>();
        let seq7 = vec!['A'; 101];
        let seq8 = vec!['-'; 100];
        let seq9 = vec!['-'; 100];
        let seq10 = vec!['.'; 100];
        let seq11:Vec<char> = vec![];
        let seq12:Vec<char> = vec![];
        // seq should be the proportion of identical characters in two strings
        let f = identity_fraction;
        assert_eq!(f(&seq1, &seq2, false).unwrap(), 1.0);
        assert_eq!(f(&seq1, &seq2, true).unwrap(), 1.0);
        assert_eq!(f(&seq1, &seq3, false).unwrap(), 0.98);
        assert_eq!(f(&seq1, &seq3, true).unwrap(), 0.98);
        assert_eq!(f(&seq1, &seq4, false).unwrap(), 0.62);
        assert_eq!(f(&seq1, &seq4, true).unwrap(), 0.62);

        // gap should be taken care
        assert_eq!(f(&seq1, &seq5, false).unwrap(), 0.5);
        assert_eq!(f(&seq1, &seq5, true).unwrap(), 0.5);
        assert_eq!(f(&seq5, &seq6, false).unwrap(), 2.0/4.0);
        assert_eq!(f(&seq5, &seq6, true).unwrap(), 2.0/3.0);

        // ValueError should be thrown if the sequences length is not the same
        assert!(f(&seq1, &seq7, false).is_err());

        // Seq id should be 0 if there are not non-gap chars in sequence
        assert_eq!(f(&seq8, &seq9, false).unwrap(), 0f64);
        assert_eq!(f(&seq8, &seq10, false).unwrap(), 0f64);
        assert_eq!(f(&seq11, &seq12, false).unwrap(), 0f64);
    }
    #[test]
    fn test_clustering() {
        let mut sequences = SequenceCollection::new();
        [
            "AAAAA",
            "AAAAB",
            "AAABB",
            "AAABC",
            "AACDE",
            "AACFG"
        ].iter()
        .enumerate()
        .map(|(i, x)|
            AnnotatedSequence::from_string(
                format!("Seq_{}", i),
                String::from(*x)
            )
        ).for_each(
            |ann| sequences.add(ann).unwrap()
        );
        let sequences = sequences.to_msa().ok().unwrap();
        let mut cl = Hobohm1::new(&sequences);
        cl.clusterize().unwrap();

        let clusters = cl.clusters().unwrap();
        assert_eq!(clusters.len(), 4);
        assert_eq!(clusters[0].members.len(), 2);
        assert_eq!(clusters[1].members.len(), 2);
        assert_eq!(clusters[2].members.len(), 1);
        assert_eq!(clusters[3].members.len(), 1);

        let expected_weight = [0.5, 0.5, 0.5, 0.5, 1.0, 1.0];
        cl.weights()
            .unwrap()
            .iter()
            .zip(expected_weight.iter())
            .for_each(|(a, b)| assert_eq!(a, b))
    }
}