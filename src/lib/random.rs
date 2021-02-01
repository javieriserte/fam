use crate::seqs::{Alignment, SequenceAccesors};
extern crate rand;
use rand::{prelude::SliceRandom, thread_rng};
pub trait RandomGen {
    fn shuffle(&mut self, fixed_gaps: bool);
}

impl RandomGen for Alignment {
    fn shuffle(&mut self, fixed_gaps: bool) {
        let mut rng = thread_rng();
        match fixed_gaps {
            true=> {
                for i in 0..self.size() {
                    let seq = self.get_mut(i).unwrap().seq_mut().unwrap();
                    let mut chars = seq
                        .iter()
                        .cloned()
                        .filter(|x| *x != '-')
                        .collect::<Vec<_>>();
                    let slice: &mut [char] = &mut chars;
                    slice.shuffle(&mut rng);
                    let mut j=0;
                    for i in 0..seq.len() {
                        if seq[i] != '-' {
                            seq[i] = slice[j];
                            j+=1;
                        }
                    }
                }
            },
            false=> {}
        };
    }
}

#[cfg(test)]
mod test {
    use crate::seqs::{Alignment, AnnotatedSequence, SequenceAccesors};

    use super::RandomGen;

    fn sample_msa() -> Alignment {
        let mut msa = Alignment::new();
        let s1 = AnnotatedSequence::from_string(
            String::from("s1"),
            String::from("ACT--AC--GA--TGGABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"),
        );
        let s2 = AnnotatedSequence::from_string(
            String::from("s2"),
            String::from("CCTG-CTG-CTG-CT-ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"),
        );
        let s3 = AnnotatedSequence::from_string(
            String::from("s3"),
            String::from("-CAGACAGACAGACA-ABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"),
        );
        msa.add(s1).unwrap();
        msa.add(s2).unwrap();
        msa.add(s3).unwrap();
        msa
    }
    #[test]
    fn test_random_msa() {
        let mut msa=sample_msa();
        let prev = msa.get(0).unwrap().clone().seq_as_string();
        let gapped_prev = prev
            .chars()
            .enumerate()
            .filter_map(
                |(i, c)| if c=='-' { None } else {Some(i)})
            .collect::<Vec<_>>();


        msa.shuffle(true);

        let post = msa.get(0).unwrap().clone().seq_as_string();
        let gapped_post = post
            .chars()
            .enumerate()
            .filter_map(
                |(i, c)| if c=='-' { None } else {Some(i)})
            .collect::<Vec<_>>();
        assert_ne!(prev, post);
        assert_eq!(gapped_prev, gapped_post);

    }
}

