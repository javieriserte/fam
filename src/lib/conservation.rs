use crate::seqs::Alignment;

pub trait Conservation {
    fn protein_entropy(&self) -> Vec<f64>;
    fn dna_entropy(&self) -> Vec<f64>;
}

fn dna_index(c:char) -> usize {
    match c {
        'A' => 0,
        'C' => 1,
        'T' => 2,
        'G' => 3,
        '-' => 4,
        _ => 5
    }
}

fn amino_index(c: char) -> usize {
    match c {
        'A' => 0,
        'C' => 1,
        'D' => 2,
        'E' => 3,
        'F' => 4,
        'G' => 5,
        'H' => 6,
        'I' => 7,
        'K' => 8,
        'L' => 9,
        'M' => 10,
        'N' => 11,
        'P' => 12,
        'Q' => 13,
        'R' => 14,
        'S' => 15,
        'T' => 16,
        'V' => 17,
        'W' => 18,
        'Y' => 19,
        '-' => 20,
        _ => 21
    }
}

fn _entropy(
        msa: &Alignment,
        alphabet_size: usize,
        f_index: &dyn Fn(char)->usize)
        -> Vec<f64> {
    let mut aas = vec![0usize; alphabet_size+2];
    msa.columns()
        .into_iter()
        .map(
            |col| {
                aas.iter_mut().for_each(|x| *x = 0);
                col.iter().for_each(|x|
                    aas[f_index(**x)]+=1
                );
                let nchars: usize = aas[0..alphabet_size].iter().sum();
                aas.iter().fold(
                    0f64,
                    |a, b| {
                        let pi = (*b as f64)/(nchars as f64);
                        if pi>0f64{
                            a - pi * pi.log(2f64)
                        } else {
                            a
                        }
                    }
                )
            })
        .collect()
}

impl Conservation for Alignment{
    fn protein_entropy(&self) -> Vec<f64> {
        _entropy(self, 20, &amino_index)
    }

    fn dna_entropy(&self) -> Vec<f64> {
        _entropy(self, 4, &dna_index)
    }
}

    #[cfg(test)]
    mod test{
        use crate::seqs::{AnnotatedSequence, SequenceAccesors, SequenceCollection};

        use super::Conservation;
        #[test]
        fn test_iter_mut() {
            let mut a = vec![1,2,3,4,5];
            a.iter_mut().for_each(|x| *x = 0);
            println!("{:?}", a);
        }
        #[test]
        fn test_entropy_proteins() {
            let mut input = SequenceCollection::new();
            let seqs = vec![
                "SHSMRYFYTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIRNTRNVKAQSQTD",
                "SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIGETRKVKAHSQTH",
                "SHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRTEPRPPWIRNTQIFKTNTQTY",
                "SHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRTEPRPPWIRNTQIFKTNTQTY",
                "SHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRTEPRAPWIRNTQIFKTNTQTY",
                "SHSMRYFHTSVSRPGRGEPRFITVGYVDDTLFVRFDSDAASPREEPRAPWIRETQICKAKAQTD",
                "SHSMRYFYTAVSRPGRGEPHFIAVGYVDDTQFVRFDSDAASPRGEPRAPWVRETQKYKRQAQTD",
                "SHSMRYFSTSVSWPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPREPWVRETQKYKRQAQAD",
                "PHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMRETQKAKGNEQSF",
                "PHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMRETQKAKGQEQWF"
        ];
        seqs.iter()
            .enumerate()
            .map(|(i, x)|
                AnnotatedSequence::from_string(
                    String::from(format!("Seq_{}", i)),
                    String::from(*x))
                )
            .for_each(|x| input.add(x).ok().unwrap());
        let expected = vec![
            0.722, 0.000, 0.000, 0.469, 0.000, 0.000, 0.000, 2.161,
            0.000, 0.971, 0.881, 0.000, 0.469, 0.000, 0.000, 0.722,
            0.469, 0.000, 0.000, 0.469, 0.722, 0.469, 1.357, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.469, 0.469, 1.157, 0.000,
            0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.722,
            0.722, 0.722, 0.000, 2.246, 0.000, 0.000, 0.000, 1.157,
            0.469, 0.000, 1.371, 0.469, 0.971, 0.000, 0.722, 1.361,
            2.246, 0.000, 1.971, 1.722, 1.971, 0.000, 1.357, 1.846,];
        let real = input.to_msa().map(|x| {
                x.protein_entropy()}).ok().unwrap();
        assert_eq!(real.len(), 64);
        assert!(real.iter().zip(expected).all(|(a, b)| (*a-b).abs()<0.001f64));

    }
}