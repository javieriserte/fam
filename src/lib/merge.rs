use crate::seqs::{ApplyBufferedSequenceCollection, BufferedSeqCollection};
use crate::seqs::{
    AnnotatedSequence, SeqError, SequenceAccesors, SequenceCollection,
};
use std::collections::HashMap;
use std::collections::HashSet;

/// Join two or more sequence collections vertically
///
/// Assumes that seqs contains at least one sequence collection.
pub fn join<T: SequenceAccesors>(
    seqs: Vec<T>,
) -> Result<impl SequenceAccesors, SeqError> {
    let mut result = SequenceCollection::new();
    for sq in seqs {
        for s in sq.iter() {
            result.add(s.clone())?
        }
    }
    Ok(result)
}

/// Concatenates two or more sequences collections horizontally
///
/// Assumes that seqs contains at least one sequence collection.
pub fn concat<T: SequenceAccesors>(
    seqs: Vec<T>,
) -> Result<impl SequenceAccesors, SeqError> {
    let order = seqs
        .first()
        .unwrap()
        .iter()
        .map(|x| x.id())
        .collect::<Vec<&str>>();
    let mut result = seqs
        .first()
        .unwrap()
        .iter()
        .map(|x| (x.id(), vec![]))
        .collect::<HashMap<&str, Vec<char>>>();
    for sa in seqs.iter() {
        let mut used = order.iter().collect::<HashSet<_>>();
        for ann_seq in sa.iter() {
            match result.get_mut(ann_seq.id()) {
                None => {
                    return Err(SeqError::NonExistenId(String::from(
                        ann_seq.id(),
                    )))
                }
                Some(x) => {
                    x.extend(ann_seq.seq().unwrap());
                    used.remove(&ann_seq.id());
                }
            }
        }
        if !used.is_empty() {
            let a = used.iter().next().unwrap();
            return Err(SeqError::MissingID(String::from(**a)));
        }
    }
    let mut r = SequenceCollection::new();
    for o in order {
        let seq = AnnotatedSequence::new(
            o.to_string(),
            result.get(o).unwrap().clone(),
        );
        r.add(seq)?;
    }
    Ok(r)
}

mod test {
    #[allow(unused_imports)]
    use super::*;
    #[allow(unused_imports)]
    use crate::seqs::AnnotatedSequence;
    #[test]
    fn join_returns_the_same_with_one_seqcol_input() {
        let mut r = SequenceCollection::new();
        r.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let r2 = join(vec![r]).unwrap();
        assert_eq!(r2.size(), 2);
        assert_eq!(r2.get(0).unwrap().id(), "S1");
        assert_eq!(r2.get(1).unwrap().id(), "S2");
    }
    #[test]
    fn join_with_two_seqcol_input() {
        let mut r1 = SequenceCollection::new();
        r1.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r1.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let mut r2 = SequenceCollection::new();
        r2.add(AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("CTCG"),
        ))
        .unwrap();
        r2.add(AnnotatedSequence::from_string(
            String::from("S4"),
            String::from("DTCG"),
        ))
        .unwrap();
        let r3 = join(vec![r1, r2]).unwrap();
        assert_eq!(r3.size(), 4);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(2).unwrap().id(), "S3");
        assert_eq!(r3.get(3).unwrap().id(), "S4");
    }
    #[test]
    fn join_fails_with_two_equal_ids() {
        let mut r1 = SequenceCollection::new();
        r1.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r1.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let mut r2 = SequenceCollection::new();
        r2.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("CTCG"),
        ))
        .unwrap();
        r2.add(AnnotatedSequence::from_string(
            String::from("S4"),
            String::from("DTCG"),
        ))
        .unwrap();
        let r3 = join(vec![r1, r2]);
        match r3 {
            Err(SeqError::DuplicatedId(x)) => assert_eq!(&x, "S1"),
            _ => panic!(
                "This should fail with \
                SeqError::DuplicatedID('S1')"
            ),
        }
    }
    #[test]
    fn concat_returns_the_same_with_one_sequence() {
        let mut r = SequenceCollection::new();
        r.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let r2 = concat(vec![r]).unwrap();
        assert_eq!(r2.size(), 2);
        assert_eq!(r2.get(0).unwrap().id(), "S1");
        assert_eq!(r2.get(1).unwrap().id(), "S2");
    }
    #[test]
    fn concat_with_two_seqcol_input() {
        let mut r1 = SequenceCollection::new();
        r1.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r1.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let mut r2 = SequenceCollection::new();
        r2.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("CTCG"),
        ))
        .unwrap();
        r2.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("DTCG"),
        ))
        .unwrap();
        let r3 = concat(vec![r1, r2]).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCGCTCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
    }
    #[test]
    fn concat_fails_with_different_ids() {
        let mut r1 = SequenceCollection::new();
        r1.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r1.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let mut r2 = SequenceCollection::new();
        r2.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("CTCG"),
        ))
        .unwrap();
        r2.add(AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("DTCG"),
        ))
        .unwrap();
        let r3 = concat(vec![r1, r2]);
        match r3 {
            Err(SeqError::NonExistenId(x)) => assert_eq!(x, "S3"),
            _ => panic!("Should throw SeqError::NonExistentId('S3')"),
        }
    }
    #[test]
    fn concat_fails_with_less_seqs() {
        let mut r1 = SequenceCollection::new();
        r1.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("ATCG"),
        ))
        .unwrap();
        r1.add(AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("BTCG"),
        ))
        .unwrap();
        let mut r2 = SequenceCollection::new();
        r2.add(AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("CTCG"),
        ))
        .unwrap();
        let r3 = concat(vec![r1, r2]);
        match r3 {
            Err(SeqError::MissingID(x)) => assert_eq!(x, "S2"),
            Err(x) => eprintln!("Error: {}", x),
            _ => panic!("Should throw SeqError::MissingID('S2')"),
        }
    }
}
