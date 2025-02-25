use crate::seqs::{
    AnnotatedSequence, SeqError, SequenceAccesors, SequenceCollection,
};
use std::collections::HashSet;

/// Combines two or more sequence collections vertically
///
/// Assumes that seqs contains at least one sequence collection.
pub fn concat<T: SequenceAccesors>(
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

/// Combines two or more sequences collections horizontally
///
/// - Ignores identifiers, identifiers returned are the ones of the
/// first sequence.
/// - Assumes that seqs contains at least one sequence collection.
/// - Assumes that seqs has the same size.
pub fn join<T: SequenceAccesors>(
    seqs: Vec<T>,
) -> Result<impl SequenceAccesors, SeqError> {
    let sequence_ids = seqs
        .get(0)
        .ok_or(SeqError::Empty)?
        .iter()
        .map(|x| x.id())
        .collect::<Vec<&str>>();
    let mut seq_accumulator = sequence_ids
        .iter()
        .map(|_| vec![])
        .collect::<Vec<Vec<char>>>();
    seqs.iter()
        .for_each(
            |sc| seq_accumulator
                    .iter_mut()
                    .zip(sc.iter())
                    .for_each(
                        |(&mut ref mut r, s)| {
                            r.extend(s.seq().unwrap_or(&vec![]));
                        }
                    )
            );
    let annotated_sequences = sequence_ids
        .iter()
        .zip(seq_accumulator.into_iter())
        .map(|(sid, vecseq)| (sid.to_string(), vecseq.to_vec()))
        .map(|(x, y)| AnnotatedSequence::new(x, y))
        .collect::<Vec<_>>();
    let mut sequence_collection = SequenceCollection::new();
    annotated_sequences.into_iter()
        .map(|x| { sequence_collection.add(x) })
        .into_iter()
        .try_for_each(|x| x)?;
    Ok(sequence_collection)
}

/// Combines two or more sequences collections horizontally by its identifier.
///
/// - Assumes that seqs contains at least one sequence collection.
/// - Tries to preserve the order of the first sequence collection.
/// - If outer is true, then the resulting sequence will contain all the
/// identifiers
/// - If outer is false, then the resulting sequence will contain only the
/// identifiers that are present in all the sequence collections.
pub fn merge<T:SequenceAccesors>(
    seqs: Vec<T>,
    outer: bool
) -> Result<impl SequenceAccesors, SeqError> {
    let order = seqs
        .get(0)
        .ok_or(SeqError::Empty)?
        .iter()
        .map(|x| x.id())
        .collect::<Vec<&str>>();
    let all_identifiers_set = seqs
        .iter()
        .map(|seq| seq.iter().map(|y| y.id()))
        .into_iter()
        .flatten()
        .collect::<HashSet<_>>();
    let final_identifiers = match outer {
        true => {
            let order_set = order.iter().collect::<HashSet<_>>();
            let extra_ids = all_identifiers_set
                .into_iter()
                .filter(|x| !order_set.contains(x))
                .collect::<Vec<_>>();
            order
                .iter()
                .map(|x| *x)
                .chain(extra_ids)
                .collect::<Vec<_>>()
        },
        false => {
            let common_ids = seqs
                .iter()
                .map(|x| x.iter().map(|y| y.id()).collect::<HashSet<_>>())
                .reduce(
                    |a, b|
                        a.intersection(&b)
                            .cloned()
                            .collect::<HashSet<_>>()
                )
                .ok_or(SeqError::Empty)?;
            order
                .iter()
                .filter(|x| common_ids.contains(*x))
                .map(|x| *x)
                .collect::<Vec<_>>()
        }
    };
    let mut seq_accumulator = final_identifiers
        .iter()
        .map(|_| vec![])
        .collect::<Vec<Vec<char>>>();
    seqs.iter()
        .for_each(
            |sc|
                final_identifiers
                    .iter()
                    .enumerate()
                    .for_each(
                        |(i, id)| {
                            let seq = sc
                                .get_by_id(*id)
                                .map(|x| x.seq())
                                .flatten()
                                .map(|x| x.to_vec())
                                .unwrap_or(vec![]);
                            seq_accumulator[i].extend(seq);
                        }
                    )
            );
    let annotated_sequences = final_identifiers
        .iter()
        .zip(seq_accumulator.into_iter())
        .map(|(sid, vecseq)| (sid.to_string(), vecseq.to_vec()))
        .map(|(x, y)| AnnotatedSequence::new(x, y))
        .collect::<Vec<_>>();
    let mut sequence_collection = SequenceCollection::new();
    annotated_sequences.into_iter()
        .map(|x| { sequence_collection.add(x) })
        .into_iter()
        .try_for_each(|x| x)?;
    Ok(sequence_collection)
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
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCGCTCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
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
            String::from("S3"),
            String::from("CTCG"),
        ))
        .unwrap();
        r2.add(AnnotatedSequence::from_string(
            String::from("S4"),
            String::from("DTCG"),
        ))
        .unwrap();
        let r3 = concat(vec![r1, r2]).unwrap();
        assert_eq!(r3.size(), 4);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(2).unwrap().id(), "S3");
        assert_eq!(r3.get(3).unwrap().id(), "S4");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCG");
        assert_eq!(r3.get(2).unwrap().seq_as_string(), "CTCG");
        assert_eq!(r3.get(3).unwrap().seq_as_string(), "DTCG");
    }
    #[test]
    fn concat_fails_with_two_equal_ids() {
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
        let r3 = concat(vec![r1, r2]);
        match r3 {
            Err(SeqError::DuplicatedId(x)) => assert_eq!(&x, "S1"),
            _ => panic!(
                "This should fail with \
                SeqError::DuplicatedID('S1')"
            ),
        }
    }
    #[test]
    fn merge_with_single_collection() {
        let r1 = vec![
            ("S1", "ATCG"),
            ("S2", "BTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r3 = merge(vec![r1.clone()], false).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        let r3 = merge(vec![r1], true).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
    }
    #[test]
    fn merge_with_two_seqcol_and_no_missing_id() {
        let r1 = vec![
            ("S1", "ATCG"),
            ("S2", "BTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r2 = vec![
            ("S1", "CTCG"),
            ("S2", "DTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r3 = merge(vec![r1.clone(), r2.clone()], false).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCGCTCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
        let r3 = merge(vec![r1, r2], false).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCGCTCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
    }
    #[test]
    fn merge_with_two_seqcol_and_missing_id() {
        let r1 = vec![
            ("S1", "ATCG"),
            ("S2", "BTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r2 = vec![
            ("S2", "DTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r3 = merge(vec![r1.clone(), r2.clone()], true).unwrap();
        assert_eq!(r3.size(), 2);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
        let r3 = merge(vec![r1.clone(), r2.clone()], false).unwrap();
        assert_eq!(r3.size(), 1);
        assert_eq!(r3.get(0).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "BTCGDTCG");
    }
    #[test]
    fn merge_with_two_seqcol_and_extra_ids() {
        let r1 = vec![
            ("S1", "ATCG"),
            ("S2", "BTCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r2 = vec![
            ("S2", "DTCG"),
            ("S3", "YMCG"),
        ]
            .into_iter()
            .collect::<SequenceCollection>();
        let r3 = merge(vec![r1.clone(), r2.clone()], true).unwrap();
        assert_eq!(r3.size(), 3);
        assert_eq!(r3.get(0).unwrap().id(), "S1");
        assert_eq!(r3.get(1).unwrap().id(), "S2");
        assert_eq!(r3.get(2).unwrap().id(), "S3");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "ATCG");
        assert_eq!(r3.get(1).unwrap().seq_as_string(), "BTCGDTCG");
        assert_eq!(r3.get(2).unwrap().seq_as_string(), "YMCG");
        let r3 = merge(vec![r1.clone(), r2.clone()], false).unwrap();
        assert_eq!(r3.size(), 1);
        assert_eq!(r3.get(0).unwrap().id(), "S2");
        assert_eq!(r3.get(0).unwrap().seq_as_string(), "BTCGDTCG");
    }
}
