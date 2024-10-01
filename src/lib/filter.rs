use crate::seqs::SequenceAccesors;
use crate::seqs::Alignment;
use crate::seqs::AnnotatedSequence;
use crate::seqs::ApplyBufferedSequenceCollection;
use crate::seqs::SequenceCollection;
use crate::seqs::BufferedSeqCollection;
use regex::RegexBuilder;

pub trait Filter<T> {
    fn filter_regex_id(
        &self,
        regex: &str,
        ignore_case: bool,
        keep: bool
    ) -> T;
}

pub struct FilterBufferedSequenceCollection{}

impl FilterBufferedSequenceCollection {
    pub fn filter_regex_id(
        bsc: Box<dyn BufferedSeqCollection>,
        ignore_case:bool,
        keep: bool,
        pattern: &str
    ) -> ApplyBufferedSequenceCollection {
        let re = RegexBuilder::new(pattern)
            .case_insensitive(ignore_case)
            .build()
            .unwrap();
        let filter_func = move |s: AnnotatedSequence| {
            if re.is_match(s.id()) ^ !keep {
                vec![s]
            } else {
                vec![]
            }
        };
        ApplyBufferedSequenceCollection::new(
            bsc,
            Box::new(filter_func)
        )
    }
}

impl Filter<SequenceCollection> for SequenceCollection {
    fn filter_regex_id(
        &self,
        regex: &str,
        ignore_case: bool,
        keep: bool
    ) -> SequenceCollection {
        let re = RegexBuilder::new(regex)
            .case_insensitive(ignore_case)
            .build()
            .unwrap();
        let mut result = SequenceCollection::new();
        self.iter()
            .for_each(
                |x| match re.is_match(x.id()) ^ !keep {
                    true => result.add(x.clone()).unwrap(),
                    false => ()
                }
            );
        result
    }
}

impl Filter<Alignment> for Alignment {
    fn filter_regex_id(
        &self,
        regex: &str,
        ignore_case: bool,
        keep: bool
    ) -> Alignment {
        self.seq_col()
            .filter_regex_id(regex, ignore_case, keep)
            .to_msa()
            .expect("This shouldn't fail, comes and go from and Alignment.")
    }
}

#[cfg(test)]
mod test {
    use crate::seqs::AnnotatedSequence;
    use crate::seqs::SequenceCollection;
    use crate::seqs::SequenceAccesors;

    use super::Filter;

    #[test]
    fn test_filter_seqcol_case_sentitive() {
        let mut sq = SequenceCollection::new();
        sq.add(
            AnnotatedSequence::from_string(
                "sequence_01".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "sequence_02".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "Sequence_03".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        let result = sq.filter_regex_id("sequence", false, true);
        assert!(result.size() == 2);
        assert!(result.iter().all(|x| x.id().starts_with("seq")))
    }
    #[test]
    fn test_filter_seqcol_case_insentitive() {
        let mut sq = SequenceCollection::new();
        sq.add(
            AnnotatedSequence::from_string(
                "sEquence_01".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "Cequence_02".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "Sequence_03".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        let result = sq.filter_regex_id("sequence", true, true);
        assert!(result.size() == 2);
    }
    #[test]
    fn test_filter_seqcol_case_insentitive_exclude() {
        let mut sq = SequenceCollection::new();
        sq.add(
            AnnotatedSequence::from_string(
                "sEquence_01".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "Cequence_02".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        sq.add(
            AnnotatedSequence::from_string(
                "Sequence_03".to_string(),
                "ACTATCGTCA".to_string()
            )
        ).unwrap();
        let result = sq.filter_regex_id("sequence", true, false);
        assert!(result.size() == 1);
    }
}
