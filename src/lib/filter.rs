use crate::seqs::Alignment;
use crate::seqs::SequenceAccesors;
use crate::seqs::SequenceCollection;
use regex::Regex;
use regex::RegexBuilder;

// i have this code in rust, 
//  but it does not work, because Self cannot be known at compile time.
//  How i can solve this problem

pub trait Filter<T> {
  fn filter_regex_id(
    &self,
    regex: &str
 ) -> T;
  fn filter_regex_id_ignore_case(
    &self,
    regex: &str
 ) -> T;
}

impl Filter<SequenceCollection> for SequenceCollection {
    fn filter_regex_id(
        &self,
        regex: &str
     ) -> SequenceCollection {
        let re = Regex::new(regex).unwrap();
        let mut result = SequenceCollection::new();
        self.iter()
            .for_each(
                |x| match re.is_match(x.id()) {
                    true => result.add(x.clone()).unwrap(),
                    false => ()
                }
            );
        return result
    }

    fn filter_regex_id_ignore_case(
        &self,
        regex: &str
     ) -> SequenceCollection {
        let re = RegexBuilder::new(regex)
            .case_insensitive(true)
            .build()
            .unwrap();
        let mut result = SequenceCollection::new();
        self.iter()
            .for_each(
                |x| match re.is_match(x.id()) {
                    true => result.add(x.clone()).unwrap(),
                    false => ()
                }
            );
        return result
    }
}

impl Filter<Alignment> for Alignment {
    fn filter_regex_id(
        &self,
        regex: &str
     ) -> Alignment {
        let re = Regex::new(regex).unwrap();
        let mut result = Alignment::new();
        self.iter()
            .for_each(
                |x| match re.is_match(x.id()) {
                    true => result.add(x.clone()).unwrap(),
                    false => ()
                }
            );
        return result
    }

    fn filter_regex_id_ignore_case(
        &self,
        regex: &str
     ) -> Alignment {
        let re = RegexBuilder::new(regex)
            .case_insensitive(true)
            .build()
            .unwrap();
        let mut result = Alignment::new();
        self.iter()
            .for_each(
                |x| match re.is_match(x.id()) {
                    true => result.add(x.clone()).unwrap(),
                    false => ()
                }
            );
        return result
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
        let result = sq.filter_regex_id("sequence");
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
        let result = sq.filter_regex_id_ignore_case("sequence");
        assert!(result.size() == 2);
    }
}
