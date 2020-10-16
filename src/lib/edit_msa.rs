use crate::edit::EditSequence;
use crate::seqs::Alignment;
use crate::seqs::SeqError;
use crate::seqs::SequenceAccesors;
use crate::seqs::AnnotatedSequence;

pub trait EditMSA {
    fn insert_columns(
        &mut self,
        at: usize,
        ncols: usize,
        ch: char,
    ) -> Result<(), SeqError>;
    fn insert_rows(
        &mut self,
        at: usize,
        rownames: Vec<&str>,
        ch: char,
    ) -> Result<(), SeqError>;
    fn insert_content(
        &mut self,
        at: usize,
        content: Vec<Vec<char>>,
        default_ch: char,
    ) -> Result<(), SeqError>;

    fn remove_columns(&mut self, at: usize, n: usize) -> Result<(), SeqError>;
    fn remove_rows(&mut self, at: usize, n: usize) -> Result<(), SeqError>;

    fn replace_content(
        &mut self,
        at: usize,
        content: Vec<Vec<char>>,
    ) -> Result<(), SeqError>;
}

impl EditMSA for Alignment {
    fn insert_columns(
        &mut self,
        at: usize,
        ncols: usize,
        ch: char,
    ) -> Result<(), SeqError> {
        for i in 0..self.size() {
            let s = self.get_mut(i).unwrap();
            s.edit_insert(vec![ch; ncols], at)?;
        }
        self.length = self.length.map(|x| x + ncols);
        Ok(())
    }
    fn insert_rows(
        &mut self,
        at: usize,
        rownames: Vec<&str>,
        ch: char,
    ) -> std::result::Result<(), SeqError> {
        
        if at <= self.size()  {
            for name in rownames {
                let seq =  AnnotatedSequence::new(
                    format!("{}", name),
                    vec![ch; self.length()]
                );
                self.insert(at, seq)?
            }
            Ok(())
        } else {
            Err(SeqError::EditError)
        }
    }
    fn insert_content(
        &mut self,
        _: usize,
        _: std::vec::Vec<std::vec::Vec<char>>,
        _: char,
    ) -> std::result::Result<(), SeqError> {
        todo!()
    }
    fn remove_columns(
        &mut self,
        _: usize,
        _: usize,
    ) -> std::result::Result<(), SeqError> {
        todo!()
    }
    fn remove_rows(
        &mut self,
        _: usize,
        _: usize,
    ) -> std::result::Result<(), SeqError> {
        todo!()
    }
    fn replace_content(
        &mut self,
        _: usize,
        _: std::vec::Vec<std::vec::Vec<char>>,
    ) -> std::result::Result<(), SeqError> {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::seqs::AnnotatedSequence;
    use crate::seqs::SequenceAccesors;
    // use crate::seqs::SequenceCollection;

    fn sample_msa() -> Alignment {
        let mut msa = Alignment::new();
        let s1 = AnnotatedSequence::from_string(
            String::from("s1"),
            String::from("ACTG"),
        );
        let s2 = AnnotatedSequence::from_string(
            String::from("s2"),
            String::from("CCTG"),
        );
        let s3 = AnnotatedSequence::from_string(
            String::from("s3"),
            String::from("ACAG"),
        );
        msa.add(s1).unwrap();
        msa.add(s2).unwrap();
        msa.add(s3).unwrap();
        msa
    }
    fn compare_sequence(msa: &Alignment, id: &str, seq: &str) {
        assert_eq!(msa.get_by_id(id).unwrap().seq_as_string(), seq)
    }
    fn check_internal_length_consistency(msa: &Alignment) {
        assert_eq!(msa.get(0).map(|x| x.len()).unwrap_or(0), msa.length())
    }

    #[test]
    fn insert_zero_columns_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        msa.insert_columns(0, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        msa.insert_columns(1, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        msa.insert_columns(4, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, 0, '-').expect_err("should fail");
    }
    #[test]
    fn insert_one_column_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_columns(0, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "-ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_columns(1, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "A-CTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_columns(4, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "ACTG-");
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, 1, '-').expect_err("Should fails");
    }
    #[test]
    fn insert_many_column_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_columns(0, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "---ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_columns(1, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "A---CTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_columns(4, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "ACTG---");
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, 3, '-').expect_err("Should fails");
    }

    #[test]
    fn insert_zero_rows_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_rows(0, 0, '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        assert!(msa.insert_rows(2, 0, '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        assert!(msa.insert_rows(4, 0, '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // Insert passing end
        msa.insert_rows(5, 0, '-').expect_err("Should fail");
    }
    #[test]
    fn insert_one_rows_at_any_position() {
        unimplemented!();
    }
    #[test]
    fn insert_many_rows_at_any_position() {
        unimplemented!();
    }
}
