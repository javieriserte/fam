use crate::edit::EditSequence;
use crate::seqs::Alignment;
use crate::seqs::SeqError;
use crate::seqs::SequenceAccesors;
use crate::seqs::AnnotatedSequence;

pub trait EditMSA {
    fn insert_empty_columns(
        &mut self,
        at: usize,
        ncols: usize,
        ch: char,
    ) -> Result<(), SeqError>;
    fn insert_empty_rows(
        &mut self,
        at: usize,
        rownames: Vec<&str>,
        ch: char,
    ) -> Result<(), SeqError>;
    fn insert_columns(
        &mut self,
        at: usize,
        content: Vec<Vec<char>>
    ) -> Result<(), SeqError>;
    fn insert_rows(
        &mut self,
        at: usize,
        rownames: Vec<&str>,
        content: Vec<Vec<char>>
    ) -> Result<(), SeqError>;

    fn remove_columns(&mut self, positions: Vec<usize>) -> Result<(), SeqError>;
    fn remove_rows(&mut self, positions: Vec<usize>) -> Result<(), SeqError>;

    fn replace_content(
        &mut self,
        at: usize,
        content: Vec<Vec<char>>,
    ) -> Result<(), SeqError>;

    fn reorder_rows(
        &mut self,
        order: Vec<usize>
    ) -> Result<(), SeqError>;
}

// Transpose a vector of vectors.
// All vectors in input shold have the same length and be non emtpy.
fn transpose_vector<T: Copy>(input: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    let w = input.len();
    let h = input.get(0)
        .and_then(|x| Some(x.len()))
        .unwrap_or(0);
    let all_equals = input
        .windows(2)
        .all(|x| x[0].len()==x[1].len());
    if w > 0 {
        if !all_equals {
            panic!("Columns of different length");
        }
        let mut r = vec![vec![]; h];
        (0..w).for_each(
            |x| (0..h).for_each(
                |y| r[y].push(input[x][y])
            )
        );
        r
    } else {
        vec![]
    }
}

impl EditMSA for Alignment {
    fn insert_empty_columns(
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
    fn insert_empty_rows(
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
    fn insert_columns(
        &mut self,
        at: usize,
        content: Vec<Vec<char>>
    ) -> Result<(), SeqError> {
        if content.iter().all(|x| x.len() == self.size()) && 
            at <= self.length() {
            let tcontent = transpose_vector(&content);
            for (x, y) in (0..self.size()).zip(tcontent) {
                self.get_mut(x)
                .unwrap()
                .edit_insert(y, at)?
            }
            self.length = self.length.map(|x| x + content.len()) ;
            Ok(())
        } else {
            Err(SeqError::EditError)
        }
    }
    fn insert_rows(
        &mut self,
        at: usize,
        rownames: Vec<&str>,
        content: Vec<Vec<char>>
    ) -> std::result::Result<(), SeqError> {
        if at <= self.size() {
            for (name, cont) in rownames.iter().zip(content) {
                let cseq = AnnotatedSequence::new(
                    String::from(*name),
                    cont
                );
                self.seqs.insert(at, cseq)?
            }
            Ok(())
        } else {
            Err(SeqError::EditError)
        }
    }
    fn remove_columns(
        &mut self,
        positions: Vec<usize>
    ) -> std::result::Result<(), SeqError> {
        let mut sorted = positions.clone();
        sorted.sort_unstable();
        sorted.reverse();
        for row_index in 0..self.size() {
            let row = self.get_mut(row_index).unwrap();
            for column_index in sorted.iter() {
                row.edit_delete(*column_index, 1)?
            }
        }
        self.length = self.length.map(|x| x - sorted.len());
        Ok(())
    }
    fn remove_rows(
        &mut self,
        positions: Vec<usize>
    ) -> std::result::Result<(), SeqError> {
        let mut sorted = positions;
        sorted.sort_unstable();
        sorted.reverse();
        Ok(
            for row_index in sorted {
                self.seqs.remove(row_index).ok_or_else(|| SeqError::EditError)?;
            }
        )
    }
    fn replace_content(
        &mut self,
        _: usize,
        _: std::vec::Vec<std::vec::Vec<char>>,
    ) -> std::result::Result<(), SeqError> {
        todo!()
    }
    fn reorder_rows(&mut self, _: Vec<usize>) -> Result<(), SeqError> {
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
    fn sample_large() -> Alignment {
        let mut msa = Alignment::new();
        let s1 = AnnotatedSequence::from_string(
            String::from("s1"), String::from("ACTG"),
        );
        let s2 = AnnotatedSequence::from_string(
            String::from("s2"), String::from("CCTG"),
        );
        let s3 = AnnotatedSequence::from_string(
            String::from("s3"), String::from("ACAG"),
        );
        let s4 = AnnotatedSequence::from_string(
            String::from("s4"), String::from("ACAG"),
        );
        let s5 = AnnotatedSequence::from_string(
            String::from("s5"), String::from("ACAG"),
        );
        msa.add(s1).unwrap();
        msa.add(s2).unwrap();
        msa.add(s3).unwrap();
        msa.add(s4).unwrap();
        msa.add(s5).unwrap();
        msa
    }
    fn compare_sequence(msa: &Alignment, id: &str, seq: &str) {
        assert_eq!(msa.get_by_id(id).unwrap().seq_as_string(), seq)
    }

    fn compare_sequence_at(msa: &Alignment, index: usize, seq: &str) {
        assert_eq!(msa.get(index).unwrap().seq_as_string(), seq)
    }
    fn check_internal_length_consistency(msa: &Alignment) {
        assert_eq!(msa.get(0).map(|x| x.len()).unwrap_or(0), msa.length())
    }

    #[test]
    fn insert_empty_zero_columns_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        msa.insert_empty_columns(0, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        msa.insert_empty_columns(1, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        msa.insert_empty_columns(4, 0, '-').unwrap();
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_empty_columns(5, 0, '-').expect_err("should fail");
    }
    #[test]
    fn insert_empty_one_column_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(0, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "-ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(1, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "A-CTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(4, 1, '-').is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "ACTG-");
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_empty_columns(5, 1, '-').expect_err("Should fails");
    }
    #[test]
    fn insert_empty_many_column_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(0, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "---ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(1, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "A---CTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_empty_columns(4, 3, '-').is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "ACTG---");
        check_internal_length_consistency(&msa);
        // Insert passing end
        let mut msa = sample_msa();
        msa.insert_empty_columns(5, 3, '-').expect_err("Should fails");
    }

    #[test]
    fn insert_empty_zero_rows_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(0, vec![], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        assert!(msa.insert_empty_rows(2, vec![], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // // Insert at end
        assert!(msa.insert_empty_rows(3, vec![], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        check_internal_length_consistency(&msa);
        // // Insert passing end
        msa.insert_empty_rows(4, vec![], '-').expect_err("Should fail");
    }
    #[test]
    fn insert_empty_one_rows_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(0, vec!["ns1"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence_at(&msa, 0, "----");
        compare_sequence_at(&msa, 1, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(1, vec!["ns1"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence_at(&msa, 1, "----");
        compare_sequence_at(&msa, 0, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(3, vec!["ns1"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence_at(&msa, 3, "----");
        compare_sequence_at(&msa, 0, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert passing end
        msa.insert_empty_rows(4, vec!["ns1"], '-').expect_err("Should fail");
        
    }
    #[test]
    fn insert_empty_many_rows_at_any_position() {
        // Insert at beggining
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(0, vec!["ns1", "ns2"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence(&msa, "ns2", "----");
        compare_sequence_at(&msa, 0, "----");
        compare_sequence_at(&msa, 1, "----");
        compare_sequence_at(&msa, 2, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at middle
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(1, vec!["ns1", "ns2"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence(&msa, "ns2", "----");
        compare_sequence_at(&msa, 1, "----");
        compare_sequence_at(&msa, 2, "----");
        compare_sequence_at(&msa, 0, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert at end
        let mut msa = sample_msa();
        assert!(msa.insert_empty_rows(3, vec!["ns1", "ns2"], '-').is_ok());
        assert_eq!(msa.length(), 4);
        compare_sequence(&msa, "s1", "ACTG");
        compare_sequence(&msa, "ns1", "----");
        compare_sequence(&msa, "ns2", "----");
        compare_sequence_at(&msa, 3, "----");
        compare_sequence_at(&msa, 4, "----");
        compare_sequence_at(&msa, 0, "ACTG");
        check_internal_length_consistency(&msa);
        // Insert passing end
        msa.insert_empty_rows(4, vec!["ns1", "ns2"], '-')
            .expect_err("Should fail");
    }
    #[test]
    fn insert_empty_many_rows_fails_with_repeated_names() {
        // New is equals to existent
        let mut msa = sample_msa();
        msa.insert_empty_rows(1, vec!["s1", "ns2"], '-')
            .expect_err("Should fail");
        // New has to equal names
        msa.insert_empty_rows(1, vec!["ns2", "ns2", "ns1"], '-')
            .expect_err("Should fail");
    }
    #[test]
    fn insert_zero_columns() {
        // At beggining
        let mut msa = sample_msa();
        assert!(msa.insert_columns(0, vec![]).is_ok());
        println!("SEQ: [{}]", msa.get(0).unwrap().seq_as_string());
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // At middle
        let mut msa = sample_msa();
        assert!(msa.insert_columns(1, vec![]).is_ok());
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // At end
        let mut msa = sample_msa();
        assert!(msa.insert_columns(4, vec![]).is_ok());
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
        // Passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, vec![]).expect_err("Should fail");
    }
    #[test]
    fn insert_one_column() {
        // At beggining
        let c1 = "AAA".chars().collect::<Vec<char>>();
        let mut msa = sample_msa();
        assert!(msa.insert_columns(0, vec![c1.clone()]).is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "AACTG");
        check_internal_length_consistency(&msa);
        // At middle
        let mut msa = sample_msa();
        assert!(msa.insert_columns(2, vec![c1.clone()]).is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "ACATG");
        check_internal_length_consistency(&msa);
        // At end
        let mut msa = sample_msa();
        assert!(msa.insert_columns(4, vec![c1.clone()]).is_ok());
        assert_eq!(msa.length(), 5);
        compare_sequence(&msa, "s1", "ACTGA");
        check_internal_length_consistency(&msa);
        // Passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, vec![c1.clone()]).expect_err("Should fail");
    }
    #[test]
    fn insert_columns_with_wrong_size() {
        let c1 = "AA".chars().collect::<Vec<char>>();
        let c2 = "AAAA".chars().collect::<Vec<char>>();
        let mut msa = sample_msa();
        msa.insert_columns(0, vec![c1]).expect_err("Should fail");
        msa.insert_columns(0, vec![c2]).expect_err("Should fail");
    }
    #[test]
    fn insert_many_columns() {
        // At beggining
        let c1 = "AAA".chars().collect::<Vec<char>>();
        let c2 = "EEE".chars().collect::<Vec<char>>();
        let c3 = "III".chars().collect::<Vec<char>>();
        let mut msa = sample_msa();
        assert!(msa.insert_columns(0,
            vec![c1.clone(), c2.clone(), c3.clone()]).is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "AEIACTG");
        check_internal_length_consistency(&msa);
        // At middle
        let mut msa = sample_msa();
        assert!(msa.insert_columns(2,
            vec![c1.clone(), c2.clone(), c3.clone()]).is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "ACAEITG");
        check_internal_length_consistency(&msa);
        // At end
        let mut msa = sample_msa();
        assert!(msa.insert_columns(4,
            vec![c1.clone(), c2.clone(), c3.clone()]).is_ok());
        assert_eq!(msa.length(), 7);
        compare_sequence(&msa, "s1", "ACTGAEI");
        // Passing end
        let mut msa = sample_msa();
        msa.insert_columns(5, vec![
            c1.clone(), c2.clone(), c3.clone()]).expect_err("Should fail");
    }    

    #[test]
    fn transpose_vector_test() {
        let a = vec![
            vec![1,2,3],
            vec![4,5,6]
        ];
        let b = transpose_vector(&a);
        assert_eq!(b.len(), 3);
        assert_eq!(b[0].len(), 2);
        assert_eq!(b[0], vec![1, 4]);
        assert_eq!(b[1], vec![2, 5]);
        assert_eq!(b[2], vec![3, 6]);
    }
    #[test]
    fn remove_zero_columns() {
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![]).is_ok());
        assert_eq!(msa.length(), 4);
        check_internal_length_consistency(&msa);
    }
    #[test]
    fn remove_one_columns() {
        // At beggining
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![0]).is_ok());
        assert_eq!(msa.length(), 3);
        check_internal_length_consistency(&msa);
        // At middle
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![2]).is_ok());
        assert_eq!(msa.length(), 3);
        check_internal_length_consistency(&msa);

        // At End
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![3]).is_ok());
        assert_eq!(msa.length(), 3);
        check_internal_length_consistency(&msa);

        // Passing end
        let mut msa = sample_msa();
        msa.remove_columns(vec![4]).expect_err("Should fail");
    }

    #[test]
    fn remove_many_columns() {
        // two columns sorted
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![0, 2]).is_ok());
        assert_eq!(msa.length(), 2);
        check_internal_length_consistency(&msa);
        compare_sequence(&msa, "s1", "CG");
        // two columns unsorted
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![2, 0]).is_ok());
        assert_eq!(msa.length(), 2);
        check_internal_length_consistency(&msa);
        compare_sequence(&msa, "s1", "CG");
        // three columns unsorted, including end
        let mut msa = sample_msa();
        assert!(msa.remove_columns(vec![3, 0, 1]).is_ok());
        assert_eq!(msa.length(), 1);
        check_internal_length_consistency(&msa);
        compare_sequence(&msa, "s1", "T");
        // two columns, passing end
        let mut msa = sample_msa();
        msa.remove_columns(vec![4, 0, 1]).expect_err("Should Fail");
    }
    #[test]
    fn remove_zero_rows() {
        let mut msa = sample_msa();
        let rows_to_remove = vec![];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 3);
    }
    #[test]
    fn remove_one_row() {
        // At beggining
        let mut msa = sample_large();
        let rows_to_remove = vec![0];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 4);

        // At middle
        let mut msa = sample_large();
        let rows_to_remove = vec![2];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 4);
        
        // At end
        let mut msa = sample_large();
        let rows_to_remove = vec![4];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 4);
        
        // Passing end
        let mut msa = sample_large();
        let rows_to_remove = vec![5];
        msa.remove_rows(rows_to_remove).expect_err("Should fail");
        assert_eq!(msa.size(), 5);
    }
    #[test]
    fn remove_many_rows() {
        let mut msa = sample_large();
        let rows_to_remove = vec![0, 1];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 3);
        assert_eq!(msa.get(0).unwrap().id(), "s3" );
        assert_eq!(msa.get(1).unwrap().id(), "s4" );
        assert_eq!(msa.get(2).unwrap().id(), "s5" );

        let mut msa = sample_large();
        let rows_to_remove = vec![0, 2];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 3);
        assert_eq!(msa.get(0).unwrap().id(), "s2" );
        assert_eq!(msa.get(1).unwrap().id(), "s4" );
        assert_eq!(msa.get(2).unwrap().id(), "s5" );
        
        let mut msa = sample_large();
        let rows_to_remove = vec![0, 4];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 3);
        assert_eq!(msa.get(0).unwrap().id(), "s2" );
        assert_eq!(msa.get(1).unwrap().id(), "s3" );
        assert_eq!(msa.get(2).unwrap().id(), "s4" );

        let mut msa = sample_large();
        let rows_to_remove = vec![3,2,1];
        assert!(msa.remove_rows(rows_to_remove).is_ok());
        assert_eq!(msa.size(), 2);
        assert_eq!(msa.get(0).unwrap().id(), "s1" );
        assert_eq!(msa.get(1).unwrap().id(), "s5" );

    }
    #[test]
    fn remove_rows_many_times() {
        let mut msa = sample_large();
        let rows_to_remove_first = vec![0,3];
        let rows_to_remove_second = vec![0,2];
        assert!(msa.remove_rows(rows_to_remove_first).is_ok());
        assert!(msa.remove_rows(rows_to_remove_second).is_ok());
        assert_eq!(msa.size(), 1);
        assert_eq!(msa.get(0).unwrap().id(), "s3");
    }
}
