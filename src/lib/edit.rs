use crate::seqs::{AnnotatedSequence, SeqError};
use std::cmp::min;

pub trait EditSequence {
    fn edit_insert(
        &mut self,
        new: Vec<char>,
        at: usize,
    ) -> Result<(), SeqError>;

    fn edit_replace(
        &mut self,
        new: Vec<char>,
        at: usize,
        count: usize,
    ) -> Result<(), SeqError>;

    fn edit_delete(&mut self, at: usize, count: usize) -> Result<(), SeqError>;
}

impl EditSequence for AnnotatedSequence {
    fn edit_insert(
        &mut self,
        new: Vec<char>,
        at: usize,
    ) -> Result<(), SeqError> {
        match self.seq_mut() {
            Some(x) => {
                if at <= x.len() {
                    x.splice(at..at, new);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            }
            None => Err(SeqError::Empty),
        }
    }

    fn edit_replace(
        &mut self,
        new: Vec<char>,
        at: usize,
        count: usize,
    ) -> Result<(), SeqError> {
        match self.seq_mut() {
            Some(x) => {
                if at + count <= x.len() {
                    x.splice(at..min(at + count, x.len()), new);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            }
            None => Err(SeqError::Empty),
        }
    }

    fn edit_delete(&mut self, at: usize, count: usize) -> Result<(), SeqError> {
        match self.seq_mut() {
            Some(x) => {
                if at < x.len() && (x.len() - at) >= count {
                    x.splice(at..min(at + count, x.len()), vec![]);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            }
            None => Err(SeqError::Empty),
        }
    }
}

#[cfg(test)]
mod test {
    use super::EditSequence;
    use super::{AnnotatedSequence, SeqError};
    fn test_sequence() -> AnnotatedSequence {
        AnnotatedSequence::from_string(
            String::from("seq1"),
            String::from("ACTG"),
        )
    }
    fn compare_sequence(seq: &AnnotatedSequence, seq_str: &str) {
        assert_eq!(&seq.seq_as_string(), seq_str);
    }
    fn should_fail_with_edit_error(r: Result<(), SeqError>) {
        match r {
            Ok(()) => panic!("This should throw SeqError::EditError"),
            Err(x) => match x {
                SeqError::EditError => {}
                _ => panic!("This should throw SeqError::EditError"),
            },
        }
    }

    #[test]
    fn insert_one_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X'], 0).unwrap();
        compare_sequence(&s1, "XACTG");
    }

    #[test]
    fn insert_one_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X'], 2).unwrap();
        compare_sequence(&s1, "ACXTG");
    }

    #[test]
    fn insert_one_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X'], 4).unwrap();
        compare_sequence(&s1, "ACTGX");
    }
    #[test]
    fn insert_one_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_insert(vec!['X'], 5));
    }

    #[test]
    fn insert_many_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X', 'Y', 'Z'], 0).unwrap();
        compare_sequence(&s1, "XYZACTG");
    }

    #[test]
    fn insert_many_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X', 'Y', 'Z'], 2).unwrap();
        compare_sequence(&s1, "ACXYZTG");
    }

    #[test]
    fn insert_many_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X', 'Y', 'Z'], 4).unwrap();
        compare_sequence(&s1, "ACTGXYZ");
    }
    #[test]
    fn insert_many_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_insert(vec!['X', 'Y', 'Z'], 5));
    }
    #[test]
    fn insert_many_char_twice_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_insert(vec!['X', 'Y', 'Z'], 4).unwrap();
        s1.edit_insert(vec!['H', 'I', 'J', 'K'], 4).unwrap();
        compare_sequence(&s1, "ACTGHIJKXYZ");
    }
    #[test]
    fn delete_zero_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_delete(0, 0).unwrap();
        compare_sequence(&s1, "ACTG");
    }

    #[test]
    fn delete_zero_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_delete(2, 0).unwrap();
        compare_sequence(&s1, "ACTG");
    }

    #[test]
    fn delete_zero_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_delete(3, 0).unwrap();
        compare_sequence(&s1, "ACTG");
    }
    #[test]
    fn delete_zero_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_delete(4, 0));
    }

    #[test]
    fn delete_one_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_delete(0, 1).unwrap();
        compare_sequence(&s1, "CTG");
    }

    #[test]
    fn delete_one_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_delete(2, 1).unwrap();
        compare_sequence(&s1, "ACG");
    }

    #[test]
    fn delete_one_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_delete(3, 1).unwrap();
        compare_sequence(&s1, "ACT");
    }
    #[test]
    fn delete_one_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_delete(4, 0));
    }

    #[test]
    fn delete_many_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_delete(0, 3).unwrap();
        compare_sequence(&s1, "G");
    }

    #[test]
    fn delete_many_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_delete(2, 2).unwrap();
        compare_sequence(&s1, "AC");
    }

    #[test]
    fn delete_many_char_at_end_without_passing_end() {
        let mut s1 = test_sequence();
        s1.edit_delete(3, 1).unwrap();
        compare_sequence(&s1, "ACT");
    }
    #[test]
    fn delete_many_char_at_end_passing_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_delete(3, 2));
    }
    #[test]
    fn delete_many_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_delete(4, 2));
    }
    #[test]
    fn delete_many_char_twice_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_delete(2, 1).unwrap();
        s1.edit_delete(0, 1).unwrap();
        compare_sequence(&s1, "CG");
    }

    #[test]
    fn replace_zero_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 0, 0).unwrap();
        compare_sequence(&s1, "ACTG");
        let mut s2 = test_sequence();
        s2.edit_replace(vec![], 0, 1).unwrap();
        compare_sequence(&s2, "CTG");
    }

    #[test]
    fn replace_zero_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 2, 0).unwrap();
        compare_sequence(&s1, "ACTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 2, 1).unwrap();
        compare_sequence(&s1, "ACG");
    }

    #[test]
    fn replace_zero_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 3, 0).unwrap();
        compare_sequence(&s1, "ACTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 3, 1).unwrap();
        compare_sequence(&s1, "ACT");
    }
    #[test]
    fn replace_zero_char_at_after_end() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec![], 4, 0).unwrap();
        compare_sequence(&s1, "ACTG");

        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec![], 4, 1));
    }

    #[test]
    fn replace_one_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 0, 0).unwrap();
        compare_sequence(&s1, "XACTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 0, 1).unwrap();
        compare_sequence(&s1, "XCTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 0, 4).unwrap();
        compare_sequence(&s1, "X");
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X'], 0, 5));
    }

    #[test]
    fn replace_one_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 2, 0).unwrap();
        compare_sequence(&s1, "ACXTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 2, 1).unwrap();
        compare_sequence(&s1, "ACXG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 2, 2).unwrap();
        compare_sequence(&s1, "ACX");
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X'], 2, 3));
    }
    #[test]
    fn replace_one_char_at_end() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 3, 1).unwrap();
        compare_sequence(&s1, "ACTX");
    }
    #[test]
    fn replace_one_char_at_after_end() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X'], 4, 0).unwrap();
        compare_sequence(&s1, "ACTGX");
    }

    #[test]
    fn replace_many_char_at_begging() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 0, 0).unwrap();
        compare_sequence(&s1, "XYZACTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 0, 2).unwrap();
        compare_sequence(&s1, "XYZTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 0, 4).unwrap();
        compare_sequence(&s1, "XYZ");
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X', 'Y', 'Z'], 0, 5));
    }

    #[test]
    fn replace_many_char_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 2, 0).unwrap();
        compare_sequence(&s1, "ACXYZTG");
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 2, 2).unwrap();
        compare_sequence(&s1, "ACXYZ");
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X', 'Y', 'Z'], 2, 4));
    }

    #[test]
    fn replace_many_char_at_end_without_passing_end() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 4, 0).unwrap();
        compare_sequence(&s1, "ACTGXYZ");
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X', 'Y', 'Z'], 4, 1));
    }
    #[test]
    fn replace_many_char_at_end_passing_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X', 'Y', 'Z'], 4, 1));
    }
    #[test]
    fn replace_many_char_at_after_end() {
        let mut s1 = test_sequence();
        should_fail_with_edit_error(s1.edit_replace(vec!['X', 'Y', 'Z'], 5, 0));
    }
    #[test]
    fn replace_many_char_twice_at_middle() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 2, 1).unwrap();
        s1.edit_replace(vec!['W', 'Q', 'R'], 3, 2).unwrap();
        compare_sequence(&s1, "ACXWQRG");
    }

    #[test]
    fn all_edit_operations() {
        let mut s1 = test_sequence();
        s1.edit_replace(vec!['X', 'Y', 'Z'], 2, 1).unwrap();
        s1.edit_insert(vec!['Q', 'W'], 0).unwrap();
        s1.edit_delete(1, 3).unwrap();
        compare_sequence(&s1, "QXYZG");
    }
}
