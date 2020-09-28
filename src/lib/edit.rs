use std::cmp::min;
use crate::seqs::{
    SequenceAccesors, 
    SequenceCollection, 
    SeqError,
    AnnotatedSequence};

pub trait EditSequence {
    fn edit_insert(
        &mut self,
        new: Vec<char>,
        at: usize
    ) -> Result<(), SeqError>;
    
    fn edit_replace(
        &mut self, 
        new: Vec<char>, 
        at: usize
    ) -> Result<(), SeqError>;

    fn edit_delete(
        &mut self,
        at: usize,
        count: usize
    ) -> Result<(), SeqError>;
}

impl EditSequence for AnnotatedSequence {
    fn edit_insert(
        &mut self,
        new: Vec<char>,
        at: usize
    ) -> Result<(), SeqError> {
        match self.seq_mut() {
            Some(x) => {
                if at <= x.len()  {
                    x.splice(at..at, new);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            },
            None => Err(SeqError::Empty)
        }
    }

    fn edit_replace(
        &mut self, 
        new: Vec<char>, 
        at: usize
    ) -> Result<(), SeqError> {
        match self.seq_mut(){
            Some(x) => {
                if at < x.len()  {
                    x.splice(at..min(at+new.len(), x.len()), new);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            },
            None => Err(SeqError::Empty)
        }
    }

    fn edit_delete(
        &mut self,
        at: usize,
        count: usize
    ) -> Result<(), SeqError> {
        match self.seq_mut(){
            Some(x) => {
                if at < x.len() {
                    x.splice(at..min(at+count, x.len()), vec![]);
                    Ok(())
                } else {
                    Err(SeqError::EditError)
                }
            },
            None => Err(SeqError::Empty)
        }
    }
}

mod test {
    use super::*;
    #[test]
    fn insert_one_char_at_begging() {
        let mut s1 = AnnotatedSequence::from_string(
            String::from("seq1"),
            String::from("ACTG")
        );
        s1.edit_insert(vec!['X'], 0).unwrap();
        assert_eq!(
            s1.seq_as_string(), String::from("XACTG"));
    }

    #[test]
    fn insert_one_char_at_middle() {
        let mut s1 = AnnotatedSequence::from_string(
            String::from("seq1"),
            String::from("ACTG")
        );
        s1.edit_insert(vec!['X'], 2).unwrap();

        assert_eq!(
            s1.seq_as_string(), String::from("ACXTG"));
    }
    
    #[test]
    fn insert_one_char_at_end() {
        let mut s1 = AnnotatedSequence::from_string(
            String::from("seq1"),
            String::from("ACTG")
        );
        s1.edit_insert(vec!['X'], 4).unwrap();
        
        assert_eq!(
            s1.seq_as_string(), String::from("ACTGX"));
    }
    
    #[test]
    fn insert_one_char_at_after_end() {
        let mut s1 = AnnotatedSequence::from_string(
            String::from("seq1"),
            String::from("ACTG")
        );
        match s1.edit_insert(vec!['X'], 5) {
            Ok(()) => panic!("This should throw SeqError::EditError"),
            Err(x) => match x {
                SeqError::EditError => {},
                _ => panic!("This should throw SeqError::EditError")
            }
        }
    }

}