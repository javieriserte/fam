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
                if at < x.len()  {
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
    #[test]
    fn name() {
        unimplemented!();
    }
}