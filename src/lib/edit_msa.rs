use crate::seqs::SeqError;

pub trait EditMSA {
    fn insert_columns(
        &mut self,
        at: usize,
        n: usize,
        ch: char,
    ) -> Result<(), SeqError>;
    fn insert_rows(
        &mut self,
        at: usize,
        n: usize,
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
