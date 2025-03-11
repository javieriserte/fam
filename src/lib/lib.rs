pub mod edit;
pub mod edit_msa;
pub mod fastaio;
pub mod combine;
pub mod random;
pub mod conservation;
pub mod plotting;
pub mod clustering;
pub mod random_voronoi;
pub mod matrices;
pub mod filter;
pub mod degap;
pub mod gapping;
pub mod trim;

pub mod seqs {
    use std::{
        cell::RefCell,
        cmp::{max, min},
        collections::{HashMap, HashSet},
        io::{BufRead, ErrorKind}, iter::FromIterator
    };
    use std::fmt::{Display, Error, Formatter};
    use std::iter::{IntoIterator, Iterator};
    use std::collections::hash_map::Entry::{Vacant, Occupied};

    use crate::fastaio::{
        reader_for, InputFormats, SequenceReader
    };

    #[derive(Debug)]
    pub enum SeqError {
        DuplicatedId(String),
        DifferentLength,
        NonExistenId(String),
        MissingID(String),
        EditError,
        Empty,
    }

    impl Display for SeqError {
        fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
            match self {
                SeqError::DuplicatedId(x) => write!(
                    f,
                    "Attempted to add a sequence with a \
                    duplicated ID ['{}'].",
                    x
                ),
                SeqError::NonExistenId(x) => {
                    write!(f, "Attempted to get a non existing ID [{}].", x)
                }
                SeqError::DifferentLength => write!(
                    f,
                    "Attempted to add a sequence with \
                    different length to a MSA."
                ),
                SeqError::MissingID(x) => {
                    write!(f, "A sequence with ID [{}] is required", x)
                }
                SeqError::EditError => {
                    write!(f, "Attempted to edit a sequence beyond its length")
                }
                SeqError::Empty => {
                    write!(f, "Attempted to access an empty sequence")
                }
            }
        }
    }

    impl From<SeqError> for std::io::Error {
        fn from(x: SeqError) -> Self {
            std::io::Error::new(
                ErrorKind::Other,
                format!("{}.\n", x)
            )
        }
    }

    pub trait SequenceAccesors {
        fn get_mut(&mut self, index: usize) -> Option<&mut AnnotatedSequence>;
        // Retrieve de number of sequences (rows) in the collection.
        fn size(&self) -> usize;
        fn contains(&self, id: &str) -> bool;
        fn get(&self, index: usize) -> Option<&AnnotatedSequence>;
        fn get_by_id(&self, id: &str) -> Option<&AnnotatedSequence>;
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(), SeqError>;
        fn insert(
            &mut self,
            index: usize,
            seq: AnnotatedSequence,
        ) -> Result<(), SeqError>;
        fn remove(&mut self, index: usize) -> Option<AnnotatedSequence>;
        fn remove_by_id(&mut self, id: &str) -> Option<AnnotatedSequence>;
        fn move_up(&mut self, id: &str) -> Result<(), SeqError>;
        fn iter(&self) -> SequenceIterable;
        fn reorder(&mut self, order: Vec<usize>) -> Result<(), SeqError>;
    }

    /// Struct to represent a single sequence of a MSA
    #[derive(Clone, PartialEq, Debug, Eq)]
    pub struct AnnotatedSequence {
        id: String,
        sequence: Option<Vec<char>>,
    }

    impl AnnotatedSequence {
        /// Creates a new AnnotatedSequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::new(
        ///     String::from("S1"),
        ///     vec!['A', 'T', 'C', 'A', 'T', 'G',
        ///          'C', 'T', 'A', 'C', 'T', 'G']);
        /// assert_eq!(a.id(), "S1");
        /// assert_eq!(a.seq(), Some(
        ///     &vec!['A', 'T', 'C', 'A', 'T', 'G',
        ///           'C', 'T', 'A', 'C', 'T', 'G']));
        /// ```
        pub fn new(id: String, sequence: Vec<char>) -> Self {
            AnnotatedSequence {
                id,
                sequence: Some(sequence),
            }
        }

        pub fn from_string<T: ToString, U: ToString>(id: T, sequence: U) -> Self {
            AnnotatedSequence {
                id: id.to_string(),
                sequence: Some(
                    sequence.to_string().chars().collect::<Vec<char>>()
                ),
            }
        }

        /// Creates a new empty AnnotatedSequence, contains only id.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(a.id() , "S1");
        /// assert_eq!(a.seq() , None);
        /// ```
        pub fn empty(id: String) -> Self {
            AnnotatedSequence { id, sequence: None }
        }

        /// Set or change sequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATC")
        /// );
        /// assert_eq!(a.seq(), Some(&vec!['A', 'T', 'C']));
        /// a.set_sequence(vec!['G', 'G', 'G']);
        /// assert_eq!(a.seq(), Some(&vec!['G', 'G', 'G']));
        ///
        /// let mut b = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(b.seq(), None);
        /// b.set_sequence(vec!['G', 'G', 'G']);
        /// assert_eq!(b.seq(), Some(&vec!['G', 'G', 'G']));
        ///
        /// ```
        pub fn set_sequence(&mut self, seq: Vec<char>) {
            self.sequence = Some(seq);
        }

        pub fn set_sequence_as_string(&mut self, seq: String) {
            if seq.is_empty() {
                self.sequence = None;
            } else {
                self.sequence = Some(seq.chars().collect());
            }
        }

        /// Moves out the sequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATC")
        /// );
        /// assert_eq!(a.seq(), Some(&vec!['A', 'T', 'C']));
        /// let b = a.take_sequence();
        /// assert_eq!(a.seq() , None);
        /// assert_eq!(b , Some(vec!['A', 'T', 'C']));
        /// ```
        pub fn take_sequence(&mut self) -> Option<Vec<char>> {
            self.sequence.take()
        }

        /// Modify the id
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(a.id(), "S1");
        /// a.set_id(String::from("S2"));
        /// assert_eq!(a.id() , "S2");
        /// ```
        pub fn set_id(&mut self, id: String) {
            self.id = id;
        }

        /// Retrieves a reference of the sequence.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATCG")
        /// );
        /// assert_eq!(a.seq() , Some(&vec!['A', 'T', 'C', 'G']));
        /// ```
        pub fn seq(&self) -> Option<&Vec<char>> {
            match &self.sequence {
                None => None,
                Some(x) => Some(x),
            }
        }

        pub fn seq_mut(&mut self) -> Option<&mut Vec<char>> {
            match &mut self.sequence {
                None => None,
                Some(x) => Some(x),
            }
        }

        pub fn seq_as_string(&self) -> String {
            match &self.sequence {
                None => String::from(""),
                Some(x) => x.iter().collect::<String>(),
            }
        }

        /// Retrieves a reference of the sequence.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(a.len(), 0);
        /// a.set_sequence_as_string(String::from("ATCG"));
        /// assert_eq!(a.len(), 4);
        /// a.set_sequence_as_string(String::from("AT"));
        /// assert_eq!(a.len(), 2);
        /// ```
        pub fn len(&self) -> usize {
            match &self.sequence {
                None => 0,
                Some(x) => x.len(),
            }
        }

        /// Checks if the sequence is empty.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::empty(String::from("S1"));
        /// assert!(a.is_empty());
        /// ```
        pub fn is_empty(&self) -> bool {
            self.len() == 0
        }

        /// Gets a copy the sequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::empty(String::from("S1"));
        /// assert!(a.seq_copy().is_err());
        /// a.set_sequence(vec!['A', 'T', 'C', 'G']);
        /// assert_eq!(a.seq_copy().unwrap(), vec!['A', 'T', 'C', 'G']);
        /// assert_eq!(a.seq(), Some(&vec!['A', 'T', 'C', 'G']));
        /// ```
        pub fn seq_copy(&self) -> Result<Vec<char>, String> {
            match &self.sequence {
                None => Err(String::from("Sequence is empty")),
                Some(x) => Ok(x.clone()),
            }
        }

        /// Gets a reference of the id
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::new(String::from("S1"), vec!['A']);
        /// assert_eq!(a.id() , "S1");
        /// ```
        pub fn id(&self) -> &str {
            &self.id
        }

        pub fn remove_positions(&mut self, positions: &Vec<usize>) {
            let seq = self.take_sequence().unwrap();
            let mut new_seq = vec![];
            let to_remove = HashSet::<usize>::from_iter(
                positions.iter().cloned()
            );
            for (i, x) in seq.iter().enumerate() {
                if !to_remove.contains(&i) {
                    new_seq.push(*x);
                }
            }
            self.set_sequence(new_seq);
        }

        /// Trims a sequence with a fixed length from the left and right.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut seq = AnnotatedSequence::from_string(
        ///     "S1",
        ///     "ABCDEFGHIJ"
        /// );
        /// seq.trim_fixed(2, 3);
        /// assert_eq!(seq.seq_as_string(), "CDEFG");
        /// let mut seq = AnnotatedSequence::from_string(
        ///     "S1",
        ///     "ABCDEFGHIJ"
        /// );
        /// seq.trim_fixed(5, 5);
        /// assert_eq!(seq.seq_as_string(), "");
        /// let mut seq = AnnotatedSequence::from_string(
        ///     "S1",
        ///     "ABCDEFGHIJ"
        /// );
        /// seq.trim_fixed(6, 5);
        /// assert_eq!(seq.seq_as_string(), "");
        /// ```
        pub fn trim_fixed(&mut self, left: usize, right:usize) {
            let seq = self.take_sequence().unwrap();
            if left + right > seq.len() {
                let new_seq = vec![];
                self.set_sequence(new_seq);
            } else {
                let new_seq = (&seq[left .. seq.len() - right]).to_vec();
                self.set_sequence(new_seq);
            }
        }
    }

    #[derive(Clone, PartialEq, Eq, Debug)]
    pub struct SequenceCollection {
        sequences: Vec<AnnotatedSequence>,
        ids: HashMap<String, usize>, // Maps id to index in sequences vector
    }

    impl SequenceCollection {
        /// Creates a new SequenceCollection
        /// ```
        /// use famlib::seqs::SequenceCollection;
        /// let a = SequenceCollection::new();
        /// ```
        pub fn new() -> Self {
            SequenceCollection {
                sequences: vec![],
                ids: HashMap::new(),
            }
        }

        /// Mute the SequenceCollection to a Alignment
        pub fn to_msa(self) -> Result<Alignment, SeqError> {
            let mut msa = Alignment::new();
            if self.size() == 0 {
                msa.seqs = self;
                Ok(msa)
            } else {
                let ref_len = self.sequences[0].len();
                if self.sequences.iter().all(|x| x.len() == ref_len) {
                    msa.seqs = self;
                    msa.length = Some(ref_len);
                    Ok(msa)
                } else {
                    Err(SeqError::DifferentLength)
                }
            }
        }

        pub fn to_buffered(self) -> BufferedSeqCollectionFromSeqCol {
            BufferedSeqCollectionFromSeqCol::new(self)
        }
    }

    impl Default for SequenceCollection {
        fn default() -> Self {
            Self::new()
        }
    }

    impl Display for SequenceCollection {
        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
/* Output Example
         1         2         3         4         5         6         7         8
.........0.........0.........0.........0.........0.........0.........0.........0
> Seq1       ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ 50
> Seq2       ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ            40
> Seq3       ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ...DEFGHIJ 60
...
> Seq5       ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ...DEFGHIJ 60
> Sequence6  ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ...DEFGHIJ 60
> My_e...e_7 ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ABCDEFGHIJ ...DEFGHIJ 60
*/
            let col_width = 10;
            let head = 3;
            let max_printed_cols=5;
            let printed_id_length=10;
            let nseqs = self.size();
            let mut widths = vec![];
            let mut max_ncol = 0;
            for i in 0..nseqs {
                if i <head || i>=nseqs-3 {
                    let c_len = self.get(i).map_or_else(||0, |x| x.len());
                    let ncol = match c_len>0 {
                        true => (c_len-1)/col_width+1,
                        false => 0
                    };
                    widths.push((i, c_len, ncol));
                    max_ncol = max(max_ncol, ncol);
                }
            }
            max_ncol = min(max_ncol, max_printed_cols);
            for (row, len, cols) in widths {
                let s = self.get(row).unwrap();
                let mut id = String::from(s.id());
                if id.len()>printed_id_length {
                    let header = (printed_id_length-2)/2;
                    let tail = printed_id_length-header-3;
                    id = format!(
                        "{}...{}",
                        &id[0..header],
                        &id[id.len()-tail..id.len()]);
                }
                let mut seq = s.seq_as_string();
                if cols > max_ncol {
                    seq = format!(
                        "{}...{}",
                        &seq[0..(max_ncol-1)*col_width],
                        &seq[seq.len()-(col_width-3)..seq.len()]);
                } else {
                    let mut sb = seq.chars().collect::<Vec<char>>();
                    sb.resize(max_ncol*col_width, ' ');
                    seq = sb.into_iter().collect();
                }
                seq = seq.chars().enumerate().flat_map(
                    |(i, x)| if i>0 && (i+1)%10==1 {
                        vec![' ', x]
                    } else {
                        vec![x]
                    }
                ).collect();
                writeln!(
                    f,
                    "> {:<w$} {} {}", id, seq, len, w=printed_id_length
                )?;
                if row==2 && self.size()>6 {writeln!(f, "...")?;}
            }
            Ok(())
        }
    }

    impl SequenceAccesors for SequenceCollection {
        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(
        ///     String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::from_string(
        ///     String::from("S3"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.insert(0, c);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S3");
        /// ```
        fn insert(
            &mut self,
            index: usize,
            seq: AnnotatedSequence,
        ) -> Result<(), SeqError> {
            if self.ids.contains_key(&seq.id().to_string()) {
                Err(SeqError::DuplicatedId(String::from(seq.id())))
            } else {
                for x in self.ids.values_mut() {
                    if *x >= index {
                        *x += 1;
                    }
                }
                self.ids.insert(seq.id().to_string(), index);
                self.sequences.insert(index, seq);
                Ok(())
            }
        }

        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(
        ///     String::from("S2"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S1");
        /// assert_eq!(seqs.get(1).unwrap().id(), "S2");
        /// ```
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(), SeqError> {
            let last_index = self.ids.len();
            let entry = self.ids.entry(seq.id().to_string());
            match entry {
                Vacant(e) => {
                    e.insert(last_index);
                    self.sequences.push(seq);
                    Ok(())
                },
                Occupied(_) => {
                    Err(SeqError::DuplicatedId(seq.id().to_string()))
                }
            }
        }

        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// assert_eq!(seqs.size(), 0);
        /// seqs.add(a);
        /// assert_eq!(seqs.size(), 1);
        /// ```
        fn size(&self) -> usize {
            self.sequences.len()
        }

        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// assert_eq!(seqs.get(0).unwrap().id(), String::from("S1"));
        /// ```
        fn get(&self, index: usize) -> Option<&AnnotatedSequence> {
            self.sequences.get(index)
        }

        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// assert_eq!(seqs.get_by_id("S2"), None);
        /// assert_eq!(
        ///     seqs.get_by_id("S1").unwrap().seq_as_string(),
        ///     String::from("ATCATGCTACTG")
        /// );
        /// ```
        fn get_by_id(&self, id: &str) -> Option<&AnnotatedSequence> {
            match self.ids.get(id) {
                None => None,
                Some(x) => self.get(*x),
            }
        }

        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     AnnotatedSequence,
        ///     SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(
        ///     String::from("S2"),
        ///     String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// let c = seqs.remove_by_id(&String::from("S1")).unwrap();
        /// assert_eq!(c.id(), "S1");
        /// assert_eq!(seqs.size(), 1);
        /// let d = seqs.remove_by_id(&String::from("S2")).unwrap();
        /// assert_eq!(d.id(), "S2");
        /// assert_eq!(seqs.size(), 0);
        /// ```
        fn remove_by_id(&mut self, id: &str) -> Option<AnnotatedSequence> {
            if self.ids.contains_key(id) {
                let index = self.ids.remove(id).unwrap();
                for x in self.ids.values_mut() {
                    if *x > index {
                        *x -= 1;
                    }
                }
                Some(self.sequences.remove(index))
            } else {
                None
            }
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::from_string(String::from("S3"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.add(c);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S1");
        /// assert_eq!(seqs.get(1).unwrap().id(), "S2");
        /// assert_eq!(seqs.get(2).unwrap().id(), "S3");
        /// seqs.move_up(&String::from("S3"));
        /// assert_eq!(seqs.get(0).unwrap().id(), "S3");
        /// assert_eq!(seqs.get(1).unwrap().id(), "S1");
        /// assert_eq!(seqs.get(2).unwrap().id(), "S2");
        /// ```
        fn move_up(&mut self, id: &str) -> Result<(), SeqError> {
            match self.remove_by_id(id) {
                None => Err(SeqError::NonExistenId(id.to_string())),
                Some(x) => self.insert(0, x),
            }
        }
        fn iter(&self) -> SequenceIterable<'_> {
            SequenceIterable {
                sequences: self,
                next: 0,
            }
        }
        fn contains(&self, id: &str) -> bool {
            self.ids.contains_key(id)
        }

        /// Remove row in Sequence collection.
        /// args: index: 0-based index of the element to be removed.
        fn remove(&mut self, index: usize) -> Option<AnnotatedSequence> {
            if !self.ids.is_empty() && index < self.size() {
                let r = self.sequences.remove(index);
                self.ids.remove(r.id());
                for x in self.ids.values_mut() {
                    if *x > index {
                        *x -= 1
                    }
                }
                Some(r)
            } else {
                None
            }
        }

        /// Get a mutable Annotated sequence from a SequenceCollection
        /// ```
        /// use famlib::seqs::{
        ///     AnnotatedSequence,
        ///     SequenceCollection,
        ///     SequenceAccesors};
        /// let mut sq = SequenceCollection::new();
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("ACTG"));
        /// sq.add(a);
        /// let mut b = sq.get_mut(0).unwrap();
        /// b.set_sequence_as_string(String::from("AT"));
        /// assert_eq!(sq.get(0).unwrap().seq_as_string(), "AT");
        /// ```
        fn get_mut(&mut self, index: usize) -> Option<&mut AnnotatedSequence> {
            self.sequences.get_mut(index)
        }

        /// Changes the order of the rows in the sequence collection
        /// order should have repeated elements and the largest
        /// element should be at most as large as the number of sequences
        /// minus 1.
        fn reorder(
                &mut self,
                order: std::vec::Vec<usize>)
                -> std::result::Result<(), SeqError> {
            let nrows = self.size();
            if !order.iter().all(|x| *x<nrows) {
                return Err(SeqError::EditError)
            };
            let mut order_counter = HashMap::<usize, usize>::new();
            for i in order.iter(){
                order_counter.insert(
                    *i,
                    match order_counter.get(i){
                        None => 1,
                        Some(x) => x+1
                    }
                );
            };
            if !order_counter.values().all(|x|*x<=1) {
                return Err(SeqError::EditError)
            }
            let mut own_sequences: Vec<_> = self
                .sequences
                .iter()
                .cloned()
                .map(Some)
                .collect();
            self.sequences = order
                .iter()
                .map(|x| own_sequences[*x].take().unwrap())
                .collect::<Vec<_>>();
            self.ids = self
                .sequences
                .iter()
                .enumerate()
                .map(|(i, x)| (String::from(x.id()), i))
                .collect::<HashMap<_,_>>();
            Ok(())
        }
    }

    /// IntoIterator implementation for SequenceCollection
    /// ```
    /// use famlib::seqs::{
    ///     AnnotatedSequence,
    ///     SequenceCollection,
    ///     SequenceAccesors
    /// };
    /// fn build_a() -> SequenceCollection {
    ///     let a = AnnotatedSequence::from_string(
    ///         String::from("S1"),
    ///         String::from("A--A--C--C--"));
    ///     let b = AnnotatedSequence::from_string(
    ///         String::from("S2"),
    ///         String::from("TAGTACGATGAC"));
    ///     let c = AnnotatedSequence::from_string(
    ///         String::from("S3"),
    ///         String::from("TAGTACGATGAC"));
    ///     let mut seqcol = SequenceCollection::new();
    ///     seqcol.add(a);
    ///     seqcol.add(b);
    ///     seqcol.add(c);
    ///     return seqcol;
    /// }
    /// let seqcol = build_a();
    /// assert_eq!(seqcol.size(), 3);
    /// let x = seqcol.into_iter().collect::<Vec<_>>();
    /// assert_eq!(x.len(), 3);
    /// ```
    impl IntoIterator for SequenceCollection {
        type Item = AnnotatedSequence;
        // type IntoIter = SequenceCollectionIntoiter;
        type IntoIter = <Vec<AnnotatedSequence> as IntoIterator>::IntoIter;
        fn into_iter(self) -> <Self as IntoIterator>::IntoIter {
            self.sequences.into_iter()
            // SequenceCollectionIntoiter{data:self}
        }
    }

    /// Implements the `FromIterator` trait for `SequenceCollection`, allowing
    /// it to be created from an iterator of `AnnotatedSequence` items.
    ///
    /// # Arguments
    ///
    /// * `iter` - An iterator that yields `AnnotatedSequence` items.
    ///
    /// # Returns
    ///
    /// A `SequenceCollection` containing all the `AnnotatedSequence` items from
    /// the iterator.
    ///
    /// # Example
    ///
    /// ```
    /// use famlib::seqs::SequenceCollection;
    /// use famlib::seqs::AnnotatedSequence;
    /// use famlib::seqs::SequenceAccesors;
    ///
    /// fn build_a() -> SequenceCollection {
    ///     let a = AnnotatedSequence::from_string(
    ///         String::from("S1"),
    ///         String::from("A--A--C--C--"));
    ///     let b = AnnotatedSequence::from_string(
    ///         String::from("S2"),
    ///         String::from("TAGTACGATGAC"));
    ///     let c = AnnotatedSequence::from_string(
    ///         String::from("S3"),
    ///         String::from("TAGTACGATGAC"));
    ///     let mut seqcol = SequenceCollection::new();
    ///     seqcol.add(a);
    ///     seqcol.add(b);
    ///     seqcol.add(c);
    ///     return seqcol;
    /// }
    /// let sequences = build_a();
    /// let collection: SequenceCollection = sequences.into_iter().collect();
    /// assert_eq!(collection.size(), 3);
    /// ```
    impl FromIterator<AnnotatedSequence> for SequenceCollection {
        fn from_iter<I: IntoIterator<Item = AnnotatedSequence>>(
            iter: I
        ) -> Self {
            let mut seqcol = SequenceCollection::new();
            for x in iter { seqcol.add(x).unwrap(); }
            seqcol
        }
    }


    /// Implements the `FromIterator` trait for `SequenceCollection`, allowing
    /// it to be created from an iterator of `(Into<String>, Into<String>)`
    /// items.  The first element of the tuple is the sequence ID and the second
    /// element is the sequence itself.
    ///
    /// # Arguments
    /// * `iter` - An iterator that yields `(Into<String>, Into<String>)` items.
    ///
    /// # Returns
    /// A `SequenceCollection` containing all the `AnnotatedSequence` items from
    /// the iterator.
    /// ```
    /// use famlib::seqs::SequenceCollection;
    /// use famlib::seqs::SequenceAccesors;
    /// use std::iter::FromIterator;
    /// let seqs = vec![
    ///     (String::from("S1"), "ATCATGCTACTG"),
    ///     (String::from("S2"), "TAGTACGATGAC")
    /// ];
    /// let seqcol = SequenceCollection::from_iter(seqs.clone());
    /// let seqcol2 = seqs.into_iter().collect::<SequenceCollection>();
    /// assert_eq!(seqcol.size(), 2);
    /// assert_eq!(seqcol2.size(), 2);
    /// ```
    impl <K, V> FromIterator<(K, V)> for SequenceCollection
        where K: Into<String>, V: Into<String> {
        fn from_iter<I: IntoIterator<Item = (K, V)>>(
            iter: I
        ) -> Self {
            let mut seqcol = SequenceCollection::new();
            for (id, seq) in iter {
                let a = AnnotatedSequence::from_string(
                    id.into(),
                    seq.into()
                );
                seqcol.add(a).unwrap();
            }
            seqcol
        }
    }


    pub trait BufferedSeqCollection {
        fn next_sequence(&self) -> Option<AnnotatedSequence>;
        fn to_sequence_collection(&self) -> SequenceCollection {
            let mut sc = SequenceCollection::new();
            loop {
                match self.next_sequence() {
                    None => break,
                    Some(x) => sc.add(x).unwrap()
                }
            }
            sc
        }
    }

    pub struct BufferedSeqCollectionFromRead {
        buffer: RefCell<Option<Box<dyn BufRead>>>,
        reader: RefCell<Box<dyn SequenceReader>>,
    }

    impl BufferedSeqCollectionFromRead {
        pub fn new(buffer: Box<dyn BufRead>, format: InputFormats) -> Self {
            let reader = reader_for(format);
            BufferedSeqCollectionFromRead {
                buffer: RefCell::new(Some(buffer)),
                reader: RefCell::new(reader)
            }
        }
        pub fn is_consumed(&self) -> bool {
            return self.reader.borrow().consumed();
        }
        pub fn next_line(&self) -> Option<String> {
            let mut line = String::new();
            let bytes_read = self
                .buffer
                .borrow_mut()
                .as_mut()?
                .read_line(&mut line)
                .ok()?;
            (bytes_read > 0)
                .then_some(line)
        }
    }

    impl BufferedSeqCollection for BufferedSeqCollectionFromRead {
        fn next_sequence(&self) -> Option<AnnotatedSequence> {
            let mut reader = self.reader.borrow_mut();
            if reader.consumed() {
                return None;
            }
            loop {
                let next_line = self.next_line();
                match next_line {
                    Some(line) => reader.add_line(line),
                    None => reader.end_input(),
                }
                let returning = reader.try_build();
                if returning.is_some() || reader.consumed() {
                    return returning;
                }
            }
            // loop {
            //     let mut returning: Option<AnnotatedSequence> = None;
            //     let next_line = self.next_line();
            //     let eof = next_line.is_none();
            //     let mut line = next_line.unwrap_or_default();
            //     if eof || line.starts_with(">") {
            //         if let Some(id) = self.take_current_id() {
            //             let ann_seq = AnnotatedSequence::from_string(
            //                 id,
            //                 self.get_current_seq(),
            //             );
            //             self.set_consumed(eof);
            //             returning = Some(ann_seq);
            //         };
            //         self.clear_current_seq();
            //         if line.starts_with(">") {
            //             line.truncate(line.trim_end().len());
            //             self.set_current_id(line.split_off(1));
            //         }
            //     } else {
            //         line.truncate(line.trim_end().len());
            //         self.add_line_to_sequence(line);
            //     }
            //     if returning.is_some() || eof {
            //         return returning
            //     }
            // }
        }
    }

    pub struct BufferedSeqCollectionFromSeqCol {
        seqcol: SequenceCollection,
        current_index: RefCell<usize>,
        consumed: RefCell<bool>
    }

    impl BufferedSeqCollectionFromSeqCol {
        pub fn new(seqcol: SequenceCollection) -> Self {
            BufferedSeqCollectionFromSeqCol {
                seqcol,
                current_index: RefCell::new(0),
                consumed: RefCell::new(false)
            }
        }
    }

    impl BufferedSeqCollection for BufferedSeqCollectionFromSeqCol {
        fn next_sequence(&self) -> Option<AnnotatedSequence> {
            if *self.consumed.borrow() {
                return None;
            }
            let index = *self.current_index.borrow();
            let seq = &self.seqcol;
            if seq.get(index).is_none() {
                self.consumed.replace(true);
                None
            } else {
                self.current_index.replace(index + 1);
                seq.get(index).cloned()
            }
        }
    }

    pub struct ApplyBufferedSequenceCollection {
        source: RefCell<Box<dyn BufferedSeqCollection>>,
        apply_fun: Box<dyn Fn(AnnotatedSequence) -> Vec<AnnotatedSequence>>,
        consumed: RefCell<bool>,
        interal_buffered_sequences: RefCell<Vec<AnnotatedSequence>>
    }

    impl ApplyBufferedSequenceCollection {
        pub fn new(
            source: Box<dyn BufferedSeqCollection>,
            apply_fun: Box<dyn Fn(AnnotatedSequence) -> Vec<AnnotatedSequence>>
        ) -> ApplyBufferedSequenceCollection {
            ApplyBufferedSequenceCollection {
                source: RefCell::new(source),
                apply_fun,
                consumed: RefCell::new(false),
                interal_buffered_sequences: RefCell::new(vec![])
            }
        }
    }

    impl BufferedSeqCollection for ApplyBufferedSequenceCollection {
        fn next_sequence(&self) -> Option<AnnotatedSequence> {
            loop {
                if !self.interal_buffered_sequences.borrow().is_empty() {
                    return self
                        .interal_buffered_sequences
                        .borrow_mut()
                        .pop();
                } else {
                    match self.source.borrow().next_sequence() {
                        None => {
                            self.consumed.replace(true);
                            return None;
                        }
                        Some(seq) => {
                            let mut ib = self
                                .interal_buffered_sequences
                                .borrow_mut();
                            for x in (self.apply_fun)(seq).into_iter() {
                                ib.push(x);
                            }
                        }
                    }
                }
            }
        }
    }

    #[derive(Clone, Debug)]
    pub struct Alignment {
        pub (crate) seqs: SequenceCollection,
        // The numner of columns of the MSA
        pub (crate) length: Option<usize>,
    }

    impl Alignment {
        pub fn new() -> Self {
            Alignment {
                seqs: SequenceCollection::new(),
                length: None,
            }
        }
        /// Make A gapstripped copy of the alignment
        /// ```
        /// use famlib::seqs::{
        ///     SequenceCollection,
        ///     Alignment,
        ///     SequenceAccesors,
        ///     AnnotatedSequence};
        /// let a = AnnotatedSequence::from_string(
        ///     String::from("S1"),
        ///     String::from("A--A--C--C--"));
        /// let b = AnnotatedSequence::from_string(
        ///     String::from("S2"),
        ///     String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::from_string(
        ///     String::from("S3"),
        ///     String::from("TAGTACGATGAC"));
        /// let mut seqs = Alignment::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.add(c);
        /// let gs = seqs.gapstrip();
        /// assert_eq!(gs.length(), 4);
        /// assert_eq!(
        ///     gs.get(0).unwrap().seq().unwrap(),
        ///     &vec!['A', 'A', 'C', 'C']);
        /// assert_eq!(
        ///     gs.get(1).unwrap().seq().unwrap(),
        ///     &vec!['T', 'T', 'G', 'G']);
        /// ```
        pub fn gapstrip(&self) -> Self {
            let reference = self.get(0).unwrap().seq().unwrap();
            let mut aln = Alignment::new();
            let new_seq_iter = self.seqs.sequences.iter().map(|x| {
                (
                    reference
                        .iter()
                        .zip(x.seq().unwrap())
                        .filter_map(
                            |(a, b)| if *a != '-' { Some(b) } else { None },
                        )
                        .collect::<String>(),
                    x.id.clone(),
                )
            });
            for (new_seq, new_id) in new_seq_iter {
                aln.add(AnnotatedSequence::from_string(new_id, new_seq))
                    .unwrap();
            }
            aln
        }

        /// Return the number of columns of the MSA
        pub fn length(&self) -> usize {
            self.length.unwrap_or(0)
        }

        pub fn column(&self, index: usize) -> Option<Vec<char>> {
            if self.length() > index {
                Some(
                    self.seqs
                        .iter()
                        .map(|x| x.seq().unwrap()[index])
                        .collect::<Vec<char>>(),
                )
            } else {
                None
            }
        }

        pub fn column_ref(&self, index: usize) -> Option<Vec<&char>> {
            if self.length() > index {
                Some(
                    self.seqs
                        .iter()
                        .map(|x| &x.seq().unwrap()[index])
                        .collect(),
                )
            } else {
                None
            }
        }

        pub fn columns(&self) -> ColumnIterable {
            ColumnIterable {
                msa: self,
                next: 0,
            }
        }

        pub fn col_gap_frq(&self) -> Option<Vec<f64>> {
            if self.size() > 0 {
                Some(
                    self.columns()
                        .map(|x| {
                            x.iter().fold(0f64, |a, b| {
                                a + if **b == '-' { 1f64 } else { 0f64 }
                            }) / self.size() as f64
                        })
                        .collect::<Vec<f64>>(),
                )
            } else {
                None
            }
        }
        pub fn seq_col(&self) -> &SequenceCollection {
            &self.seqs
        }
        pub fn seq_col_owned(self) -> SequenceCollection {
            self.seqs
        }

        pub fn remove_all_gap_columns(&mut self) {
            let mut to_remove = vec![];
            for i in 0..self.length() {
                match self.column(i) {
                    None => (),
                    Some(column) => {
                        if column.iter().all(|x| *x == '-' || *x == '.') {
                            to_remove.push(i);
                        }
                    }
                }
            }
            for ann_seq in self.seqs.sequences.iter_mut() {
                ann_seq.remove_positions(&to_remove);
            }
        }

        pub fn remove_frq_gap_columns(&mut self, threshold: f64) {
            let mut to_remove = vec![];
            if let Some(frqs) = self.col_gap_frq() {
                for (i, frq) in frqs.iter().enumerate() {
                    if *frq > threshold {
                        to_remove.push(i);
                    }
                }
            }
            for ann_seq in self.seqs.sequences.iter_mut() {
                ann_seq.remove_positions(&to_remove);
            }
        }
    }

    impl Default for Alignment {
        fn default() -> Self {
            Self::new()
        }
    }

    impl SequenceAccesors for Alignment {
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(), SeqError> {
            match self.length {
                Some(x) => {
                    if x == seq.len() {
                        self.seqs.add(seq)
                    } else {
                        Err(SeqError::DifferentLength)
                    }
                }
                None => {
                    self.length = Some(seq.len());
                    self.seqs.add(seq)
                }
            }
        }
        // Retrieves the number of rows in the alignment.
        fn size(&self) -> usize {
            self.seqs.sequences.len()
        }
        fn get(&self, index: usize) -> Option<&AnnotatedSequence> {
            self.seqs.get(index)
        }

        fn get_by_id(&self, id: &str) -> Option<&AnnotatedSequence> {
            self.seqs.get_by_id(id)
        }

        fn remove(&mut self, index: usize) -> Option<AnnotatedSequence> {
            self.seqs.remove(index)
        }
        fn insert(
            &mut self,
            index: usize,
            sequence: AnnotatedSequence,
        ) -> Result<(), SeqError> {
            self.seqs.insert(index, sequence)
        }
        fn move_up(&mut self, id: &str) -> Result<(), SeqError> {
            self.seqs.move_up(id)
        }
        fn iter(&self) -> SequenceIterable<'_> {
            SequenceIterable {
                sequences: &self.seqs,
                next: 0,
            }
        }
        fn contains(&self, id: &str) -> bool {
            self.seqs.ids.contains_key(id)
        }
        fn remove_by_id(&mut self, id: &str) -> Option<AnnotatedSequence> {
            self.seqs.remove_by_id(id)
        }
        fn get_mut(&mut self, index: usize) -> Option<&mut AnnotatedSequence> {
            self.seqs.get_mut(index)
        }
        fn reorder(
            &mut self,
            order: std::vec::Vec<usize>) -> std::result::Result<(), SeqError> {
            self.seqs.reorder(order)
        }
    }

    impl Display for Alignment {
        fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
            self.seqs.fmt(f)
        }
    }

    pub struct SequenceIterable<'seqcol> {
        sequences: &'seqcol SequenceCollection,
        next: usize,
    }

    impl<'a> Iterator for SequenceIterable<'a> {
        type Item = &'a AnnotatedSequence;
        fn next(&mut self) -> Option<<Self as Iterator>::Item> {
            if self.sequences.size() > 0 {
                if self.next < self.sequences.size() {
                    let current = &self.sequences.sequences[self.next];
                    self.next += 1;
                    Some(current)
                } else {
                    None
                }
            } else {
                None
            }
        }
    }
    #[test]
    fn test_sequence_iterator() {
        let a = AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("A--A--C--C--"),
        );
        let b = AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("TAGTACGATGAC"),
        );
        let c = AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("TAGTACGATGAC"),
        );
        let mut seqs = Alignment::new();
        seqs.add(a).unwrap();
        seqs.add(b).unwrap();
        seqs.add(c).unwrap();
        for x in seqs.iter() {
            assert_eq!(x.len(), 12);
        }
    }

    pub struct ColumnIterable<'msa> {
        msa: &'msa Alignment,
        next: usize,
    }

    impl<'msa> Iterator for ColumnIterable<'msa> {
        type Item = Vec<&'msa char>;
        fn next(&mut self) -> Option<<Self as Iterator>::Item> {
            if self.msa.length() > 0 {
                if self.next < self.msa.length() {
                    let current = self.msa.column_ref(self.next);
                    self.next += 1;
                    current
                } else {
                    None
                }
            } else {
                None
            }
        }
    }

    impl IntoIterator for Alignment {
        type Item = AnnotatedSequence;
        type IntoIter = <Vec<AnnotatedSequence> as IntoIterator>::IntoIter;
        fn into_iter(self) -> <Self as IntoIterator>::IntoIter {
            self.seqs.into_iter()
        }
    }

    #[test]
    fn test_column_iterator() {
        let a = AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("A--A--C--C--"),
        );
        let b = AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("TAGTACGATGAC"),
        );
        let c = AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("TAGTACGATGAC"),
        );
        let mut seqs = Alignment::new();
        seqs.add(a).unwrap();
        seqs.add(b).unwrap();
        seqs.add(c).unwrap();
        let clms = seqs.columns().collect::<Vec<Vec<&char>>>();
        assert_eq!(clms[0], vec![&'A', &'T', &'T']);
        assert_eq!(clms[1], vec![&'-', &'A', &'A']);
        assert_eq!(clms[11], vec![&'-', &'C', &'C']);
    }
}

#[cfg(test)]
mod test{
    use crate::seqs::BufferedSeqCollection;
    use crate::seqs::SequenceCollection;
    use crate::seqs::Alignment;
    use crate::seqs::AnnotatedSequence;
    use crate::seqs::SequenceAccesors;
    use crate::seqs::SeqError;
    fn sample_sc() -> Result<SequenceCollection, SeqError> {
        let a = AnnotatedSequence::from_string(
            String::from("S1"),
            String::from("A"),
        );
        let b = AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("AA"),
        );
        let c = AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("AAA"),
        );
        let d = AnnotatedSequence::from_string(
            String::from("S4"),
            String::from("AAAA"),
        );
        let mut sq = SequenceCollection::new();
        sq.add(a)?;
        sq.add(b)?;
        sq.add(c)?;
        sq.add(d)?;
        Ok(sq)
    }
    fn sample_display_sc() -> SequenceCollection {
        let mut sq = SequenceCollection::new();
        for i in 0..10 {
            let s = match i==0 {
                true => String::from(""),
                false => {
                    (0..10*(i-1))
                        .into_iter()
                        .map(|_| 'A')
                        .collect::<String>()
                        + &String::from("AAAAABBBBB")
                }
            };
            if i < 5 {
                let ca = AnnotatedSequence::from_string(
                    format!("Seq_{}", i), s);
                sq.add(ca).unwrap();
            }
            else {
                let ca = AnnotatedSequence::from_string(
                    format!("Sequence_00{}", i), s);
                sq.add(ca).unwrap();
            }
        }
        sq
    }
    fn sample_display_msa() -> Alignment {
        let mut sq = SequenceCollection::new();
        for i in 0..10 {
            let ca = AnnotatedSequence::from_string(
                format!("Seq_{}", i),
                String::from("ATGKKKGTGCATTAA")
            );
            sq.add(ca).unwrap();
        }
        sq.to_msa().ok().unwrap()
    }
    fn sample_gapped_msa() -> Alignment {
        let mut sq = SequenceCollection::new();
        let seqs = vec![
            "-----TT-CA",
            "AT---TT-CA",
            "ATG-----CA",
            "ATGC----CA",
        ];
        for (i, s) in seqs.into_iter().enumerate() {
            let ca = AnnotatedSequence::from_string(
                format!("Seq_{}", i),
                String::from(s)
            );
            sq.add(ca).unwrap();
        }
        sq.to_msa().ok().unwrap()
    }
    #[test]
    fn reorder_equal_order() {
        let mut msa = sample_sc().unwrap();
        assert!(msa.reorder(vec![0,1,2,3]).is_ok());
        assert_eq!(msa.get(0).unwrap().id(), "S1");
        assert_eq!(msa.get(1).unwrap().id(), "S2");
        assert_eq!(msa.get(2).unwrap().id(), "S3");
        assert_eq!(msa.get(3).unwrap().id(), "S4");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "A");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "AAAA");
    }
    #[test]
    fn reorder_inverse_order() {
        let mut msa = sample_sc().unwrap();
        assert!(msa.reorder(vec![3,2,1,0]).is_ok());
        assert_eq!(msa.get(0).unwrap().id(), "S4");
        assert_eq!(msa.get(1).unwrap().id(), "S3");
        assert_eq!(msa.get(2).unwrap().id(), "S2");
        assert_eq!(msa.get(3).unwrap().id(), "S1");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "AAAA");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "AA");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "A");
    }
    #[test]
    fn reorder_random_order() {
        let mut msa = sample_sc().unwrap();
        assert!(msa.reorder(vec![2,3,0,1]).is_ok());
        assert_eq!(msa.get(0).unwrap().id(), "S3");
        assert_eq!(msa.get(1).unwrap().id(), "S4");
        assert_eq!(msa.get(2).unwrap().id(), "S1");
        assert_eq!(msa.get(3).unwrap().id(), "S2");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AAAA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "A");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "AA");
    }
    #[test]
    fn reorder_with_big_indexes() {
        let mut msa = sample_sc().unwrap();
        msa.reorder(vec![2,5,0,1]).expect_err("Should fail");
        assert_eq!(msa.get(0).unwrap().id(), "S1");
        assert_eq!(msa.get(1).unwrap().id(), "S2");
        assert_eq!(msa.get(2).unwrap().id(), "S3");
        assert_eq!(msa.get(3).unwrap().id(), "S4");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "A");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "AAAA");
    }
    #[test]
    fn reorder_with_repeated_indexes() {
        let mut msa = sample_sc().unwrap();
        msa.reorder(vec![2,0,0,1]).expect_err("Should fail");
        assert_eq!(msa.get(0).unwrap().id(), "S1");
        assert_eq!(msa.get(1).unwrap().id(), "S2");
        assert_eq!(msa.get(2).unwrap().id(), "S3");
        assert_eq!(msa.get(3).unwrap().id(), "S4");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "A");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "AAAA");
    }
    #[test]
    fn reorder_with_few_indexes() {
        let mut msa = sample_sc().unwrap();
        assert!(msa.reorder(vec![2,0,1]).is_ok());
        assert_eq!(msa.get(0).unwrap().id(), "S3");
        assert_eq!(msa.get(1).unwrap().id(), "S1");
        assert_eq!(msa.get(2).unwrap().id(), "S2");
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "AAA");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "A");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "AA");
    }

    #[test]
    fn test_display_seq_col() {
        let expected = "\
> Seq_0                                                             0
> Seq_1      AAAAABBBBB                                             10
> Seq_2      AAAAAAAAAA AAAAABBBBB                                  20
...
> Sequ...007 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA ...AABBBBB 70
> Sequ...008 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA ...AABBBBB 80
> Sequ...009 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA ...AABBBBB 90
";
        let scol = sample_display_sc();
        assert_eq!(format!("{}", scol), expected);
    }
    #[test]
    fn test_display_msa() {
        let expected = "\
> Seq_0      ATGKKKGTGC ATTAA      15
> Seq_1      ATGKKKGTGC ATTAA      15
> Seq_2      ATGKKKGTGC ATTAA      15
...
> Seq_7      ATGKKKGTGC ATTAA      15
> Seq_8      ATGKKKGTGC ATTAA      15
> Seq_9      ATGKKKGTGC ATTAA      15
";
        let msa = sample_display_msa();
        let display_test = format!("{}", msa);
        assert_eq!(display_test, expected);
    }
    #[test]
    fn test_buffered_seq_collection_from_read() {
        use std::io::Cursor;
        use crate::fastaio::InputFormats;
        let buffer = Cursor::new(
            ">S1\nATC\n>S2\nATG\n>S3\nATT\n");
        let bsc = crate::seqs::BufferedSeqCollectionFromRead::new(
            Box::new(buffer),
            InputFormats::Fasta
        );
        let seq1 = bsc.next_sequence().unwrap();
        assert_eq!(seq1.id(), "S1");
        assert_eq!(seq1.seq_as_string(), "ATC");
        let seq2 = bsc.next_sequence().unwrap();
        assert_eq!(seq2.id(), "S2");
        assert_eq!(seq2.seq_as_string(), "ATG");
        let seq3 = bsc.next_sequence().unwrap();
        assert_eq!(seq3.id(), "S3");
        assert_eq!(seq3.seq_as_string(), "ATT");
    }
    #[test]
    fn test_remove_all_gap_columns() {
        let mut msa = sample_gapped_msa();
        // input:
        // -----TT-CA
        // AT---TT-CA
        // ATG-----CA
        // ATGC----CA
        // expected output:
        // ----TTCA
        // AT--TTCA
        // ATG---CA
        // ATGC--CA
        msa.remove_all_gap_columns();
        assert_eq!(msa.get(0).unwrap().seq_as_string(), "----TTCA");
        assert_eq!(msa.get(1).unwrap().seq_as_string(), "AT--TTCA");
        assert_eq!(msa.get(2).unwrap().seq_as_string(), "ATG---CA");
        assert_eq!(msa.get(3).unwrap().seq_as_string(), "ATGC--CA");
    }
}
