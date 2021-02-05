pub mod edit;
pub mod edit_msa;
pub mod fastaio;
pub mod merge;
pub mod random;

pub mod seqs {
    use std::{collections::HashMap, io::ErrorKind};
    use std::fmt::{Display, Error, Formatter};
    use std::iter::{IntoIterator, Iterator};

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
            return std::io::Error::new(
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
    #[derive(Clone, PartialEq, Debug)]
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

        pub fn from_string(id: String, sequence: String) -> Self {
            AnnotatedSequence {
                id,
                sequence: Some(sequence.chars().collect::<Vec<char>>()),
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
        /// let mut a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATC"));
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
        pub fn set_sequence(&mut self, seq: Vec<char>) -> () {
            self.sequence = Some(seq);
        }

        pub fn set_sequence_as_string(&mut self, seq: String) -> () {
            if seq.len() == 0 {
                self.sequence = None;
            } else {
                self.sequence = Some(seq.chars().collect());
            }
        }

        /// Moves out the sequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATC"));
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
        pub fn set_id(&mut self, id: String) -> () {
            self.id = id;
        }

        /// Retrieves a reference of the sequence.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCG"));
        /// assert_eq!(a.seq() , Some(&vec!['A', 'T', 'C', 'G']));
        /// ```
        pub fn seq(&self) -> Option<&Vec<char>> {
            match &self.sequence {
                None => None,
                Some(x) => Some(&x),
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
    }

    #[derive(Clone, PartialEq)]
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
        pub fn to_msa(self) -> Result<Alignment, Self> {
            let mut msa = Alignment::new();
            if self.size() == 0 {
                msa.seqs = self;
                Ok(msa)
            } else {
                let ref_len = self.sequences[0].len();
                if self.sequences.iter().all(|x| x.len() == ref_len) {
                    msa.seqs = self;
                    Ok(msa)
                } else {
                    Err(self)
                }
            }
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
                        *x = *x + 1;
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
            if self.ids.contains_key(&seq.id().to_string()) {
                Err(SeqError::DuplicatedId(seq.id().to_string()))
            } else {
                self.ids.insert(seq.id().to_string(), self.ids.len());
                self.sequences.push(seq);
                Ok(())
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
                        *x = *x - 1;
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
                sequences: &self,
                next: 0,
            }
        }
        fn contains(&self, id: &str) -> bool {
            self.ids.contains_key(id)
        }

        /// Remove row in Sequence collection.
        /// args: index: 0-based index of the element to be removed.
        fn remove(&mut self, index: usize) -> Option<AnnotatedSequence> {
            if self.ids.len() > 0 && index < self.size() {
                let r = self.sequences.remove(index);
                self.ids.remove(r.id());
                for x in self.ids.values_mut() {
                    if *x > index {
                        *x = *x - 1
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
                .to_owned()
                .into_iter()
                .map(|x| Some(x))
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
            match self.length {
                None => 0,
                Some(x) => x,
            }
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
                msa: &self,
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
    use crate::seqs::SequenceCollection;
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
}
