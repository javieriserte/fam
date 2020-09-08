pub mod fastaio;
pub mod merge;

pub mod seqs {
    use std::fmt::{Display, Error, Formatter};
    use std::collections::HashSet;
    use std::iter::{IntoIterator, Iterator};

    #[derive(Debug)]
    pub enum SeqError{
        DuplicatedId(String),
        DifferentLength,
        NonExistenId(String),
        MissingID(String)
    }
    impl Display for SeqError{
        fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
            match self{
                SeqError::DuplicatedId(x) => write!(f,
                    "Attempted to add a sequence with a \
                    duplicated ID ['{}'].", x),
                SeqError::NonExistenId(x) => write!(f,
                    "Attempted to get a non existing ID [{}].", x),
                SeqError::DifferentLength => write!(f,
                    "Attempted to add a sequence with \
                    different length to a MSA."),
                SeqError::MissingID(x) => write!(f,
                    "A sequence with ID [{}] is required", x)
            }
        }
    }

    pub trait SequenceAccesors {
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(),SeqError> ;
        fn get(&self, index: usize) -> Option<&AnnotatedSequence>;
        fn get_mut(&mut self, index: usize) -> Option<&mut AnnotatedSequence>;
        fn remove(&mut self, id: &String) -> Option<AnnotatedSequence>;
        fn add_on_top(&mut self, seq: AnnotatedSequence) -> Result<(),SeqError>;
        fn size(&self) -> usize;
        fn move_up(&mut self, id: &String) -> Result<(), SeqError>;
        fn iter(&self) -> SequenceIterable;
        fn contains(&self, id: &str) -> bool;
        fn take_first(&mut self) -> Option<AnnotatedSequence>;
    }

    /// Struct to represent a single sequence of a MSA
    #[derive(Clone, PartialEq, Debug)]
    pub struct AnnotatedSequence {
        id: String ,
        sequence: Option<Vec<char>>
    }

    impl AnnotatedSequence {
        /// Creates a new AnnotatedSequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::new(String::from("S1"), vec!['A', 'T', 'C', 'A', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G']);
        /// assert_eq!(a.id(), "S1");
        /// assert_eq!(a.seq(), Some(&vec!['A', 'T', 'C', 'A', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G']));
        /// ```
        pub fn new(id: String, sequence: Vec<char>) -> Self {
            AnnotatedSequence{
                id,
                sequence:Some(sequence)
            }
        }

        pub fn from_string(id: String, sequence: String) -> Self {
            AnnotatedSequence{
                id,
                sequence:Some(sequence.chars().collect::<Vec<char>>())
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
            AnnotatedSequence{
                id,
                sequence: None
            }
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
                Some(x) => Some(&x)
            }
        }

        pub fn seq_as_string(&self) -> String {
            match &self.sequence {
                None => String::from(""), 
                Some(x) => x.iter().collect::<String>()
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
                Some(x) => x.len()
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
                Some(x) => Ok(x.clone())
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
        ids: HashSet<String> // TODO: replace with a hashmap to the
                             // sequence indexes.
    }

    impl SequenceCollection {
        /// Creates a new SequenceCollection
        /// ```
        /// use famlib::seqs::SequenceCollection;
        /// let a = SequenceCollection::new();
        /// ```
        pub fn new() -> Self {
            SequenceCollection{
                sequences: vec![],
                ids: HashSet::new()
            }
        }

        /// Mute the SequenceCollection to a Alignment
        pub fn to_msa(self) -> Result<Alignment, Self> {
            let mut msa = Alignment::new();
            if self.size()==0 {
                msa.seqs = self;
                Ok(msa)
            } else {
                let ref_len = self.sequences[0].len();
                if self.sequences.iter().all(|x| x.len()==ref_len) {
                    msa.seqs = self;
                    Ok(msa)
                } else {
                    Err(self)
                }
            }
        }
    }

    impl SequenceAccesors for SequenceCollection {
        // fn contains(&self, id: &str) -> bool{
        //     self.ids.contains(id)
        // }
        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::from_string(String::from("S3"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.add_on_top(c);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S3");
        /// ```
        fn add_on_top(&mut self, seq: AnnotatedSequence) -> Result<(), SeqError> {
            if self.ids.contains(&seq.id().to_string()) {
                Err(SeqError::DuplicatedId(String::from(seq.id())))
                // Err(String::from("Repeated ID."))
            } else {
                self.ids.insert(seq.id().to_string());
                self.sequences.insert(0, seq);
                Ok(())
            }
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::from_string(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S1");
        /// assert_eq!(seqs.get(1).unwrap().id(), "S2");
        /// ```
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(), SeqError> {
            if self.ids.contains(&seq.id().to_string()) {
                Err(SeqError::DuplicatedId(seq.id().to_string()))
            } else {
                self.ids.insert(seq.id().to_string());
                self.sequences.push(seq);
                Ok(())
            }
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// assert_eq!(seqs.size(), 0);
        /// seqs.add(a);
        /// assert_eq!(seqs.size(), 1);
        /// ```
        fn size(&self) -> usize {
            self.sequences.len()
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// assert_eq!(seqs.get(0).unwrap().id(), String::from("S1"));
        /// ```
        fn get(&self, index: usize) -> Option<&AnnotatedSequence> {
            self.sequences.get(index)
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::from_string(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// let b = seqs.remove(&String::from("S1")).unwrap();
        /// assert_eq!(b.id(), "S1");
        /// assert_eq!(seqs.size(), 0);
        /// ```
        fn remove(&mut self, id: &String) -> Option<AnnotatedSequence> {
            if self.ids.contains(id) {
                self.ids.remove(id);
                let found:Vec<usize> = self.sequences.iter().enumerate().filter_map(|(i, annseq)| match annseq.id() == id {
                    true => Some(i),
                    false => None
                }).collect();
                if found.len() == 1 {
                    Some(self.sequences.remove(found[0]))
                } else {
                    panic!("Weird stuff happened")
                }
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
        fn move_up(&mut self, id: &String) -> Result<(), SeqError>{
            match self.remove(id) {
                None => Err(SeqError::NonExistenId(id.to_string())),
                Some(x) => self.add_on_top(x)
            }
        }    
        fn iter(&self) -> SequenceIterable<'_> { 
            SequenceIterable{
                sequences: &self,
                next:0
            }
        }
        fn contains(&self, id: &str) -> bool {
            self.ids.contains(id)
        }

        // fn take(&mut self, id: &str) -> Option<AnnotatedSequence> {
        //     if self.contains(id) {
        //         for (i, s) in self.sequences.iter().enumerate() {
        //             if s.id() == id {
        //                 return Some(self.sequences.remove(i));
        //             }
        //         };
        //         None
        //     } else {
        //         None
        //     }
        // }

        // Consumes the first element of the sequence Collection
        fn take_first(&mut self) -> Option<AnnotatedSequence> {
            if self.ids.len()>0 {
                let r = self.sequences.remove(0);
                self.ids.remove(r.id());
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
        seqs: SequenceCollection,
        // The numner of columns of the MSA
        length: Option<usize>
    }

    impl Alignment {
        pub fn new() -> Self {
            Alignment{
                seqs: SequenceCollection::new(),
                length: None
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
            let new_seq_iter = self.seqs.sequences
                .iter()
                .map(|x| (
                    reference.iter()
                            .zip(x.seq().unwrap())
                            .filter_map(
                                |(a, b)| if *a != '-' {Some(b)} else {None})
                            .collect::<String>(),
                    x.id.clone()));
            for (new_seq, new_id) in new_seq_iter{
                aln.add(AnnotatedSequence::from_string(new_id, new_seq)).unwrap();
            };
            aln
        }

        /// Return the number of columns of the MSA
        pub fn length(&self) -> usize {
            match self.length {
                None => 0,
                Some(x) => x
            }
        }

        pub fn column(&self, index: usize) -> Option<Vec<char>> {
            if self.length() > index {
                Some(self.seqs.iter().map(|x| x.seq().unwrap()[index]).collect::<Vec<char>>())
            } else {
                None
            }
        }

        pub fn column_ref(&self, index: usize) -> Option<Vec<&char>> {
            if self.length() > index {
                Some(self.seqs.iter().map(|x| &x.seq().unwrap()[index]).collect())
            } else {
                None
            }
        }

        pub fn columns(&self) -> ColumnIterable {
            ColumnIterable{
                msa: &self,
                next:0
            }
        }

        pub fn col_gap_frq(&self) -> Option<Vec<f64>> {
            if self.size()>0 {
            Some(self.columns()
                .map(|x| {
                    x.iter()
                     .fold(0f64, |a, b| a + if **b == '-' {1f64} else {0f64}) / self.size() as f64
                }).collect::<Vec<f64>>())
            } else {
                None
            }
        }
    }

    impl SequenceAccesors for Alignment {
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(),SeqError> {
            match self.length {
                Some(x) => {
                    if x == seq.len() {
                        self.seqs.add(seq)
                    } else {
                        Err(SeqError::DifferentLength)
                    }
                },
                None => {
                    self.length = Some(seq.len());
                    self.seqs.add(seq)
                }
            }
        }
        fn size(&self) -> usize {
            self.seqs.sequences.len()
        }
        fn get(&self, index : usize) -> Option<&AnnotatedSequence> {
            self.seqs.get(index)
        }
        fn remove(&mut self, id: &String) -> Option<AnnotatedSequence> {
            self.seqs.remove(id)
        }
        fn add_on_top(&mut self, sequence: AnnotatedSequence) -> Result<(), SeqError> {
            self.seqs.add_on_top(sequence)
        }
        fn move_up(&mut self, id: &String) -> Result<(), SeqError> {
            self.seqs.move_up(id)
        }
        fn iter(&self) -> SequenceIterable<'_> { 
            SequenceIterable{
                sequences: &self.seqs,
                next:0
            }
        }
        fn contains(&self, id: &str) -> bool {
            self.seqs.ids.contains(id)
        }
        fn take_first(&mut self) -> Option<AnnotatedSequence> {
            self.seqs.take_first()
        }
        fn get_mut(&mut self, index: usize) -> Option<&mut AnnotatedSequence> {
            self.seqs.get_mut(index)
        }
    }

    pub struct SequenceIterable<'seqcol>{
        sequences: &'seqcol SequenceCollection,
        next: usize,
    }

    impl<'a> Iterator for SequenceIterable<'a> {
        type Item = &'a AnnotatedSequence;
        fn next(&mut self) -> Option<<Self as Iterator>::Item> { 
            if self.sequences.size()>0 {
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
            String::from("A--A--C--C--"));
        let b = AnnotatedSequence::from_string(
            String::from("S2"),
            String::from("TAGTACGATGAC"));
        let c = AnnotatedSequence::from_string(
            String::from("S3"),
            String::from("TAGTACGATGAC"));
        let mut seqs = Alignment::new();
        seqs.add(a).unwrap();
        seqs.add(b).unwrap();
        seqs.add(c).unwrap();
        for x in seqs.iter(){
            assert_eq!(x.len(), 12);
        }
    }

    pub struct ColumnIterable<'msa> {
        msa: &'msa Alignment,
        next: usize
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
        let a = AnnotatedSequence::from_string(String::from("S1"), String::from("A--A--C--C--"));
        let b = AnnotatedSequence::from_string(String::from("S2"), String::from("TAGTACGATGAC"));
        let c = AnnotatedSequence::from_string(String::from("S3"), String::from("TAGTACGATGAC"));
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