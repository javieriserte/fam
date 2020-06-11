pub mod fastaio;

pub mod seqs {
    use std::collections::HashSet;

    pub trait SequenceAccesors {
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(),String> ;
        fn get(&self, index: usize) -> Option<&AnnotatedSequence>;
        fn remove(&mut self, id: &String) -> Option<AnnotatedSequence>;
        fn add_on_top(&mut self, seq: AnnotatedSequence) -> Result<(),String>;
        fn size(&self) -> usize;
        fn move_up(&mut self, id: &String) -> Result<(), String>;
        fn iter(&self) -> SequenceIterable;
    }


    /// Struct to represent a single sequence of a MSA
    #[derive(Clone, PartialEq, Debug)]
    pub struct AnnotatedSequence {
        id: String ,
        sequence: Option<String>
    }

    impl AnnotatedSequence {
        /// Creates a new AnnotatedSequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// assert_eq!(a.id() , "S1");
        /// assert_eq!(a.seq() , Some("ATCATGCTACTG"));
        /// ```
        pub fn new(id: String, sequence: String) -> Self {
            AnnotatedSequence{
                id,
                sequence:Some(sequence)
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
        /// let mut a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// assert_eq!(a.seq() , Some("ATCATGCTACTG"));
        /// a.set_sequence(String::from("GGGG"));
        /// assert_eq!(a.seq() , Some("GGGG"));
        /// 
        /// let mut b = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(b.seq() , None);
        /// b.set_sequence(String::from("GGGG"));
        /// assert_eq!(b.seq() , Some("GGGG"));
        /// 
        /// ```
        pub fn set_sequence(&mut self, seq: String) -> () {
            self.sequence = Some(seq);
        }

        /// Moves out the sequence
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// assert_eq!(a.seq(), Some("ATCATGCTACTG"));
        /// let b = a.take_sequence();
        /// assert_eq!(a.seq() , None);
        /// assert_eq!(b , Some(String::from("ATCATGCTACTG")));
        /// ```
        pub fn take_sequence(&mut self) -> Option<String> {
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
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCG"));
        /// assert_eq!(a.seq() , Some("ATCG"));
        /// ```
        pub fn seq(&self) -> Option<&str> {
            match &self.sequence {
                None => None, 
                Some(x) => Some(&x)
            }
        }

        /// Retrieves a reference of the sequence.
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let mut a = AnnotatedSequence::empty(String::from("S1"));
        /// assert_eq!(a.len(), 0);
        /// a.set_sequence(String::from("ATCG"));
        /// assert_eq!(a.len(), 4);
        /// a.set_sequence(String::from("AT"));
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
        /// a.set_sequence(String::from("ATCG"));
        /// assert_eq!(a.seq_copy().unwrap(), String::from("ATCG"));
        /// assert_eq!(a.seq(), Some("ATCG"));
        /// ```
        pub fn seq_copy(&self) -> Result<String, String> {
            match &self.sequence {
                None => Err(String::from("Sequence is empty")),
                Some(x) => Ok(x.clone())
            }
        }
        
        /// Gets a reference of the id
        /// ```
        /// use famlib::seqs::AnnotatedSequence;
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// assert_eq!(a.id() , "S1");
        /// ```
        pub fn id(&self) -> &str {
            &self.id
        }

    }

    #[derive(Clone, PartialEq)]
    pub struct SequenceCollection {
        sequences: Vec<AnnotatedSequence>,
        ids: HashSet<String>
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
        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::new(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::new(String::from("S3"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.add_on_top(c);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S3");
        /// ```
        fn add_on_top(&mut self, seq: AnnotatedSequence) -> Result<(), String> {
            if self.ids.contains(&seq.id().to_string()) {
                Err(String::from("Repeated ID."))
            } else {
                self.ids.insert(seq.id().to_string());
                self.sequences.insert(0, seq);
                Ok(())
            }
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::new(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// assert_eq!(seqs.get(0).unwrap().id(), "S1");
        /// assert_eq!(seqs.get(1).unwrap().id(), "S2");
        /// ```
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(), String> {
            if self.ids.contains(&seq.id().to_string()) {
                Err(String::from("Repeated ID."))
            } else {
                self.ids.insert(seq.id().to_string());
                self.sequences.push(seq);
                Ok(())
            }
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
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
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let mut seqs = SequenceCollection::new();
        /// seqs.add(a);
        /// assert_eq!(seqs.get(0).unwrap().id(), String::from("S1"));
        /// ```
        fn get(&self, index: usize) -> Option<&AnnotatedSequence> {
            self.sequences.get(index)
        }

        /// ```
        /// use famlib::seqs::{SequenceCollection, AnnotatedSequence, SequenceAccesors};
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
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
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("ATCATGCTACTG"));
        /// let b = AnnotatedSequence::new(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::new(String::from("S3"), String::from("TAGTACGATGAC"));
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
        fn move_up(&mut self, id: &String) -> Result<(), String>{
            match self.remove(id) {
                None => Err(String::from("ID not exists")),
                Some(x) => {
                    self.add_on_top(x)
                }
            }
        }    
        fn iter(&self) -> SequenceIterable<'_> { 
            SequenceIterable{
                sequences: &self,
                next:0
            }
        }
    }

    pub struct Alignment {
        seqs: SequenceCollection,
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
        /// use famlib::seqs::{SequenceCollection, Alignment, SequenceAccesors, AnnotatedSequence};
        /// let a = AnnotatedSequence::new(String::from("S1"), String::from("A--A--C--C--"));
        /// let b = AnnotatedSequence::new(String::from("S2"), String::from("TAGTACGATGAC"));
        /// let c = AnnotatedSequence::new(String::from("S3"), String::from("TAGTACGATGAC"));
        /// let mut seqs = Alignment::new();
        /// seqs.add(a);
        /// seqs.add(b);
        /// seqs.add(c);
        /// let gs = seqs.gapstrip();
        /// assert_eq!(gs.length(), 4);
        /// assert_eq!(gs.get(0).unwrap().seq().unwrap(), "AACC");
        /// assert_eq!(gs.get(1).unwrap().seq().unwrap(), "TTGG");
        /// ```
        pub fn gapstrip(&self) -> Self {
            let reference = self.get(0).unwrap().seq().unwrap();
            let mut aln = Alignment::new();
            let new_seq_iter = self.seqs.sequences
                .iter()
                .map(|x| (
                    reference.chars()
                            .zip(x.seq().unwrap().chars())
                            .filter_map(|(a, b)| if a != '-' {Some(b)} else {None})
                            .collect::<String>(),
                    x.id.clone()));
            for (new_seq, new_id) in new_seq_iter{
                aln.add(AnnotatedSequence::new(new_id, new_seq)).unwrap();
            };
            aln
        }
        pub fn length(&self) -> usize {
            match self.length {
                None => 0,
                Some(x) => x
            }
        }
    }

    impl SequenceAccesors for Alignment {
        fn add(&mut self, seq: AnnotatedSequence) -> Result<(),String> {
            match self.length {
                Some(x) => {
                    if x == seq.len() {
                        self.seqs.add(seq)
                    } else {
                        Err(format!("Sequence {} has no the same length than the previous sequence", seq.id()))
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
        fn add_on_top(&mut self, sequence: AnnotatedSequence) -> Result<(), String> {
            self.seqs.add_on_top(sequence)
        }
        fn move_up(&mut self, id: &String) -> Result<(), String> {
            self.seqs.move_up(id)
        }
        fn iter(&self) -> SequenceIterable<'_> { 
            SequenceIterable{
                sequences: &self.seqs,
                next:0
            }
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
    fn teste_sequence_iterator() {
        let a = AnnotatedSequence::new(String::from("S1"), String::from("A--A--C--C--"));
        let b = AnnotatedSequence::new(String::from("S2"), String::from("TAGTACGATGAC"));
        let c = AnnotatedSequence::new(String::from("S3"), String::from("TAGTACGATGAC"));
        let mut seqs = Alignment::new();
        seqs.add(a).unwrap();
        seqs.add(b).unwrap();
        seqs.add(c).unwrap();
        for x in seqs.iter(){
            assert_eq!(x.len(), 12);
        }
    }
}