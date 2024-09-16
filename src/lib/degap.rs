use crate::seqs::{
  Alignment,
  AnnotatedSequence,
  ApplyBufferedSequenceCollection,
  BufferedSeqCollection,
  SequenceAccesors,
  SequenceCollection
};

pub trait Degap<T> {
    fn degap(&self) -> T;
}

pub struct DegapBufferedSequenceCollection{}

impl Degap<AnnotatedSequence> for AnnotatedSequence {
    fn degap(&self) -> AnnotatedSequence {
        let seq = self
          .seq()
          .map(
            |x| x.clone().into_iter().filter(|y| *y != '-').collect()
          )
          .unwrap_or(vec![]);
        AnnotatedSequence::new(self.id().to_string(), seq)
    }
}

impl DegapBufferedSequenceCollection {
    pub fn degap(
        bsc: Box<dyn BufferedSeqCollection>
    ) -> ApplyBufferedSequenceCollection {
        let filter_func = move |s: AnnotatedSequence| { vec![s.degap()] };
        ApplyBufferedSequenceCollection::new(
            bsc,
            Box::new(filter_func)
        )
    }
}

impl Degap<SequenceCollection> for SequenceCollection {
    fn degap(&self) -> SequenceCollection {
        let mut result = SequenceCollection::new();
        self.iter()
            .for_each(
                |x| result.add(x.degap()).unwrap()
            );
        result
    }
}

impl Degap<Alignment> for Alignment {
    fn degap(&self) -> Alignment {
        let mut result = Alignment::new();
        self.iter()
            .for_each(
                |x| result.add(x.degap()).unwrap()
            );
        result
    }
}

#[cfg(test)]
mod test {
  use crate::seqs::AnnotatedSequence;
  use crate::seqs::SequenceAccesors;
  use crate::seqs::SequenceCollection;
  use super::Degap;

  #[test]
  fn test_degap_annotated_sequence() {
    use crate::seqs::AnnotatedSequence;
    let seq = AnnotatedSequence::new("id".to_string(), "A--T".chars().collect());
    let degapped = seq.degap();
    assert_eq!(degapped.seq().unwrap(), &vec!['A', 'T']);
  }

  #[test]
  fn test_degap_sequence_collection() {
    let mut sq = SequenceCollection::new();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_01".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_02".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "Sequence_03".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    let result = sq.degap();
    assert!(result.size() == 3);
    assert!(result.iter().all(|x| x.seq_as_string().eq("ACA")))
  }
  #[test]
  fn test_degap_alignment() {
    let mut sq = SequenceCollection::new();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_01".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_02".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "Sequence_03".to_string(),
            "AC-------A".to_string()
        )
    ).unwrap();
    let result = sq.to_msa().unwrap().degap();
    assert!(result.size() == 3);
    assert!(result.iter().all(|x| x.seq_as_string().eq("ACA")))
  }
}