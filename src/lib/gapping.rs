use crate::seqs::{
  Alignment,
  AnnotatedSequence,
  ApplyBufferedSequenceCollection,
  BufferedSeqCollection,
  SequenceAccesors,
  SequenceCollection
};

pub trait PadWithGaps<T> {
  fn pad_with_gaps(&self, max_length: usize) -> T;
  fn pad_with_gaps_to_max_length(&self) -> T;
}
pub trait PadWithGapsOnce<T> {
  fn pad_with_gaps(self, max_length: usize) -> T;
  fn pad_with_gaps_to_max_length(&self) -> T;
}

pub struct PadWithGapsBufferedSequenceCollection{}

impl PadWithGaps<AnnotatedSequence> for AnnotatedSequence {
  fn pad_with_gaps(&self, max_length: usize) -> AnnotatedSequence {
    let seq = self
      .seq()
      .map(
        |x| {
          let mut result = x.clone();
          while result.len() < max_length {
            result.push('-');
          }
          result
        }
      )
      .unwrap_or(vec![]);
    AnnotatedSequence::new(self.id().to_string(), seq)
  }
  fn pad_with_gaps_to_max_length(&self) -> AnnotatedSequence {
    unimplemented!();
  }
}

impl PadWithGapsOnce<ApplyBufferedSequenceCollection>
    for Box<dyn BufferedSeqCollection> {
  fn pad_with_gaps(self, max_length: usize) -> ApplyBufferedSequenceCollection {
    let filter_func = move |s: AnnotatedSequence| {
      vec![s.pad_with_gaps(max_length)]
    };
    ApplyBufferedSequenceCollection::new(
      self,
      Box::new(filter_func)
    )
  }
  fn pad_with_gaps_to_max_length(&self) -> ApplyBufferedSequenceCollection {
    unimplemented!();
  }
}

impl PadWithGaps<SequenceCollection> for SequenceCollection {
  fn pad_with_gaps(&self, max_length: usize) -> SequenceCollection {
    self.iter()
      .map(|x| x.pad_with_gaps(max_length))
      .collect::<SequenceCollection>()
  }
  fn pad_with_gaps_to_max_length(&self) -> SequenceCollection {
    let max_length = self.iter().map(|x| x.seq().unwrap().len()).max();
    self.pad_with_gaps(max_length.unwrap())
  }
}

impl PadWithGaps<Alignment> for SequenceCollection {
  fn pad_with_gaps(&self, max_length: usize) -> Alignment {
    let padded: SequenceCollection = self.pad_with_gaps(max_length);
    padded.to_msa().unwrap()
  }
  fn pad_with_gaps_to_max_length(&self) -> Alignment {
    let max_length = self.iter().map(|x| x.seq().unwrap().len()).max();
    self.pad_with_gaps(max_length.unwrap())
  }
}

#[cfg(test)]
mod test {
  use crate::seqs::AnnotatedSequence;
  use crate::seqs::SequenceAccesors;
  use crate::seqs::SequenceCollection;
  use crate::seqs::Alignment;
  use super::PadWithGaps;

  #[test]
  fn test_pad_annotated_sequence() {
    use crate::seqs::AnnotatedSequence;
    let seq = AnnotatedSequence::new("id".to_string(), vec!['A', 'C']);
    let padded = seq.pad_with_gaps(5);
    assert_eq!(padded.seq_as_string(), "AC---");
  }
  #[test]
  fn test_pad_sequence_collection() {
    let seq1 = AnnotatedSequence::new("id1".to_string(), vec!['A', 'C']);
    let seq2 = AnnotatedSequence::new("id2".to_string(), vec!['A', 'C', 'G']);
    let mut seqs = SequenceCollection::new();
    seqs.add(seq1).ok();
    seqs.add(seq2).ok();
    let padded: SequenceCollection = seqs.pad_with_gaps(5);
    assert_eq!(padded.get(0).unwrap().seq_as_string(), "AC---");
    assert_eq!(padded.get(1).unwrap().seq_as_string(), "ACG--");
  }
  #[test]
  fn test_pad_alignment() {
    let seq1 = AnnotatedSequence::new("id1".to_string(), vec!['A', 'C']);
    let seq2 = AnnotatedSequence::new("id2".to_string(), vec!['A', 'C', 'G']);
    let mut seqs = SequenceCollection::new();
    seqs.add(seq1).ok();
    seqs.add(seq2).ok();
    let padded: Alignment = seqs.pad_with_gaps(5);
    assert_eq!(padded.get(0).unwrap().seq_as_string(), "AC---");
    assert_eq!(padded.get(1).unwrap().seq_as_string(), "ACG--");
  }
}