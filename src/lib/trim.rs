use crate::edit_msa::EditMSA;
use crate::gapping::PadWithGaps;
use crate::seqs::SequenceAccesors;
use crate::seqs::Alignment;
use crate::seqs::AnnotatedSequence;
use crate::seqs::ApplyBufferedSequenceCollection;
use crate::seqs::SequenceCollection;
use crate::seqs::BufferedSeqCollection;

pub trait Trim<T> {
  fn trim_fixed(
    &self,
    right: usize,
    left: usize
  ) -> SequenceCollection;
  fn trim_by_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment;
  fn trim_by_terminal_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment;
}

pub struct TrimBufferedSequenceCollection{}

impl TrimBufferedSequenceCollection {
  pub fn trim_fixed(
    bsc: Box<dyn BufferedSeqCollection>,
    left: usize,
    right: usize
  ) -> ApplyBufferedSequenceCollection {
    let filter_func = move |mut s: AnnotatedSequence| {
      loop {
        s.trim_fixed(left, right);
        return vec![s]
      }
    };
    ApplyBufferedSequenceCollection::new(
      bsc,
      Box::new(filter_func)
    )
  }
}

impl Trim<SequenceCollection> for SequenceCollection {
  fn trim_fixed(
    &self,
    left: usize,
    right: usize
  ) -> SequenceCollection {
    self
      .iter()
      .cloned()
      .map(
        |mut s| {
          s.trim_fixed(left, right);
          s
        }
      )
      .collect::<SequenceCollection>()
  }

  fn trim_by_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment {
    let padded: SequenceCollection = self.pad_with_gaps_to_max_length();
    padded.to_msa().unwrap().trim_by_gaps(right, left)
  }

  fn trim_by_terminal_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment {
    let padded: SequenceCollection = self.pad_with_gaps_to_max_length();
    padded.to_msa().unwrap().trim_by_terminal_gaps(right, left)
  }
}

impl Trim<Alignment> for Alignment {
  fn trim_fixed(
    &self,
    right: usize,
    left: usize
  ) -> SequenceCollection {
    self
      .seqs
      .iter()
      .cloned()
      .map(
        |mut s| {
          s.trim_fixed(left, right);
          s
        }
      )
      .collect::<SequenceCollection>()
  }

  fn trim_by_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment {
    let ncols = self.length();
    let seq_len = self.length();
    let mut gapped_cols_left = vec![false; ncols];
    let mut gapped_cols_right = vec![false; ncols];
    let gapped_columns = self
      .columns()
      .map(
        |col| col.contains(&&'-') || col.contains(&&'.')
      )
      .collect::<Vec<_>>();
    if left {
      for (i, gap) in (0..seq_len).zip(gapped_columns.iter()) {
        match (i, gap) {
          (0, x) => gapped_cols_left[0] = *x,
          (_, false) => (),
          (i, true) => gapped_cols_left[i] = gapped_cols_left[i-1],
        };
      }
    }
    if right {
      for (i, gap) in (0..seq_len).zip(gapped_columns.iter()).rev() {
        match(i, gap) {
          (j, x) if j+1 == seq_len => gapped_cols_right[j] = *x,
          (_, false) => (),
          (j, true) => gapped_cols_right[j] = gapped_cols_right[j+1],
        }
      }
    }
    let cols_to_drop = gapped_cols_left
      .iter()
      .zip(gapped_cols_right)
      .map(|(a, b)| *a || b)
      .enumerate()
      .filter(|(_, x)| *x)
      .map(|(i, _)| i)
      .collect::<Vec<_>>();
    let mut new_aln = self
      .iter()
      .cloned()
      .collect::<SequenceCollection>()
      .to_msa()
      .unwrap();
    new_aln.remove_columns(cols_to_drop).unwrap();
    new_aln
  }

  fn trim_by_terminal_gaps(
    &self,
    right: bool,
    left: bool
  ) -> Alignment {
    let mut terminal_gaps_left: usize = 0;
    let mut terminal_gaps_right: usize = 0;
    let length = self.length();
    for seq in self.iter() {
      let seq = seq.seq().unwrap();
      let mut left_gaps = seq.len();
      let mut right_gaps = seq.len();
      for (i, c) in seq.iter().enumerate() {
        match c {
          '-' | '.' => {}
          _ => {
            left_gaps = i;
            break
          },
        }
      }
      terminal_gaps_left = std::cmp::max(terminal_gaps_left, left_gaps);
      for (i, c) in seq.iter().rev().enumerate() {
        match c {
          '-' | '.' => {}
          _ => {
            right_gaps = i;
            break
          },
        }
      }
      terminal_gaps_right = std::cmp::max(terminal_gaps_right, right_gaps);
    }
    let terminal_gaps = match (left, right) {
      (true, true) | (false, false) => {
        let mut terminal_gaps = vec![];
        for i in 0..terminal_gaps_left {
          terminal_gaps.push(i);
        }
        for i in (length-terminal_gaps_right)..length {
          terminal_gaps.push(i);
        }
        terminal_gaps
      },
      (true, false) => {
        let mut terminal_gaps = vec![];
        for i in 0..terminal_gaps_left {
          terminal_gaps.push(i);
        }
        terminal_gaps
      },
      (false, true) => {
        let mut terminal_gaps = vec![];
        for i in (length-terminal_gaps_right)..length {
          terminal_gaps.push(i);
        }
        terminal_gaps
      }
    };
    let mut new_aln = self.clone();
    new_aln
      .remove_columns(terminal_gaps)
      .unwrap();
    new_aln
  }
}

#[cfg(test)]
mod test {
  use crate::seqs::AnnotatedSequence;
  use crate::seqs::SequenceCollection;
  use crate::seqs::SequenceAccesors;

  use super::Trim;

  #[test]
  fn test_trim_seqcol_fixed() {
    let mut sq = SequenceCollection::new();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_01".to_string(),
            "ACTATCGTCA".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "sequence_02".to_string(),
            "ACTATCGT".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "Sequence_03".to_string(),
            "ACTATCG".to_string()
        )
    ).unwrap();
    sq.add(
        AnnotatedSequence::from_string(
            "Sequence_04".to_string(),
            "ACTA".to_string()
        )
    ).unwrap();
    let result = sq.trim_fixed(2, 3);
    let seqs = result
      .into_iter()
      .map(|x| x.seq_as_string())
      .collect::<Vec<_>>();
    assert_eq!(
      seqs,
      vec!["TATCG", "TAT", "TA", ""]
    );
  }
  #[test]
  fn test_trim_aln_by_gaps_no_gaps() {
    let mut sq = SequenceCollection::new();
      sq.add(
          AnnotatedSequence::from_string(
              "sequence_01".to_string(),
              "ACTATCGTCA".to_string()
          )
      ).unwrap();
      sq.add(
          AnnotatedSequence::from_string(
              "sequence_02".to_string(),
              "ACTATCGTCA".to_string()
          )
      ).unwrap();
      sq.add(
          AnnotatedSequence::from_string(
              "Sequence_03".to_string(),
              "ACTATCGTCA".to_string()
          )
      ).unwrap();
      sq.add(
          AnnotatedSequence::from_string(
              "Sequence_04".to_string(),
              "ACTATCGTCA".to_string()
          )
      ).unwrap();
    let msa = sq.clone().to_msa().unwrap();
    let new_msa = msa.trim_by_gaps(true, true);
    assert_eq!(
      new_msa
        .seqs
        .iter()
        .map(|x| x.seq_as_string().len())
        .collect::<Vec<_>>(),
        vec![10, 10, 10, 10]
    );
    assert!(
      !new_msa
        .seqs
        .iter()
        .map(|x| x.seq().unwrap().iter())
        .flatten()
        .any(|c| *c == '-' || *c == '.')
    )
  }
  #[test]
  fn test_trim_aln_by_gaps_with_gaps() {
    let msa = vec![
        ("1", "ABCDEFG"),
        ("2", "-BCDE--"),
        ("3", "--CDEF-"),
      ]
      .into_iter()
      .collect::<SequenceCollection>()
      .to_msa()
      .unwrap();
    let new_aln = msa.trim_by_gaps(true, true);
    assert_eq!(
      new_aln
        .seqs
        .iter()
        .map(|x| x.seq_as_string().len())
        .collect::<Vec<_>>(),
        vec![3, 3, 3]
    );
    assert!(
      !new_aln
        .seqs
        .iter()
        .map(|x| x.seq().unwrap().iter())
        .flatten()
        .any(|c| *c == '-' || *c == '.')
    )
  }
  #[test]
  fn test_trim_aln_by_gaps_with_gaps_in_the_middle() {
    let msa = vec![
        ("1", "ABC---G"),
        ("2", "ABCD---"),
        ("3", "ABCDEF-"),
      ]
      .into_iter()
      .collect::<SequenceCollection>()
      .to_msa()
      .unwrap();
    let new_aln = msa.trim_by_gaps(true, true);
    assert_eq!(
      new_aln
        .seqs
        .iter()
        .map(|x| x.seq_as_string().len())
        .collect::<Vec<_>>(),
        vec![3, 3, 3]
    );
  }
  #[test]
  fn test_trim_aln_by_terminal_gaps_with_gaps_in_the_middle() {
    let msa = vec![
        ("1", "ABC---G--"),
        ("2", "ABCD-----"),
        ("3", "-BCDEF---"),
      ]
      .into_iter()
      .collect::<SequenceCollection>()
      .to_msa()
      .unwrap();
    let new_aln = msa.trim_by_terminal_gaps(true, true);

    assert_eq!(
      new_aln
        .seqs
        .iter()
        .map(|x| x.seq_as_string().len())
        .collect::<Vec<_>>(),
        vec![3, 3, 3]
    );
  }
}
