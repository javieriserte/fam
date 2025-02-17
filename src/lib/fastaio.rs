use crate::seqs::BufferedSeqCollectionFromRead;
use crate::seqs::{
  AnnotatedSequence,
  SequenceAccesors,
  SequenceCollection,
  BufferedSeqCollection
};
use std::io::Error;
use std::io::Write;

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;

use std::io::{self};

pub struct FastaReaderFromLines {
  local_lines: Vec<String>,
  end_of_input: bool
}

impl FastaReaderFromLines {
  pub fn new() -> FastaReaderFromLines {
    FastaReaderFromLines {
      local_lines: vec![],
      end_of_input: false
    }
  }

  pub fn add_line(&mut self, line: Option<String>) {
    if let Some(l) = line {
      self.local_lines.push(l);
    } else {
      self.end_of_input = true;
    }
  }

  fn build(&mut self, end:usize) -> Option<AnnotatedSequence> {
    if end <= 1 {
      return None;
    }
    let id = self.local_lines[0].split_off(1);
    let seq = self.local_lines[1..end].join("");
    self.local_lines.drain(0..end);
    Some(
      AnnotatedSequence::from_string(id, seq)
    )
  }

  pub fn first_id_index(&self) -> usize {
    self
      .local_lines[1..]
      .iter()
      .position(|x| x.starts_with(">"))
      .unwrap_or(self.local_lines.len()-1)
  }

  pub fn try_build(&mut self) -> Option<AnnotatedSequence> {
    match self.local_lines.last() {
      Some(line) if line.starts_with( ">") => {
        let end = self.first_id_index();
        self.build(end+1)
      }
      Some(_) if self.end_of_input => {
        let end = self.first_id_index();
        self.build(end+1)
      }
      Some(_) => None,
      None => None
    }
  }

  pub fn consumed(&self) -> bool {
    self.local_lines.is_empty() && self.end_of_input
  }
}

pub fn sequence_collection_from_bufread<T: BufRead>(
  mut reader: T,
) -> Result<SequenceCollection, Error> {
  let mut msa = SequenceCollection::new();
  let mut current_sequence: Vec<String> = vec![];
  let mut current_id: Option<String> = None;
  loop {
    let mut line = String::new();
    let len = reader.read_line(&mut line)?;
    if len == 0 || line.starts_with(">") {
      if let Some(id) = &current_id {
        let ann_seq = AnnotatedSequence::from_string(
          id.clone(),
          current_sequence.join(""),
        );
        msa.add(ann_seq)?
      };
      current_sequence.clear();
      if line.starts_with(">") {
        line.truncate(line.trim_end().len());
        current_id = Some(line.split_off(1));
      }
    } else {
      line.truncate(line.trim_end().len());
      current_sequence.push(line);
    }
    if len == 0 {
      break;
    }
  }
  Ok(msa)
}

pub fn sequence_collection_from_file(
    path: &Path,
) -> Result<SequenceCollection, Error> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    sequence_collection_from_bufread(reader)
}

pub fn sequence_collection_from_stdin() -> Result<SequenceCollection, Error> {
    sequence_collection_from_bufread(io::stdin().lock())
}

pub fn buffered_sequence_collection_from_stdin()
  -> Result<BufferedSeqCollectionFromRead, Error> {
    let buffer = io::stdin().lock();
    Ok(BufferedSeqCollectionFromRead::new(Box::new(buffer)))
}

pub fn buffered_sequence_collection_from_file(
  path: &Path
) -> Result<BufferedSeqCollectionFromRead, Error> {
    let f = File::open(path)?;
    let reader = BufReader::new(f);
    Ok(BufferedSeqCollectionFromRead::new(Box::new(reader)))
}

pub fn write_sequence_collection<T1: SequenceAccesors, T2: Write>(
    seqs: &T1,
    writer: T2,
) -> Result<(), io::Error> {
    let mut bw = BufWriter::new(writer);
    for annseq in seqs.iter() {
        bw.write_fmt(format_args!(">{}\n", annseq.id()))?;
        bw.write_fmt(format_args!("{}\n", annseq.seq_as_string()))?;
    }
    Ok(())
}

pub fn write_buffered_sequence_collection<T1: BufferedSeqCollection, T2: Write>(
    seqs: &T1,
    writer: T2,
) -> Result<(), io::Error> {
    let mut bw = BufWriter::new(writer);
    loop {
        let s = seqs.next_sequence();
        match s {
            Some(s) => {
                bw.write_fmt(format_args!(">{}\n", s.id()))?;
                bw.write_fmt(format_args!("{}\n", s.seq_as_string()))?;
            }
            None => break,
        }
    }
    Ok(())
}

mod test {
    #[allow(unused_imports)]
    use crate::fastaio::sequence_collection_from_bufread;
    #[allow(unused_imports)]
    use crate::seqs::{SequenceAccesors, SequenceCollection};
    #[test]
    fn test_sequence_collection_from_bufread_when_ok() {
        let input = ">S1\r\nATCTCG\n>S2\nTCT\nCGA\r\n>S3\nATG\r\nTAG";
        let b = input.as_bytes();
        let msa = sequence_collection_from_bufread(b).unwrap();
        assert_eq!(msa.get(0).unwrap().id(), "S1");
        assert_eq!(
            *msa.get(0).unwrap().seq().unwrap(),
            vec!['A', 'T', 'C', 'T', 'C', 'G']
        );
        assert_eq!(msa.get(1).unwrap().id(), "S2");
        assert_eq!(
            *msa.get(1).unwrap().seq().unwrap(),
            vec!['T', 'C', 'T', 'C', 'G', 'A']
        );
        assert_eq!(msa.get(2).unwrap().id(), "S3");
        assert_eq!(
            *msa.get(2).unwrap().seq().unwrap(),
            vec!['A', 'T', 'G', 'T', 'A', 'G']
        );
    }
    #[test]
    fn test_fastareader_from_lines_all_together() {
        // FIXME: Make tests
        use crate::fastaio::FastaReaderFromLines;
        let mut fr = FastaReaderFromLines::new();
        fr.add_line(Some(">S1".to_string()));
        fr.add_line(Some("ATCTCG".to_string()));
        fr.add_line(Some(">S2".to_string()));
        fr.add_line(Some("TCT".to_string()));
        fr.add_line(Some("CGA".to_string()));
        fr.add_line(Some(">S3".to_string()));
        fr.add_line(Some("ATG".to_string()));
        fr.add_line(Some("TAG".to_string()));
        fr.add_line(None);
        assert_eq!(fr.try_build().unwrap().id(), "S1");
        assert_eq!(fr.try_build().unwrap().id(), "S2");
        assert_eq!(fr.try_build().unwrap().id(), "S3");
        assert_eq!(fr.try_build(), None);
    }#[test]
    fn test_fastareader_from_lines_one_by_one() {
        // FIXME: Make tests
        use crate::fastaio::FastaReaderFromLines;
        let mut fr = FastaReaderFromLines::new();
        fr.add_line(Some(">S1".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some("ATCTCG".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some(">S2".to_string()));
        assert_eq!(fr.try_build().unwrap().id(), "S1");
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some("TCT".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some("CGA".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some(">S3".to_string()));
        assert_eq!(fr.try_build().unwrap().id(), "S2");
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some("ATG".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(Some("TAG".to_string()));
        assert_eq!(fr.try_build(), None);
        assert_eq!(fr.consumed(), false);
        fr.add_line(None);
        assert_eq!(fr.try_build().unwrap().id(), "S3");
        assert_eq!(fr.consumed(), true);
    }
}
