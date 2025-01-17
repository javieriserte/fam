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
}
