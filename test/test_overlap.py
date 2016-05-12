from unittest import TestCase

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from overlap import find_overlaps


class OverlapTestCase(TestCase):
  def test_empty_seqs(self):
    self.assertEqual(find_overlaps(Seq(''), Seq('')), [])
    self.assertEqual(find_overlaps(Seq(''), Seq('A')), [])
    self.assertEqual(find_overlaps(Seq('G'), Seq('')), [])

  def test_nonoverlapping(self):
    self.assertEqual(find_overlaps(Seq('AGCT', generic_dna), Seq('GAGCT', generic_dna)), [])
    self.assertEqual(find_overlaps(Seq('AGCT', generic_dna), Seq('G', generic_dna)), [])
    self.assertEqual(find_overlaps(Seq('AGCT', generic_dna), Seq('GTCA', generic_dna)), [])

  def test_overlapping_single(self):
    self.assertEqual(
      map(str, find_overlaps(Seq('AGCT', generic_dna), Seq('T', generic_dna))),
      ['T'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCT', generic_dna), Seq('TG', generic_dna))),
      ['T'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCT', generic_dna), Seq('TGC', generic_dna))),
      ['T'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCT', generic_dna), Seq('CT', generic_dna))),
      ['CT'],
    )

  def test_overlapping_multiple(self):
    self.assertEqual(
      map(str, find_overlaps(Seq('AGCTCTAA', generic_dna), Seq('ACTCTGAG', generic_dna))),
      ['A'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCTCTAA', generic_dna), Seq('AACTCTGAG', generic_dna))),
      ['A', 'AA'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCTCTAAA', generic_dna), Seq('AACTCTGAG', generic_dna))),
      ['A', 'AA'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCTCTAAA', generic_dna), Seq('AAACTCTGAG', generic_dna))),
      ['A', 'AA', 'AAA'],
    )

    self.assertEqual(
      map(str, find_overlaps(Seq('AGCTCT', generic_dna), Seq('CTCTGAG', generic_dna))),
      ['CT', 'CTCT'],
    )
