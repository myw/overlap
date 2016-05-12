from pprint import pprint

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def get_seqs():
  seqs = {
    'lambda': ['TATCACCGCAAGGGATA', 'AATAC', 'TAACACCGTGCGTGTTG', 'ACTATTT', 'TACCTCTGGCGGTGATA'],
    '434': ['ACAAGAAAACTGT', 'ATTTGCAA', 'ACAAGATACATTGT', 'ATGAAAT', 'ACAAGAAAGTTTGT'],
  }

  for phage_name, seq_strings in seqs.iteritems():
    seqs[phage_name] = [Seq(seq_string, generic_dna) for seq_string in seq_strings]

  return seqs


def find_overlaps(seq1, seq2):
  overlaps = []

  for ix in xrange(1, len(seq2) + 1):
    subseq = seq2[0:ix]

    if seq1.endswith(subseq):
      overlaps.append(subseq)

  return tuple(overlaps)


def pull_out_motifs(seqs):
  """ The motifs are in between filler sequences; we usually don't care about those """
  return seqs[0::2]


def calculate_potential_strands(seq_map):
  motifs = {
    phage_name: pull_out_motifs(seqs) for phage_name, seqs in seq_map.iteritems()
  }

  # find the complement strand to the lambda sequences
  motifs['lambda'] = [seq.complement() for seq in motifs['lambda']]

  # can start with either lambda or 434 sites upstream; order must be preserved
  overlaps = {
    'lambda': [],
    '434': [],
  }

  for ix, motif_lambda in enumerate(motifs['lambda']):
    motif_434 = motifs['434'][ix]

    # starting overlap: lambda site upstream, 434 site downstream
    overlap_lambda = [find_overlaps(motif_lambda, motif_434)]

    # ending overlap: 434 site upstream, next lambda site downstream
    if ix + 1 < len(motifs['lambda']):
      overlap_lambda.append(find_overlaps(motif_434, motifs['lambda'][ix + 1]))
    else:
      overlap_lambda.append(())

    overlaps['lambda'].append(overlap_lambda)

    # starting overlap: 434 site upstream, lambda site downstream
    overlap_434 = [find_overlaps(motif_434, motif_lambda)]

    # ending overlap: lambda site upstream, next 434 site downstream
    if ix + 1 < len(motifs['434']):
      overlap_434.append(find_overlaps(motif_lambda, motifs['434'][ix + 1]))
    else:
      overlap_434.append(())

    overlaps['434'].append(overlap_434)

  return overlaps


def main():
  pprint(calculate_potential_strands(get_seqs()))


if __name__ == '__main__':
  main()
