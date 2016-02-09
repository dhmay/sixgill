#!/usr/bin/env python
""" Utilities for DNA manipulation """

import logging

__author__ = "Damon May"
__copyright__ = "Copyright (c) 2016 Damon May"
__license__ = "Apache 2.0"
__version__ = "0.1"

log = logging.getLogger(__name__)

COMPLEMENT_DICT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

# this is copied from http://www.petercollingridge.co.uk/python-bioinformatics-tools/codon-table
DNA_BASES = ['T', 'C', 'A', 'G']
CODONS = [a + b + c for a in DNA_BASES for b in DNA_BASES for c in DNA_BASES]
AMINO_ACIDS_FORTRANSLATE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS_FORTRANSLATE))


def forward_translate_dna_firstframe(dna_sequence, stop_before_stop=False):
    """forward-translate the DNA sequence in Frame 1.
    If stop_before_stop, end with the last AA before the first in-frame stop codon
    :param dna_sequence:
    :param stop_before_stop:
    :return:
    """
    result = ""
    pos = 0
    while pos <= len(dna_sequence) - 3:
        aa = CODON_TABLE[dna_sequence[pos:pos + 3]]
        if stop_before_stop and aa == '*':
            break
        result = result + CODON_TABLE[dna_sequence[pos:pos + 3]]
        pos += 3
    return result


def forward_translate_dna_oneframe(dna_sequence, frame, stop_before_stop=False):
    """
    Frames are 0-based, 0-5
    :param dna_sequence:
    :param frame:
    :param stop_before_stop:
    :return:
    """
    if frame < 0 or frame > 5:
        raise ValueError('forward_translate_dna_oneframe with invalid frame %d' % frame)
    seq_to_translate = dna_sequence
    if frame > 2:
        seq_to_translate = reverse_complement(dna_sequence)
    frame %= 3
    if frame > 0:
        seq_to_translate = seq_to_translate[frame:]
    return forward_translate_dna_firstframe(seq_to_translate, stop_before_stop=stop_before_stop)


def reverse_complement(dna_sequence):
    return complement(reverse(dna_sequence))


def reverse(dna_sequence):
    """simple string reversal"""
    return dna_sequence[::-1]


def complement(dna_sequence):
    result_list = list()
    for nt in dna_sequence:
        if nt == '-':
            result_list.append(nt)
        else:
            result_list.append(COMPLEMENT_DICT.get(nt.upper(), ''))

    return "".join(result_list)

