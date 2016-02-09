#!/usr/bin/env python
"""
The central module for building and interacting with proteinlet databases.

Proteinlet databases are bgzipped tab-separated value files. This module does not handle
the compression! As far as this module is concerned it is reading and writing .tsv files.
"""

import logging
import re
import csv
import sys
from sixgill import dna

csv.field_size_limit(sys.maxsize)


__author__ = "Damon May"
__copyright__ = "Copyright (c) 2016 Damon May"
__license__ = "Apache 2.0"
__version__ = "0.1"

logger = logging.getLogger(__name__)

# the names of the columns in a tide-index read-peptides file
TIDEINDEX_COLUMNS = ['sequence', 'protein id']

# minimum length of a proteinlet to include
DEFAULT_MIN_PROTEINLET_AALENGTH = 10

# this regex will match fully-tryptic peptides
TRYPTIC_CLEAVAGE_REGEX = re.compile(".(?:(?<![KR](?!P)).)*")

# defaults for various parameters

# default minimum basecall quality score over any amino acid in the sequence
DEFAULT_MIN_AA_QUALSCORE = 30
# default minimum partial-ORF length (as bounded by stop codons)
DEFAULT_MIN_ORF_LENGTH = 40
# default minimum AA length of longest tryptic peptide
DEFAULT_MIN_LONGEST_PEPTIDE_LENGTH = 7
# default number of reads to keep a proteinlet
DEFAULT_MIN_READ_COUNT = 2

PROTEINLETDB_COLUMNS = ['sequence', 'length', 'min_qualscore', 'n_reads',
                        'partial_orf_length', 'nt_sequence']


def load_peptideproteins_from_tide_index(read_tide_index_file, peptides_to_keep):
    """
    Load a set of all proteins from a tide index that have a peptide in peptides_to_keep
    :param read_tide_index_file:
    :param peptides_to_keep:
    :return:
    """
    result = set()
    for row in csv.DictReader(read_tide_index_file, delimiter='\t'):
        peptide = row['sequence']
        if peptide in peptides_to_keep:
            row_proteins = [protstr.strip() for protstr in row['protein id'].split(';')]
            result.update(row_proteins)
    return result


def load_peptide_protein_map_from_tide_index(read_tide_index_file, peptides_to_keep=None):
    """
    load a map from peptide to a list of all proteins containing that peptide, from a Tide read-tide-index output
    file. Optionally, only keep peptides in a given collection
    :param read_tide_index_file:
    :param peptides_to_keep:
    :return:
    """
    result = {}
    for row in csv.DictReader(read_tide_index_file, delimiter='\t'):
        peptide = row['sequence']
        # if we have a list of peptides to keep, and this isn't on it, ignore
        if peptides_to_keep and peptide not in peptides_to_keep:
            continue
        row_proteins = [protstr.strip() for protstr in row['protein id'].split(';')]
        result[peptide] = row_proteins
    return result


def extract_read_proteinlets(ntseq, phred_qualscores, min_proteinlet_length, min_basecall_qual,
                             min_orf_length, min_longest_tryppep_length):
    """
    Extract all proteinlets from a read that match the filtering criteria. Don't check for duplicates, return
    everything. This method does all the heavy lifting.

    Return value is a list of all the proteinlets and then a tuple of accounting information for use
    in charting.

    :param ntseq: read nucleotide sequence
    :param phred_qualscores: quality scores for each base
    :param min_proteinlet_length: minimum length of a proteinlet to return
    :param min_basecall_qual: minimum minimum basecall quality of a proteinlet to return
    :param min_orf_length: minimum length of the ORF to return
    :param min_longest_tryppep_length:
    :return: a tuple of [result proteinlets], (n_discarded_tooshort, n_discarded_minqualscore,
                                n_discarded_stopcodon, n_discarded_longestpep_tooshort,
                                n_discarded_toofew_trypticsites))
    """
    ntseq_revcomp = dna.reverse_complement(ntseq)

    # if we're filtering on quality, then we need to keep track of the qual scores in both
    # forward and reverse order.
    if min_basecall_qual > 0:
        phred_qualscores_reverse = phred_qualscores[::-1]
    else:
        phred_qualscores_reverse = None

    # debug logging
    logger.debug("read: %s" % ntseq)
    logger.debug("  quals: %s" % phred_qualscores)

    # accounting
    n_discarded_tooshort = 0
    n_discarded_minqualscore = 0
    n_discarded_stopcodon = 0
    n_discarded_longestpep_tooshort = 0
    n_discarded_toofew_trypticsites = 0

    # iterate over all the frames and find the proteinlets.
    result_proteinlets = []
    for frame in xrange(0, 6):
        is_reverse = frame > 2
        nt_shift = frame % 3

        # figure out whether we're translating, and assessing quality, for the forward
        # sequence or the reverse
        seq_to_translate = ntseq
        read_qualscores_list = phred_qualscores
        # if this frame is a reverse frame, reverse ntseq and score list
        if is_reverse:
            seq_to_translate = ntseq_revcomp
            read_qualscores_list = phred_qualscores_reverse

        # calculate amino acid sequence
        aa_seq = dna.forward_translate_dna_oneframe(seq_to_translate, nt_shift, stop_before_stop=False)

        # if too short, discard
        if len(aa_seq) < min_proteinlet_length:
            n_discarded_tooshort += 1
            continue

        # no stop codons allowed
        if '*' in aa_seq:
            n_discarded_stopcodon += 1
            continue

        # get the NT sequence that covers the full read AA seq of this frame
        nt_seq_for_full_aaseq = seq_to_translate[nt_shift:nt_shift + 3 * len(aa_seq)]

        # if debug mode, check translation against aa seq
        if logger.isEnabledFor(logging.DEBUG):
            assert(dna.forward_translate_dna_firstframe(nt_seq_for_full_aaseq) == aa_seq)

        # find tryptic peptides
        tryptic_peptides = TRYPTIC_CLEAVAGE_REGEX.findall(aa_seq)
        logger.debug("        Tryptic peps: %s" % ",".join(tryptic_peptides))

        # if we don't have two internal tryptic ends, discard
        if len(tryptic_peptides) < 3:
            n_discarded_toofew_trypticsites += 1
            continue

        # move in from the left to the first tryptic end
        # move in from the right to the last tryptic end
        tryppeps_in_proteinletseq = tryptic_peptides[1:len(tryptic_peptides) - 1]
        aa_startindex = len(tryptic_peptides[0])

        # construct the proteinlet AA sequence
        proteinlet_seq = ''.join(tryppeps_in_proteinletseq)
        logger.debug("        proteinlet_seq: %s" % proteinlet_seq)

        # if too short, discard
        if len(proteinlet_seq) < min_proteinlet_length:
            n_discarded_tooshort += 1
            continue

        # check to see if any tryptic peptides are long enough
        max_tryppep_length = max([len(tryppep) for tryppep in tryppeps_in_proteinletseq])
        if max_tryppep_length < min_longest_tryppep_length:
            # if there are no tryptic peptides long enough, return None
            n_discarded_longestpep_tooshort += 1
            logger.debug("          no tryptic peps_long enough")
            continue

        # big savings if we don't have to deal with qual scores
        min_aaqual = 0
        if min_basecall_qual > 0:
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("        qualscores: %s" %
                             ','.join([str(x) for x in read_qualscores_list[aa_startindex:aa_startindex + len(proteinlet_seq)]]))
            aa_minqualscore_list = [0] * len(proteinlet_seq)
            should_stop_lowqual = False
            for proteinletseq_idx in xrange(0, len(proteinlet_seq)):
                aa_idx = aa_startindex + proteinletseq_idx
                first_codonpos = aa_idx * 3 + nt_shift
                aa_minqual = min(read_qualscores_list[first_codonpos:first_codonpos + 3])
                if aa_minqual < min_basecall_qual:
                    should_stop_lowqual = True
                    break
                aa_minqualscore_list[proteinletseq_idx] = aa_minqual
            if should_stop_lowqual:
                logger.debug("        low quality.")
                n_discarded_minqualscore += 1
                continue
            min_aaqual = min(aa_minqualscore_list)

        logger.debug("        keeping proteinlet!")

        # calculate the nucleotide sequence of this proteinlet
        proteinlet_ntseq_startpos = 3 * aa_startindex
        proteinlet_nt_seq = nt_seq_for_full_aaseq[proteinlet_ntseq_startpos:proteinlet_ntseq_startpos + 3 * len(proteinlet_seq)]

        # debug logging
        logger.debug("** %d  %d" % (aa_startindex, proteinlet_ntseq_startpos))
        logger.debug("        nt seq = %s" % proteinlet_nt_seq)
        if logger.isEnabledFor(logging.DEBUG):
            assert(len(proteinlet_seq) * 3 == len(proteinlet_nt_seq))
            if dna.forward_translate_dna_firstframe(proteinlet_nt_seq) != proteinlet_seq:
                raise ValueError("forward translation of nt seq:\n %s, does not match aa seq: \n %s" %
                                 (dna.forward_translate_dna_firstframe(proteinlet_nt_seq),
                                  proteinlet_seq))

        # add the proteinlet to the list for this read
        result_proteinlets.append(Proteinlet(proteinlet_seq, min_aaqual,
                                             len(nt_seq_for_full_aaseq), [proteinlet_nt_seq], 1))

    # return the proteinlets and accounting information for this read
    return result_proteinlets, (n_discarded_tooshort, n_discarded_minqualscore,
                                n_discarded_stopcodon, n_discarded_longestpep_tooshort,
                                n_discarded_toofew_trypticsites)


def read_proteinlets(proteinlet_file):
    """
    Read proteinlets from file.
    :param proteinlet_database_file:
    :param countsfile:
    :return:
    """
    metadata_reader = csv.DictReader(proteinlet_file, delimiter='\t')
    for row in metadata_reader:
        proteinlet_seq = row['sequence']
        min_qualscore = int(row['min_qualscore'])
        partial_orf_length = int(row['partial_orf_length'])
        nt_seqs = row['nt_sequence'].split(',')
        n_reads = int(row['n_reads'])

        proteinlet = Proteinlet(proteinlet_seq, min_qualscore,
                                partial_orf_length, nt_seqs, n_reads)
        yield proteinlet


def write_proteinlets(proteinlets, out_file):
    """
    Write proteinlets to a file
    :param proteinlets: iterator over proteinlets
    :param out_file: output file
    :return: number of proteinlets written
    """

    out_file.write('\t'.join(PROTEINLETDB_COLUMNS) + "\n")
    n_written = 0
    for proteinlet in proteinlets:
        out_file.write(proteinlet.make_output_line() + '\n')
        n_written += 1
    return n_written


def filter_proteinlets(proteinlet_generator, min_orf_len, min_aa_seq_len,
                       min_readcount, min_qualscore, min_longest_tryp_pep_len,
                       max_proteinlets=None):
    """
    Apply filtering criteria to proteinlets, yielding the ones that pass as they're
    encountered.
    :param proteinlet_generator:
    :param min_orf_len:
    :param min_aa_seq_len:
    :param min_readcount:
    :param min_qualscore:
    :param min_longest_tryp_pep_len:
    :param max_proteinlets: maximum number of proteinlets to return (for testing)
    :return:
    """
    n_yielded = 0
    for proteinlet in proteinlet_generator:
        if not proteinlet.passes_filter(min_aa_seq_len,
                                        min_orf_len,
                                        min_readcount, min_qualscore,
                                        min_longest_tryp_pep_len):
            continue
        yield proteinlet
        n_yielded += 1
        if max_proteinlets is not None and n_yielded >= max_proteinlets:
            logger.debug("filter_proteinlets stopping early, after %d." % n_yielded)
            break


class Proteinlet(object):
    """
    A class representing a proteinlet, carrying around all the information we might need later. Does not
    track individual read IDs
    """
    def __init__(self, sequence, min_qualscore, partial_orf_len, nt_seqs,
                 n_reads):
        self.sequence = sequence
        self.min_qualscore = min_qualscore
        self.partial_orf_len = partial_orf_len
        self.nt_seqs = nt_seqs
        self.n_reads = n_reads

    def calc_longest_tryppep_length(self):
        """
        Calculate the length of the longest tryptic peptide in this proteinlet
        :return:
        """
        tryp_peps = TRYPTIC_CLEAVAGE_REGEX.findall(self.sequence)
        return max([len(tryppep) for tryppep in tryp_peps])

    def passes_filter(self, min_aa_seq_length,
                      min_orf_length,
                      min_reads, min_qualscore,
                      min_longest_peptide_length):
        """
        Does this proteinlet pass a set of filtering criteria?
        :param min_aa_seq_length:
        :param min_orf_length:
        :param min_reads:
        :param min_qualscore:
        :param min_longest_peptide_length:
        :return:
        """
        if len(self.sequence) < min_aa_seq_length:
            return False
        if self.partial_orf_len < min_orf_length:
            return False
        if self.n_reads < min_reads:
            return False
        if self.min_qualscore < min_qualscore:
            return False
        if min_longest_peptide_length > 1 and self.calc_longest_tryppep_length() < min_longest_peptide_length:
            return False
        return True

    def make_output_line(self):
        """
        make a line for this proteinlet in a database file.
        does not include \n
        :return:
        """
        nt_sequence_field = ','.join(self.nt_seqs)
        return '\t'.join([self.sequence, str(len(self.sequence)),
                          str(self.min_qualscore),
                          str(self.n_reads),
                          str(self.partial_orf_len),
                          nt_sequence_field])


