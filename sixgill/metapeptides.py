#!/usr/bin/env python
"""
The central module for building and interacting with metapeptide databases.

Metapeptide databases are bgzipped tab-separated value files. This module does not handle
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

# minimum length of a metapeptide to include
DEFAULT_MIN_METAPEPTIDE_AALENGTH = 10

# this regex will match fully-tryptic peptides
TRYPTIC_CLEAVAGE_REGEX = re.compile(".(?:(?<![KR](?!P)).)*")

# defaults for various parameters

# default minimum basecall quality score over any amino acid in the sequence
DEFAULT_MIN_AA_QUALSCORE = 30
# default minimum partial-ORF length (as bounded by stop codons)
DEFAULT_MIN_ORF_LENGTH = 40
# default minimum AA length of longest tryptic peptide
DEFAULT_MIN_LONGEST_PEPTIDE_LENGTH = 7
# default number of reads to keep a metapeptide
DEFAULT_MIN_READ_COUNT = 2

METAPEPTIDEDB_COLUMNS = ['sequence', 'length', 'min_qualscore',
                        'partial_orf_length', 'metagene_score', 'read_ids']

METAPEPTIDE_STATUS_OK = 0
METAPEPTIDE_STATUS_BAD_TOOSHORT = 1
METAPEPTIDE_STATUS_BAD_MINQUALSCORE = 2
METAPEPTIDE_STATUS_BAD_STOPCODON = 3
METAPEPTIDE_STATUS_BAD_LONGESTPEP_TOOSHORT = 4
METAPEPTIDE_STATUS_BAD_TOOFEW_TRYPTICSITES = 5
METAPEPTIDE_STATUS_BAD_AMBIGUOUSDNA = 6

# sentinel value that indicates that a metapeptide does not have a MetaGene score
METAGENE_SCORE_MISSING = -1.0


def extract_read_metapeptides(ntseq, phred_qualscores, min_metapeptide_length, min_basecall_qual,
                             min_orf_length, min_longest_tryppep_length, read_id):
    """
    Extract all metapeptides from a read that match the filtering criteria. Don't check for duplicates, return
    everything. This method does all the heavy lifting.

    Return value is a list of all the metapeptides and then a tuple of accounting information for use
    in charting.

    :param ntseq: read nucleotide sequence
    :param phred_qualscores: quality scores for each base
    :param min_metapeptide_length: minimum length of a metapeptide to return
    :param min_basecall_qual: minimum minimum basecall quality of a metapeptide to return
    :param min_orf_length: minimum length of the ORF to return
    :param min_longest_tryppep_length:
    :param read_id: identifier of read
    :return: a tuple of [result metapeptides], (n_discarded_tooshort, n_discarded_minqualscore,
                                n_discarded_stopcodon, n_discarded_longestpep_tooshort,
                                n_discarded_toofew_trypticsites, n_discarded_ambiguousdna))
    """
    # debug logging
    logger.debug("read: %s" % ntseq)
    logger.debug("  quals: %s" % phred_qualscores)

    # accounting
    n_discarded_tooshort = 0
    n_discarded_minqualscore = 0
    n_discarded_stopcodon = 0
    n_discarded_longestpep_tooshort = 0
    n_discarded_toofew_trypticsites = 0
    n_discarded_ambiguousdna = 0

    # iterate over all the frames and find the metapeptides.
    result_metapeptides = []
    for frame in xrange(0, 6):
        is_reverse = frame > 2
        nt_shift = frame % 3
        # add the metapeptide to the list for this read
        metapeptide, status = extract_frame_metapeptide(ntseq, phred_qualscores, min_metapeptide_length,
                                                        min_basecall_qual, min_orf_length,
                                                        min_longest_tryppep_length, is_reverse,
                                                        nt_shift, 0, len(ntseq), read_id)
        if status == METAPEPTIDE_STATUS_OK:
            result_metapeptides.append(metapeptide)
        elif status == METAPEPTIDE_STATUS_BAD_TOOSHORT:
            n_discarded_tooshort += 1
        elif status == METAPEPTIDE_STATUS_BAD_MINQUALSCORE:
            n_discarded_minqualscore += 1
        elif status == METAPEPTIDE_STATUS_BAD_STOPCODON:
            n_discarded_stopcodon += 1
        elif status == METAPEPTIDE_STATUS_BAD_LONGESTPEP_TOOSHORT:
            n_discarded_longestpep_tooshort += 1
        elif status == METAPEPTIDE_STATUS_BAD_TOOFEW_TRYPTICSITES:
            n_discarded_toofew_trypticsites += 1
        elif status == METAPEPTIDE_STATUS_BAD_AMBIGUOUSDNA:
            n_discarded_ambiguousdna += 1

    # return the metapeptides and accounting information for this read
    return result_metapeptides, (n_discarded_tooshort, n_discarded_minqualscore,
                                 n_discarded_stopcodon, n_discarded_longestpep_tooshort,
                                 n_discarded_toofew_trypticsites,
                                 n_discarded_ambiguousdna)


def extract_frame_metapeptide(ntseq, phred_qualscores, min_metapeptide_length, min_basecall_qual,
                              min_orf_length, min_longest_tryppep_length,
                              is_minus_strand, frame, startpos, endpos, read_id,
                              should_keep_with_cterm_stop=False,
                              metagene_score=METAGENE_SCORE_MISSING):
    """
    Try to extract a metapeptide from this frame.
    :param ntseq:
    :param phred_qualscores:
    :param min_metapeptide_length:
    :param min_basecall_qual:
    :param min_longest_tryppep_length:
    :param is_minus_strand:
    :param frame:
    :param startpos:
    :param endpos:
    :param read_id: identifier of the read the frame came from
    :return: metapeptide, reason.
    """

    # figure out whether we're translating, and assessing quality, for the forward
    # sequence or the reverse
    seq_to_translate = ntseq[startpos:endpos+1]
    if 'N' in seq_to_translate:
        return None, METAPEPTIDE_STATUS_BAD_AMBIGUOUSDNA
    read_qualscores_list = phred_qualscores[startpos:endpos+1]
    # if this frame is a reverse frame, reverse ntseq and score list
    if is_minus_strand:
        seq_to_translate = dna.reverse_complement(seq_to_translate)
        read_qualscores_list = list(reversed(read_qualscores_list))

    # calculate amino acid sequence

    aa_seq = dna.forward_translate_dna_oneframe(seq_to_translate, frame, stop_before_stop=False)

    # if too short, discard
    if len(aa_seq) < min_metapeptide_length:
        return None, METAPEPTIDE_STATUS_BAD_TOOSHORT

    # for MetaGene output that ends with a *
    aaseq_ends_with_stop = False
    if should_keep_with_cterm_stop:
        if aa_seq.endswith('*'):
            aa_seq = aa_seq[:-1]
            aaseq_ends_with_stop = True

    # no stop codons allowed
    if '*' in aa_seq:
        return None, METAPEPTIDE_STATUS_BAD_STOPCODON

    # get the NT sequence that covers the full read AA seq of this frame
    nt_seq_for_full_aaseq = seq_to_translate[frame:frame + 3 * len(aa_seq)]

    # if debug mode, check translation against aa seq
    if logger.isEnabledFor(logging.DEBUG):
        assert(dna.forward_translate_dna_firstframe(nt_seq_for_full_aaseq) == aa_seq)

    # find tryptic peptides
    tryptic_peptides = TRYPTIC_CLEAVAGE_REGEX.findall(aa_seq)
    logger.debug("        Tryptic peps: %s" % ",".join(tryptic_peptides))

    # by default, we require two tryptic sites. But the rightmost one is waived if
    # a stop codon occurs after the coding sequence, and the leftmost one is waived if
    # we didn't start at position 0
    n_required_tryptic_peptides = 3
    first_tryppep_idx = 1
    last_tryppep_idx_offset = 1
    if aaseq_ends_with_stop or (is_minus_strand and startpos > 0) or (not is_minus_strand and endpos < len(ntseq) - 1):
        # if this gene sequence doesn't *end* at one end or the other of the read
        n_required_tryptic_peptides -= 1
        last_tryppep_idx_offset = 0
    if (not is_minus_strand and startpos > 0) or (is_minus_strand and endpos < len(ntseq) - 1):
        # if this gene sequence doesn't *begin* at one end or the other of the read
        n_required_tryptic_peptides -= 1
        first_tryppep_idx = 0
    # if we don't have the necessary number of internal tryptic ends, discard
    if len(tryptic_peptides) < n_required_tryptic_peptides:
        return None, METAPEPTIDE_STATUS_BAD_TOOFEW_TRYPTICSITES

    # move in from the left to the first tryptic end
    # move in from the right to the last tryptic end
    tryppeps_in_metapeptideseq = tryptic_peptides[first_tryppep_idx:len(tryptic_peptides) - last_tryppep_idx_offset]
    aa_startindex = 0
    if first_tryppep_idx == 1:
        aa_startindex = len(tryptic_peptides[0])

    # construct the metapeptide AA sequence
    metapeptide_seq = ''.join(tryppeps_in_metapeptideseq)
    logger.debug("        metapeptide_seq: %s" % metapeptide_seq)

    # if too short, discard
    if len(metapeptide_seq) < min_metapeptide_length:
        return None, METAPEPTIDE_STATUS_BAD_TOOSHORT

    # check to see if any tryptic peptides are long enough
    max_tryppep_length = max([len(tryppep) for tryppep in tryppeps_in_metapeptideseq])
    if max_tryppep_length < min_longest_tryppep_length:
        # if there are no tryptic peptides long enough, return None
        logger.debug("          no tryptic peps_long enough")
        return None, METAPEPTIDE_STATUS_BAD_LONGESTPEP_TOOSHORT
    #print("@@@ %d %d %d %d" % (aa_startindex, len(metapeptide_seq), len(read_qualscores_list), len(ntseq)))
    # big savings if we don't have to deal with qual scores
    min_aaqual = 0
    if min_basecall_qual > 0:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("        qualscores: %s" %
                         ','.join([str(x) for x in read_qualscores_list[aa_startindex:aa_startindex + len(metapeptide_seq)]]))
        aa_minqualscore_list = [0] * len(metapeptide_seq)
        should_stop_lowqual = False
        for metapeptideseq_idx in xrange(0, len(metapeptide_seq)):
            aa_idx = aa_startindex + metapeptideseq_idx
            first_codonpos = aa_idx * 3 + frame
            aa_minqual = min(read_qualscores_list[first_codonpos:first_codonpos + 3])
            if aa_minqual < min_basecall_qual:
                should_stop_lowqual = True
                break
            aa_minqualscore_list[metapeptideseq_idx] = aa_minqual
        if should_stop_lowqual:
            logger.debug("        low quality.")
            return None, METAPEPTIDE_STATUS_BAD_MINQUALSCORE
        min_aaqual = min(aa_minqualscore_list)

    logger.debug("        keeping metapeptide!")

    # calculate the nucleotide sequence of this metapeptide
    metapeptide_ntseq_startpos = 3 * aa_startindex
    metapeptide_nt_seq = nt_seq_for_full_aaseq[metapeptide_ntseq_startpos:metapeptide_ntseq_startpos + 3 * len(metapeptide_seq)]

    # debug logging
    logger.debug("** %d  %d" % (aa_startindex, metapeptide_ntseq_startpos))
    logger.debug("        nt seq = %s" % metapeptide_nt_seq)
    if logger.isEnabledFor(logging.DEBUG):
        assert(len(metapeptide_seq) * 3 == len(metapeptide_nt_seq))
        if dna.forward_translate_dna_firstframe(metapeptide_nt_seq) != metapeptide_seq:
            raise ValueError("forward translation of nt seq:\n %s, does not match aa seq: \n %s" %
                             (dna.forward_translate_dna_firstframe(metapeptide_nt_seq),
                              metapeptide_seq))
    metapeptide = Metapeptide(metapeptide_seq, min_aaqual, len(nt_seq_for_full_aaseq), [read_id],
                              metagene_score)
    return metapeptide, METAPEPTIDE_STATUS_OK


def read_metapeptides(metapeptide_file):
    """
    Read metapeptides from file.
    :param metapeptide_database_file:
    :param countsfile:
    :return:
    """
    metadata_reader = csv.DictReader(metapeptide_file, delimiter='\t')
    for row in metadata_reader:
        metapeptide_seq = row['sequence']
        min_qualscore = int(row['min_qualscore'])
        partial_orf_length = int(row['partial_orf_length'])
        read_ids = row['read_ids'].split(',')
        metagene_score = METAGENE_SCORE_MISSING
        if 'metagene_score' in row.keys():
            metagene_score = float(row['metagene_score'])

        metapeptide = Metapeptide(metapeptide_seq, min_qualscore,
                                  partial_orf_length, read_ids, metagene_score)
        yield metapeptide


def write_metapeptides(metapeptides, out_file):
    """
    Write metapeptides to a file
    :param metapeptides: iterator over metapeptides
    :param out_file: output file
    :return: number of metapeptides written
    """

    out_file.write('\t'.join(METAPEPTIDEDB_COLUMNS) + "\n")
    n_written = 0
    for metapeptide in metapeptides:
        out_file.write(metapeptide.make_output_line() + '\n')
        n_written += 1
    return n_written


def filter_metapeptides(metapeptide_generator, min_orf_len, min_aa_seq_len,
                        min_readcount, min_qualscore, min_longest_tryp_pep_len,
                        min_metagene_score,
                        max_metapeptides=None):
    """
    Apply filtering criteria to metapeptides, yielding the ones that pass as they're
    encountered.
    :param metapeptide_generator:
    :param min_orf_len:
    :param min_aa_seq_len:
    :param min_readcount:
    :param min_qualscore:
    :param min_longest_tryp_pep_len:
    :param min_metagene_score
    :param max_metapeptides: maximum number of metapeptides to return (for testing)
    :return:
    """
    n_yielded = 0
    for metapeptide in metapeptide_generator:
        if not metapeptide.passes_filter(min_aa_seq_len,
                                         min_orf_len,
                                         min_readcount, min_qualscore,
                                         min_longest_tryp_pep_len,
                                         min_metagene_score):
            continue
        yield metapeptide
        n_yielded += 1
        if max_metapeptides is not None and n_yielded >= max_metapeptides:
            logger.debug("filter_metapeptides stopping early, after %d." % n_yielded)
            break


class Metapeptide(object):
    """
    A class representing a metapeptide, carrying around all the information we might need later.
    """
    def __init__(self, sequence, min_qualscore, partial_orf_len,
                 read_ids, metagene_score):
        self.sequence = sequence
        self.min_qualscore = min_qualscore
        self.partial_orf_len = partial_orf_len
        self.read_ids = read_ids
        self.metagene_score = metagene_score

    def get_readcount(self):
        """
        Count the reads
        """
        return len(self.read_ids)

    def calc_longest_tryppep_length(self):
        """
        Calculate the length of the longest tryptic peptide in this metapeptide
        :return:
        """
        tryp_peps = TRYPTIC_CLEAVAGE_REGEX.findall(self.sequence)
        return max([len(tryppep) for tryppep in tryp_peps])

    def passes_filter(self, min_aa_seq_length,
                      min_orf_length,
                      min_reads, min_qualscore,
                      min_longest_peptide_length,
                      min_metagene_score):
        """
        Does this metapeptide pass a set of filtering criteria?
        :param min_aa_seq_length:
        :param min_orf_length:
        :param min_reads:
        :param min_qualscore:
        :param min_longest_peptide_length:
        :param min_metagene_score:
        :return:
        """
        if len(self.sequence) < min_aa_seq_length:
            return False
        if self.partial_orf_len < min_orf_length:
            return False
        if self.get_readcount() < min_reads:
            return False
        if self.min_qualscore < min_qualscore:
            return False
        if min_longest_peptide_length > 1 and self.calc_longest_tryppep_length() < min_longest_peptide_length:
            return False
        if min_metagene_score > METAGENE_SCORE_MISSING and self.metagene_score < min_metagene_score:
            return False
        return True

    def make_output_line(self):
        """
        make a line for this metapeptide in a database file.
        does not include \n
        :return:
        """
        read_ids_field = ','.join(self.read_ids)
        return '\t'.join([self.sequence, str(len(self.sequence)),
                          str(self.min_qualscore),
                          str(self.partial_orf_len),
                          str(self.metagene_score),
                          read_ids_field])





