from nose.tools import assert_equal
import sixgill
import os
from subprocess import call
import shutil
import hashlib

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_DATAFILE_DIRNAME = os.path.join(THIS_DIR, 'datafiles')
TEST_OUTPUT_DIRNAME = os.path.join(THIS_DIR, 'test_output')

TEST_DB_NOMETAGENE_FILENAME = os.path.join(TEST_OUTPUT_DIRNAME, 'testdb_nometagene.metapeptides.gz')
TEST_DB_METAGENE_FILENAME = os.path.join(TEST_OUTPUT_DIRNAME, 'testdb_metagene.metapeptides.gz')
TEST_DATAFILE_SMALLFASTQGZ_FILENAME = os.path.join(TEST_DATAFILE_DIRNAME, 'small.fq.gz')
TEST_METAGENE_FILENAME = os.path.join(TEST_DATAFILE_DIRNAME, 'metagene_output.txt')
TEST_DB_NOMETAGENE_METAGENE_MERGED_FILENAME = os.path.join(TEST_OUTPUT_DIRNAME,
                                                           'metagene_nometagene_merged.metapeptides.gz')
TEST_DB_NOMETAGENE_METAGENE_MERGED_MIN2READS_FILENAME = os.path.join(TEST_OUTPUT_DIRNAME,
                                                                     'metagene_nometagene_merged.min2reads.metapeptides.gz')
TEST_AA_FASTA_FILENAME = os.path.join(TEST_OUTPUT_DIRNAME,
                                      'testdb_metagene.metapeptides.fasta')

MD5_BUILD_NOMETAGENE = 'b715d9bb3ab938787e035f716be25b0c'
MD5_BUILD_METAGENE = '774eba0107c2565c62ba9af4cf24912b'
MD5_MERGE = '88707bd1c66d8bb1209c43a7862c8bed'
MD5_FILTER = '93efd321f18aaff426036f34a1ccd383'
MD5_MAKEFASTA_AA = '9740ad1737954d2bf22b65b4b88b0381'


def setup_module():
    print "SETUP."
    print("THIS_DIR=%s" % THIS_DIR)
    print("TEST_DATAFILE_DIR=%s" % TEST_DATAFILE_DIRNAME)
    try:
        os.makedirs(TEST_OUTPUT_DIRNAME)
    except OSError:
        None
    print("created test output dir %s" % TEST_OUTPUT_DIRNAME)
    print("SETUP DONE.")


def teardown_module():
    print "TEAR DOWN."
    for the_file in os.listdir(TEST_OUTPUT_DIRNAME):
        file_path = os.path.join(TEST_OUTPUT_DIRNAME, the_file)
        os.unlink(file_path)
    print("TEAR DOWN DONE.")


def test_build_nometagene():
    print "Build test without metagene..."
    command = "sixgill_build --minreadcount=1 --out=%s %s" % (TEST_DB_NOMETAGENE_FILENAME,
                                                              TEST_DATAFILE_SMALLFASTQGZ_FILENAME)
    print("Running build command:")
    print(command)
    call(command.split(' '))
    assert_equal(md5(TEST_DB_NOMETAGENE_FILENAME), MD5_BUILD_NOMETAGENE)
    print "Build test without metagene complete."


def test_build_metagene():
    print "Build test with metagene..."
    command = "sixgill_build --minreadcount=1 --metagenefile=%s --out=%s %s" % (TEST_METAGENE_FILENAME,
                                                               TEST_DB_METAGENE_FILENAME,
                                                               TEST_DATAFILE_SMALLFASTQGZ_FILENAME)
    print("Running build command:")
    print(command)
    call(command.split(' '))
    assert_equal(md5(TEST_DB_METAGENE_FILENAME), MD5_BUILD_METAGENE)
    print "Build test with metagene complete."


def test_merge():
    print "Merge test"
    command = "sixgill_merge --out=%s %s %s" % (TEST_DB_NOMETAGENE_METAGENE_MERGED_FILENAME,
                                                TEST_DB_NOMETAGENE_FILENAME,
                                                TEST_DB_METAGENE_FILENAME)
    print("Running merge command:")
    print(command)
    call(command.split(' '))
    assert_equal(md5(TEST_DB_NOMETAGENE_METAGENE_MERGED_FILENAME), MD5_MERGE)
    print "Merge test complete."


def test_filter():
    print "Filter test"
    command = "sixgill_filter --out=%s --minreadcount=2 %s" % (TEST_DB_NOMETAGENE_METAGENE_MERGED_MIN2READS_FILENAME,
                                                               TEST_DB_NOMETAGENE_METAGENE_MERGED_FILENAME)
    print("Running filter command:")
    print(command)
    call(command.split(' '))
    assert_equal(md5(TEST_DB_NOMETAGENE_METAGENE_MERGED_MIN2READS_FILENAME), MD5_FILTER)
    print "Filter test complete."


def test_makefasta_aa():
    print("Make fasta test")
    command = "sixgill_makefasta --type=aa --out=%s %s" % (TEST_AA_FASTA_FILENAME,
                                                           TEST_DB_METAGENE_FILENAME)
    print("Running makefasta command:")
    print(command)
    call(command.split(' '))
    assert_equal(md5(TEST_AA_FASTA_FILENAME), MD5_MAKEFASTA_AA)
    print "Filter test complete."


def md5(fname):
    # calculate an MD5 sum on a file
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
