# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import subprocess
import os
import sys
import tqdm
from tempfile import NamedTemporaryFile as ntf


def check_args(dbfile, query):
    """Checks that diamond database and query file exist"""
    if not os.path.exists(dbfile):
        sys.stderr.write("ERROR: Diamond database {} does not exist\n".format(dbfile))
        sys.exit()
    if not os.path.exists(query):
        sys.stderr.write("ERROR: Query file {} does not exist\n".format(query))
        sys.exit()


def check_dirs(outfile, tmpdir):
    """Sets up output and temporary directory if they do not exist"""
    outdir = os.path.dirname(outfile)
    if outdir == "":
        outdir = "."
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not tmpdir:
        tmpdir = outdir
    else:
        tmpdir = tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    return tmpdir


def filter_seqs_by_len(infile, outfile, minlen):
    """Filters sequences by length

    Parameters
    ----------
    infile: str
        Sequence file in fasta format
    outfile: str
        File in fasta format containing sequences longer than minlen
    minlen: int
        Minimum size of sequences to keep
    """

    from Bio.SeqIO import parse, write
    i = 0
    for record in tqdm.tqdm(parse(infile, 'fasta'), unit=" sequences", ncols=100, desc="Filtering sequences"):
        if len(record) >= minlen:
            write(record, outfile, "fasta")
            i += 1
    sys.stderr.write("{} sequences longer than {} written to {}\n".format(i, minlen, outfile))


def diamond(query, outfile, dbfile, mode="blastx", cpus=1, evalue=0.001, top=10, blocksize=2.0, chunks=4,
            tmpdir=False, minlen=False):
    """Runs diamond blast with query file

    This is a wrapper function for running diamond aligner in either blastx or blastp mode (default is blastx). Some
    initial checks on input is performed as well as optional filtering of input sequences by length.

    Parameters
    ----------
    query: str
        Fasta file with sequences to query against database
    outfile: str
        Path to output file with diamond results
    dbfile: str
        Path to diamond database to use for searching
    mode: str
        Mode to use for diamond, either blastx or blastp
    cpus: int
        Number of cpus to use for diamond
    evalue: float
        Maximum allowed e-value to report hits for
    top: int
        Keep hits within top percent of best scoring hit
    blocksize: float
        Sequence block size in billions of letters (default=2.0). Set to 20 on clusters.
    chunks: int
        Number of chunks for index processing (default=4). Setting to one will shorten runtime but increase memory
        requirements.
    tmpdir: str
        Temporary directory for output
    minlen: int
        Minimum length for input sequences.
    """

    # Make sure that diamond database and query file exist
    check_args(dbfile, query)
    # Make sure tmpdir exists if specified
    tmpdir = check_dirs(outfile, tmpdir)
    # Set cpus to minimum allowed
    cpus = max(1, cpus)
    # Filter input sequences if minlen is specified
    if minlen:
        q = ntf(mode="w", delete=False)
        filter_seqs_by_len(query, q.name, minlen)
        query = q.name
    outfile = str(outfile).replace(".gz", "")
    # Use blast tabular output with taxonomy id in last column
    outfmt = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids"
    status = subprocess.run('diamond {m} -q {q} -p {p} -f 6 {f} --top {top} -e {e} -b {b} -c {c} --tmpdir {tmpdir} \
                   --more-sensitive --compress 1 -d {db} -o {out}'.format(m=mode,
                                                                          q=query,
                                                                          p=cpus,
                                                                          f=outfmt,
                                                                          top=top,
                                                                          e=evalue,
                                                                          b=blocksize,
                                                                          c=chunks,
                                                                          out=outfile,
                                                                          db=dbfile,
                                                                          tmpdir=tmpdir),
                            shell=True)
    if minlen:
        sys.stderr.write("Removing temporary file {}\n".format(q.name))
        q.close()
        os.remove(q.name)
    return status.returncode
