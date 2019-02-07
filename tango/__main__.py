#!/usr/bin/env python
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

from argparse import ArgumentParser
from tango import prepare
from tango import search
from tango import assign
from time import time
import sys
import os


def download(args):
    """Handles downloading (fasta and NCBI taxonomy)"""
    if args.db == "taxonomy":
        prepare.download_ncbi_taxonomy(args)
        prepare.init_sqlite_taxdb(args.taxdir)
    elif args.db == "idmap":
        prepare.download_nr_idmap(args)
    else:
        prepare.download_fasta(args)


def reformat(args):
    """Reformats a database to comply with diamond

    diamond makedb has some requirements when adding taxonomic info
    1. protein seqids cannot be longer than 14 characters
    2. nodes.dmp and names.dmp must be supplied
    3. a prot.accession2taxid.gz file mapping protein ids to taxonomy ids must be supplied
    """
    prepare.format_fasta(args)


def update(args):
    """Updates the protein -> taxid idmap

    If protein accessions were too long these were remapped during the format stage. This requires a new
    map file to be created with the new (shorter) protein ids.
    """
    prepare.update_idmap(args)


def build(args):
    """Builds the diamond database from downloaded fasta and taxonomy files"""
    prepare.build_diamond_db(args)


def run_diamond(args):
    """Runs diamond"""
    start_time = time()
    sys.stderr.write("Running diamond blastx with {} threads\n".format(args.threads))
    search.diamond(args)
    end_time = time()
    run_time = round(end_time - start_time, 1)
    sys.stderr.write("Total time: {}s\n".format(run_time))


def assign_taxonomy(args):
    """Parses diamond results and assigns taxonomy"""
    # Check outfile
    if os.path.isdir(args.outfile):
        sys.exit("\nERROR: Outfile {} is a directory\n".format(args.outfile))
    outdir = os.path.dirname(args.outfile)
    if outdir != "" and not os.path.exists(outdir):
        os.makedirs(outdir)
    start_time = time()
    sys.stderr.write("Assigning taxonomy with {} threads\n".format(args.threads))
    # Parse diamond file
    assign.parse_hits(args)
    end_time = time()
    run_time = round(end_time - start_time, 1)
    sys.stderr.write("Total time: {}s\n".format(run_time))


def get_version():
    from tango import __version__
    return '%(prog)s {version}'.format(version=__version__)


def usage(args):
    if args.version:
        import tango
        print(tango.version)
    else:
        print("""
        To print help message: tango -h
        """)


def main():
    parser = ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version=get_version())
    parser.set_defaults(func=usage)
    subparser = parser.add_subparsers(title="subcommands",
                                      description="valid subcommands")
    # Download parser
    download_parser = subparser.add_parser("download", help="Download fasta file and NCBI taxonomy files")
    download_parser.add_argument("db", choices=["uniref100", "uniref90", "uniref50", "nr", "taxonomy", "idmap"],
                                 default="uniref100", help="Database to download. Defaults to 'uniref100'")
    download_parser.add_argument("-d", "--dldir",
                                 help="Write files to this directory. Defaults to db name in current directory. \
                                 Will be created if missing.")
    download_parser.add_argument("--tmpdir", type=str,
                                 help="Temporary directory for downloading files")
    download_parser.add_argument("-t", "--taxdir", default="./taxonomy",
                                 help="Directory to store NCBI taxdump files. \
                                Defaults to 'taxonomy/' in current directory")
    download_parser.add_argument("-f", "--force", action="store_true",
                                 help="Overwrite downloaded files")
    download_parser.add_argument("--skip_check", action="store_true",
                                 help="Skip check of downloaded fasta file. Default: False")
    download_parser.set_defaults(func=download)  # Call download function with arguments
    # Format parser
    format_parser = subparser.add_parser("format", help="Format fasta file for diamond and create protein2taxid map")
    format_parser.add_argument("fastafile", type=str,
                               help="Specify protein fasta to reformat")
    format_parser.add_argument("reformatted", type=str,
                               help="Path to reformatted fastafile")
    format_parser.add_argument("-f", "--force", action="store_true",
                               help="Force overwrite of existing reformatted fastafile")
    format_parser.add_argument("--forceidmap", action="store_true",
                               help="Force overwrite of existing accession2taxid mapfile")
    format_parser.add_argument("-m", "--taxidmap", type=str,
                               help="Protein accession to taxid mapfile. For UniRef this file is created from \
                                    information in the fasta headers and stored in a file named prot.accession2taxid.gz\
                                    in the same directory as the reformatted fasta file. Specify another path here.")
    format_parser.add_argument("--maxidlen", type=int, default=14,
                               help="""Maximum allowed length of sequence ids. Defaults to 14 (required by diamond \
                                     for adding taxonomy info to database). Ids longer than this are written to \
                                     a file with the original id""")
    format_parser.add_argument("--tmpdir", type=str,
                               help="Temporary directory for writing fasta files")
    format_parser.set_defaults(func=reformat)  # Call format function with arguments
    # Update parser
    update_parser = subparser.add_parser("update", help="Update protein to taxid map file with new sequence ids")
    update_parser.add_argument("taxonmap", type=str,
                               help="Existing prot.accession2taxid.gz file")
    update_parser.add_argument("idfile", type=str,
                               help="File mapping long sequence ids to new ids")
    update_parser.add_argument("newfile", type=str,
                               help="Updated mapfile")
    update_parser.set_defaults(func=update)  # Call update function with arguments
    # Build parser
    build_parser = subparser.add_parser("build", help="Build diamond database from downloaded files")
    build_parser.add_argument("fastafile", help="Specify (reformatted) fasta file")
    build_parser.add_argument("taxonmap", help="Protein accession to taxid mapfile (must be gzipped)")
    build_parser.add_argument("taxonnodes", help="nodes.dmp file from NCBI taxonomy database")
    build_parser.add_argument("-d", "--dbfile", help="Name of diamond database file. Defaults to diamond.dmnd in \
                              same directory as the protein fasta file")
    build_parser.add_argument("-p", "--threads", type=int, default=1,
                              help="Number of threads to use when building (defaults to 1)")
    build_parser.set_defaults(func=build)  # Call build function with arguments
    # Search parser
    search_parser = subparser.add_parser("search", help="Run diamond blastx with nucleotide fasta file")
    search_parser.add_argument("query", type=str,
                               help="Query contig nucleotide file")
    search_parser.add_argument("dbfile", type=str,
                               help="Diamond database file")
    search_parser.add_argument("outfile", type=str,
                               help="Diamond output file")
    search_parser.add_argument("-p", "--threads", default=1, type=int,
                               help="Threads to allocate for diamond")
    search_parser.add_argument("-b", "--blocksize", type=float, default=2.0,
                               help="Sequence block size in billions of letters (default=2.0). Set to 20 on clusters")
    search_parser.add_argument("-c", "--chunks", type=int, default=4,
                               help="Number of chunks for index processing (default=4)")
    search_parser.add_argument("--top", type=int, default=10,
                               help="report alignments within this percentage range of top alignment score (default=10)")
    search_parser.add_argument("-e", "--evalue", default=0.001, type=float,
                               help="maximum e-value to report alignments (default=0.001)")
    search_parser.add_argument("-t", "--tmpdir", type=str,
                               help="directory for temporary files")
    search_parser.set_defaults(func=run_diamond)
    # Assign parser
    assign_parser = subparser.add_parser("assign", help="Assigns taxonomy from diamond output")
    assign_parser.add_argument("diamond_results", type=str,
                               help="Diamond blastx results")
    assign_parser.add_argument("outfile", type=str,
                               help="Output file")
    assign_parser.add_argument("-m", "--mode", type=str, default="rank_lca", choices=['rank_lca', 'rank_vote', 'score'],
                               help="Mode to use for parsing taxonomy: 'rank_lca' (default), 'rank_vote' or 'score'")
    assign_parser.add_argument("-t", "--taxdir", type=str, default="./taxonomy",
                               help="Directory specified during 'tango download taxonomy'. Defaults to taxonomy/.")
    assign_parser.add_argument("-T", "--top", type=int, default=10,
                               help="Top percent of best score to consider hits for")
    assign_parser.add_argument("-r", "--ranks", nargs="*", default=["superkingdom", "phylum", "class", "order",
                                                                    "family", "genus", "species"])
    assign_parser.add_argument("-p", "--threads", type=int, default=1,
                               help="Number of threads to use. Defaults to 1.")
    assign_parser.add_argument("-c", "--chunksize", type=int, default=1,
                               help="Size of chunks sent to process pool. For large input files using a large chunksize\
                                can make the job complete much faster than using the default value of 1.")
    assign_parser.add_argument("--rank_thresholds", nargs="*", default=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
    assign_parser.add_argument("--vote_threshold", default=0.5, type=float,
                               help="Minimum fraction required when voting on rank assignments.")
    assign_parser.add_argument("--blobout", type=str,
                               help="Output hits.tsv table compatible with blobtools")
    assign_parser.set_defaults(func=assign_taxonomy)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
