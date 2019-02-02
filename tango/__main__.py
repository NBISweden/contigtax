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
from prepare import download_fasta, download_ncbi_taxonomy, download_nr_idmap
from prepare import format_fasta, update_idmap, build_diamond_db
from search import diamond

__version__ = '0.0.1'


def download(args):
    """Handles downloading (fasta and NCBI taxonomy)"""
    if args.db == "taxonomy":
        download_ncbi_taxonomy(args)
    elif args.db == "idmap":
        download_nr_idmap(args)
    else:
        download_fasta(args)


def reformat(args):
    """Reformats a database to comply with diamond

    diamond makedb has some requirements when adding taxonomic info
    1. protein seqids cannot be longer than 14 characters
    2. nodes.dmp and names.dmp must be supplied
    3. a prot.accession2taxid.gz file mapping protein ids to taxonomy ids must be supplied
    """
    format_fasta(args)


def update(args):
    update_idmap(args)


def build(args):
    """Builds the diamond database from downloaded fasta and taxonomy files"""
    build_diamond_db(args)


def run_diamond(args):
    """Runs diamond"""
    diamond(args)


def usage(args):
    if args.version:
        print(__version__)
    else:
        print("""
        To print help message: tango -h
        """)


def main():
    parser = ArgumentParser()
    parser.add_argument("-v", "--version", action="store_true")
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
    assign_parser.add_argument("diamond_results")
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
