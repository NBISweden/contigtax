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

import os
import sys
import logging
from argparse import ArgumentParser
from tango_prepare import download_fasta, format_fasta, download_ncbi_taxonomy, download_nr_idmap, build_diamond_db

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def download(args):
    """Handles downloading (fasta and NCBI taxonomy)"""
    if args.db == "taxonomy":
        download_ncbi_taxonomy(args)
    elif args.db == "idmap":
        download_nr_idmap(args)
    else:
        download_fasta(args)


def format(args):
    """Reformats a database to comply with diamond

    diamond makedb has some requirements when adding taxonomic info
    1. protein seqids cannot be longer than 14 characters
    2. nodes.dmp and names.dmp must be supplied
    3. a prot.accession2taxid.gz file mapping protein ids to taxonomy ids must be supplied
    """
    idmap = format_fasta(args)


def build(args):
    """Builds the diamond database from downloaded fasta and taxonomy files"""
    build_diamond_db(args)


def usage(args):
    args.print_help()


def main():
    parser = ArgumentParser()
    subparser = parser.add_subparsers(title="subcommands",
                                      description="valid subcommands")
    parser.set_defaults(func=usage)
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
    format_parser.set_defaults(func=format)  # Call format function with arguments
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
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
