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
from tango_prepare import download
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


def prepare(args):
    #logging.info("Downloading {db} database to {outdir}/".format(db=args.db, outdir=args.dldir))
    download(args.dldir, args.db)



def main():
    parser = ArgumentParser()
    subparser = parser.add_subparsers(title="subcommands",
                                      description="valid subcommands")
    # Download parser
    download_parser = subparser.add_parser("prepare", help="Download protein fasta file")
    download_parser.add_argument("dldir",
                                 help="Write files to this directory. Will be created if missing.")
    download_parser.add_argument("db", choices=["uniref100","uniref90","uniref50"],
                                 help="Database to download.")
    download_parser.set_defaults(func=prepare)
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
