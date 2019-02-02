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


def check_args(args):
    if not os.path.exists(args.dbfile):
        sys.stderr.write("ERROR: Diamond database {} does not exist\n".format(args.dbfile))
        sys.exit()
    if not os.path.exists(args.query):
        sys.stderr.write("ERROR: Query file {} does not exist\n".format(args.query))


def check_dirs(args):
    outdir = os.path.dirname(args.outfile)
    if outdir == "":
        outdir = "."
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not args.tmpdir:
        tmpdir = outdir
    else:
        tmpdir = args.tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    return tmpdir


def diamond(args):
    """Runs diamond blastx"""
    check_args(args)
    tmpdir = check_dirs(args)
    outfile = str(args.outfile).replace(".gz", "")
    outfmt = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids slen"
    status = subprocess.run('diamond blastx -q {q} -p {p} -f 6 {f} --top {top} -e {e} -b {b} -c {c} --tmpdir {tmpdir} \
                   --more-sensitive --compress 1 -d {db} -o {out}'.format(q=args.query, p=args.threads, f=outfmt,
                                                                          top=args.top, e=args.evalue, b=args.blocksize,
                                                                          c=args.chunks, out=outfile, db=args.dbfile,
                                                                          tmpdir=tmpdir),
                            shell=True)
    return status.returncode
