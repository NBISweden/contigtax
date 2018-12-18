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
import urllib.request
import sys


def download(dldir, db):
    def progress(count, blocksize, totalsize):
        percent = round(count * blocksize * 100 / totalsize, 2)
        mb = round(count*blocksize/1024**2, 2)
        total_mb = round(totalsize*blocksize/1024**2)
        sys.stderr.write("Downloading: {p}% complete ({t}/{tot}) Mb)\r".format(p=percent, t=mb, tot=total_mb))
        sys.stdout.flush()

    if not os.path.exists(dldir):
        os.makedirs(dldir)
    url = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/{db}/{db}.fasta.gz".format(db=db)
    sys.stderr.write("Downloading {} to {}\n".format(url, dldir))
    urllib.request.urlretrieve(url, os.path.join(dldir,"{db}.fasta.gz".format(db=db)), reporthook=progress)
