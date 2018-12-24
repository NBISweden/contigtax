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
import ftplib
import hashlib
import sys
from Bio.SeqIO import parse
import gzip as gz
import shutil
from collections import defaultdict
import subprocess
import tarfile
import time


def download_ncbi_taxonomy(args):
    """Downloads the taxdump.tar.gz file from NCBI and unpacks it"""
    if not os.path.exists(args.taxdir):
        os.makedirs(args.taxdir)
    if os.path.exists("{}/nodes.dmp".format(args.taxdir)) and os.path.exists("{}/names.dmp".format(args.taxdir)):
        if args.force:
            sys.stderr.write("Re-downloading because force=True\n")
            pass
        else:
            return 0
    url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    local = "{}/taxdump.tar.gz".format(args.taxdir)
    sys.stderr.write("Downloading NCBI taxdump.tar.gz\n")
    urllib.request.urlretrieve(url, local, reporthook=progress)
    tar = tarfile.open(local, "r:gz")
    members = [tar.getmember(x) for x in ["nodes.dmp", "names.dmp"]]
    sys.stderr.write("Extracting nodes.dmp and names.dmp to {}\n".format(args.taxdir))
    tar.extractall(path=args.taxdir, members=members)
    return 0


def get_dl_status(fasta, force, skip_check):
    """Checks whether to perform download of fasta

    First the function checks if fastafile exists, then runs check_gzip to verify integrity.
    If either of these checks fails the function returns True signalling that download should proceed
    """
    if os.path.exists(fasta):
        if force:
            sys.stderr.write("Downloading because force=True\n")
            return True
        else:
            if skip_check:
                return False
            else:
                sys.stderr.write("Checking integrity of download\n")
                if check_gzip(fasta):
                    sys.stderr.write("Downloaded file is OK. Skipping download\n")
                    return False
                else:
                    sys.stderr.write("Downloaded file NOT OK. Re-downloading\n")
                    return True
    else:
        return True


def check_gzip(f):
    """Checks if a local gzipped file is ok

    Runs gunzip -t on the file using subprocess. Returns True on returncode 0.
    """
    status = subprocess.run(["gzip", "-t", f], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if status.returncode == 0:
        return True
    return False


def progress(count, blocksize, totalsize):
    """Prints progress on download to terminal"""
    # Calculate percent and do some string formatting for constant output length
    percent = round(count * blocksize * 100 / totalsize, 2)
    l, d = len(str(percent).split(".")[0]), len(str(percent).split(".")[1])
    percent = "{}{}{}".format((3-l)*" ", percent, (2-d)*"0")
    # Calculate Mb downloaded so far and format string output
    mb = round(count*blocksize/(1024**2), 2)
    total_mb = round(totalsize / 1024)
    t = len(str(total_mb))
    l, d = len(str(mb).split(".")[0]), len(str(mb).split(".")[1])
    mb = "{}{}".format(mb, (2-d)*"0")
    mb_l = len(mb)
    msg = "({mb}/{tot} Mb){trail}".format(mb=mb, tot=total_mb, trail=" "*(t+3-mb_l))
    sys.stderr.write("Downloading: {p}% complete {msg}\r".format(p=percent, mb=mb, msg=msg))
    sys.stdout.flush()


def setup_download_dirs(args):
    if not args.dldir:
        dldir = "./{db}".format(db=args.db)
    else:
        dldir = os.path.expandvars(args.dldir)
    if args.tmpdir:
        tmpdir = os.path.expandvars(args.tmpdir)
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    else:
        tmpdir = dldir
    if not os.path.exists(dldir):
        os.makedirs(dldir)
    return dldir, tmpdir


def download_fasta(args):
    """Downloads a Uniref fasta file

    If db is any of the 'UniRef' versions this function downloads a UniRef fasta file from ftp.uniprot.org/.

    If the specified db is 'nr' the latest available release of nr.gz file is downloaded from ftp.ncbi.nlm.nih.gov/.
    In this case, the protein accession mapfile is also downloaded.
    """
    dldir, tmpdir = setup_download_dirs(args)
    # Local path to target
    fasta = "{}/{}.fasta.gz".format(dldir, args.db)
    tmp_fasta = "{}/{}.fasta.gz".format(os.path.expandvars(tmpdir), args.db)
    release_note = "{}/{}.release_note".format(dldir, args.db)
    if "uniref" in args.db:
        url_base = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref"
        url = "{base}/{db}/{db}.fasta.gz".format(base=url_base, db=args.db)
        url_release_note = "{base}/{db}/{db}.release_note".format(base=url_base, db=args.db)
        urllib.request.urlretrieve(url_release_note, release_note)
    elif args.db == "nr":
        url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
        with open(release_note, 'w') as fhout:
            timestamp = time.asctime(time.gmtime())
            fhout.write("Downloaded {}\n".format(timestamp))
    # Perform checks to see whether fasta file should be downloaded
    dl = get_dl_status(fasta, args.force, args.skip_check)
    if dl:
        sys.stderr.write("Downloading {} to {}\n".format(url, tmp_fasta))
        urllib.request.urlretrieve(url, tmp_fasta, reporthook=progress)  # Download the fasta
        move(tmp_fasta, fasta)  # Move fasta from temporary directory
        sys.stderr.write("Download complete\n")
        # Check integrity
        sys.stderr.write("Checking integrity of fasta file\n")
        if check_gzip(fasta):
            sys.stderr.write("Downloaded file is OK\n")
            return 0
        else:
            sys.stderr.write("Downloaded file NOT OK!\n")
            return 1
    else:
        return 0


def move(src, target):
    if src == target:
        return
    shutil.move(src, target)


def read_relase_notes(f):
    with open(f, 'r') as fhin:
        for line in fhin:
            line = line.rstrip().lstrip()
            if line.split(":")[0] == "Number of clusters":
                n = int(line.split(":")[-1].replace(",",""))
                return n
    return -1


def zip_file(src, target):
    """Compresses an existing file"""
    with open(src, 'rb') as f_in:
        with gz.open(target, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def progress_msg(i, n):
    if i % 1000 == 0:
        if n > 0:
            p = (float(i) / n) * 100
            msg = "{}% ({}/{}) records parsed\r".format(round(p, 2), i, n)
        else:
            msg = "{} records parsed\r".format(i)
        sys.stderr.write(msg)
        sys.stdout.flush()


def parse_seqid(record, db):
    if "uniref" in db:
        newid = (record.id).split("_")[-1]  # Remove the 'UniRefxx_' string
        taxid = [x.split("=")[1] for x in (record.description).split(" ") if "TaxID" in x][0]
        return newid, taxid
    else:
        return record.id, None


def format_fasta(args):
    """Reformats a protein fasta file by stripping 'UniRef_' prefix and outputs taxidmap

    Because diamond makedb requires protein accessions to be at most 14 characters
    this function strips the 'UniRefxx_' prefix from all protein ids.

    Taxonomy Ids are also parsed from (UniRef) fasta headers and written to idmap.gz in the format used by NCBI
    for the prot.accession2taxid.gz file.
    """
    fastadir = ""
    if not os.path.exists(dldir):
        os.makedirs(dldir)
    if not args.tmpdir:
        tmpdir = dldir
    tmpdir = os.path.expandvars(tmpdir)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    release_notes = "{}/{}.release_note".format(dldir, args.db)
    if os.path.exists(release_notes):
        n = read_relase_notes(release_notes)
    else:
        n = 0
    # Fasta and reformatted fasta files
    fasta = "{}/{}.fasta.gz".format(dldir, args.db)
    fasta_r = "{}/{}.reformat.fasta.gz".format(dldir, args.db)
    _fasta_r = "{}/{}.reformat.fasta.gz".format(tmpdir, args.db)
    # Mapping file with protein accessions to taxids
    mapfile = "{}/prot.accession2taxid.gz".format(dldir)
    _mapfile = "{}/prot.accession2taxid.gz".format(tmpdir)
    # Counts file with number of proteins per taxid
    countsfile = "{}/taxoncounts.tsv.gz".format(dldir)
    _countsfile = "{}/taxoncounts.tsv.gz".format(tmpdir)
    taxoncounts = defaultdict(int)
    # IDmap dictionary for ids that are longer than 14 characters
    idmap = {}
    _idmapfile = "{}/idmap.tsv.gz".format(tmpdir)
    idmapfile = "{}/idmap.tsv.gz".format(dldir)
    j = 1
    sys.stderr.write("Reformatting {}\n".format(fasta))
    with gz.open(fasta, 'rt') as fhin, gz.open(_mapfile, 'wb') as fhidmap, gz.open(_fasta_r, 'wb') as fhreformat:
        idmap_string = "{}\t{}\t{}\t{}\n".format("accession","accession_version","taxid","gi")
        fhidmap.write(idmap_string.encode())
        for i, record in enumerate(parse(fhin, "fasta"), start=0):
            progress_msg(i, n)  # Handle progress report output
            newid, taxid = parse_seqid(record, args.db)
            if "uniref" in args.db:
                taxoncounts[taxid]+=1  # Count taxids
                idmap_string = "{}\t{}\t{}\t{}\n".format(newid, newid, taxid, newid)  # Write newid to taxid mapping
                fhidmap.write(idmap_string.encode())
            if len(newid) > args.maxidlen:
                newid = "id{j}".format(j=j)
                idmap[record.id] = newid
                j += 1
            fasta_string = ">{}\n{}\n".format(newid, str(record.seq))  # Write record with new id
            fhreformat.write(fasta_string.encode())
    if "uniref" in args.db:
        # Write taxon counts
        with gz.open(_countsfile, 'wb') as fhout:
            for taxid, taxcount in taxoncounts.items():
                s = "{}\t{}\n".format(taxid, taxcount)
                fhout.write(s.encode())
        move(_countsfile, countsfile)
        move(_mapfile, mapfile)
    move(_fasta_r, fasta_r)
    if len(idmap.keys()) > 0:
        with gz.open(_idmapfile, 'wb') as fhout:
            for old, new in idmap.items():
                line = "{}\t{}\n".format(old, new)
                fhout.write(line.encode())
        move(_idmapfile, idmapfile)
    return 0
