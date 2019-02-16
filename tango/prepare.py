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
from ete3 import NCBITaxa
import urllib.request
import sys
from Bio.SeqIO import parse
import gzip as gz
import shutil
from pathlib import Path
import subprocess
import tarfile
import time
import tqdm
from glob import glob


def progress(count, blocksize, totalsize):
    """Prints progress on download to terminal"""
    # Calculate percent and do some string formatting for constant output length
    percent = round(count * blocksize * 100 / totalsize, 2)
    l, d = len(str(percent).split(".")[0]), len(str(percent).split(".")[1])
    percent = "{}{}{}".format((3-l)*" ", percent, (2-d)*"0")
    # Calculate Mb downloaded so far and format string output
    mb = round(count*blocksize/(1024**2), 2)
    total_mb = round(totalsize / (1024**2), 2)
    t = len(str(total_mb))
    l, d = len(str(mb).split(".")[0]), len(str(mb).split(".")[1])
    mb = "{}{}".format(mb, (2-d)*"0")
    mb_l = len(mb)
    msg = "({mb}/{tot} Mb){trail}".format(mb=mb, tot=total_mb, trail=" "*(t+3-mb_l))
    sys.stderr.write("Downloading: {p}% complete {msg}\r".format(p=percent, mb=mb, msg=msg))
    sys.stdout.flush()

def my_hook(t):
    """Wraps tqdm instance.
    https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    Don't forget to close() or __exit__()
    the tqdm instance once you're done with it (easiest using `with` syntax).
    -------
    """
    last_b = [0]

    def update_to(b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b
    return update_to


def setup_download_dirs(args):
    """Creates directories for downloading"""
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


def build_diamond_db(args):
    """Runs diamond makedb

    Builds a diamond database with taxonomy id mappings"""
    fastadir = os.path.dirname(args.fastafile)
    if not args.dbfile:
        dbfile = "{}/diamond.dmnd".format(fastadir)
    else:
        dbfile = args.dbfile
    subprocess.run('gunzip -c {fastafile} | diamond makedb --taxonmap {taxonmap} --taxonnodes {taxonnodes} \
                          -d {dbfile} -p {threads}'.format(fastafile=args.fastafile, taxonmap=args.taxonmap,
                                                           taxonnodes=args.taxonnodes, dbfile=dbfile,
                                                           threads=args.threads), shell=True)
    return 0


def download_nr_idmap(args):
    """Downloads the nr protein to taxonomy id map"""
    if not args.dldir:
        dldir = "./nr"
    else:
        dldir = args.dldir
    if not os.path.exists(dldir):
        os.makedirs(dldir)
    if not args.tmpdir:
        tmpdir = dldir
    else:
        tmpdir = args.tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    if os.path.exists(os.path.join(dldir,"prot.accession2taxid.gz")):
        if args.force:
            sys.stderr.write("Re-downloading idmap because force=True\n")
            pass
        else:
            sys.stderr.write("{}/prot.accession2taxid.gz already exists. Force download with '-f'\n".format(dldir))
            return 0
    url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
    _local = "{}/prot.accession2taxid.gz".format(tmpdir)
    local = "{}/prot.accession2taxid.gz".format(dldir)
    sys.stderr.write("Downloading prot.accession2taxid.gz to {}\n".format(dldir))
    with tqdm.tqdm(ncols=100, unit="bytes") as t:
        reporthook = my_hook(t)
        urllib.request.urlretrieve(url, _local, reporthook=reporthook)
    move(_local, local)
    return 0


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
    with tqdm.tqdm(ncols=100, unit=" bytes") as t:
        reporthook = my_hook(t)
        urllib.request.urlretrieve(url, local, reporthook=reporthook)
    tar = tarfile.open(local, "r:gz")
    members = [tar.getmember(x) for x in ["nodes.dmp", "names.dmp"]]
    sys.stderr.write("Extracting nodes.dmp and names.dmp to {}\n".format(args.taxdir))
    tar.extractall(path=args.taxdir, members=members)
    return 0


def init_sqlite_taxdb(taxdir, dbname):
    """Creates ete3 sqlite database"""
    taxdump_file = os.path.join(taxdir,"taxdump.tar.gz")
    dbfile = os.path.join(taxdir, dbname)
    if not os.path.exists(taxdump_file) or os.path.exists(dbfile):
        taxdump_file = None
    if not os.path.exists(taxdir):
        os.makedirs(taxdir)
    Path(dbfile).touch()
    ncbi_taxa = NCBITaxa(dbfile=dbfile, taxdump_file=taxdump_file)
    return ncbi_taxa


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
        idmap_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        with open(release_note, 'w') as fhout:
            timestamp = time.asctime(time.gmtime())
            fhout.write("Downloaded {}\n".format(timestamp))
    # Perform checks to see whether fasta file should be downloaded
    dl = get_dl_status(fasta, args.force, args.skip_check)
    if dl:
        sys.stderr.write("Downloading {} to {}\n".format(url, tmp_fasta))
        with tqdm.tqdm(ncols=100, unit=" bytes") as t:
            reporthook = my_hook(t)
            urllib.request.urlretrieve(url, tmp_fasta, reporthook=reporthook)  # Download the fasta
        if args.db == "nr" and not args.skip_idmap:
            _idmap = "{}/prot.accession2taxid.gz".format(os.path.expandvars(tmpdir))
            idmap = "{}/prot.accession2taxid.gz".format(dldir)
            with tqdm.tqdm(ncols=100, unit=" bytes") as t:
                reporthook = my_hook(t)
                urllib.request.urlretrieve(idmap_url, _idmap, reporthook=reporthook)  # Download the nr idmap
            move(_idmap, idmap)
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
    """Moves file if source and target are different"""
    if src == target:
        return
    shutil.move(src, target)


def read_relase_notes(f):
    """Attempts to extract total number of sequences from relase notes"""
    with open(f, 'r') as fhin:
        for line in fhin:
            line = line.rstrip().lstrip()
            if line.split(":")[0] == "Number of clusters":
                n = int(line.split(":")[-1].replace(",",""))
                return n
    return -1


def parse_seqid(record):
    """Parses sequence and taxids from fasta"""
    taxid = None
    if "UniRef" in record.id:
        db_type = "uniref"
        newid = (record.id).split("_")[-1]  # Remove the 'UniRefxx_' string
        try:
            taxid = [x.split("=")[1] for x in (record.description).split(" ") if "TaxID" in x][0]
        except IndexError:
            taxid = None
    else:
        db_type = "non_uniref"
        newid = record.id
    return newid, taxid, db_type


def setup_format_dirs(args):
    """Checks what directories and files to create for formatting"""
    fastadir = os.path.dirname(args.fastafile)
    formatdir = os.path.dirname(args.reformatted)
    if not os.path.exists(formatdir) and formatdir != "":
        os.makedirs(formatdir)
    if not args.tmpdir:
        tmpdir = formatdir
    else:
        tmpdir = args.tmpdir
    tmpdir = os.path.expandvars(tmpdir)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    rn_glob = glob("{}/*.release_note".format(fastadir))
    if len(rn_glob) > 0:
        n = read_relase_notes(rn_glob[0])
    else:
        n = 0
    return formatdir, tmpdir, n


def update_idmap(args):
    """Updates taxonmap ids

    If some sequence ids were too long (e.g. longer than 14 characters) the taxonmap file needs to be updated
    with the shorter sequence ids.
    """
    idmap = {}
    sys.stderr.write("Reading idmap from {}\n".format(args.idfile))
    with gz.open(args.idfile, 'rt') as fhin:
        for line in fhin:
            old, new = line.rstrip().split("\t")
            idmap[old] = new
    sys.stderr.write("Creating updated taxonmap {} from {}\n".format(args.newfile, args.taxonmap))
    with gz.open(args.taxonmap, 'rt') as fhin, gz.open(args.newfile, 'wb') as fhout:
        s = "{}\t{}\t{}\t{}\n".format("accession", "accession.version", "taxid", "gi")
        fhout.write(s.encode())
        for i, line in enumerate(tqdm.tqdm(fhin, ncols=100, total=None, unit=" lines")):
            if i == 0:
                continue
            line = line.rstrip()
            acc, accv, taxid, gi = line.split("\t")
            try:
                new = idmap[acc]
                s = "{}\t{}\t{}\t{}\n".format(new, new, taxid, new)
            except KeyError:
                s = "{}\n".format(line)
            fhout.write(s.encode())
    return 0


def write_idmap(s, fhidmap):
    if fhidmap:
        fhidmap.write(s.encode())
    return


def format_fasta(args):
    """Reformats a protein fasta file and outputs taxidmap

    Because diamond makedb requires protein accessions to be at most 14 characters
    this function forces sequence ids to conform to this length.

    Taxonomy Ids are also parsed from (UniRef) fasta headers and written to idmap.gz in the format used by NCBI
    for the prot.accession2taxid.gz file.
    """
    if os.path.exists(args.reformatted) and not args.force:
        sys.stderr.write("{} already exists. Specify '-f' to force overwrite\n".format(args.reformatted))
        return 0
    formatdir, tmpdir, n = setup_format_dirs(args)
    if n < 0:
        N = None
    else:
        N = n
    # Fasta and reformatted fasta files
    _fasta_r = "{}/{}".format(tmpdir, os.path.basename(args.reformatted))
    # Mapping file with protein accessions to taxids
    if not args.taxidmap:
        mapfile = "{}/prot.accession2taxid.gz".format(formatdir)
        _mapfile = "{}/prot.accession2taxid.gz".format(tmpdir)
    else:
        mapfile = args.taxidmap
        _mapfile = "{}/{}".format(tmpdir, os.path.basename(args.taxidmap))
    if args.forceidmap or not os.path.exists(mapfile):
        fhidmap = gz.open(_mapfile, 'wb')
    else:
        sys.stderr.write("{} already exists.\n".format(mapfile))
        fhidmap = False
    # IDmap dictionary for ids that are longer than --maxidlen characters
    idmap = {}
    idmapfile = "{}/idmap.tsv.gz".format(formatdir)
    _idmapfile = "{}/idmap.tsv.gz".format(tmpdir)
    j = 1
    uniref = False
    sys.stderr.write("Reformatting {}\n".format(args.fastafile))
    with gz.open(args.fastafile, 'rt') as fhin, gz.open(_fasta_r, 'wb') as fhreformat:
        idmap_string = "{}\t{}\t{}\t{}\n".format("accession", "accession.version", "taxid", "gi")
        write_idmap(idmap_string, fhidmap)
        for i, record in enumerate(tqdm.tqdm(parse(fhin, "fasta"), unit=" records", ncols=100, total=N)):
            newid, taxid, format = parse_seqid(record)
            if len(newid) > args.maxidlen:
                newid = "id{j}".format(j=j)
                idmap[record.id] = newid
                j += 1
            # Skip UniRef entries with no taxid
            if format == "uniref":
                if not taxid:
                    continue
                idmap_string = "{}\t{}\t{}\t{}\n".format(newid, newid, taxid, newid)  # Write newid to taxid mapping
                write_idmap(idmap_string, fhidmap)
            fasta_string = ">{}\n{}\n".format(newid, str(record.seq))  # Write record with new id
            fhreformat.write(fasta_string.encode())
    sys.stderr.write("{}/{} records parsed\n".format(i, i))
    # Close idmap file handle
    fhidmap.close()
    # If uniref database the mapfile has been created during reformatting
    if uniref:
        move(_mapfile, mapfile)
    # Move reformmated fasta
    move(_fasta_r, args.reformatted)
    # Check if there are remapped sequences
    if len(idmap.keys()) > 0:
        with gz.open(_idmapfile, 'wb') as fhout:
            for old, new in idmap.items():
                line = "{}\t{}\n".format(old, new)
                fhout.write(line.encode())
        move(_idmapfile, idmapfile)
        sys.stderr.write("{}/{} sequence ids longer than {} characters written to {}\n".format(len(idmap), i,
                                                                                               args.maxidlen, idmapfile))
        if not uniref:
            sys.stderr.write(
                ">>>Make sure to update your prot.accession2taxid.gz file with updated ids. See 'tango update'<<<\n")

    return 0
