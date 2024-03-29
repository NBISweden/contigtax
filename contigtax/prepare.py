from os import makedirs
from os.path import exists, basename, expandvars, dirname, join as opj
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


def setup_download_dirs(dldir, db, tmpdir):
    """Creates directories for downloading"""
    if not dldir:
        dldir = "./{db}".format(db=db)
    else:
        dldir = expandvars(dldir)
    if tmpdir:
        tmpdir = expandvars(tmpdir)
        if not exists(tmpdir):
            makedirs(tmpdir)
    else:
        tmpdir = dldir
    if not exists(dldir):
        makedirs(dldir)
    return dldir, tmpdir


def get_dl_status(fasta, force, skip_check):
    """Checks whether to perform download of fasta

    First the function checks if fastafile exists, then runs check_gzip to
    verify integrity. If either of these checks fails the function returns
    True signalling that download should proceed
    """

    if exists(fasta):
        if force:
            sys.stderr.write("Downloading because force=True\n")
            return True
        else:
            if skip_check:
                return False
            else:
                sys.stderr.write("Checking integrity of download\n")
                if check_gzip(fasta):
                    sys.stderr.write("Downloaded file is OK. Skipping "
                                     "download\n")
                    return False
                else:
                    sys.stderr.write("Downloaded file NOT OK. "
                                     "Re-downloading\n")
                    return True
    else:
        return True


def check_gzip(f):
    """Checks if a local gzipped file is ok

    Runs gunzip -t on the file using subprocess. Returns True on returncode 0.
    """

    status = subprocess.run(["gzip", "-t", f], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    if status.returncode == 0:
        return True
    return False


def build_diamond_db(fastafile, taxonmap, taxonnodes, dbfile, cpus=1):
    """Creates a diamond database

    The diamond database is built with taxonomic information taken from the
    taxonmap file in the format of the NCBI prot.accession2taxid.gz file:

    accession	accession.version	taxid	gi
    Q6GZX4	Q6GZX4	1	Q6GZX4
    Q6GZX3	Q6GZX3	10486	Q6GZX3

    Parameters
    ----------
    fastafile: str
        Input fastafile
    taxonmap: str
        File mapping protein ids to taxids
    taxonnodes: str
        nodes.dmp file from NCBI taxonomy ftp (download with contigtax download
        taxonomy
    dbfile: str
        diamond database to create
    cpus: int
        number of cpus to use
    """
    from contigtax import diamond_legacy
    if diamond_legacy():
        tmap_string = ""
    else:
        tmap_string = "--taxonmap {} --taxonnodes {}".format(taxonmap,
                                                             taxonnodes)
    fastadir = dirname(fastafile)
    if not dbfile:
        dbfile = "{}/diamond.dmnd".format(fastadir)
    else:
        dbfile = dbfile
    p = subprocess.run("diamond makedb --in {fastafile} {tmap}"
                       " -d {dbfile} -p {cpus}".format(fastafile=fastafile,
                                                       tmap=tmap_string,
                                                       dbfile=dbfile,
                                                       cpus=cpus), shell=True)
    p.check_returncode()


def download_nr_idmap(dldir="./nr", tmpdir=False, force=False):
    """Downloads the nr protein to taxonomy id map"""
    if not dldir:
        dldir = "./nr"
    else:
        dldir = dldir
    if not exists(dldir):
        makedirs(dldir)
    if not tmpdir:
        tmpdir = dldir
    else:
        tmpdir = tmpdir
        if not exists(tmpdir):
            makedirs(tmpdir)
    if exists(opj(dldir, "prot.accession2taxid.gz")):
        if force:
            sys.stderr.write("Re-downloading idmap because force=True\n")
            pass
        else:
            sys.stderr.write(
                "{}/prot.accession2taxid.gz already exists. Force download "
                "with '-f'\n".format(dldir))
            return 0
    url_base = "ftp://ftp.ncbi.nih.gov/pub/taxonomy"
    url = "{}/accession2taxid/prot.accession2taxid.gz".format(url_base)
    _local = "{}/prot.accession2taxid.gz".format(tmpdir)
    local = "{}/prot.accession2taxid.gz".format(dldir)
    sys.stderr.write(
        "Downloading prot.accession2taxid.gz to {}\n".format(dldir))
    with tqdm.tqdm(ncols=100, unit="bytes") as t:
        reporthook = my_hook(t)
        urllib.request.urlretrieve(url, _local, reporthook=reporthook)
    move(_local, local)
    return 0


def download_ncbi_taxonomy(taxdir="taxonomy", force=False):
    """
    Downloads the taxdump.tar.gz file from NCBI and unpacks it in the
    specified directory
    """
    if not exists(taxdir):
        makedirs(taxdir)
    if exists("{}/nodes.dmp".format(taxdir)) and exists(
            "{}/names.dmp".format(taxdir)):
        if force:
            sys.stderr.write("Re-downloading because force=True\n")
            pass
        else:
            return 0
    url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    local = "{}/taxdump.tar.gz".format(taxdir)
    sys.stderr.write("Downloading NCBI taxdump.tar.gz\n")
    with tqdm.tqdm(ncols=100, unit=" bytes") as t:
        reporthook = my_hook(t)
        urllib.request.urlretrieve(url, local, reporthook=reporthook)
    tar = tarfile.open(local, "r:gz")
    members = [tar.getmember(x) for x in ["nodes.dmp", "names.dmp"]]
    sys.stderr.write(
        "Extracting nodes.dmp and names.dmp to {}\n".format(taxdir))
    tar.extractall(path=taxdir, members=members)
    return 0


def init_sqlite_taxdb(taxdir, dbname, force=False):
    """Creates ete3 sqlite database"""
    taxdump_file = opj(taxdir, "taxdump.tar.gz")
    dbfile = opj(taxdir, dbname)
    if not exists(taxdump_file) or exists(dbfile) and not force:
        taxdump_file = None
    if not exists(taxdir):
        makedirs(taxdir)
    Path(dbfile).touch()
    ncbi_taxa = NCBITaxa(dbfile=dbfile, taxdump_file=taxdump_file)
    return ncbi_taxa


def download_fasta(dldir, db, tmpdir=False, force=False, skip_check=False,
                   skip_idmap=False):
    """Downloads a protein fasta file

    If db is any of the 'UniRef' versions this function downloads a UniRef
    fasta file from ftp.uniprot.org/.

    If the specified db is 'nr' the latest available release of nr.gz file is
    downloaded from ftp.ncbi.nlm.nih.gov/.
    In this case, the protein accession mapfile is also downloaded.

    Relase notes are also created in the same directory used to store the
    download

    Parameters
    ----------
    dldir: str
        Directory to use for download
    db: str
        Database to download. One of uniref100, uniref90, uniref50, nr.
    tmpdir: str
        Temporary directory to use for download.
    force: bool
        Should existing files be overwritten?
    skip_check: bool
        Skip checking integrity of download
    skip_idmap: bool
        Skip download of idmap file (only applies to nr)

    Returns
    -------
    Exit code 0 or 1
    """

    dldir, tmpdir = setup_download_dirs(dldir, db, tmpdir)
    # Local path to target
    fasta = "{}/{}.fasta.gz".format(dldir, db)
    tmp_fasta = "{}/{}.fasta.gz".format(expandvars(tmpdir), db)
    release_note = "{}/{}.release_note".format(dldir, db)
    _ = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid"
    idmap_url = "{}/prot.accession2taxid.gz".format(_)
    if "uniref" in db:
        url_base = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref"
        url = "{base}/{db}/{db}.fasta.gz".format(base=url_base, db=db)
        url_release_note = "{base}/{db}/{db}." \
                           "release_note".format(base=url_base, db=db)
        urllib.request.urlretrieve(url_release_note, release_note)
    elif db == "nr":
        url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
        with open(release_note, 'w') as fhout:
            timestamp = time.asctime(time.gmtime())
            fhout.write("Downloaded {}\n".format(timestamp))
    # Perform checks to see whether fasta file should be downloaded
    dl = get_dl_status(fasta, force, skip_check)
    if dl:
        sys.stderr.write("Downloading {} to {}\n".format(url, tmp_fasta))
        with tqdm.tqdm(ncols=100, unit=" bytes") as t:
            reporthook = my_hook(t)
            # Download the fasta
            urllib.request.urlretrieve(url, tmp_fasta, reporthook=reporthook)
        if db == "nr" and not skip_idmap:
            _idmap = "{}/prot.accession2taxid.gz".format(
                expandvars(tmpdir))
            idmap = "{}/prot.accession2taxid.gz".format(dldir)
            sys.stderr.write("Downloading taxidmap "
                             "from {}\n".format(idmap_url))
            with tqdm.tqdm(ncols=100, unit=" bytes") as t:
                reporthook = my_hook(t)
                # Download the nr idmap
                urllib.request.urlretrieve(idmap_url, _idmap,
                                           reporthook=reporthook)
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
                n = int(line.split(":")[-1].replace(",", ""))
                return n
    return -1


def parse_seqid(r):
    """Parses sequence and taxids from fasta

    UniRef records have the format 'UniRefxx_Q6GZX3' where xx is one of 100,
    90 or 50. The fasta description for these entries contain the lowest
    common taxonomic id for all sequences belonging to each uniref id:

    Example:
    >UniRef50_Q6GZX3 protein n=51 Tax=Iridoviridae TaxID=10486 RepID=002L_FRG3G

    This function trims the 'UniRefxx' part and extracts the taxonomy id from
    the fasta description.

    If the format doesn't correspond to UniRef (e.g. when the records are from
    'nr' instead) the original sequence id
    is returned instead.

    Parameters
    ----------
    r: Bio.SeqRecord.SeqRecord
        Sequence record as parsed using Bio.SeqIO.parse

    Returns
    -------
    newid: str
        Sequence identifier string
    taxid: str
        Taxonomy identifier string
    db_type: str
        Text string showing whether sequence is from a UniRef fasta file or not
    """
    taxid = None
    if "UniRef" in r.id:
        db_type = "uniref"
        newid = (r.id).split("_")[-1]  # Remove the 'UniRefxx_' string
        # Try to extract the taxid
        try:
            items = (r.description).split(" ")
            taxid = [x.split("=")[1] for x in items if "TaxID" in x][0]
        except IndexError:
            taxid = None
    else:
        db_type = "non_uniref"
        newid = r.id
    return newid, taxid, db_type


def setup_format_dirs(fasta, reformat, tmpdir):
    """Checks what directories and files to create for formatting"""
    fastadir = dirname(fasta) if dirname(fasta) else "."
    formatdir = dirname(reformat) if dirname(reformat) else "."
    if not exists(formatdir):
        makedirs(formatdir)
    if not tmpdir:
        tmpdir = formatdir
    else:
        tmpdir = tmpdir
    tmpdir = expandvars(tmpdir)
    if not exists(tmpdir) and tmpdir != "":
        makedirs(tmpdir)
    rn_glob = glob("{}/*.release_note".format(fastadir))
    if len(rn_glob) > 0:
        n = read_relase_notes(rn_glob[0])
    else:
        n = 0
    return formatdir, tmpdir, n


def update_idmap(idfile, taxonmap, newfile):
    """Updates taxonmap ids

    If some sequence ids were too long (e.g. longer than 14 characters) the
    taxonmap file needs to be updated
    with the shorter sequence ids.
    """
    idmap = {}
    sys.stderr.write("Reading idmap from {}\n".format(idfile))
    with gz.open(idfile, 'rt') as fhin:
        for line in fhin:
            old, new = line.rstrip().split("\t")
            idmap[old] = new
    sys.stderr.write(
        "Creating updated taxonmap {} from {}\n".format(newfile, taxonmap))
    with gz.open(taxonmap, 'rt') as fhin, gz.open(newfile, 'wb') as fhout:
        s = "{}\t{}\t{}\t{}\n".format("accession", "accession.version",
                                      "taxid", "gi")
        fhout.write(s.encode())
        for i, line in enumerate(
                tqdm.tqdm(fhin, ncols=100, total=None, unit=" lines")):
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
    """Writes output to taxidmap file"""
    if fhidmap:
        fhidmap.write(s.encode())
    return


def format_fasta(fastafile, reformatted, tmpdir=False, force=False,
                 taxidmap=False, forceidmap=False, maxidlen=14):
    """Reformats a protein fasta file and outputs taxidmap

    The reformatting part includes stripping the 'UniRefxx' prefix from
    UniRef sequence ids. Taxonomy Ids are also parsed from (UniRef) fasta
    headers and written to idmap.gz in the format used by NCBI for the
    prot.accession2taxid.gz file.

    Because diamond makedb requires protein accessions to be at most 14
    characters this function forces sequence ids to conform to this length by
    default.

    Parameters
    ----------
    fastafile: str
        Path to original fastafile
    reformatted: str
        Path to reformatted output
    force: bool
        Should existing reformatted output be overwritten?
    tmpdir: str
        Temporary directory to write reformatted fasta to
    taxidmap: str
        Path to store sequence->taxid map information. Defaults to
        'prot.accession2taxid.gz' in same directory as
        reformatted output.
    forceidmap: bool
        Should existing taxidmap be overwritten?
    maxidlen: int
        Maximum length of sequence ids.

    Returns
    -------
    Exit code 0 upon successful completion
    """

    if exists(reformatted) and not force:
        sys.stderr.write("{} already exists. Specify '-f' to "
                         "force overwrite\n".format(reformatted))
        return 0
    formatdir, tmpdir, n = setup_format_dirs(fastafile, reformatted, tmpdir)
    if n < 0:
        N = None
    else:
        N = n
    # Fasta and reformatted fasta files
    _fasta_r = "{}/{}".format(tmpdir, basename(reformatted))
    # Mapping file with protein accessions to taxids
    if not taxidmap:
        mapfile = "{}/prot.accession2taxid.gz".format(formatdir)
        _mapfile = "{}/prot.accession2taxid.gz".format(tmpdir)
    else:
        mapfile = taxidmap
        _mapfile = "{}/{}".format(tmpdir, basename(taxidmap))
    if forceidmap or not exists(mapfile):
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
    sys.stderr.write("Reformatting {}\n".format(fastafile))
    with gz.open(fastafile, 'rt') as fhin, gz.open(_fasta_r,
                                                   'wb') as fhreformat:
        idmap_string = "{}\t{}\t{}\t{}\n".format("accession",
                                                 "accession.version", "taxid",
                                                 "gi")
        write_idmap(idmap_string, fhidmap)
        for i, record in enumerate(
                tqdm.tqdm(parse(fhin, "fasta"), unit=" records", ncols=100,
                          total=N), start=1):
            newid, taxid, format = parse_seqid(record)
            if taxid is not None:
                # Check that the taxid is valid
                # Skip the entry if it is not
                try:
                    taxid = int(taxid)
                except ValueError:
                    continue

            if len(newid) > maxidlen:
                newid = "id{j}".format(j=j)
                idmap[record.id] = newid
                j += 1
            # Skip UniRef entries with no taxid
            if format == "uniref":
                if not uniref:
                    uniref = True
                if not taxid:
                    continue
                # Write newid to taxid mapping
                idmap_string = "{}\t{}\t{}\t{}\n".format(newid, newid, taxid,
                                                         newid)
                write_idmap(idmap_string, fhidmap)
            fasta_string = ">{}\n{}\n".format(newid, str(
                record.seq))  # Write record with new id
            fhreformat.write(fasta_string.encode())
    sys.stderr.write("{}/{} records parsed\n".format(i, i))
    # Close idmap file handle if open
    if fhidmap:
        fhidmap.close()
    # If uniref database the mapfile has been created during reformatting
    if uniref:
        move(_mapfile, mapfile)
    # Move reformmated fasta
    move(_fasta_r, reformatted)
    # Check if there are remapped sequences
    if len(idmap.keys()) > 0:
        with gz.open(_idmapfile, 'wb') as fhout:
            for old, new in idmap.items():
                line = "{}\t{}\n".format(old, new)
                fhout.write(line.encode())
        move(_idmapfile, idmapfile)
        sys.stderr.write(
            "{}/{} sequence ids longer than {} characters written "
            "to {}\n".format(len(idmap), i, maxidlen, idmapfile))
        if not uniref:
            sys.stderr.write(
                ">>>Make sure to update your prot.accession2taxid.gz file "
                "with updated ids. See 'contigtax update'<<<\n")
    return 0
