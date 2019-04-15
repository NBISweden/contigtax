# TANGO
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/tango/README.html) ![dependencies](https://img.shields.io/conda/pn/bioconda/tango.svg) ![GNU General Public License](https://img.shields.io/aur/license/tango.svg) ![update history](https://img.shields.io/github/last-commit/johnne/tango/master.svg)

TANGO (**T**axonomic cl**A**ssificatio**N** of meta**G**enomic c**O**ntigs)
is a tool that assigns taxonomy to metagenomic contigs by querying contig
nucleotide sequences against a protein database using `diamond blastx`
and parses hits using rank-specific thresholds. The use of rank-specific
 thresholds was first introduced by [Luo et al 2014](https://academic.oup.com/nar/article/42/8/e73/1076763)
 and used with some modification as explained in [Alneberg et al 2018](https://www.nature.com/articles/sdata2018146).

## Usage

1. Download fasta file
```
tango download uniref100
```

2. Download NCBI taxonomy
```
tango download taxonomy
```

3. Reformat fasta file and create taxonmap
```
tango format uniref100/uniref100.fasta.gz uniref100/uniref100.reformat.fasta.gz
```

4. Build diamond database
```
tango build uniref100/uniref100.reformat.fasta.gz uniref100/prot.accession2taxid.gz taxonomy/nodes.dmp
```

5. Search
```
tango search -p 4 assembly.fa uniref100/diamond.dmnd assembly.tsv.gz
```

6. Assign
```
tango assign -p 4 assembly.tsv.gz assembly.taxonomy.tsv
```
