# contigtax
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/contigtax/README.html) 
![dependencies](https://img.shields.io/conda/pn/bioconda/contigtax.svg) 
![Docker](https://img.shields.io/docker/pulls/nbisweden/contigtax)
![update history](https://img.shields.io/github/last-commit/NBISweden/contigtax/main.svg) 
![CI](https://github.com/NBISweden/contigtax/workflows/CI/badge.svg?branch=main)

contigtax is a tool that assigns taxonomy to metagenomic contigs by querying 
contig nucleotide sequences against a protein database using `diamond blastx` 
and parses hits using rank-specific thresholds. The use of rank-specific 
thresholds was first introduced by [Luo et al 2014](https://academic.oup.com/nar/article/42/8/e73/1076763)
and used with some modification as explained in [Alneberg et al 2018](https://www.nature.com/articles/sdata2018146).

## Install
Simplest way to install `contigtax` is via conda:
```bash
conda install -c bioconda contigtax
```

Alternatively, pull the docker image:
```bash
docker pull nbisweden/contigtax
```

## Usage

1. Download fasta file
```
contigtax download uniref100
```

2. Download NCBI taxonomy
```
contigtax download taxonomy
```

3. Reformat fasta file and create taxonmap
```
contigtax format uniref100/uniref100.fasta.gz uniref100/uniref100.reformat.fasta.gz
```

4. Build diamond database
```
contigtax build uniref100/uniref100.reformat.fasta.gz uniref100/prot.accession2taxid.gz taxonomy/nodes.dmp
```

5. Search (here assembled contigs are in file `assembly.fa`)
```
contigtax search -p 4 assembly.fa uniref100/diamond.dmnd assembly.tsv.gz
```

6. Assign (here output from the `contigtax search` step are in file `assembly.tsv.gz`)
```
contigtax assign -p 4 assembly.tsv.gz assembly.taxonomy.tsv
```

### Running contigtax with Docker

To run contigtax with docker simply substitute contigtax in the commands above with
`docker run --rm -v $(pwd):/work nbisweden/contigtax`, _e.g._:

1. Download fasta file 
```
docker run --rm -v $(pwd):/work nbisweden/contigtax download uniref100
```

2. Download NCBI taxonomy
```
docker run --rm -v $(pwd):/work nbisweden/contigtax download taxonomy
```

3. Reformat fasta file and create taxonmap
```
docker run --rm -v $(pwd):/work nbisweden/contigtax format uniref100/uniref100.fasta.gz uniref100/uniref100.reformat.fasta.gz
```

4. Build diamond database
```
docker run --rm -v $(pwd):/work nbisweden/contigtax build uniref100/uniref100.reformat.fasta.gz uniref100/prot.accession2taxid.gz taxonomy/nodes.dmp
```

5. Search (here assembled contigs are in file `assembly.fa`)
```
docker run --rm -v $(pwd):/work nbisweden/contigtax search -p 4 assembly.fa uniref100/diamond.dmnd assembly.tsv.gz
```

6. Assign (here output from the `contigtax search` step are in file `assembly.tsv.gz`)
```
docker run --rm -v $(pwd):/work nbisweden/contigtax assign -p 4 assembly.tsv.gz assembly.taxonomy.tsv
```
