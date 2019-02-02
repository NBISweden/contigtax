# TANGO

**T**axonomic cl**A**ssificatio**N** of meta**G**enomic c**O**ntigs

## Usage

1. Download fasta file
```
tango download uniref100
```

2. Download NCBI taxdump
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
tango search assembly.fa uniref100/diamond.dmnd
```