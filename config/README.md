# Configuration
In order to use this workflow, you need to customise the given example `config/config.yaml` file for your purposes. Please see the example here:
```yaml
samples_directory: resources/samples # This could be customized to any directory in which the sequencing samples are stored.
experiments:
  reference:
    genome_sequence: resources/reference_genome_sequence.fasta     # The fasta and gff files need to be provided (TriTrypDB sequence files can be automatically downloaded).
    genome_annotations: resources/reference_genome_annotations.gff

  cellline1:
    parent: reference
    samples:
      - SRR0000000      # samples can alternatively be SRA accession numbers, which will be automatically downloaded.
      - cellline1_run2

  cellline2:
    parent: reference
    samples:
      - cellline2_run1  # This is assumed to refer to the paired end fastq files <samples_directory>/cellline2/cellline2_run1_1.fq.gz and <samples_directory>/cellline2/cellline2_run1_2.fq.gz
      - cellline2_run2

  cellline1_clone1:
    parent: cellline1
    samples:
      - cellline1_clone1_run1
      - cellline1_clone1_run2

  cellline1_clone2:
    parent: cellline1
    samples:
      - cellline1_clone2_run1
      - cellline1_clone2_run2

  cellline2_clone1:
    parent: cellline2
    samples:
      - cellline2_clone1_run1
      - cellline2_clone1_run2

  cellline2_clone2:
    parent: cellline2
    samples:
      - cellline2_clone2_run1
      - cellline2_clone2_run2
```

## Details
* `samples_directory`: This points to the location of the sequencing samples fastq files. It could be anywhere, but is usually in `resources/samples`.
* `experiments`: This is a list of cell lines describing the relationship tree we want to analyse. Every entry underneath is the name of a cell line.
* `genome_sequence`: For reference genomes, this is the location of the `.fasta` file holding the genome sequence. It's not needed for any other entry in the `experiments` list, only for reference genomes. In special cases, this can be downloaded automatically (see below).
* `genome_annotations`: Similar to `genome_sequence`, this is the location of the `.gff` file holding the genome annotation features. 
* `parent`: This gives the name of the cell line entry of the parent of the current cell line. If it is not given, the current entry is assumed to be a reference genome.
* `samples`: Holds a list of sequencing samples associated with this cell line entry. Any entry in this list can be in the form of `SRR123456789` when it is assumed to be an SRA accession number and will be automatically downloaded. Otherwise, if given as `cellline1_clone1_run1`, it is assumed to refer to the paired-end fastq files in the location `<samples_directory>/cellline1_clone1_run1_1.fq.gz` and `<samples_directory>/cellline1_clone1_run1_2.fq.gz`.

## Automatic downloading from TriTrypDB
If the `genome_sequence` entry is of the form `TriTrypDB-<version>_<organism>_Genome.fasta`, it will be automatically downloaded from TriTrypDB.