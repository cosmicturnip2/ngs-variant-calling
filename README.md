# NGS Variant Calling

A reproducible, dockerised pipeline for end-to-end NGS variant calling and
annotation: from raw FASTQ reads through quality control, trimming,
alignment, duplicate marking, base quality score recalibration, variant
calling, and functional annotation. Built with conda/bioconda inside a
Miniforge-based image, set up to run on a MacBook without a native Linux
machine.

## Tools included

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — read quality control
- [Trimmomatic](https://github.com/usadellab/Trimmomatic) — adapter/quality trimming
- [BWA-MEM](https://github.com/lh3/bwa) — short-read alignment
- [HISAT2](http://daehwankimlab.github.io/hisat2/) — spliced aligner (RNA-seq)
- [GATK4](https://gatk.broadinstitute.org/) — variant calling, BQSR, and annotation (Funcotator)
- [samtools](https://www.htslib.org/) — SAM/BAM file manipulation

## Prerequisites

- [Docker Desktop](https://docs.docker.com/desktop/) installed and running

## Build the image

```bash
git clone https://github.com/<your-username>/ngs-variant-calling.git
cd ngs-variant-calling
docker build -t ngs-tools .
```

## Run it

```bash
mkdir -p ~/ngs-data
docker run -it --rm -v ~/ngs-data:/data ngs-tools bash
```

Any files placed in `~/ngs-data` on your host machine are available inside
the container at `/data`.

## Running the pipeline

See `pipeline/run.sh` for the full pipeline, covering:

1. Quality control (FastQC)
2. Paired-end adapter/quality trimming (Trimmomatic)
3. Reference indexing (BWA, samtools, GATK)
4. Alignment with read-group tagging (BWA-MEM)
5. SAM → BAM conversion (samtools)
6. Mate-information fixing and duplicate marking (samtools)
7. Base quality score recalibration (GATK4)
8. Variant calling (GATK4 HaplotypeCaller)
9. dbSNP annotation (GATK4 VariantAnnotator)
10. SNP/INDEL splitting for review (GATK4 SelectVariants)
11. Functional annotation (GATK4 Funcotator)

The result is a functionally annotated VCF file listing candidate genetic
variants for the sample, along with their predicted gene/protein effects.

### Reference files needed

- A reference genome FASTA (e.g. `hg38.fa`)
- A known-sites VCF (e.g. dbSNP) for base quality score recalibration and
  variant annotation
- A Funcotator data sources folder for functional annotation, downloaded
  separately with `gatk FuncotatorDataSourceDownloader`

## Project structure

```
ngs-variant-calling/
├── Dockerfile
├── README.md
├── LICENSE
├── .dockerignore
└── pipeline/
    └── run.sh
```

## License

MIT — see [LICENSE](LICENSE).
