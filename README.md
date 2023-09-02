
# bulkatac2fragments


This software takes paired-end ATAC-seq files and generes 10x-like fragment files.
Thus they can be analyzed using single-cell workflows (such as signac or archr)

## Requirements

The samtools and bedtools commands are used internally and must be on the path.

## Usage

The software processes one bam file at a time, which must be paired end. Set up a loop or SLURM script
to process all samples. They will be deduplicated as part of the process. Currently the number of read
copies is not stored in the output file. Just specifying the input.bam is usually good enough. The
samples can be renamed later anyway. The BAM file must be position sorted.

java -jar bulkatac2fragments.jar input.bam [output.tsv [barcode_to_write]]

Once all files processed, they can be merged:

bedtools merge -i *.atac_fragments.tsv > atac_fragments.tsv

Note that bgzip must be used for compression:

bgzip -@ 8 -i atac_fragments.tsv

These files are then best indexed right away, unless your downstream framework takes care of it for you:

tabix -p vcf atac_fragments.tsv.gz



## Downstream analysis

It is recommended that the data is aligned with a CellRanger-ATAC reference genome, as all tutorials for
e.g. Signac or Archr assume this. Thus it would minimize the surprises. CellRanger internally uses STAR
as the aligner, and thus it would make sense to stick to it overall.
