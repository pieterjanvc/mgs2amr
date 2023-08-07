# MGS2AMR revisions

The files in this folder are **not** part of the MGS2AMR core pipeline, 
but are part of the accompanying manuscript's supplemental materials.

The pipeline described below was suggest by reviewers as an alternative
approach to the described MGS2AMR pipeline, and was created to compare 
performance, both in accuracy and speed.

1) Metagenome assembly is performed by using metaSPAdes
2) Any regions containing ARG are detected using DIAMOND
3) Custom R scripts analyse the DIAMOND and metaSPAdes output to decide on 
specific ARG presence and select relevant scaffolds with regions up to 5000 bp 
around any ARG (similar to step 2 in the MGS2AMR pipeline)
4) BLASTn is used to align the scaffolds to the nucleotide database to discover 
relevant bacteria (similar to step 3 in the MGS2AMR pipeline
with identical settings)
5) Custom R scripts analyse the output and annotate ARG with specific 
bacteria (similar to step4 in the MGS2AMR)

The following dependencies need to be installed and accessible through `$PATH`
(Tool versions later than the described ones can be used, but the pipeline was
created and tested with the provided version numbers)

- spades/3.15.5
- python3
- sratoolkit/3.0.3
- R/4.2.3
- bbmap/38.96.0
- sqlite3/3.33.0
- pigz/2.6.0
- blast/2.13.0
- diamond/2.1.8
- [SEQ2MGS](https://github.com/pieterjanvc/seq2mgs)

The pipeline can be run using the provided `revPipeline.sh` script using the
following parameters

- -i isolate SRA ID to use for mix-in by SEQ2MGS
- -b background metagenome SRA ID to use by SEQ2MGS
- -r relative abundance of isolate in background (between 0 - 1)

Alternatively, you can continue from a previous pipelineId using
- -p continue from existing pipelineId

*NOTE: No support will be provided in case of issues with this pipeline, 
as it is not part of the MGS2AMR release*
