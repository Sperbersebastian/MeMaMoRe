Inputs
- SRA/assemblies/spades/<SAMPLE>/contigs.fasta
- SRA/assemblies/spades/<SAMPLE>/assembly_graph_with_scaffolds.gfa
- SRA/assemblies/spades/<SAMPLE>/map/<SAMPLE>.sorted.bam

Run
- Called by bin/modules/binning.sh → comebin/exec.sh
- Uses envs: comebin_env (COMEBin), binny_env (jgi_summarize_bam_contig_depths)

Outputs
- SRA/binning/comebin/<SAMPLE>/bins/*.fa
- SRA/binning/comebin/<SAMPLE>/contigs_to_bin.with_set.tsv
- SRA/binning/comebin/<SAMPLE>/run.log
