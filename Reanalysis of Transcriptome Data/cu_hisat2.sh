
hisat2-build GCA_000698785.1_ASM69878v1_genomic.fna ./hisat2/Genome_Index/nvie


hisat2  -q -p 16 --rna-strandness R --summary-file ./hisat2/Sample_hisat2_summary -x ./hisat2/Genome_Index/nvie -U Unaligned_Transcript_File -S ./hisat2/Aligned_File.sam
