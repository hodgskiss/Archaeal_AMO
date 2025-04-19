

mafft --reorder Selected_Sequences > Selected_Sequences_msa


hmmbuild --amino Subunit_hmm Selected_Sequences_msa
