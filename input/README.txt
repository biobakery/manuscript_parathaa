The following inputs are used at some point in the R scripts that make the figures. 
These notes describe where to get or create each input file.

silva.seed_v138_1.tax : see https://mothur.org/blog/2021/SILVA-v138_1-reference-files/

taxmap_slv_ssu_ref_138.1.txt : https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz

20231215.silva.seed_v138_1.ng.dada.fasta,
20231215_silva.seed_v138_1.ng.dada.sp.fasta : run src/create.seedDB.dada2.R <- Note that on 12/15/23 Meg noticed that dada2 dbs were only Bacteria,
and changed the script to include Archaea as well. Performance of dada2 changes slightly, so figs/tables should be updated.

SILVAsubsample_SeedGenera_V1V2.pcr.fasta,
SILVAsubsample_SeedGenera_V4V5.pcr.fasta: run Synthetic.Reads.R, and run the additional commands that are commented out from Synthetic.Reads.R (the latter should be put in a script)

silva.seed_v138_1.ng.fasta : run degap.seqs() on the file from: https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz (needs a script)

silva_species_assignment_v138.1.fa : https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1

Mock (Figure 2) Files: 
SRR3225703.fasta : accession numbers in manuscript, from SRA
SRR3225701.unique.fasta : Use accession numbers in manuscript to get SRR3225701.fasta from SRA; run "unique.seqs" in mothur (needs script; also outputs SRR3225701.count_table )
SRR3225701.count_table : (See SRR3225701.unique.fasta)

ASD (Figure 3) Files:
all_samples_taxonomy_closed_reference.tsv : From Kelsey
all_samples_taxonomy_closed_reference_withseqs.tsv : From Kelsey
amplicons_fromKelsey.fasta : From Kelsey