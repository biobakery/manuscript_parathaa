import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Analysis Template"     #Update the description as needed
    ) 



# Add workflow arguments
workflow.add_argument(
    name="paraDir",
    desc="PATH to parathaa github repo"
)

workflow.add_argument(
    name="threads",
    desc="Number of threads",
    default=1
)


workflow.add_argument(
    name="benchonly",
    desc="Only run benchmarks",
    action="store_true"
)


workflow.add_argument(
    name="sensitive",
    desc="run parathaa in sensitive mode",
    action="store_true"
)

workflow.add_argument(
    name="bench_FL",
    desc="Set this to run full length benchmarks on all datasets along with V1V2 and V4V5 subregions",
    action="store_true"

)



# Parsing the workflow arguments
args = workflow.parse_args()

#add target files
silva_seed_db="input/silva.seed_v138_1.align"
silva_seed_tax="input/silva.seed_v138_1.tax"
silva_seed_ng="input/silva.seed_v138_1.ng.fasta"
dada2_spec_assign="input/silva_species_assignment_v138.1.fa" #what is this used for????
silva_taxonomy_file="input/taxmap_slv_ssu_ref_138.1.txt"
dada2_seed_db="input/20231215.silva.seed_v138_1.ng.dada.fasta"
dada2_seed_db_sp="input/20231215_silva.seed_v138_1.ng.dada.sp.fasta"
dada2_seed_db_FL="input/20231215.silva.seed_v138_1.ng.dada_FL.fasta"

syn_IDs="input/subsampleIDs_SeedGenera.txt"
silva_full_db="input/SILVA_138.1_SSURef_tax_silva.fasta"
silva_full_db_DNA="input/SILVA_138.1_SSURef_tax_silva.DNA.fasta"
FL_syn_reads="input/SILVAsubsample_SeedGenera.fasta"
V4V5_syn_reads="input/SILVAsubsample_SeedGenera_V4V5.pcr.fasta"
V1V2_syn_reads="input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta"
V1V2_db="input/SILVA_V1V2"
V4V5_db="input/SILVA_V4V5"
FL_db="input/SILVA_FL"
V1V2_assignments="output/V1V2_syn/"
V4V5_assignments="output/V4V5_syn/"
V1V2_para_taxa="output/V1V2_syn/taxonomic_assignments.tsv"
V4V5_para_taxa="output/V4V5_syn/taxonomic_assignments.tsv"

V4V5_mock_fasta="input/Mock_data/SRR3225703.fasta"
V4V5_mock_out="output/V4V5_Mock"
V4V5_mock_assignment="output/V4V5_Mock/taxonomic_assignments.tsv"
V1V2_mock_fasta_uni="input/Mock_data/SRR3225701.unique.fasta"
V1V2_mock_out="output/V1V2_Mock"
V1V2_mock_assignment="output/V1V2_Mock/taxonomic_assignments.tsv"
V1V2_mock_fasta="input/Mock_data/SRR3225701.fasta"
V1V2_mock_counts="input/Mock_data/SRR3225701.count_table"
Fig2="output/Figures/Fig_2A_MockHeatmap.png"


ASD_reads="input/ASD_data/amplicons_fromKelsey.fasta"
ASD_out="output/V4V5_ASD"
ASD_assignments="output/V4V5_ASD/taxonomic_assignments.tsv"
ASD_counts="input/ASD_data/all_samples_taxonomy_closed_reference.tsv"
ASD_countSeq="input/ASD_data/all_samples_taxonomy_closed_reference_withseqs.tsv"
ASD_figure="output/BrayCurtisPCoA.png"

even_genus_IDs="input/even_queryIDs.txt"
novel_genus_IDs="input/novel_query_ids.txt"
hold1_genus_IDs="input/holdout_query1.txt"
hold2_genus_IDs="input/holdout_query2.txt"

FL_even_genus_reads="input/FL_even_genus.fasta"
FL_novel_genus_reads="input/FL_novel_genus.fasta"
FL_holdout1_reads="input/FL_holdout1.fasta"
FL_holdout2_reads="input/FL_holdout2.fasta"


V4V5_novel_reads="input/V4V5_novel_genus.pcr.fasta"
V1V2_novel_reads="input/V1V2_novel_genus.pcr.fasta"

V4V5_even_reads="input/V4V5_even_genus.pcr.fasta"
V1V2_even_reads="input/V1V2_even_genus.pcr.fasta"

V4V5_holdout1_reads="input/V4V5_holdout1.pcr.fasta"
V1V2_holdout1_reads="input/V1V2_holdout1.pcr.fasta"

V4V5_holdout2_reads="input/V4V5_holdout2.pcr.fasta"
V1V2_holdout2_reads="input/V1V2_holdout2.pcr.fasta"

V1V2_even_assignments="output/V1V2_even/"
V1V2_even_tax="output/V1V2_even/taxonomic_assignments.tsv"
V4V5_even_assignments="output/V4V5_even/"
V4V5_even_tax="output/V4V5_even/taxonomic_assignments.tsv"


V1V2_novel_assignments="output/V1V2_novel/"
V1V2_novel_tax="output/V1V2_novel/taxonomic_assignments.tsv"
V4V5_novel_assignments="output/V4V5_novel/"
V4V5_novel_tax="output/V4V5_novel/taxonomic_assignments.tsv"

V1V2_holdout1_assignments="output/V1V2_holdout1/"
V1V2_holdout1_tax="output/V1V2_holdout1/taxonomic_assignments.tsv"
V4V5_holdout1_assignments="output/V4V5_holdout1/"
V4V5_holdout1_tax="output/V4V5_holdout1/taxonomic_assignments.tsv"

V1V2_holdout2_assignments="output/V1V2_holdout2/"
V1V2_holdout2_tax="output/V1V2_holdout2/taxonomic_assignments.tsv"
V4V5_holdout2_assignments="output/V4V5_holdout2/"
V4V5_holdout2_tax="output/V4V5_holdout2/taxonomic_assignments.tsv"

original_bench_out="output/benchmarks/original/"
original_V1V2_bench="output/benchmarks/original/Figures/synth_mult_arc/V1V2_full_comparisons.RData"
original_V4V5_bench="output/benchmarks/original/Figures/synth_mult_arc/V4V5_full_comparisons.RData"


even_bench_out="output/benchmarks/even/"
even_V1V2_bench="output/benchmarks/even/Figures/synth_mult_arc/V1V2_full_comparisons.RData"
even_V4V5_bench="output/benchmarks/even/Figures/synth_mult_arc/V4V5_full_comparisons.RData"

novel_bench_out="output/benchmarks/novel/"
novel_V1V2_bench="output/benchmarks/novel/Figures/synth_mult_arc/V1V2_full_comparisons.RData"
novel_V4V5_bench="output/benchmarks/novel/Figures/synth_mult_arc/V4V5_full_comparisons.RData"

holdout1_bench_out="output/benchmarks/holdout1/"
holdout1_V1V2_bench="output/benchmarks/holdout1/Figures/synth_mult_arc/V1V2_full_comparisons.RData"
holdout1_V4V5_bench="output/benchmarks/holdout1/Figures/synth_mult_arc/V4V5_full_comparisons.RData"

holdout2_bench_out="output/benchmarks/holdout2/"
holdout2_V1V2_bench="output/benchmarks/holdout2/Figures/synth_mult_arc/V1V2_full_comparisons.RData"
holdout2_V4V5_bench="output/benchmarks/holdout2/Figures/synth_mult_arc/V4V5_full_comparisons.RData"

#Full length target files
FL_original_tax="output/FL_original/taxonomic_assignments.tsv"
FL_original_assignments="output/FL_original/"

FL_even_tax="output/FL_even/taxonomic_assignments.tsv"
FL_even_assignments="output/FL_even/"

FL_novel_tax="output/FL_novel/taxonomic_assignments.tsv"
FL_novel_assignments="output/FL_novel/"

FL_holdout1_tax="output/FL_holdout1/taxonomic_assignments.tsv"
FL_holdout1_assignment="output/FL_holdout1/"

FL_holdout2_tax="output/FL_holdout2/taxonomic_assignments.tsv"
FL_holdout2_assignment="output/FL_holdout2/"

FL_original_bench_out="output/benchmarks/original/FL/"
FL_original_bench="output/benchmarks/original/FL/FL_full_comparisons.RData"

FL_even_bench_out="output/benchmarks/even/FL/"
FL_even_bench="output/benchmarks/even/FL/FL_full_comparisons.RData"

FL_novel_bench_out="output/benchmarks/novel/FL/"
FL_novel_bench="output/benchmarks/novel/FL/FL_full_comparisons.RData"

FL_novel_holdout1_out="output/benchmarks/holdout1/FL/"
FL_novel_holdout1="output/benchmarks/holdout1/FL/FL_full_comparisons.RData"

FL_novel_holdout2_out="output/benchmarks/holdout2/FL/"
FL_novel_holdout2="output/benchmarks/holdout2/FL/FL_full_comparisons.RData"


oligos_v4v5=os.path.join(args.paraDir, "input/primers/V4V5.oligos")
oligos_v1v2=os.path.join(args.paraDir, "input/primers/V1V2.oligos")


#### Prep required files ####
workflow.add_task(
    "gunzip input/Mock_data/SRR3225701.fasta.gz",
    targets=V4V5_mock_fasta,
    names="unzipping mock V4V5 fasta"
)

## downlaod silva files
workflow.add_task(
    "wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_1.tgz -P input/; tar -xf input/silva.seed_v138_1.tgz -C input/",
    targets=[silva_seed_db, silva_seed_tax],
    name="download silva data"
)

workflow.add_task(
    "wget https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 -P input/; mv input/silva_species_assignment_v138.1.fa.gz?download=1 input/silva_species_assignment_v138.1.fa.gz; gunzip input/silva_species_assignment_v138.1.fa.gz",
    targets=dada2_spec_assign,
    name="Download dada2 info"
)

workflow.add_task(
    "wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz -P input/; gunzip input/taxmap_slv_ssu_ref_138.1.txt.gz",
    targets=silva_taxonomy_file,
    name="download full silva taxonomy file"
)

workflow.add_task(
   "wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_tax_silva.fasta.gz -P input/; gunzip input/SILVA_138.1_SSURef_tax_silva.fasta.gz",
   targets=silva_full_db,
   name="download full silva database"
)

#download pre-computed Databases

workflow.add_task(
    "wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V1V2.tar.gz -P input/; tar -xf input/SILVA_V1V2.tar.gz -C input/; rm input/SILVA_V1V2.tar.gz",
    targets=V1V2_db,
    name="download pre-computed V1V2 database"
)

workflow.add_task(
    "wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V4V5.tar.gz -P input/; tar -xf input/SILVA_V4V5.tar.gz -C input/; rm input/SILVA_V4V5.tar.gz",
    targets=V4V5_db,
    name="download pre-computed V4V5 database"
)

workflow.add_task(
    "wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_FL.tar.gz -P input/; tar -xf input/SILVA_FL.tar.gz -C input/; rm input/SILVA_FL.tar.gz",
    targets=V4V5_db,
    name="download pre-computed V4V5 database"
)

#Generate DADA2 files

workflow.add_task(
    "mothur '#degap.seqs(fasta=[args[0]])'",
    args=silva_seed_db,
    depends=silva_seed_db,
    name="degapping input sequences",
    targets=silva_seed_ng
)

            
workflow.add_task(
    "Rscript src/create.seedDB.dada2.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
    args=args.paraDir,
    depends=[silva_seed_ng, silva_taxonomy_file],
    targets=[dada2_seed_db, dada2_seed_db_sp, dada2_seed_db_FL],
    name="generating DADA2 silva.seed DB"
)

# Generate Synthetic reads for original bench

workflow.add_task(
   "Rscript src/Synthetic.Reads.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
   args=args.paraDir,
   depends=[silva_seed_tax, silva_taxonomy_file],
   targets=[syn_IDs],
   name="generating synthetic read IDs"
)

# Generate Synthetic reads for even genus bench
workflow.add_task(
    "Rscript src/Generate_even_genera_query.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
    args=args.paraDir,
    depends=[silva_seed_tax, silva_taxonomy_file],
    targets=[even_genus_IDs],
    name="generate even genus read IDs"
)

# Generate Synthetic reads for novel genus bench
workflow.add_task(
   "Rscript src/Generate_novel_genera_query_dev.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
   args=args.paraDir,
   depends=[silva_seed_tax, silva_taxonomy_file],
   targets=[novel_genus_IDs],
   name="generate novel genus read IDs"
)

# Generate Synthetic reads for hold out 1 and 2
workflow.add_task(
    "Rscript src/Generate_more_holdout_data.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
    args=args.paraDir,
    depends=[silva_seed_tax, silva_taxonomy_file],
    targets=[hold1_genus_IDs, hold2_genus_IDs],
    name="generate holdout 1 and 2 read IDs"
)



#convert silva db from RNA to DNA
workflow.add_task(
   "sed '/^[^>]/s/U/T/g' [depends[0]] > [targets[0]]",
   depends=[silva_full_db],
   targets=silva_full_db_DNA,
   name="Converting silva RNA DB to DNA" 
)

#generate the full length benchmark read files

#original
workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, syn_IDs],
   targets=FL_syn_reads,
   name="generating full length synthetic read file"
)

#even genus
workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, even_genus_IDs],
   targets=FL_even_genus_reads,
   name="generating full length even genus reads"
)

#novel genus
workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, novel_genus_IDs],
   targets=FL_novel_genus_reads,
   name="generating full length novel genus reads"
)

workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, hold1_genus_IDs],
   targets=FL_holdout1_reads,
   name="generating full length holdout1 reads"
)

workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, hold2_genus_IDs],
   targets=FL_holdout2_reads,
   name="generating full length holdout2 reads"
)


#PCR syn reads
workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/SILVAsubsample_SeedGenera.pcr.fasta [targets[0]]",
   depends=FL_syn_reads,
   args=oligos_v4v5,
   targets=V4V5_syn_reads,
   name="generating v4v5 synthetic reads"
)

workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/SILVAsubsample_SeedGenera.pcr.fasta [targets[0]]",
   depends=FL_syn_reads,
   args=oligos_v1v2,
   targets=V1V2_syn_reads,
   name="generating v1v2 synthetic reads"
)

## Filter the reads to minimum length of 200
workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V1V2_syn_reads
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V4V5_syn_reads
)

#PCR novel genus reads
workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_novel_genus.pcr.fasta [targets[0]]",
   depends=FL_novel_genus_reads,
   args=oligos_v4v5,
   targets=V4V5_novel_reads,
   name="generating v4v5 novel reads"
)

workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_novel_genus.pcr.fasta [targets[0]]",
   depends=FL_novel_genus_reads,
   args=oligos_v1v2,
   targets=V1V2_novel_reads,
   name="generating v1v2 novel reads"
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V1V2_novel_reads
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V4V5_novel_reads
)


#PCR even genus reads
workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_even_genus.pcr.fasta [targets[0]]",
   depends=FL_even_genus,
   args=oligos_v4v5,
   targets=V4V5_even_reads,
   name="generating v4v5 even reads"
)

workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_even_genus.pcr.fasta [targets[0]]",
   depends=FL_even_genus,
   args=oligos_v1v2,
   targets=V1V2_even_reads,
   name="generating v1v2 even reads"
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V1V2_even_reads
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V4V5_even_reads
)


#PCR Holdout1
workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_holdout1.pcr.fasta [targets[0]]",
   depends=FL_holdout1_reads,
   args=oligos_v4v5,
   targets=V4V5_holdout1_reads,
   name="generating v4v5 holdout1 reads"
)

workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_holdout1.pcr.fasta [targets[0]]",
   depends=FL_holdout1_reads,
   args=oligos_v1v2,
   targets=V1V2_holdout1_reads,
   name="generating v1v2 holdout1 reads"
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V1V2_holdout1_reads
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V4V5_holdout1_reads
)


#PCR Holdout2
workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_holdout2.pcr.fasta [targets[0]]",
   depends=FL_holdout2_reads,
   args=oligos_v4v5,
   targets=V4V5_holdout2_reads,
   name="generating v4v5 holdout2 reads"
)

workflow.add_task(
   "mothur '#pcr.seqs(fasta=[depends[0]], oligos=[args[0]], pdiffs=0, rdiffs=0)'; mv input/FL_holdout2.pcr.fasta [targets[0]]",
   depends=FL_holdout2_reads,
   args=oligos_v1v2,
   targets=V1V2_holdout2_reads,
   name="generating v1v2 holdout2 reads"
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V1V2_holdout2_reads
)

workflow.add_task(
    "seqkit seq --min-len 200 [depends[0]] > [depends[0]]",
    depends=V4V5_holdout2_reads
)



#run parathaa on V1V2 synthetic using default specific mode
if(not args.sensitive):
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V1V2_db, V1V2_syn_reads],
        targets=V1V2_para_taxa,
        args=[args.threads,V1V2_assignments],
        name="Assigning taxonomy to V1V2 synthetic reads"
    )

    #run parathaa on V4V5 synthetic using defualt specific mode

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V4V5_db, V4V5_syn_reads],
        targets=V4V5_para_taxa,
        args=[args.threads, V4V5_assignments],
        name="Assigning taxonomy to V4V5 synthetic reads"
    )


    # Run parathaa on V1V2 and V4V5 even genus
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V1V2_db, V1V2_even_reads],
        targets=V1V2_even_tax,
        args=[args.threads,V1V2_even_assignments],
        name="Assigning taxonomy to V1V2 even genus reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V4V5_db, V4V5_even_reads],
        targets=V4V5_even_tax,
        args=[args.threads, V4V5_even_assignments],
        name="Assigning taxonomy to V4V5 even genus reads"
    )


    # Run parathaa on V1V2 and V4V5 novel genus
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V1V2_db, V1V2_novel_reads],
        targets=V1V2_novel_tax,
        args=[args.threads,V1V2_novel_assignments],
        name="Assigning taxonomy to V1V2 novel genus reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V4V5_db, V4V5_novel_reads],
        targets=V4V5_novel_tax,
        args=[args.threads, V4V5_novel_assignments],
        name="Assigning taxonomy to V4V5 novel genus reads"
    )

    # Run parathaa on V1V2 and V4V5 holdout1 
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V1V2_db, V1V2_holdout1_reads],
        targets=V1V2_holdout1_tax,
        args=[args.threads,V1V2_holdout1_assignments],
        name="Assigning taxonomy to V1V2 holdout 1 reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V4V5_db, V4V5_holdout1_reads],
        targets=V4V5_holdout1_tax,
        args=[args.threads, V4V5_holdout1_assignments],
        name="Assigning taxonomy to V4V5 holdout 1 reads"
    )

    # Run parathaa on V1V2 and V4V5 holdout2
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V1V2_db, V1V2_holdout2_reads],
        targets=V1V2_holdout2_tax,
        args=[args.threads,V1V2_holdout2_assignments],
        name="Assigning taxonomy to V1V2 holdout 2 reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
        depends=[V4V5_db, V4V5_holdout2_reads],
        targets=V4V5_holdout2_tax,
        args=[args.threads, V4V5_holdout2_assignments],
        name="Assigning taxonomy to V4V5 holdout 2 reads"
    )
    if(bench_FL):
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_syn_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_even_genus_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_novel_genus_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_holdout1_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_holdout2_reads
        )

        # original
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
            depends=[FL_db, FL_syn_reads],
            targets=FL_original_tax,
            args=[args.threads, FL_original_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Even
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
            depends=[FL_db, FL_even_genus_reads],
            targets=FL_even_tax,
            args=[args.threads, FL_even_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Novel
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
            depends=[FL_db, FL_novel_genus_reads],
            targets=FL_novel_tax,
            args=[args.threads, FL_novel_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Holdout1
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
            depends=[FL_db, FL_holdout1_reads],
            targets=FL_holdout1_tax,
            args=[args.threads, FL_holdout1_assignment],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Holdout2
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]]",
            depends=[FL_db, FL_holdout2_reads],
            targets=FL_holdout2_tax,
            args=[args.threads, FL_holdout2_assignment],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
elif(args.sensitive):
    
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V1V2_db, V1V2_syn_reads],
        targets=V1V2_para_taxa,
        args=[args.threads,V1V2_assignments],
        name="Assigning taxonomy to V1V2 synthetic reads"
    )

    #run parathaa on V4V5 synthetic using defualt specific mode

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V4V5_db, V4V5_syn_reads],
        targets=V4V5_para_taxa,
        args=[args.threads, V4V5_assignments],
        name="Assigning taxonomy to V4V5 synthetic reads"
    )


    # Run parathaa on V1V2 and V4V5 even genus
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V1V2_db, V1V2_even_reads],
        targets=V1V2_even_tax,
        args=[args.threads,V1V2_even_assignments],
        name="Assigning taxonomy to V1V2 even genus reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V4V5_db, V4V5_even_reads],
        targets=V4V5_even_tax,
        args=[args.threads, V4V5_even_assignments],
        name="Assigning taxonomy to V4V5 even genus reads"
    )


    # Run parathaa on V1V2 and V4V5 novel genus
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V1V2_db, V1V2_novel_reads],
        targets=V1V2_novel_tax,
        args=[args.threads,V1V2_novel_assignments],
        name="Assigning taxonomy to V1V2 novel genus reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V4V5_db, V4V5_novel_reads],
        targets=V4V5_novel_tax,
        args=[args.threads, V4V5_novel_assignments],
        name="Assigning taxonomy to V4V5 novel genus reads"
    )

    # Run parathaa on V1V2 and V4V5 holdout1 
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V1V2_db, V1V2_holdout1_reads],
        targets=V1V2_holdout1_tax,
        args=[args.threads,V1V2_holdout1_assignments],
        name="Assigning taxonomy to V1V2 holdout 1 reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V4V5_db, V4V5_holdout1_reads],
        targets=V4V5_holdout1_tax,
        args=[args.threads, V4V5_holdout1_assignments],
        name="Assigning taxonomy to V4V5 holdout 1 reads"
    )

    # Run parathaa on V1V2 and V4V5 holdout2
    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V1V2_db, V1V2_holdout2_reads],
        targets=V1V2_holdout2_tax,
        args=[args.threads,V1V2_holdout2_assignments],
        name="Assigning taxonomy to V1V2 holdout 2 reads"
    )

    workflow.add_task(
        "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
        depends=[V4V5_db, V4V5_holdout2_reads],
        targets=V4V5_holdout2_tax,
        args=[args.threads, V4V5_holdout2_assignments],
        name="Assigning taxonomy to V4V5 holdout 2 reads"
    )
    #run full length bench mark assignments if its set
    if(bench_FL):
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_syn_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_even_genus_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_novel_genus_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_holdout1_reads
        )
        workflow.add_task(
            "seqkit seq --min-len 1000 [depends[0]] > [depends[0]]",
            depends=FL_holdout2_reads
        )
        # original
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
            depends=[FL_db, FL_syn_reads],
            targets=FL_original_tax,
            args=[args.threads, FL_original_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Even
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
            depends=[FL_db, FL_even_genus_reads],
            targets=FL_even_tax,
            args=[args.threads, FL_even_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Novel
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
            depends=[FL_db, FL_novel_genus_reads],
            targets=FL_novel_tax,
            args=[args.threads, FL_novel_assignments],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Holdout1
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
            depends=[FL_db, FL_holdout1_reads],
            targets=FL_holdout1_tax,
            args=[args.threads, FL_holdout1_assignment],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )
        #Holdout2
        workflow.add_task(
            "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [args[1]] --threads [args[0]] --sensitive",
            depends=[FL_db, FL_holdout2_reads],
            targets=FL_holdout2_tax,
            args=[args.threads, FL_holdout2_assignment],
            name="Assigning taxonomy to V4V5 holdout 2 reads"
        )



# Run script to benchmark original V1V2 and V4V5 synthetic data

workflow.add_task(
    "Rscript src/Plots.Table.1_dev.R -p [args[0]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --queryV4V5 [depends[4]] --queryV1V2 [depends[5]] -t [depends[6]] -o [args[1]] -s [depends[7]]",
    depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_para_taxa, V1V2_para_taxa, V4V5_syn_reads, V1V2_syn_reads, silva_taxonomy_file, silva_seed_tax],
    args=[args.paraDir, original_bench_out],
    targets=[original_V1V2_bench, original_V4V5_bench],
    name="Benchmarking synthetic assignments"
)

#run script to benchmark even genus data
workflow.add_task(
    "Rscript src/Plots.Table.1_dev.R -p [args[0]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --queryV4V5 [depends[4]] --queryV1V2 [depends[5]] -t [depends[6]] -o [args[1]] -s [depends[7]]",
    depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_even_tax, V1V2_even_tax, V4V5_even_reads, V1V2_even_reads, silva_taxonomy_file, silva_seed_tax],
    args=[args.paraDir, even_bench_out],
    targets=[even_V1V2_bench, even_V4V5_bench],
    name="Benchmarking even genus assignments"

)

#run script to benchmark novel genus data
workflow.add_task(
    "Rscript src/Plots.Table.1_dev.R -p [args[0]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --queryV4V5 [depends[4]] --queryV1V2 [depends[5]] -t [depends[6]] -o [args[1]] -s [depends[7]]",
    depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_novel_tax, V1V2_novel_tax, V4V5_novel_reads, V1V2_novel_reads, silva_taxonomy_file, silva_seed_tax],
    args=[args.paraDir, novel_bench_out],
    targets=[novel_V1V2_bench, novel_V4V5_bench],
    name="Benchmarking novel genus assignments"
)

#run script to benchmark holdout 1
workflow.add_task(
    "Rscript src/Plots.Table.1_dev.R -p [args[0]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --queryV4V5 [depends[4]] --queryV1V2 [depends[5]] -t [depends[6]] -o [args[1]] -s [depends[7]]",
    depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_holdout1_tax, V1V2_holdout1_tax, V4V5_holdout1_reads, V1V2_holdout1_reads, silva_taxonomy_file, silva_seed_tax],
    args=[args.paraDir, holdout1_bench_out],
    targets=[holdout1_V1V2_bench, holdout1_V4V5_bench],
    name="Benchmarking holout 1 assignments"
)

workflow.add_task(
    "Rscript src/Plots.Table.1_dev.R -p [args[0]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --queryV4V5 [depends[4]] --queryV1V2 [depends[5]] -t [depends[6]] -o [args[1]] -s [depends[7]]",
    depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_holdout2_tax, V1V2_holdout2_tax, V4V5_holdout2_reads, V1V2_holdout2_reads, silva_taxonomy_file, silva_seed_tax],
    args=[args.paraDir, holdout2_bench_out],
    targets=[holdout2_V1V2_bench, holdout2_V4V5_bench],
    name="Benchmarking holout 1 assignments"
)

if(bench_FL):
    #Bench original FL
    workflow.add_task(
        "Rscript src/full_lenfth_bench.R -p [args[0]] --dada_db_FL [depends[0]] -t [depends[1]] -o [args[1]] --paraAssign [depends[2]] --query [depends[3]] -s [depends[4]]",
        depends=[dada2_seed_db_FL, silva_taxonomy_file, FL_original_tax, FL_syn_reads, silva_seed_tax],
        args=[args.paraDir, FL_original_bench_out],
        targets=[FL_original_bench],
        name="Benchmarking Full length original dataset"

    )
    #Bench even FL
    workflow.add_task(
        "Rscript src/full_lenfth_bench.R -p [args[0]] --dada_db_FL [depends[0]] -t [depends[1]] -o [args[1]] --paraAssign [depends[2]] --query [depends[3]] -s [depends[4]]",
        depends=[dada2_seed_db_FL, silva_taxonomy_file, FL_even_tax, FL_even_genus_reads, silva_seed_tax],
        args=[args.paraDir, FL_even_bench_out],
        targets=[FL_even_bench],
        name="Benchmarking Full length original dataset"
    )
    #Bench novel FL
    workflow.add_task(
        "Rscript src/full_lenfth_bench.R -p [args[0]] --dada_db_FL [depends[0]] -t [depends[1]] -o [args[1]] --paraAssign [depends[2]] --query [depends[3]] -s [depends[4]]",
        depends=[dada2_seed_db_FL, silva_taxonomy_file, FL_novel_tax, FL_novel_genus_reads, silva_seed_tax],
        args=[args.paraDir, FL_novel_bench_out],
        targets=[FL_novel_bench],
        name="Benchmarking Full length original dataset"
    )
    #Bench holdout1 FL
    workflow.add_task(
        "Rscript src/full_lenfth_bench.R -p [args[0]] --dada_db_FL [depends[0]] -t [depends[1]] -o [args[1]] --paraAssign [depends[2]] --query [depends[3]] -s [depends[4]]",
        depends=[dada2_seed_db_FL, silva_taxonomy_file, FL_holdout1_tax, FL_holdout1_reads, silva_seed_tax],
        args=[args.paraDir, FL_holdout1_bench_out],
        targets=[FL_holdout1_bench],
        name="Benchmarking Full length original dataset"
    )
    #Bench holdout2 FL
    workflow.add_task(
        "Rscript src/full_lenfth_bench.R -p [args[0]] --dada_db_FL [depends[0]] -t [depends[1]] -o [args[1]] --paraAssign [depends[2]] --query [depends[3]] -s [depends[4]]",
        depends=[dada2_seed_db_FL, silva_taxonomy_file, FL_holdout2_tax, FL_holdout2_reads, silva_seed_tax],
        args=[args.paraDir, FL_holdout2_bench_out],
        targets=[FL_holdout2_bench],
        name="Benchmarking Full length original dataset"
    )

########### End of benchmarking ###############


#only run if benchmark only flag is false
if(not args.benchonly):
# Run parathaa on mock V4V5
    workflow.add_task(
        "mkdir output/V4V5_Mock; parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [args[0]] --output [args[1]] --threads [args[2]]",
        depends=[V4V5_db],
        args=[V4V5_mock_fasta, V4V5_mock_out, args.threads],
        targets=V4V5_mock_assignment,
        name="Assigning parathaa taxonomy to V4V5 mock reads"
    )

    # Run parathaa on mock V1V2

    workflow.add_task(
        "mkdir output/V1V2_Mock; parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [args[0]] --output [args[1]] --threads [args[2]]",
        depends=[V1V2_db],
        args=[V1V2_mock_fasta_uni, V1V2_mock_out, args.threads],
        targets=V1V2_mock_assignment,
        name="Assigning parathaa taxonomy to V1V2 mock reads"
    )

    # Run script to generate mock figures

    workflow.add_task(
        "Rscript src/Plots.Figure.2.R --dada_db [depends[0]] --dada_db_sp [depends[1]] -o [args[0]] --paraAssignV4V5 [depends[2]] --paraAssignV1V2 [depends[3]] --fastaV4V5 [args[1]] --fastaV1V2 [args[2]] --fastaV1V2uni [args[3]] --V1V2Counts [args[4]]",
        depends=[dada2_seed_db, dada2_seed_db_sp, V4V5_mock_assignment, V1V2_mock_assignment],
        args=[args.output, V4V5_mock_fasta, V1V2_mock_fasta, V1V2_mock_fasta_uni, V1V2_mock_counts],
        targets=Fig2,
        name="Generate mock data figure (figure2)"
    )

    # Run parathaa on ASD data

    workflow.add_task(
        "mkdir output/V4V5_ASD; parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [args[0]] --output [args[1]] --threads [args[2]]",
        depends=[V4V5_db],
        args=[ASD_reads, ASD_out, args.threads],
        targets=ASD_assignments,
        name="Generating parathaa assignments for ASD dataset"
    )


    # Run script to generate ASD figures

    workflow.add_task(
        "Rscript src/Plots.Figure.3.R --counts [args[0]] --countSeq [args[1]] --fasta [args[2]] --dada_db [depends[0]] --dada_db_sp [depends[1]] --paraAssign [depends[2]] -o [args[3]]",
        depends=[dada2_seed_db, dada2_seed_db_sp, ASD_assignments],
        args=[ASD_counts, ASD_countSeq, ASD_reads, args.output],
        targets=ASD_figure,
        name="Generating ASD Figures"
    )



#done
workflow.go()
