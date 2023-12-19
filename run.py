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
syn_IDs="input/subsampleIDs_SeedGenera.txt"
silva_full_db="input/SILVA_138.1_SSURef_tax_silva.fasta"
silva_full_db_DNA="input/SILVA_138.1_SSURef_tax_silva.DNA.fasta"
FL_syn_reads="input/SILVAsubsample_SeedGenera.fasta"
V4V5_syn_reads="input/SILVAsubsample_SeedGenera_V4V5.pcr.fasta"
V1V2_syn_reads="input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta"
V1V2_db="input/SILVA_V1V2"
V4V5_db="input/SILVA_V4V5"
V1V2_assignments="output/V1V2_syn/"
V4V5_assignments="output/V4V5_syn/"



oligos_v4v5=os.path.join(args.paraDir, "input/primers/V4V5.oligos")
oligos_v1v2=os.path.join(args.paraDir, "input/primers/V1V2.oligos")
# Write workflow to download files

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
    targets=[dada2_seed_db, dada2_seed_db_sp],
    name="generating DADA2 silva.seed DB"
)

# Generate Synthetic reads

workflow.add_task(
   "Rscript src/Synthetic.Reads.R -p [args[0]] -s [depends[0]] -t [depends[1]]",
   args=args.paraDir,
   depends=[silva_seed_tax, silva_taxonomy_file],
   targets=[syn_IDs],
   name="generating synthetic read IDs"
)

#convert silva db from RNA to DNA
workflow.add_task(
   "sed '/^[^>]/s/U/T/g' [depends[0]] > [targets[0]]",
   depends=[silva_full_db],
   targets=silva_full_db_DNA,
   name="Converting silva RNA DB to DNA" 
)

workflow.add_task(
   "faSomeRecords [depends[0]] [depends[1]] [targets[0]]",
   depends=[silva_full_db_DNA, syn_IDs],
   targets=FL_syn_reads,
   name="generating full length synthetic read file"
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

#run parathaa on V1V2

workflow.add_task(
    "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [targets[0]] --threads [args[0]]",
    depends=[V1V2_db, V1V2_syn_reads],
    targets=V1V2_assignments,
    args=args.threads,
    name="Assigning taxonomy to V1V2 synthetic reads"
)

#run parathaa on V4V5

workflow.add_task(
    "parathaa_run_taxa_assignment --treeFiles [depends[0]] --query [depends[1]] --output [targets[0]] --threads [args[0]]",
    depends=[V4V5_db, V4V5_syn_reads],
    targets=V4V5_assignments,
    args=args.threads,
    name="Assigning taxonomy to V4V5 synthetic reads"
)


# Run the workflow
workflow.go()
