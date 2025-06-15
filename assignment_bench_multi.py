#anadama workflow that will run parathaa at multiple different speces multi thresholds
#it then calls a script that runs dada2 naive bayes with 0 minBoot for species and creates an ROC
#note for exact species matching this is not possible

import os
import numpy as np
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

import subprocess

workflow = Workflow(
    version="0.1.0",                    #Update the version as needed
    description="Benchmarking parathaa taxonomic assignment"     #Update the description as needed
) 


workflow.add_argument(
    name="maxMulti",
    desc="Maximum value for Multi",
    default="3"
)

workflow.add_argument(
   name="minMulti",
   desc="Minimum value for Multi",
   default="0"

)


workflow.add_argument(
    name="step",
    desc="incremental value ot increase multipler for each run",
    default="0.1"
)

workflow.add_argument(
   name="treeFilesA",
   desc="Directory with V4V5 parathaa DB"
)

workflow.add_argument(
   name="queryA",
   desc="query V4V5 sequences to classify"
)

workflow.add_argument(
   name="treeFilesB",
   desc="Directory with parathaa V1V2 DB"
)

workflow.add_argument(
   name="queryB",
   desc="query sequences to classify"
)

workflow.add_argument(
    name="FLtreeFile",
    desc="Directory with FL parathaa DB"
)

workflow.add_argument(
        name="FLquery",
        desc="query FL sequences to classify"
)

workflow.add_argument(
   name="threads",
   desc="number of threads"
)

workflow.add_argument(
   name="paraDir",
   desc="Parathaa Dir"
)

#workflow.add_argument(
#   name="dadaDb",
#   desc="dada database"
#)

#workflow.add_argument(
#   name="dadaDbSp",
#   desc="dada species database"
#)

#workflow.add_argument(
#   name="taxonomy",
#   desc="taxonomy file"
#)

args= workflow.parse_args()


## first we run an inital run that makes the trees etc. once that is done we can just call the assignment directly with the output from this first run.

def run_initial(task):
    init_dir=os.path.join(task.args[0], "init")
    os.mkdir(init_dir)
    V1V2_init_dir=os.path.join(task.args[0], "init", "V1V2")
    V4V5_init_dir=os.path.join(task.args[0], "init", "V4V5")
    FL_init_dir=os.path.join(task.args[0], "init", "FL")
    os.mkdir(V1V2_init_dir)
    os.mkdir(V4V5_init_dir)
    os.mkdir(FL_init_dir)

    callV1V2 = "parathaa_run_taxa_assignment --treeFiles " + task.args[1] + " --query " + task.args[2] + " --output " + V1V2_init_dir + " --threads " + task.args[3]
    callV4V5 = "parathaa_run_taxa_assignment --treeFiles " + task.args[4] + " --query " + task.args[5] + " --output " + V4V5_init_dir + " --threads " + task.args[3]
    callFL = "parathaa_run_taxa_assignment --treeFiles " + task.args[6] + " --query " + task.args[7] + " --output " + FL_init_dir + " --threads " + task.args[3]
    print('Running initial runs to create phylo trees')
    print(callV1V2)
    print(callV4V5)
    print(callFL)
    subprocess.call(callV1V2, shell=True)
    subprocess.call(callV4V5, shell=True)
    subprocess.call(callFL, shell=True)
    


### now we have the initial we can call the assignment script directly with the parameters that we want.
def run_bench(task):
    #generate array of numbers to test
    
    test_multi=np.arange(float(task.args[0]), float(task.args[1]), float(task.args[2]))
    for i in test_multi:
        par_dir=os.path.join(task.args[3], "multi_"+str(round(i, 5)))
        #os.mkdir(par_dir)

        V1V2_dir=os.path.join(task.args[3], "multi_"+str(round(i, 5)), "V1V2")
        V4V5_dir=os.path.join(task.args[3], "multi_"+str(round(i, 5)), "V4V5")
        FL_dir=os.path.join(task.args[3], "multi_"+str(round(i, 5)), "FL")
        
        #os.mkdir(V1V2_dir)
        #os.mkdir(V4V5_dir)
        #os.mkdir(FL_dir)

        ##call assignment script
        assignment_script=os.path.join(task.args[4], "parathaa", "utility", "tax.assign_parallel.R")
        #need to fix the args here.
        util1=os.path.join(task.args[4], "parathaa", "utility", "nearest_neighbours_parallel.R")
        V1V2_jplace=os.path.join(task.args[3], "init", "V1V2", "merged_sub.jplace")
        V1V2_tree=os.path.join(task.args[5], "resultTree_bestThresholds.RData")
        V1V2_scores=os.path.join(task.args[5], "optimal_scores.RData")

        V4V5_jplace=os.path.join(task.args[3], "init", "V4V5", "merged_sub.jplace")
        V4V5_tree=os.path.join(task.args[7], "resultTree_bestThresholds.RData")
        V4V5_scores=os.path.join(task.args[7], "optimal_scores.RData")


        FL_jplace=os.path.join(task.args[3], "init", "FL", "merged_sub.jplace")
        FL_tree=os.path.join(task.args[8], "resultTree_bestThresholds.RData")
        FL_scores=os.path.join(task.args[8], "optimal_scores.RData")

        callV1V2 = "Rscript " + assignment_script + " -j " + V1V2_jplace + " -o " + V1V2_dir + " -t " + V1V2_tree + " -s " + V1V2_scores + " --threads " + task.args[6] + " --util1 " + util1 + " -m " + str(i) 
        callV4V5 = "Rscript " + assignment_script + " -j " + V4V5_jplace + " -o " + V4V5_dir + " -t " + V4V5_tree + " -s " + V4V5_scores + " --threads " + task.args[6] + " --util1 " + util1 + " -m " + str(i)
        callFL = "Rscript " + assignment_script + " -j " + FL_jplace + " -o " + FL_dir + " -t " + FL_tree + " -s " + FL_scores + " --threads " + task.args[6] + " --util1 " + util1 + " -m " + str(i)
        print(callV1V2)
        print(callV4V5)
        print(callFL)
        #subprocess.call(callV1V2, shell=True)
        #subprocess.call(callV4V5, shell=True)
        subprocess.call(callFL, shell=True)

#workflow.add_task(
#    run_initial,
#    args=[args.output, args.treeFilesB, args.queryB, args.threads, args.treeFilesA, args.queryA, args.FLtreeFile, args.FLquery]
#)

workflow.add_task(
        run_bench,
        args=[args.minMulti, args.maxMulti, args.step, args.output, args.paraDir,  args.treeFilesB, args.threads, args.treeFilesA, args.FLtreeFile]
)


workflow.go()
