#!/bin/bash
#
#BSUB -P geantvmpi                      # project code
#BSUB -J geantvmpi                      # job name
#BSUB -W 00:10                          # wall-clock time (hrs:mins)
#BSUB -n 1                              # request one slot for this job
#BSUB -q 1nh                            # queue
#BSUB -e errors.%J.openmp               # error file name in which %J is replaced by the job ID
#BSUB -o output.%J.openmp               # output file name in which %J is replaced by the job ID
 
/usr/lib64/mpich/bin/mpirun bin/CMSApp -e 10 -f fstate_FTFP_BERT_G496p02_1mev.root -g cms2015.root -b 1 -t 2 -x xsec_FTFP_BERT_G496p02_1mev.root -p 1
