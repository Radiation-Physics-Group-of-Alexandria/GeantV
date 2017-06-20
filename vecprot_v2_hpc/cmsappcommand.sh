RES=`pwd`/resources
mpirun -n 4 ./geant/build/bin/CMSApp -E $RES/eventFiles.txt -f $RES/fstate_FTFP_BERT_G496p02_1mev.root -g $RES/cms2015.root -x $RES/xsec_FTFP_BERT_G496p02_1mev.root -p 1 -b 2 -t 2 -v 1
