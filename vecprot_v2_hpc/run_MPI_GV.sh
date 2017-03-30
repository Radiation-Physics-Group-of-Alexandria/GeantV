scl enable devtoolset-3 bash
export AR=/opt/rh/devtoolset-3/root/usr/bin/ar
export LD=/opt/rh/devtoolset-3/root/usr/bin/ld
export PATH=/usr/lib64/mpich/bin/:$PATH
cd ../builds/bin && mpirun -n 2 ./CMSApp -E /afs/cern.ch/work/g/geant/GeantV-master/builds/pp14TeVminbias.root -e 5 -f /afs/cern.ch/work/g/geant/GeantV-master/builds/fstate_FTFP_BERT_G496p02_1mev.root -g /afs/cern.ch/work/g/geant/GeantV-master/builds/cms2015.root -b 1 -t 2 -x /afs/cern.ch/work/g/geant/GeantV-master/builds/xsec_FTFP_BERT_G496p02_1mev.root -p 1
