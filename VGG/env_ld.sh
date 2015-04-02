#!/bin/csh
setenv WORKDIR /work/halla/solid/yez/dvcs/VGG
setenv LD_LIBRARY_PATH ${WORKDIR}/envr:${LD_LIBRARY_PATH}

${WORKDIR}/dvcsjlab 
