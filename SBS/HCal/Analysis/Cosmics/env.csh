## Apparently system got updated to 7.5, but does not have the
## latest compiled versions of ROOT for that release, so use the
## old one from 7.5
#setenv ROOTSYS /data/jlab_software/2.2/Linux_RedHat7.5-x86_64-gcc5.3.0/root/6.12.06
#setenv EVIO /site/12gev_phys/2.2/Linux_RedHat7.5-x86_64-gcc4.8.5/evio/5.1
#setenv EVIO_INCDIR /data/jlab_software/2.0/evio/5.1/src/libsrc
#source $ROOTSYS/bin/thisroot.csh
source /site/12gev_phys/softenv.csh 2.3
setenv ANALYZER ${HOME}/analyzer/analyzer-1.6.6
setenv LD_LIBRARY_PATH ${ANALYZER}:${LD_LIBRARY_PATH}
setenv PATH ${ANALYZER}:${PATH}
#setenv DB_DIR ${HOME}/test_fadc
setenv EVIO ${ANALYZER}/hana_decode
setenv EVIO_INCDIR ${ANALYZER}/hana_decode
setenv EVIO_LIBDIR ${ANALYZER}/hana_decode
setenv HCAL_ROOTFILES rootfiles
setenv HCAL_DATA /home/daq/data
