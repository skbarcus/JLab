# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
if [ "$PS1" ]; then
    PS1="[\u@\h \W]\\$ "
fi

alias ls='ls -h --color'
alias mv='/bin/mv -iv'
alias cp='/bin/cp -iv'
alias rm='/bin/rm -Iv'

if [ -z "$CODA" ]; then
    # CODA 3 using coda_scripts : SESSION=hcal 03feb2021
    export SESSION=hcal
    export EXPID=hcal

    source ${HOME}/coda_scripts/setupCODA3.bash

    # LINUXVME for hcal_pulser
    export LINUXVME=${HOME}/linuxvme
    export LINUXVME_INC=${HOME}/linuxvme/include
    export LINUXVME_LIB=${HOME}/linuxvme/${MACHINE}/lib
    export LINUXVME_BIN=${HOME}/linuxvme/${MACHINE}/bin

    export FADC250_PARAMS=${HOME}/cfg

    export LD_LIBRARY_PATH=.:${LINUXVME_LIB}:${LD_LIBRARY_PATH}
    export PATH=.:${LINUXVME_BIN}:${HOME}/bin:${PATH}
fi
