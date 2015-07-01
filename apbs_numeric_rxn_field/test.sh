#! /bin/bash

clear
if [ -f apbs_numeric_rxn_field ] ; then rm apbs_numeric_rxn_field ; fi
make -j 7
if [ -f apbs_numeric_rxn_field  ]; then
    ./apbs_numeric_rxn_field test/state234.pqr test/state234.78.txt test/state234.1.txt $1
    fi
