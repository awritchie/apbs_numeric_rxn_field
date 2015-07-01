#! /bin/bash

#if [ -f apbs_numeric_rxn_field ] ; then rm apbs_numeric_rxn_field ; fi
make -j 7
if [ -f apbs_numeric_rxn_field  ]; then
    ./apbs_numeric_rxn_field ../test/junk.0.pqr ../test/*.78.txt ../test/*.1.txt $1
    fi
