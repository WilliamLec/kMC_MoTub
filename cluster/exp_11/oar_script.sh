#!/bin/bash
#OAR -n efficiency_EXP11|@L@|@ite@|
#OAR -l /nodes=1/core=1,walltime=1000:20:30
#OAR -t besteffort
#OAR --stdout @L@_@ite@.out
#OAR --stderr @L@_@ite@.err
python3 kMC_main.py @seed@ @L@ @ite@
