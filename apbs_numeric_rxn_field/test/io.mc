##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Aug 22 09:53:02 2014
##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path state234.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing glen...
PBEparm_parseToken:  trying glen...
MGparm_parseToken:  trying glen...
NOsh_parseMG:  Parsing gcent...
PBEparm_parseToken:  trying gcent...
MGparm_parseToken:  trying gcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 0.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing glen...
PBEparm_parseToken:  trying glen...
MGparm_parseToken:  trying glen...
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing gcent...
PBEparm_parseToken:  trying gcent...
MGparm_parseToken:  trying gcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 0.
NOsh:  nlev = 5, dime = (193, 193, 193)
NOsh: Done parsing ELEC section (nelec = 2)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 1613 atoms
Valist_getStatistics:  Max atom coordinate:  (23.45, 17.74, 15)
Valist_getStatistics:  Min atom coordinate:  (-21.43, -18.9, -15.57)
Valist_getStatistics:  Molecule center:  (1.01, -0.58, -0.285)
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 0 (1)
NOsh_setupCalc:  Mapping ELEC statement 1 (2) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 27.8942
Vpbe_ctor2:  solute dimensions = 47.606 x 39.275 x 32.795
Vpbe_ctor2:  solute charge = -3.9989
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 1.12452
Vpbe_ctor2:  Debye length = 0.889267
Vpbe_ctor2:  zkappa2 = 1.26455
Vpbe_ctor2:  zmagic = 6999.55
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 65 table
Vclist_ctor2:  Using 75 x 75 x 65 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (57.66, 49.42, 43.35)
Vclist_setupGrid:  Grid lower corner = (-27.82, -25.29, -21.96)
Vclist_assignAtoms:  Have 3990496 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.158476
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.458749e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.595990e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 8.330160e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.517979e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 7.260479e-02
Vprtstp: contraction number = 7.260479e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.284228e-02
Vprtstp: contraction number = 1.768793e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.733182e-03
Vprtstp: contraction number = 2.128268e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 6.588940e-04
Vprtstp: contraction number = 2.410721e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.686051e-04
Vprtstp: contraction number = 2.558911e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.408724e-05
Vprtstp: contraction number = 2.614823e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.172565e-05
Vprtstp: contraction number = 2.659647e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.166485e-06
Vprtstp: contraction number = 2.700477e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 8.708630e-07
Vprtstp: contraction number = 2.750252e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 2.871130e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 3.908760e+00
Vpmg_setPart:  lower corner = (-118.99, -120.58, -120.285)
Vpmg_setPart:  upper corner = (121.01, 119.42, 119.715)
Vpmg_setPart:  actual minima = (-118.99, -120.58, -120.285)
Vpmg_setPart:  actual maxima = (121.01, 119.42, 119.715)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 3.318044271481E+03 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 4.312000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.900000e-05
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 27.8942
Vpbe_ctor2:  solute dimensions = 47.606 x 39.275 x 32.795
Vpbe_ctor2:  solute charge = -3.9989
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 1.12452
Vpbe_ctor2:  Debye length = 0.889267
Vpbe_ctor2:  zkappa2 = 1.26455
Vpbe_ctor2:  zmagic = 6999.55
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 65 table
Vclist_ctor2:  Using 75 x 75 x 65 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (57.66, 49.42, 43.35)
Vclist_setupGrid:  Grid lower corner = (-27.82, -25.29, -21.96)
Vclist_assignAtoms:  Have 3990496 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 11.505, 0.675, -1.115
VPMG::focusFillBound -- New mesh maxs = 21.505, 10.675, 8.885
VPMG::focusFillBound -- Old mesh mins = -118.99, -120.58, -120.285
VPMG::focusFillBound -- Old mesh maxs = 121.01, 119.42, 119.715
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (11.505, 0.675, -1.115)
Vpmg_setPart:  upper corner = (21.505, 10.675, 8.885)
Vpmg_setPart:  actual minima = (-118.99, -120.58, -120.285)
Vpmg_setPart:  actual maxima = (121.01, 119.42, 119.715)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (11.505, 0.675, -1.115)
VPMG::extEnergy    Disj part upper corner = (21.505, 10.675, 8.885)
VPMG::extEnergy    Old lower corner = (-118.99, -120.58, -120.285)
VPMG::extEnergy    Old upper corner = (121.01, 119.42, 119.715)
Vpmg_qmEnergy:  Calculating linear energy
VPMG::extEnergy: extQmEnergy = 42.3113 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 3230.5 kT
VPMG::extEnergy: extDiEnergy = 1563.31 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.161661
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.083945e+01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (193, 193, 193)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.320242e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (097, 097, 097)
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.988801e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 3.489840e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.512770e-01
Vprtstp: contraction number = 1.512770e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.543146e-02
Vprtstp: contraction number = 1.020080e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 1.611038e-03
Vprtstp: contraction number = 1.043996e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.731936e-04
Vprtstp: contraction number = 1.075044e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.861367e-05
Vprtstp: contraction number = 1.074732e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 2.034602e-06
Vprtstp: contraction number = 1.093068e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 2.202974e-07
Vprtstp: contraction number = 1.082754e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.825876e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 2.691070e+01
Vpmg_setPart:  lower corner = (11.505, 0.675, -1.115)
Vpmg_setPart:  upper corner = (21.505, 10.675, 8.885)
Vpmg_setPart:  actual minima = (11.505, 0.675, -1.115)
Vpmg_setPart:  actual maxima = (21.505, 10.675, 8.885)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 5.775026439731E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 3.231000e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-05
routines:  Opening virtual socket...
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 5.348939e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Fri Aug 22 09:53:56 2014
##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path state234.pqr
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing glen...
PBEparm_parseToken:  trying glen...
MGparm_parseToken:  trying glen...
NOsh_parseMG:  Parsing gcent...
PBEparm_parseToken:  trying gcent...
MGparm_parseToken:  trying gcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 0.
NOsh:  nlev = 4, dime = (97, 97, 97)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing glen...
PBEparm_parseToken:  trying glen...
MGparm_parseToken:  trying glen...
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing gcent...
PBEparm_parseToken:  trying gcent...
MGparm_parseToken:  trying gcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing ion...
PBEparm_parseToken:  trying ion...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 0.
NOsh:  nlev = 5, dime = (193, 193, 193)
NOsh: Done parsing ELEC section (nelec = 2)
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 1613 atoms
Valist_getStatistics:  Max atom coordinate:  (23.45, 17.74, 15)
Valist_getStatistics:  Min atom coordinate:  (-21.43, -18.9, -15.57)
Valist_getStatistics:  Molecule center:  (1.01, -0.58, -0.285)
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 0 (1)
NOsh_setupCalc:  Mapping ELEC statement 1 (2) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 27.8942
Vpbe_ctor2:  solute dimensions = 47.606 x 39.275 x 32.795
Vpbe_ctor2:  solute charge = -3.9989
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127327
Vpbe_ctor2:  Debye length = 7.85379
Vpbe_ctor2:  zkappa2 = 1.26455
Vpbe_ctor2:  zmagic = 6999.55
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 65 table
Vclist_ctor2:  Using 75 x 75 x 65 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (57.66, 49.42, 43.35)
Vclist_setupGrid:  Grid lower corner = (-27.82, -25.29, -21.96)
Vclist_assignAtoms:  Have 3990496 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.162848
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.462445e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 097)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.532750e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 8.368620e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 2.520680e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 7.988244e-02
Vprtstp: contraction number = 7.988244e-02
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.350987e-02
Vprtstp: contraction number = 1.691220e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 3.151785e-03
Vprtstp: contraction number = 2.332949e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 1.000954e-03
Vprtstp: contraction number = 3.175832e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 4.376529e-04
Vprtstp: contraction number = 4.372357e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 2.400681e-04
Vprtstp: contraction number = 5.485354e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.474340e-04
Vprtstp: contraction number = 6.141340e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 9.465042e-05
Vprtstp: contraction number = 6.419851e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 6.172368e-05
Vprtstp: contraction number = 6.521227e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 4.042831e-05
Vprtstp: contraction number = 6.549887e-01
Vprtstp: iteration = 11
Vprtstp: relative residual = 2.650322e-05
Vprtstp: contraction number = 6.555608e-01
Vprtstp: iteration = 12
Vprtstp: relative residual = 1.737687e-05
Vprtstp: contraction number = 6.556514e-01
Vprtstp: iteration = 13
Vprtstp: relative residual = 1.139431e-05
Vprtstp: contraction number = 6.557172e-01
Vprtstp: iteration = 14
Vprtstp: relative residual = 7.471506e-06
Vprtstp: contraction number = 6.557224e-01
Vprtstp: iteration = 15
Vprtstp: relative residual = 4.899228e-06
Vprtstp: contraction number = 6.557217e-01
Vprtstp: iteration = 16
Vprtstp: relative residual = 3.212552e-06
Vprtstp: contraction number = 6.557261e-01
Vprtstp: iteration = 17
Vprtstp: relative residual = 2.106549e-06
Vprtstp: contraction number = 6.557244e-01
Vprtstp: iteration = 18
Vprtstp: relative residual = 1.381314e-06
Vprtstp: contraction number = 6.557237e-01
Vprtstp: iteration = 19
Vprtstp: relative residual = 9.057600e-07
Vprtstp: contraction number = 6.557235e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 6.065339e+00
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.102559e+00
Vpmg_setPart:  lower corner = (-118.99, -120.58, -120.285)
Vpmg_setPart:  upper corner = (121.01, 119.42, 119.715)
Vpmg_setPart:  actual minima = (-118.99, -120.58, -120.285)
Vpmg_setPart:  actual maxima = (121.01, 119.42, 119.715)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 7.889894501085E+02 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 4.409000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.700000e-05
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 27.8942
Vpbe_ctor2:  solute dimensions = 47.606 x 39.275 x 32.795
Vpbe_ctor2:  solute charge = -3.9989
Vpbe_ctor2:  bulk ionic strength = 0.15
Vpbe_ctor2:  xkappa = 0.127327
Vpbe_ctor2:  Debye length = 7.85379
Vpbe_ctor2:  zkappa2 = 1.26455
Vpbe_ctor2:  zmagic = 6999.55
Vpbe_ctor2:  Constructing Vclist with 75 x 75 x 65 table
Vclist_ctor2:  Using 75 x 75 x 65 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 2.5 max radius
Vclist_setupGrid:  Grid lengths = (57.66, 49.42, 43.35)
Vclist_setupGrid:  Grid lower corner = (-27.82, -25.29, -21.96)
Vclist_assignAtoms:  Have 3990496 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 254.469
Vacc_storeParms:  Using 2584-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = 11.505, 0.675, -1.115
VPMG::focusFillBound -- New mesh maxs = 21.505, 10.675, 8.885
VPMG::focusFillBound -- Old mesh mins = -118.99, -120.58, -120.285
VPMG::focusFillBound -- Old mesh maxs = 121.01, 119.42, 119.715
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (11.505, 0.675, -1.115)
Vpmg_setPart:  upper corner = (21.505, 10.675, 8.885)
Vpmg_setPart:  actual minima = (-118.99, -120.58, -120.285)
Vpmg_setPart:  actual maxima = (121.01, 119.42, 119.715)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (11.505, 0.675, -1.115)
VPMG::extEnergy    Disj part upper corner = (21.505, 10.675, 8.885)
VPMG::extEnergy    Old lower corner = (-118.99, -120.58, -120.285)
VPMG::extEnergy    Old upper corner = (121.01, 119.42, 119.715)
Vpmg_qmEnergy:  Calculating linear energy
VPMG::extEnergy: extQmEnergy = 1.49959 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 779.377 kT
VPMG::extEnergy: extDiEnergy = 386.646 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 1.158250
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 2.084671e+01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (193, 193, 193)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.337501e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (097, 097, 097)
Vbuildops: Galer: (049, 049, 049)
Vbuildops: Galer: (025, 025, 025)
Vbuildops: Galer: (013, 013, 013)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 6.934938e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 3.806435e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.314264e-01
Vprtstp: contraction number = 1.314264e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.689042e-02
Vprtstp: contraction number = 1.285162e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.935155e-03
Vprtstp: contraction number = 1.737764e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 7.959560e-04
Vprtstp: contraction number = 2.711802e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 2.726126e-04
Vprtstp: contraction number = 3.424971e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 1.005141e-04
Vprtstp: contraction number = 3.687066e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 3.784660e-05
Vprtstp: contraction number = 3.765304e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 1.436400e-05
Vprtstp: contraction number = 3.795320e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 5.469193e-06
Vprtstp: contraction number = 3.807571e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 2.090923e-06
Vprtstp: contraction number = 3.823093e-01
Vprtstp: iteration = 11
Vprtstp: relative residual = 8.012924e-07
Vprtstp: contraction number = 3.832242e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 2.905344e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 3.766697e+01
Vpmg_setPart:  lower corner = (11.505, 0.675, -1.115)
Vpmg_setPart:  upper corner = (21.505, 10.675, 8.885)
Vpmg_setPart:  actual minima = (11.505, 0.675, -1.115)
Vpmg_setPart:  actual maxima = (21.505, 10.675, 8.885)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 5.526242816952E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 3.522200e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 9.000000e-06
routines:  Opening virtual socket...
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 6.745524e+01
