# Some usefull commands that demostrate the use of chapeau tools

Convert old formats

    ./chapconv ~/rs_unc/Colaboraciones/Luca/cineca/IscrC_MDKCOP/OTFP/2_Cryo/1_DEKA_K+1/OTFP/REPLICA1/namd/RUN22.ch0 
    mv converted.ch cryo.ch
    ./chapconv ~/rs_unc/Colaboraciones/Luca/cineca/IscrC_MDKCOP/OTFP/1_Model/1_DEKA_K+1/otfp/REPLICA1/namd/RUN17.ch0 
    mv converted.ch model.ch

Extend the CV space area (original 4 17 4 17 in each case) and create a offset between chapeaus

    ./chapcrop cryo.ch 4 18.40076 4 18.40076
    mv croped.ch croped_cryo.ch
    ./chapcrop model.ch 2.59924 17 2.59924 17
    mv croped.ch croped_model.ch

Attempt to add

    ./chapadd croped_*ch

attempt fail, need to adjust a little the size of one chapeau to have the same size (223 knots)

    ./chapcrop model.ch 2.59 17 2.59 17
    ./chapcrop model.ch 2.58 17 2.58 17
    ./chapcrop model.ch 2.55 17 2.55 17
    mv croped.ch croped_model.ch

Add chapeaus objects

    ./chapadd croped_*ch

