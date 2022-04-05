#!/bin/bash


#DATADIR=/home/Jnorce/Nikita/nbernier/Cad/128Cd/sept2018/128Cd_on
DATADIR=.

SCOPE=./do_it11

##CAL=../Cal04496.cal

#CT=ct_correction.cal

#CUTS=bananas2.cuts


for (( RUN=4496; RUN<=4504; RUN++ )) do
        for (( SUBRUN=0; SUBRUN<=999; SUBRUN++ )) do     
                if [ $SUBRUN -lt 10 ]; then 
                #       echo $DATADIR/run0$RUN-00$SUBRUN.mid
                        if [ -e $DATADIR/analysis0${RUN}_00${SUBRUN}.root ]
                        then
                                $SCOPE analysis0${RUN}_00${SUBRUN}.root Cal0${RUN}.cal
                        fi
                elif [ $SUBRUN -lt 100 ]; then
                #       echo $DATADIR/run0$RUN-0$SUBRUN.mid
                        if [ -e $DATADIR/analysis0${RUN}_0${SUBRUN}.root ]
                        then
                                $SCOPE analysis0${RUN}_0${SUBRUN}.root Cal0${RUN}.cal
                        fi
                else
                        if [ -e $DATADIR/analysis0${RUN}_${SUBRUN}.root ]
                        then
                                $SCOPE analysis0${RUN}_${SUBRUN}.root Cal0${RUN}.cal
                        fi
                fi

        done
done
