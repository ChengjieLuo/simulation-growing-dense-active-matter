./cmp.sh

uc=0.0005
alphamax=0.01
dt=0.005

seed=666
seedalpha=343

tmax=100000

uc=0.0005
configfile=config_${uc}_${seed}_${seedalpha}_${tmax}
logfile=log_${uc}_${seed}_${seedalpha}_${tmax}
./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile  &

uc=0.0001
configfile=config_${uc}_${seed}_${seedalpha}_${tmax}
logfile=log_${uc}_${seed}_${seedalpha}_${tmax}
./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile  &

uc=0.001
configfile=config_${uc}_${seed}_${seedalpha}_${tmax}
logfile=log_${uc}_${seed}_${seedalpha}_${tmax}
./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile  &



#tmax=20000
#configfile=config_${seed}_${seedalpha}_${tmax}
#logfile=log_${seed}_${seedalpha}_${tmax}
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile &
#
#tmax=100000
#configfile=config_${seed}_${seedalpha}_${tmax}
#logfile=log_${seed}_${seedalpha}_${tmax}
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile &

#
#tmax=6000
#configfile=config_${seed}_${seedalpha}_${tmax}
#logfile=log_${seed}_${seedalpha}_${tmax}
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile &
#
#tmax=7000
#configfile=config_${seed}_${seedalpha}_${tmax}
#logfile=log_${seed}_${seedalpha}_${tmax}
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile &
#
#tmax=8000
#configfile=config_${seed}_${seedalpha}_${tmax}
#logfile=log_${seed}_${seedalpha}_${tmax}
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile} $seed $seedalpha > $logfile &


#seed=464
#seedalpha=789
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile}_${seed}_${seedalpha} $seed $seedalpha > log_${seed}_${seedalpha}
#
#
#seed=555
#seedalpha=888
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile}_${seed}_${seedalpha} $seed $seedalpha > log_${seed}_${seedalpha}
#
#seed=390
#seedalpha=469
#./cell_simulation.exe $uc $alphamax $dt $tmax ${configfile}_${seed}_${seedalpha} $seed $seedalpha > log_${seed}_${seedalpha}


