#!/bin/bash

#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1     

#### Modules ####
#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB
module load amber

##### FUNCTIONS ####
#Function adapted from 'submit_duck_smd_gpu.csh' of the DUck std pipeline
prepare_duck_and_launch(){
   nustart=$1
   nuend=$2
   temp=$3
   nu=$nustart
   while (($nu <= $nuend)); do
      if [ "$temp" == '300K' ]; then
         dir=DUCK_${nu}
         mkdir $dir
         cd $dir
         if [ "$nu" == "0" ]; then
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi
         cd ..  
      elif [ "$temp" == '325K' ]; then
         dir=DUCK_325K_${nu}
         mkdir $dir
         cd $dir
	 if [ "$nu" == "0" ]; then
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../HMR_system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi 
         cd ..
      fi
      nu=$((nu+1))
   done

}

# Function to check if WQB is lower than 0.1 using getWqbValues.py
# getWqbValues.py is a script from Maciej modified.
# We use this one instead of the R version, as R is not available in the IQTC

check_WQB(){
   wqb_limit=$1
   lowest_wqb=$(python getWqbValues.py)
   echo $lowest_wqb > wqb.log
   are_we_above_the_limit=$(echo "$lowest_wqb < $wqb_limit" | bc )
   if [ "$are_we_above_the_limit" == "1" ]; then
      echo "Wqb lower than ${wqb_limit}, stoping DUck run"
      cp -r ./* $LIG_TARGET/
      exit 
   fi
}



#### PARAMS ####
replicas=10
min_wqb=6
    
#### Runing Duck ####
# Minimization
pmemd.cuda -O -i 1_min.in -o min.out -p HMR_system_complex.prmtop -c system_complex.inpcrd -r min.rst -ref system_complex.inpcrd

# Centering box to allow iwrap, not being centered to 0,0,0 bricks the system when using iwrap and nmr restraints together
echo -e "trajin min.rst\ncenter\ntrajout min.rst\nrun\nquit" | cpptraj -p HMR_system_complex.prmtop

# Heating & Equilibration
pmemd.cuda -O -i 2_heating150.in -o 2_heating150.out -p HMR_system_complex.prmtop -c min.rst -r  2_heating150.rst -x 2_heating150.nc -ref min.rst
pmemd.cuda -O -i 2_heating200.in -o 2_heating200.out -p HMR_system_complex.prmtop -c 2_heating150.rst -r 2_heating200.rst -x 2_heating200.nc -ref 2_heating150.rst
pmemd.cuda -O -i 2_heating250.in -o 2_heating250.out -p HMR_system_complex.prmtop -c 2_heating200.rst -r 2_heating250.rst -x 2_heating250.nc -ref 2_heating200.rst
pmemd.cuda -O -i 2_heating300.in -o 2_heating300.out -p HMR_system_complex.prmtop -c 2_heating250.rst -r 2_heating300.rst -x 2_heating300.nc -ref 2_heating250.rst
pmemd.cuda -O -i 3_eq.in -o 3_eq.out -p HMR_system_complex.prmtop -c 2_heating300.rst -r 3_eq.rst -x 3_eq.nc -ref 2_heating300.rst -e 3_eq.ene

#Launch DUck 0 and check wqb
prepare_duck_and_launch 0 0 300K
check_WQB $min_wqb

#Launch DUck_325K 0 and check wqb
prepare_duck_and_launch 0 0 325K
check_WQB $min_wqb

#For each replica wanted do: MD, prepare SMD & launch SMD
for ((i=1;i<=$replicas;++i)); do
   if [ "$i" == "1" ]; then
      pmemd.cuda -O -i md.in -o md1.out -p HMR_system_complex.prmtop -c 3_eq.rst -r md1.rst -x md1.nc -ref 3_eq.rst
   else
      pmemd.cuda -O -i md.in -o md${i}.out -p HMR_system_complex.prmtop -c md$((i-1)).rst -r md${i}.rst -x md${i}.nc -ref 3_eq.rst
   fi

   prepare_duck_and_launch $i $i 300K
   prepare_duck_and_launch $i $i 325K
   check_WQB $min_wqb

done
    
find -type f ! -regex ".*\(dat\|prmtop\|inpcrd\|in\|dist.*rst\|yaml)$" -delete


exit
