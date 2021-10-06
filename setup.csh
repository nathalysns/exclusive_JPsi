
export LANG=C

source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh -n  
#export ROOT_INCLUDE_PATH=/work/eic3/users/nathaly/software/test/macros/common:$ROOT_INCLUDE_PATH 
#export ROOT_INCLUDE_PATH=/work/eic/users/nathaly/software/test/macros/common:$ROOT_INCLUDE_PATH

export ECCE=/work/eic/users/nathaly/analysis/exclusive/dst
export MYINSTALL=$ECCE/install
export LD_LIBRARY_PATH=$MYINSTALL/lib:$LD_LIBRARY_PATH

echo $MYINSTALL

#source /cvmfs/eic.opensciencegrid.org/default/opt/fun4all/core/bin//setup_local.sh  $MYINSTALL
source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/setup_local.sh  $MYINSTALL
#export ROOT_INCLUDE_PATH=/work/eic/users/nathaly/test/software/macros/common:$ROOT_INCLUDE_PATH

export ROOT_INCLUDE_PATH=$MYINSTALL/include:$ROOT_INCLUDE_PATH

