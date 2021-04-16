if [ $4 == "4" ] 
then
    mpirun -np 4 bin/crout_mpi $1 $2 $3
else
    bin/crout $1 $2 $3 $4
fi