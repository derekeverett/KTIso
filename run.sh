#export OMP_NUM_THREADS=4 #set this to the number of threads on CPU
#set OMP_NUM_THREADS 4
# Run this script only in the conda environment with the gstools library.
cd random_fields
python generate_random_fields.py $1 $2 $3 $4
open generated_fields.png
cd ..
./RunWrapper
#cd validations
#python compare_ktiso_to_fs__after_1_5_fmc.py
#open /Users/dananjayaliyanage/git/KTIso/validations/compare_ktiso_to_fs__after_1_5_fmc/recent_diff.png
#open /Users/dananjayaliyanage/git/KTIso/validations/compare_ktiso_to_fs__after_1_5_fmc/recent_ratio.png
