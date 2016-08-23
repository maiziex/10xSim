# 10xSim

"doc/hist_reads_start_diff_abs.txt":  corrected reads start per fragment.

"src/process_bam4_3.py": plot fragments distribution and fit curve, plot corrected reads start distribution. 
For the input "molecule3.pickle" of this script, check the google drive shared folder "SharewithEric_10xSim"


For the "fragments_upto200kb.png", it fits lognorm, with below parameters:
loc  = 0
sigma = shape = 2.3515642347859629
mu = np.log(scale) = np.log(2.3352378273006389) = 0.84811373916785049

