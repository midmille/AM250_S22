"""
This is a python file for running the latency study for different sized messages as well as reading in and
plotting the results.
"""

import subprocess
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt


def Run_Latency(Ns, in_file, out_file):
    """
    This function does the following:
        1) Loops over Ns
        2) Edits the input file
        3) Executes mpirun latency
        4) Reads in a stores the result.

    """

    ## [The resulting time array for each run.]
    time = np.zeros(len(Ns))

    for k, N in enumerate(Ns): 
        
        ## [Make the input file.]
        f = open(in_file, 'wt')
        f.write(f"{N}")
        f.close()

        ## [Run the fortran latency mpi program.]
        out = subprocess.run(['mpirun', '-np', '2', 'latency'])
        
        ## [Read the output from file]
        f = open(out_file, 'r')
        val = f.read()

        ## [If the return value has E then it is scientific notation.]
        if 'E' in val: 
            ## [Removing the scientific notation and making flaot.]
            vals = val.split('E')
            val_flt = float(vals[0])*10**(float(vals[1]))
            f.close()
        else: 
            val_flt = float(val)

        ## [Store the result in the time array.]
        time[k] = val_flt

        ## [Delete the output file from this N.]
        os.remove(out_file)


    return time


def Plot_Latency(Ns, save_file, N_loop): 
    """
    This loads the resulting time array from pickle, solves for the linear fit, and plots the 
    result.

    """

    ## [load the time array.]
    time = pickle.load(open(save_file, 'rb'))
    ## [The time per message is the total time divided by 2*N_loop.]
    time_msg = time / (2*N_loop)

    ## [Linear rgression.]
    m,b = np.polyfit(Ns, time_msg, 1)

    ## [plotting the result.]
    fig, ax = plt.subplots()
    ## [plot the best fit line.]
    x = np.linspace(0, Ns[-1], 100)
    ax.plot(x, m*x +b, label = f'Best Fit Line: y = {m}x + {b}')
    ax.plot(Ns, time_msg, 'ro', fillstyle='none', label = "MPI_WTIME Data")
    ax.set_ylabel("Message Time [s]")
    ax.set_xlabel("Message Size [N]")
    ax.set_title("Latency Study")
    ax.grid()
    ax.legend()

    fig.show()
    
    
if __name__ == '__main__':

    ## [Problem parameters.]
    ## ---------------------
    ## [The number of loops for the ping ponging]
    N_loop = 10000 
    ## [The name of the file to save the time array into.]
    save_file = 'save_time.p'
    ## [The name of the input file.]
    in_file = 'input.txt'
    ## [The name of the output file.]
    out_file = 'output.txt'
    ## [The array of messages sizes to loop over.]
    Ns = np.arange(1, 110000, 10000)
    ## [Set True to run mpi latency program.]
    run = False
    ## [Set True to plot figure from pickle file.]
    plot = True


    if run: 
        ## [Run the latency program for each message size and return the result.]
        time = Run_Latency(Ns, in_file, out_file)
    
        ## [Save the result.]
        pickle.dump(time, open(save_file, 'wb'))

    if plot: 
        ## [Plot the resulting linear regression.]
        Plot_Latency(Ns, save_file, N_loop)
