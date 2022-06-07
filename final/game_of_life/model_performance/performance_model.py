"""
This python file is for the modeling and its visualization of the performance of the game of life. 
Both experimental and theoretical anaylsis are performed. 
See the reprot for all details.

This simple analysis is performed under the adssusmptions that N = Nx = Ny.

First compile the 'main.exe' executable. 
    
    For the study it is important that the following parameters are set to the following in 
    the ../src/main.f90 file.
        
        ## [This ensures that actuall game of life data is not written during the performance analysis.]
        woflag = .FALSE.

        ## [Set this flag to the paralleization method the performance is done for.]
        pflag = [Options: 'serial', 'cols', 'rows', 'tile']

Then ensure that all parameters in this file's "if __name__ == '__main__'" match the 
the parameters in the ../src/main.f90 file.

Then run the performance model. 

The output will be save to pickle files in this working directory.

TO VISUALIZE WITHOUT RUNNING: 

Simply set run flag to False, and
set the plot flag to True
"""

import subprocess
import numpy as np
import os
import pickle
import matplotlib.pyplot as plt


def Run_Performance(Ns, Np, Nxin_file, Nyin_file, out_file):
    """
    This function does the following:
        1) Loops over Ns
        2) Edits the input files for Nx and Ny
        3) Executes mpi run for the desired number of procs.
        4) Reads the time and stored result. 
    
    Parameters
    ----------
    Ns: 1-D Array, Integer
        The list of dom,ain sizes to run.
    Np: Integer
        The number of procs to use.
    Nxin_file: String
        The name and location of the Nx input file.
    Nyin_file: String
        The name and location of the Ny input file.
    out_file: String
        The name of the elapsed time file to be read in.

    """

    ## [The resulting time array for each run.]
    time = np.zeros(len(Ns))

    for k, N in enumerate(Ns): 
        
        ## [Make the input file.]
        f = open(Nxin_file, 'wt')
        f.write(f"{N}")
        f.close()

        ## [Make the input file.]
        f = open(Nyin_file, 'wt')
        f.write(f"{N}")
        f.close()

        ## [Get current dir,]
        cwd = os.getcwd()
        ## [cd to main.exe workdir.] 
        os.chdir('..')
        ## [Run the game of life program.]
        out = subprocess.run(['mpirun', '-np', f'{Np}', 'main.exe'])
        ## [cd to back.] 
        os.chdir(cwd)
        
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

def Get_tc(Nt):
    """
    This algorithm get the t_c variable from serial run. It calculates the t_c for each N=Nx=Ny in the run. Then it averages over them all.
    """

    ## [Download the serial run.]
    Ns_time = pickle.load(open('save_time_serial.p', 'rb'))

    ## [Split it.]
    Ns = Ns_time[:,0]
    time = Ns_time[:,1]

    ## [Divide the time by the number of iterations.]
    tpi = time/Nt

    tc = tpi / (Ns**2)

    tc_avg = np.mean(tc)

    return tc_avg

def Analytical_Model()

if __name__ == '__main__':

    ## ---------------------
    ## [Problem parameters.]
    ## ---------------------
    """
    Please ensure that these parameters match those you set in ../src/main.f90
    """
    Nt = 80 
    pflag = 'serial'
    ## [The name of the input file.]
    Nxin_file = '../input/Nx_in.dat'
    Nyin_file = '../input/Ny_in.dat'
    ## [The name of the time output file.]
    out_file = '../output/time.dat'

    ## [These are performance model specific parameters.]
    ## [The name of the file to save the time array into.]
    if pflag == 'serial': 
        save_file = f'save_time_{pflag}.p'
    else: 
        save_file = f'save_time_{pflag}_Np{Np}.p'
    ## [The array of messages sizes to loop over.]
    Ns = np.arange(1, 1000, 100)
    ## [The number of processors to use.]
    Np = 2
    ## [Set True to run mpi latency program.]
    run = False
    ## [Set True to plot figure from pickle file.]
    plot = True

    if run: 
        ## [Run the latency program for each message size and return the result.]
        time = Run_Performance(Ns, Np, Nxin_file, Nyin_file, out_file)
        Ns_time = np.zeros((len(Ns), 2))
        Ns_time[:,0] = Ns
        Ns_time[:,1] = time
    
        ## [Save the result.]
        pickle.dump(Ns_time, open(save_file, 'wb'))

    if plot: 
        if run == False: 
            Ns_time = pickle.load(open(save_file, 'rb'))

        ## [Plot the resulting linear regression.]
        Plot_Performance(Ns, time, Nt)