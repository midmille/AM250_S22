"""
AM250_S22/final/game_of_life/plot/plot_life.py
Author : Miles D. Miller, Univeristy of California Santa Cruz
Created : 11:53 PM, June 3, 2022
About: This python file is for the visualization of the results of the game of life. The file reads in the 
       game of life output and then plots the output. 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def Read_Mat(Nx, Ny, filename):
    """
    This function reads in the matrix from a given data file.

    Parameters
    ----------
    Nx: Integer
        The number of columns.
    Ny: Integer
        The number of rows.
    filename: String
        The name of the file to read in. 

    Returns
    -------
    A: 2-D Array, [Ny, Nx]
        The data array for the time step.
    """

    ## [Init the output array.]
    A = np.zeros((Ny, Nx))

    ## [Open the file for reeading.]
    with open(filename, 'r') as f:
        ## [Loop the rows.]
        for k, line in enumerate(f): 
            ## [Loop over the columns by splitting with white spaces.]
            for j,val in enumerate(line.split()):
                ## [Save the val to array.]
                A[k,j] = val
    return A


def Read_Output(Nx, Ny, Nt, Nw, savefile_head, outdir): 
    """
    This function reads all the output arrays is the designated output directory and 
    then store the result into a 3-D Array where the 3rd dimension is for the 
    time coordiante.

    Parameters
    ----------

    Returns
    -------
    """

    ## [Make the array of time steps that were actually written, 
    ##  the Nt+Nw is because arange doesn't include the endpoint.]
    tstps = np.arange(0, Nt+Nw, Nw)
    ## [Fortran is 1-indexed based.]
    tstps[0] = 1

    ## [Init the 3-D array.]
    life_arr = np.zeros((Ny, Nx, len(tstps)))

    ## [Loop over the tstps.]
    for k, tstp in enumerate(tstps): 
        filename = f"{outdir}{savefile_head}{100000+tstp}.dat"
        ## [Read in the array for the above file.]
        life_arr[:,:,k] = Read_Mat(Nx, Ny, filename)
    
    return life_arr, tstps


def Plot_Life(Nx, Ny, life_array, tindx):
    """
    This file plots a single frame of the life output evolution.

    Parameters
    ----------

    Returns
    -------
    """

    fig, ax = plt.subplots()
    
    ax.imshow(life_array[:,:,tindx], extent=[0, Nx, 0, Ny], cmap='Blues')
    ax.set_xticks(np.arange(0, Nx), [])
    ax.set_yticks(np.arange(0, Ny), [])
    ax.grid(True, c='black')
    
    fig.show()

    return 


def Animate_Life(Nx, Ny, life_array, tstps):
    """
    This files animates every time step of life.
    """
    frames = len(tstps)

    def animate(i): 
        ax.clear() 
        ax.imshow(life_array[:,:,i], extent=[0, Nx, 0, Ny], cmap='Blues')
        ax.set_title(f't = {tstps[i]}')
        ax.set_xticks(np.arange(0, Nx), [])
        ax.set_yticks(np.arange(0, Ny), [])
        ax.grid(True, c='black')
        return 

    fig, ax = plt.subplots()
    
    ani = FuncAnimation(fig, animate, frames=frames, interval=250, repeat=True)
    
    plt.show()

    return


if __name__ == '__main__': 
    
    ## [Parameters from main.f90]
    Nx = 20
    Ny = 20
    Nt = 80
    Nw = 1
    ## [The parallelkization flag.]
    pflag = 'serial'
    savefile_head = "life_out"
    ## [The out directory for the given parralleization.]
    if pflag == 'serial': 
        outdir = '../output/serial_out/'
    if pflag == 'cols':
        outdir = '../output/cols_out/'

    ## [Read the output.]
    life_arr, tstps = Read_Output(Nx, Ny, Nt, Nw, savefile_head, outdir)

    ## [Animate fig.]
    Animate_Life(Nx, Ny, life_arr, tstps)



