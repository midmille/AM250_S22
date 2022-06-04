"""
This file is for the functions that plot the results of the FORTRAN routines
for the final project of AM213A Problem a. 
"""

import numpy as np
import matplotlib.pyplot as plt

def Load_Image(rows, cols, filename, read_dims=True): 
    """
    This function loads the image into an array.
    
    Parameters
    ----------
    rows: Integer
        The number of rows(lines) in the file. Pass as None if read_dims == True. 
    cols: Integer
        The number of columns(entries per line) in file. Pass as None if read_dims == True. 

    filename: String
        The name of the file. 
    read_dims: ooptional, Boolean, default=True
        This flag is whether the data file being read had the matrix dimensions
        in the first line.

    Returns 
    -------
    im_dat: Float, Array, [rows, cols]
        The numpy array of the data. 
    """
    ## [The dimensions are wrritten on the first line of the file.]
    if read_dims == True:

        with open(filename, 'r') as f: 
            ## [Get the dimensions from the first line.]
            nrows, ncols = f.readline().split()
            ## [String to integer and subtracting one for fortran->python]
            nrows = int(nrows)
            ncols = int(ncols)

        ## [Make empty array with the newly found dimensions from file.]
        im_dat = np.zeros((nrows,ncols))

        ## [Read in from the given data file.]
        with open(filename, 'r') as f: 
            ## [loop over the rows.]
            for k, line in enumerate(f): 
                ## [This moves k back one so that the lines of the file matches with the 
                ## enumeration given by k, ie k=0, line=1.]
                k=k-1
                ## [skip the first row.]
                if k==-1:
                    continue 
                ## [Loop over the Split via white spaces.]
                for j,val in enumerate(line.split()):
                    ## [Save to array.]
                    im_dat[k,j] = val

    elif read_dims == False:

        ## [The empty data array.]
        im_dat = np.zeros((rows,cols))

        ## [Read in from the given data file.]
        with open(filename, 'r') as f: 
            ## [loop over the rows.]
            for k, line in enumerate(f): 
                ## [Loop over the Split via white spaces.]
                for j,val in enumerate(line.split()):
                    ## [Save to array.]
                    im_dat[k,j] = val

    ## [Return the data array.]
    return im_dat


def Show_Image(im_dat): 
    """
    This function plots the given image file using matplotlib

    Parameters
    ----------
    im_dat: Float, Array
        The numpy array of the image data to be plotted. 

    Returns
    -------
    matplotlib figure

    """

    fig, ax = plt.subplots()

    ax.imshow(im_dat , cmap='gray')

    fig.show()

    return 

def Show_All_Compressed_Images(k_list, filehead, origim_file, rows, cols): 
    """
    This function loads the data for all of the files resulting from the compressions of
    the given kth singular values in the given k_list. It also plots the original
    image. 

    Parameters
    ----------
    k_list: list
        The list of k values. The k-values corresponds to the kth singular value compression stored to file. 
    file_head: String 
        The name of the file without the k-value in it. 
    origim_file: String
        The name of the file with the original image init. 
    rows: Integer
        The number of rows of the original matrix A. 
    cols: Integer
        The number of cols of the original matrix A. 

    Returns
    -------
    None, plots the all images as subplots in figure.
    """

    fig, axes = plt.subplots(ncols=3, nrows=3)
    axes = axes.flatten()
    ## [Load and plot the original image.]
    im_dat = Load_Image(rows, cols, origim_file, read_dims=False)
    axes[0].imshow(im_dat, cmap='gray')
    axes[0].set_title('Original Image')
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    ## [Loop over the kvalue list.]
    for i, kval in enumerate(k_list): 
        ## [Formatting the kval into the file name]
        filename = f'{filehead}{100000+kval}.dat'

        ## [Load and plot the kval imag.]
        im_dat = Load_Image(None, None, filename, read_dims=True)
        axes[i+1].imshow(im_dat, cmap='gray')
        axes[i+1].set_title(f'k={kval} Image')
        axes[i+1].set_xticks([])
        axes[i+1].set_yticks([])

    fig.tight_layout()
    fig.show()


    return 


def Calc_Compression_Error(k_list, filehead, origim_file, rows, cols):
    """
    This file calculates the differences of each kval compressed image from 
    the original image using the frobenius norm. 

    Parameters
    ----------
    k_list: list
        The list of k values. The k-values corresponds to the kth singular value compression stored to file. 
    file_head: String 
        The name of the file without the k-value in it. 
    origim_file: String
        The name of the file with the original image init. 
    rows: Integer
        The number of rows of the original matrix A. 
    cols: Integer
        The number of cols of the original matrix A. 

    Returns
    -------
    err_arr: Float, 1-D Array
        This is the array of error corresponding to each kval. 
    """

    ## [Load original image.]
    origim_dat = Load_Image(rows, cols, origim_file, read_dims=False)

    ## [Init the err_arr.]
    err_arr = np.zeros(len(k_list))

    ## [Loop over the kvalue list.]
    for i, kval in enumerate(k_list): 
        ## [Formatting the kval into the file name]
        filename = f'{filehead}{100000+kval}.dat'

        ## [Load and plot the kval imag.]
        im_dat = Load_Image(None, None, filename, read_dims=True)
        
        ## [Calculate the frobenius norm of the difference of the matrices.]
        ## ['fro' stands for the frobenius norm.]
        f_norm = np.linalg.norm(origim_dat - im_dat, 'fro')
        ## [Divide the f_norm by the product of rows and columns.]
        err_arr[i] = f_norm / (rows*cols)

    return err_arr


def Plot_Compression_Error(k_list, err_arr): 
    """
    Plot the resulting error array calculated by Calc_Compression_Error function. 

    Parameters
    ----------
    k_list: Integer, list, [N]
        The list of k values.
    err_arr: Float, 1-D Array, [N]
        The err_arr output for the SVD compressions. 

    Returns
    -------
    None, just plots the error as a function of k value. 

    """

    fig, ax = plt.subplots()

    ax.plot(k_list, err_arr)
    ax.set_title('Error as Function of Number of \n Singular Values Used in Image Compression')
    ax.set_ylabel(r'$E_k = \frac{|| \mathbf{A} - \mathbf{A}_{\sigma k} ||_F}{mn}$')
    ax.set_xlabel('Number of Singular Values k')
    ax.grid()

    fig.show()

if __name__ == '__main__': 
    
    rows = 3355
    cols = 5295

    ## [Low res file name.]
    origim_file= 'dog_bw_data.dat'
    filehead='Images/Image_appn_'

    ## [The list of kvals.]
    k_list = [20,40,80,160,320,640,1280,2560]

    #Show_All_Compressed_Images(k_list, filehead, origim_file, rows, cols)

    err_arr = Calc_Compression_Error(k_list, filehead, origim_file, rows, cols)

    Plot_Compression_Error(k_list,  err_arr)
