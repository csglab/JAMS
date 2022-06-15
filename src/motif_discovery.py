import numpy as np
from numba import njit

@njit
def shift_per_row_py_numpy( shift_vector, region_df, \
                            region_length, new_region_df ):
    """
    """
    shift_vector = shift_vector - 1
    counter = 0

    for line in region_df:
    
        new_line = line[ shift_vector[counter] : \
                         ( shift_vector[counter]+region_length ) ]
        
        new_region_df[counter,] = new_line
        counter += 1

    return new_region_df


def shift_per_row_py( shift_vector, region_df, region_length ):

    try:
        region_df = np.array( region_df, dtype = float )
        shift_vector = np.array( shift_vector, dtype = int )
        region_length = int(region_length)
    except:
        print("Error: shift_per_row_py parameter are not the correct type")
        exit(2)
        
    new_region_df = np.zeros( ( region_df.shape[0], region_length ), \
                              dtype = float )
    
    
    if shift_vector.shape[0] != region_df.shape[0]:
        print("positions vector and dataframe have different n rows")
        exit(2)

    region_df = np.array( region_df, dtype = float )
    shift_vector = np.array( shift_vector, dtype = int )
    region_length = int(region_length)
    new_region_df = np.zeros( ( region_df.shape[0], region_length ), \
                              dtype = float )
    
    new_matrix = shift_per_row_py_numpy( shift_vector, \
                                         region_df, \
                                         region_length, \
                                         new_region_df )
    return( new_matrix )
