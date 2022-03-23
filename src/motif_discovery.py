import numpy as np

def shift_per_row_py( shift_vector, region_df, region_length ):
    """
    """
    region_length = int( region_length )
    shift_vector = np.array( shift_vector )
    region_df = np.array( region_df )

    if shift_vector.shape[0] != region_df.shape[0]:
        exit(2)
    
    all_list = []
    counter = 0

    for peak in region_df:

        new_peak = peak[ int(shift_vector[counter])-1: ( int(shift_vector[counter])+region_length-1 ) ]
        # a = int( region_length - new_peak.shape[0] )
        # if a < 0:
        #     a = 0 
        # lst = ["0"] * a
        # new_list = list( list(new_peak) + list(lst) )
        # new_list = [item for sublist in new_list for item in sublist]
        # print(new_list)
        all_list.append( new_peak )
        counter += 1
    
    return(all_list)
