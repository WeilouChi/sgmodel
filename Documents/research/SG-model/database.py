import pandas as pd
import numpy as np

def bondLengthDatabase():
    '''to do: 鍵長的database'''
    bondDF = {
        "bond" : ['CH', 'CC1', 'CC2', 'CC3', 'OH', 'CO1', 'CO2', 'CO1.5', 'CN1.5'],
        "distance" : [1.09, 1.52, 1.33, 1.21, 0.97, 1.42, 1.21, 1.32, 1.34]
    }
    return pd.DataFrame(bondDF)

    
    
def segmentVolumeDatabase():
    volumeDF = {
        'segment' : ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
        'volume' :[0.444, 0.444, 4.33, 10.69, 4.26, 7.32, 5.13, 3.48]
    }
    return pd.DataFrame(volumeDF)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    