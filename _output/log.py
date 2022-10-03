# IMPORT PACKAGES -------------------------------------------------------------
import os
import logging
from datetime import datetime

# IMPORT MODULES --------------------------------------------------------------

# ------------------------------------------------------------- method - FOOTER
def footer(total_time,logg):
    
    logg.info('\n')
    logg.info(' '.ljust(10)+'   CPU TIME: {} seconds'.format(total_time))
    logg.info('\n')
    logg.info('~'*35+'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    logg.info('\n')

# ------------------------------------------------------------- method - HEADER
def header(logg):

    os.system('cls')
    now = datetime.now()
    now = now.strftime('%Y-%m-%d %H:%M:%S')

    logg.info('\n')
    logg.info(' '.ljust(35)+'ParamID - Parameter Identification Software')
    logg.info('~'*78)
    logg.info(' '.ljust(78-len(now))+now)
    logg.info('\n')

# ------------------------------------------------------------- method - LOGGER
def logger():

    path = 'Output/ParamID.log'

    try:
        os.remove(path)
    except:
        pass

    logg = logging.getLogger()
    logg.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO) # or any other level
    logg.addHandler(ch)
    fh = logging.FileHandler(path)
    fh.setLevel(logging.DEBUG)
    logg.addHandler(fh)

    return logg