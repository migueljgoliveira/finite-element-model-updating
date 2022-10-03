"""
  PROGRAM: 		ParamID

  PURPOSE: 		Material Parameter IDentification

  DEVELOPER: 	M9'Calibre
"""

# IMPORT PACKAGES
from time import process_time

# IMPORT MODULES
import _input
import _output
import _algorithms
import _simulation

# ---------------------------------------------------------------------- main()
def main():

    # create logger
    logg = _output.logger()

    start_time = process_time()

    # print header to command line and log file
    _output.header(logg)

    # ------------------------------------------------------------ INPUT MODULE
    options,tests = _input.main(logg)
    
    # ------------------------------------------------------- SIMULATION MODULE
    if options.info.calibration == False:
        out = _simulation.main(options,tests)

        # OUTPUT
        run = -1
        _output.main(out,options,tests,run)

    # ----------------------------------------------------- OPTIMISATION MODULE
    if options.info.calibration == True:
        _algorithms.main(options,tests,logg)

    total_time = process_time() - start_time

    # print footer to command line and log file
    _output.footer(total_time,logg)

    return

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()