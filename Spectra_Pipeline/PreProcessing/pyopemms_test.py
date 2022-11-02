import numpy as np
import pandas as pd
from pyopenms import *


#m = RawQuant.RawFileReader.open_raw_file("PNACIC_UnkLip_BTLE_P2_QE_A_POS_13Oct20_Lola-WCSH315112.raw")

inp = MSExperiment()
MzMLFile().load("PNACIC_UnkLip_BTLE_P2_QE_A_POS_13Oct20_Lola-WCSH315112.raw", inp)
scan_nrs = [0, 2, 5, 7]


e = MSExperiment()
for k, s in enumerate(inp):
  if k in scan_nrs and s.getMSLevel() == 1:
    e.addSpectrum(s)


MzMLFile().store("test_filtered.mzML", e)
