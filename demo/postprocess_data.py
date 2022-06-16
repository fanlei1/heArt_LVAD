import sys, pdb
from dolfin import *

sys.path.append("/mnt/home/fanlei1")
from heArtLVAD.src.sim_protocols.run_BiV_ClosedLoop import run_BiV_ClosedLoop as run_BiV_ClosedLoop
from heArtLVAD.src.postprocessing.postprocessdataBiV2 import postprocessdata as postprocessdata

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {#"casename" : "biv_HF",
             "outputfolder" : './outputs_LVAD_homo/',
             "caseID" : 'rpm12k'}

Circparam = {     "stop_iter" : 16,
                  };

SimDetails = {
                  "HeartBeatLength": 800.0,
                  "closedloopparam": Circparam,
                 }



# Postprocessing
postprocessdata(IODet=IODetails, SimDet=SimDetails)

#  - - - - - - - - - - -
