import sys, pdb
from dolfin import * 

sys.path.append("/mnt/home/fanlei1")
from heArtLVAD.src.sim_protocols.run_BiV_ClosedLoop_Ees import run_BiV_ClosedLoop as run_BiV_ClosedLoop
from heArtLVAD.src.postprocessing.postprocessdataBiV2 import postprocessdata as postprocessdata
from FittingLVAD.LVADFn import LVAD as LVAD

HeartMate = LVAD("./FittingLVAD/HeartMate.npz")
print HeartMate

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
IODetails = {"casename" : "biv_HF",
             "directory_me" : '../BiVMesh/',
             "directory_ep" : '../BiVMesh/', 
             "outputfolder" : './outputs_LVAD_homo/',
             "folderName" : '',
             "caseID" : 'rpm10kpre1'}

contRactility = 100e3*0.85

GuccioneParams = {"ParamsSpecified" : True, 
     		  "Passive model": {"Name": "Guccione"},
     		  "Passive params": {"Cparam": Constant(100.0*25), 
		                     "bff"  : Constant(29.0*1.1),
		                     "bfx"  : Constant(13.3*1.1),
                  		     "bxx"  : Constant(26.6*1.1),
				   },
	         "Active model": {"Name": "Time-varying"},
     	         "Active params": {"tau" : 25, "t_trans" : 300, "B" : 4.75,  "t0" : 275,  "l0" : 1.58, \
				   "Tmax" : Constant(contRactility), "Ca0" : 4.35, "Ca0max" : 4.35, "lr" : 1.85},
                 "HomogenousActivation": True,
		 "deg" : 4, 
                 "Kappa": 1e5,
                 "incompressible" : True,
                 }

		  # Systemic
Circparam = {     "Ees_la": 350,
        	  "A_la": 58.67,
        	  "B_la": 0.049,
        	  "V0_la": 10,
        	  "Tmax_la": 120,
        	  "tau_la": 75,
		  "tdelay_la": 160,
		  "Csa": 0.0052,
		  "Cad": 0.0330,
    		  "Csv" : 0.3,
    		  "Vsa0" : 700*0.8, 
    		  "Vsv0" : 2500.0*0.8*0.9, #2500.0,-preload
		  "Vad0" : 40,
    		  "Rav" : 500.0,
    		  "Rsv" : 100.0,
    		  "Rsa" : 18000,
		  "Rad" : 106000*2.0,
    		  "Rmv" : 1500.0,
		  # Pulmonary 
		  "Ees_ra": 81.33,
		  "A_ra": 466.6,
		  "B_ra": 0.033,
		  "V0_ra": 20,
		  "Tmax_ra": 120,
		  "tau_ra": 25,
		  "tdelay_ra": 100,
		  "Cpa": 0.0125/5, 
    		  "Cpv" : 0.9,
    		  "Vpa0" : 360,
    		  "Vpv0" : 15*10,
    		  "Rpv" : 500.0,
    		  "Rtv" : 400.0,
    		  "Rpa" : 10000.08*2.0,
    		  "Rpvv" : 400, #400, 
		  # flow rate
		  "Q_lvad" : 0.0,
		  "Q_sv" : 0.0659387173129,
		  "Q_av" : 0.0,
		  "Q_sa" : 0.0113526347222,
		  "Q_ad" : 0.0797504135138,
		  "Q_mv" : 0.0,
		  "Q_tv" : 0.0,
		  "Q_pa" : 0.001213376269,
		  "Q_pv" : 0.0646105314407,
		  "Q_pvv" : 0.0,
		  # volumes

                  "V_sv" : 2461.74158826,
                  "V_LV" : 203.128641682,
                  "V_sa" : 621.337347089,
                  "V_ad" : 422.501367564,
                  "V_LA" : 23.9463825761,
                  "V_pv" : 3363.41063636,
                  "V_RV" : 170.010754459,
                  "V_pa" : 369.292168663,
                  "V_RA" : 61.4069030347,

                  "stop_iter" : 15, 
		  # LVAD
		  'Q_lvad_characteristic': HeartMate,
		  'Q_lvad_rpm' : 10000,
		  'Q_lvad_scale' : 1.0
		  };


SimDetails = {    
                  "HeartBeatLength": 800.0,
                  "dt": 1.0,
                  "writeStep": 10.0,
                  "GiccioneParams" : GuccioneParams, 
                  "nLoadSteps": 4,
                  "DTI_EP": False,
                  "DTI_ME": False,
                  "d_iso": 1.5*0.01, 
                  "d_ani_factor": 4.0, 
                  "ploc": [[-0.083, 5.6,-1.16, 2.0, 1.0]],
                  "pacing_timing": [[4.0, 20.0]],
		  "closedloopparam": Circparam,
                  "Ischemia": False,
             	  "isLV" : False,
                  "topid" : 4,
                  "LVendoid" : 2,
                  "RVendoid" : 3,
                  "epiid" : 1,
		  "abs_tol": 1e-9,
		  "rel_tol": 1e-9,
                 }


# Run Simulation
run_BiV_ClosedLoop(IODet=IODetails, SimDet=SimDetails)
# Postprocessing
#postprocessdata(IODet=IODetails, SimDet=SimDetails)

#  - - - - - - - - - - - -- - - - - - - - - - - - - - - -- - - - - - - 
    
