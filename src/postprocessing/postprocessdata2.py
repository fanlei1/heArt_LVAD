from postprocessdatalib2 import *

def postprocessdata(IODet, SimDet):

    directory = IODet["outputfolder"] + '/' 
    casename =  IODet["caseID"] 
    BCL = SimDet["HeartBeatLength"]
    cycle = SimDet["closedloopparam"]["stop_iter"]


    for ncycle in range(cycle-1,cycle):
    
         filename = directory+casename+"/"+"BiV_PV.txt"
         homo_tptt, homo_LVP, homo_LVV, homo_Qmv = extract_PV(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_Q.txt"
         homo_tptt, homo_Qao, homo_Qmv, homo_Qper, homo_Qla, homo_Qlad, homo_Qlcx, Q_lvad = extract_Q(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_P.txt"
         homo_tptt, homo_Pven, homo_LVPP, homo_Part, homo_PLA = extract_P(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_IMP_InC.txt"
         homo_tpt_IMP, homo_IMP = extract_probe(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_fiberStrain.txt"
         homo_tpt_Eff, homo_Eff = extract_probe(filename, BCL, ncycle)
    
         filename = directory+casename+"/"+"BiV_fiberStress.txt"
         homo_tpt_Sff, homo_Sff = extract_probe(filename, BCL, ncycle)
    
         ESP, ESV, ESind = extractESP(homo_LVP, homo_LVV)
    
         EDP, EDV, EDind = extractEDP(homo_LVP, homo_LVV)
    
         SBP = max(homo_Part)*0.0075
         DBP = min(homo_Part)*0.0075
    
         print filename
         print "EF = ", (max(homo_LVV) - min(homo_LVV))/max(homo_LVV), " EDV = ", max(homo_LVV), " ESV = ", min(homo_LVV),\
         		" EDP = ", EDP, " SBP = ", SBP, " DBP = ", DBP, "ESt = ", ESind, "EDt = ", EDind
    
         print "Peak LV pressure = ", max(homo_LVP)
    
         homo_directory = directory+casename+"/"
    
         tpt_array = readtpt(homo_directory + "tpt.txt")
         ind = np.where((tpt_array>(ncycle)*BCL)*(tpt_array<(ncycle+1)*BCL))
         tpt = tpt_array[ind]
         
         ## Get Point cloud for probing
         ptcloud, radialpos, vtkradialpos = getpointclouds(homo_directory, clipoffset=5e-1, npts=10000)
         vtk_py.writeXMLPData(vtkradialpos, casename+".vtp")
    
         # Get transmural variation of IMP
         index = find_nearest(tpt, homo_tptt[ np.argmax(homo_LVPP)]) # Find ID correspond to peak LV pressure
         imp = probeqty(homo_directory,  'ME/imp_constraint', ptcloud, ind, index)
         imp = imp*0.0075
         
    
         ## Get transmural variation of WD
         Sff = probetimeseries(homo_directory, "ME/fstress", ptcloud, ind, "DG", 0)
         Eff = probetimeseries(homo_directory, "ME/Eff", ptcloud, ind, "DG", 0)
         WD = np.array([-1.0*np.trapz(Sff[:,i]*0.0075, Eff[:,i]) for i in range(0,len(Sff[1,:]))])
    
         # Convert to vtp flie
         for i in range(0,len(Sff[:,1])):
         	pdata = vtk.vtkPolyData()
         	pdata.DeepCopy(vtkradialpos)
         	Sff_VTK_data = numpy_support.numpy_to_vtk(num_array=0.0075*Sff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	Sff_VTK_data.SetName("fstress_")
         	pdata.GetPointData().AddArray(Sff_VTK_data)
         	Eff_VTK_data = numpy_support.numpy_to_vtk(num_array=Eff[i,:].ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	Eff_VTK_data.SetName("Eff_")
         	pdata.GetPointData().AddArray(Eff_VTK_data)
         	WD_VTK_data = numpy_support.numpy_to_vtk(num_array=WD.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
         	WD_VTK_data.SetName("WD_")
         	pdata.GetPointData().AddArray(WD_VTK_data)
         	#vtk_py.writeXMLPData(pdata, casename+"fstress"+str(i)+".vtp")
    
    
         ## Get Ecc
         Ecc = probetimeseries(homo_directory, "ME/Ecc", ptcloud, ind, "DG", 0)
         peakEcc = np.max(np.abs(np.mean(Ecc, axis=1)*100))
         print "Peak Ecc = ", peakEcc
         
         # Get Ell
         Ell = probetimeseries(homo_directory, "ME/Ell", ptcloud, ind, "DG", 0)
         peakEll = np.max(np.abs(np.mean(Ell, axis=1)*100))
         print "Peak Ell = ", peakEll
    
    
         np.savez(directory+casename+"/"+casename+".npz", \
         	 homo_tptt    = homo_tptt,\
         	 homo_LVP     = homo_LVP,\
         	 homo_LVV     = homo_LVV,\
         	 homo_Qmv     = homo_Qmv,\
         	 homo_Qao     = homo_Qao,\
         	 homo_Qper    = homo_Qper,\
         	 homo_Qla     = homo_Qla,\
         	 homo_Qlad    = homo_Qlad,\
         	 homo_Pven    = 0.0075*homo_Pven,\
         	 homo_LVPP    = 0.0075*homo_LVPP,\
         	 homo_Part    = 0.0075*homo_Part,\
         	 homo_PLA     = 0.0075*homo_PLA,\
                  homo_tpt_IMP = 0.0075*homo_tpt_IMP,\
         	 homo_IMP     = homo_IMP,\
         	 homo_tpt_Eff = homo_tpt_Eff,\
         	 homo_Eff     = homo_Eff,\
                  homo_tpt_Sff = homo_tpt_Sff,\
         	 homo_Sff     = homo_Sff,\
         	 ESP = ESP, ESV = ESV, EDP = EDP, EDV = EDV, SBP = SBP, DBP = DBP, ESind = ESind, EDind = EDind\
         	 #Qtotal       = Qtotal,\
         	 imp          =imp,\
         	 radialpos    =radialpos,\
         	 Eff	      =Eff, \
         	 Sff	      =Sff, \
                  WD           =WD,\
         	 Ecc          =Ecc,\
         	 Ell          =Ell,\
         	 BCL	      =BCL,\
                  tpt          =tpt,\
         	 ncycle       =ncycle
         	)
    
	
	
		
