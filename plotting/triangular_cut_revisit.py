from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import patches
import os, time
from math import sqrt
from samples import *
from utils import *
import itertools

getHistos = False
lumi = 138000.0

indir = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/LeptonChecks/QCDBp_scenarios/"
outdir = os.getcwd()+'/plots_triangularCut/'
if not os.path.exists(outdir): os.system('mkdir -p ' + outdir)

samples = {"Bprime800":"Bprime800_scenarios.root",
           "Bprime1400":"Bprime1400_scenarios.root",
           "Bprime2000":"Bprime2000_scenarios.root",                    
           "QCD200":"QCD200_scenarios.root",   
           "QCD300":"QCD300_scenarios.root",
           "QCD500":"QCD500_scenarios.root",
           "QCD700":"QCD700_scenarios.root",
           "QCD1000":"QCD1000_scenarios.root",
           "QCD1500":"QCD1500_scenarios.root",
           "QCD2000":"QCD2000_scenarios.root"
}

nRun = {'Bprime800':99800.,
        'Bprime1400':99600.,
        'Bprime2000':99600.,  # fixme, need a sum-of-gen-weights calculated...
        'QCD200':61542214.,
        'QCD300':56214199.,
        'QCD500':61097673.,
        'QCD700':47314826.,
        'QCD1000':15230975.,
        'QCD1500':11887406.,
        'QCD2000':5710430.
}
xsec = {'Bprime800':1.0,'Bprime1400':1.0,'Bprime2000':1.0, #signals all the same at 1pb for now, predictions vary 
        'QCD200':1712000,
        'QCD300':347700,
        'QCD500':32100,
        'QCD700':6831,
        'QCD1000':1207,
        'QCD1500':119.9,
        'QCD2000':25.24
}

dphi ='''
using namespace ROOT::VecOps;
float dphi(const float& phi1, const float& phi2) {
    return DeltaPhi(phi1, phi2);
}
'''
gInterpreter.Declare(dphi)
gInterpreter.Declare("""
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}
""")

if(getHistos):
    histfile_name = "triangular_cut_histo.root"
    histfile = TFile.Open(histfile_name, "RECREATE")
    for sample in samples:
        filename1 = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/LeptonChecks/QCDBp_short/"+sample+"_MET.root"
        tfile1 = TFile.Open(filename1)

        filename2 = indir+samples[sample]
        tfile2 = TFile.Open(filename2)
        ftree2 = tfile2.Get("Events")
        ftree2.AddFriend("Events", tfile1)

        LepIsoC = RDataFrame(ftree2)
        ElIsoC = LepIsoC.Filter("isEl").Define("MET_Lep_DeltaPhi","dphi(LepIsoC_phi, MET_phi)").Define("weights","weights(Generator_weight,{},{},{})".format(lumi,xsec[sample],nRun[sample]))
    
        ElIsoC.Snapshot("Events", "triangular_cut_{}.root".format(sample))
        histfile.cd()
        histo = ElIsoC.Histo2D(("PhiPt_histo","PhiPt_histo", 25,0,3.15,25,50,400),"MET_Lep_DeltaPhi","MET_pt", "weights")
        histo.Write(sample)
    histfile.Close()
        
histfile = TFile.Open(histfile_name, "READ")

chain_bkg = TChain("Events")
for sample in samples:
    if "QCD" in sample:
        chain_bkg.Add("triangular_cut_{}.root".format(sample))
        
Events_bkg = RDataFrame("Events", chain_bkg)
Events_Bp800 = RDataFrame("Events", "triangular_cut_Bprime800.root")

kList = np.arange(100, 205, 5).tolist()
ptList = np.arange(50,75,5).tolist()
permutations = itertools.product(ptList, kList)

#for k, pt in permutations:
    #MET_ptThreshold = (k/1.5) * El_DelPhi - k
