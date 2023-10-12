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
getEffPlots = False
getWMTHistos = True
plotWMTHistos = True
lumi = 138000.0

indir = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/LeptonChecks/QCDBp_scenarios/"
outdir = os.getcwd()+'/plots_WMT/'
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
gInterpreter.Declare(
"""
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}
"""
)
gInterpreter.Declare(
"""
float MET_phi_threshold(float MET_Lep_DeltaPhi, float k){
return (k/1.5) * MET_Lep_DeltaPhi - k;
}
"""
)

histfile_name = "triangular_cut_histo.root"
if(getHistos):
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


# set up files and choices
histfile = TFile.Open(histfile_name, "UPDATE")

k_choices = np.arange(50, 255, 5)
pt_choices = np.arange(50,75,5) # 50, 55, 60, 65, 70                                                                                                                                             

nPermute = len(k_choices)
print(nPermute)

N_before = {"Bprime800":np.zeros(nPermute),
            "Bprime1400":np.zeros(nPermute),
            "Bprime2000":np.zeros(nPermute),
            "QCD":np.zeros(nPermute),
}
N_after = {"Bprime800":np.zeros(nPermute),
           "Bprime1400":np.zeros(nPermute),
           "Bprime2000":np.zeros(nPermute),
           "QCD":np.zeros(nPermute),
}


# Load events                                                                                                                                                                                           
chain_bkg = TChain("Events")
chain_bkg.Add("triangular_cut_QCD300.root")
chain_bkg.Add("triangular_cut_QCD500.root")
chain_bkg.Add("triangular_cut_QCD700.root")
chain_bkg.Add("triangular_cut_QCD1000.root")
chain_bkg.Add("triangular_cut_QCD1500.root")
chain_bkg.Add("triangular_cut_QCD2000.root")

Bprime800 = RDataFrame("Events", "triangular_cut_Bprime800.root")
Bprime1400 = RDataFrame("Events", "triangular_cut_Bprime1400.root")
Bprime2000 = RDataFrame("Events", "triangular_cut_Bprime2000.root")
QCD = RDataFrame(chain_bkg)
        
def getCounts(sample, Events, MET_cut, i):
    print("Processing {} for MET_pt>{}".format(sample, MET_cut))

    count_before = Events.Count()
    count_after = Events.Filter("MET_pt>{}".format(MET_cut)).Define("MET_phi_threshold", "MET_phi_threshold(MET_Lep_DeltaPhi, {})".format(k_choices[i])).Filter("MET_pt>MET_phi_threshold").Count()

    N_before[sample][i] = count_before.GetValue()
    N_after[sample][i]= count_after.GetValue()

def getWMT(sample, Events, MET_cut, k_choices):
    for k in k_choices:
        WMT_after = Events.Filter("MET_pt>{}".format(MET_cut)).Define("MET_phi_threshold", "MET_phi_threshold(MET_Lep_DeltaPhi, {})".format(k)).Filter("MET_pt>MET_phi_threshold").Histo1D(("", "", 100, 0, 2500), "LepIsoC_pt", "weights")
    
        histfile.cd()
        WMT_after.Write("lepton_pt_{}_{}_{}".format(sample, MET_cut, k))

def getWMTHistos(MET_pt):
    print("Getting lepton_pt histos for MET_pt>{}".format(MET_pt))
    getWMT("Bprime800", Bprime800, MET_pt, k_choices)
    getWMT("Bprime1400", Bprime1400, MET_pt, k_choices)
    getWMT("Bprime2000", Bprime2000, MET_pt, k_choices)
    getWMT("QCD", QCD, MET_pt, k_choices)

def plotWMTHistos(MET_pt, k_choices):
    for k in k_choices:
        print("Creating histogram plots for MET_pt>{}, k={}".format(MET_pt, k))
        Bprime800_hist = histfile.Get("lepton_pt_{}_{}_{}".format("Bprime800", MET_cut, k))
        Bprime1400_hist = histfile.Get("lepton_pt_{}_{}_{}".format("Bprime1400", MET_cut, k))
        Bprime2000_hist = histfile.Get("lepton_pt_{}_{}_{}".format("Bprime2000", MET_cut, k))
        QCD_hist = histfile.Get("lepton_pt_{}_{}_{}".format("QCD", MET_cut, k))

        c1 = TCanvas("c1", "c1", 600, 400)
        
        Legend.AddEntry(Bprime800_hist, "Bprime800", "l")
        Legend.AddEntry(Bprime1400_hist, "Bprime1400", "l")
        Legend.AddEntry(Bprime2000_hist, "Bprime2000", "l")
        Legend.AddEntry(QCD_hist, "QCD", "f")
        
        Bprime800_hist.Draw()
        Bprime1400_hist.Draw()
        Bprime2000_hist.Draw()

        rightmax = 1.1*histo_sum_bkg.GetMaximum()
        scale = gPad.GetUymax()/rightmax
        QCD_hist = QCD_hist(scale)
        
        QCD_hist.Draw("HIST")
        Bprime800_hist.Draw("HIST")
        Bprime1400_hist.Draw("HIST")
        Bprime2000_hist.Draw("HIST")

        QCD_hist.GetXaxis().SetTitle("lepton_pt")
        QCD_hist.GetYaxis().SetTitle("Events/bin")
        c1.Modified()
        c1.SaveAs(outdir+"lepton_pt_{}_{}.png".format(MET_pt, k))
        c1.Close()
        

def getPlots(MET_cut):
    for i in range(nPermute):
        getCounts("Bprime800",  Bprime800, MET_cut, i)
        getCounts("Bprime1400", Bprime1400, MET_cut, i)
        getCounts("Bprime2000", Bprime2000, MET_cut, i)
        getCounts("QCD", QCD, MET_cut, i)

    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_size_inches(12, 6)
    fig.suptitle('Optimization scan for triangular cut with MET_pt>{}'.format(MET_cut))

    ax1.plot(k_choices, N_after["Bprime800"]/np.sqrt(N_after["QCD"]), label="Bprime800")
    ax1.plot(k_choices, N_after["Bprime1400"]/np.sqrt(N_after["QCD"]), label = "Bprime1400")
    ax1.plot(k_choices, N_after["Bprime2000"]/np.sqrt(N_after["QCD"]), label = "Bprime2000")

    ax2.plot(k_choices, N_after["Bprime800"]/N_before["Bprime800"], label="Bprime800")
    ax2.plot(k_choices, N_after["Bprime1400"]/N_before["Bprime1400"], label = "Bprime1400")
    ax2.plot(k_choices, N_after["Bprime2000"]/N_before["Bprime2000"], label = "Bprime2000")
    
    ax1.set_xlabel('k [GeV]')
    ax2.set_xlabel('k [GeV]')
    ax1.set_ylabel('S/$\sqrt{B}$')
    ax2.set_ylabel('Signal Efficiency')
    
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    plt.show()
    fig.savefig('triangular_cut_pt_{}.png'.format(MET_cut))

if(getEffPlots):
    getPlots(pt_choices[0])
    getPlots(pt_choices[1])
    getPlots(pt_choices[2])
    getPlots(pt_choices[3])
    getPlots(pt_choices[4])

if(getWMTHistos):
    getWMTHistos(pt_choices[0])
    getWMTHistos(pt_choices[1])
    getWMTHistos(pt_choices[2])
    getWMTHistos(pt_choices[3])
    getWMTHistos(pt_choices[4])

if(plotWMTHistos):
    plotWMTHistos(pt_choices[0], k_choices)
    plotWMTHistos(pt_choices[1], k_choices)
    plotWMTHistos(pt_choices[2], k_choices)
    plotWMTHistos(pt_choices[3], k_choices)
    plotWMTHistos(pt_choices[4], k_choices)


