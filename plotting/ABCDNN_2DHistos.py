import os, sys, time
import itertools as it
import numpy as np
from math import sqrt
from ROOT import *
from samples import *
from utils import *

# Options
case = "_BdecayCase1"
getHistos = True
getPlots = True
lumi = 138000.0

# Setup input dir
step1dir = 'root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Sep2023_2018/'

# Setup output dir
outdir = os.getcwd()+'/plots_bkgsig2D/'
subdir = case[1:]+"/"
if not os.path.exists(outdir + subdir): os.system('mkdir -p ' + outdir + subdir)

# Define branches
branches = {
#            "lepton_pt":["lepton_pt", 50, 0, 1000, "[GeV]"],
#            "lepton_eta":["lepton_eta", 40, -4, 4, ""],
#            "lepton_miniIso":["lepton_miniIso", 50, 0, 0.2, ""],
            "NJets_central":["NJets_central", 10, 0, 10, ""],
            "NJets_DeepFlavL":["NJets_DeepFlavL", 10, 0, 10, ""],
            "NJets_forward":["NJets_forward", 10, 0, 10, ""],
            "NFatJets":["NFatJets", 10, 0, 10, ""],
            "NOS_gcJets_central":["NOS_gcJets_central", 10, 0, 10, ""],
            "NSS_gcJets_central":["NSS_gcJets_central", 10, 0, 10, ""],
            "NOS_gcJets_DeepFlavL":["NOS_gcJets_DeepFlavL", 10, 0, 10, ""],
            "NSS_gcJets_DeepFlavL":["NSS_gcJets_DeepFlavL", 10, 0, 10, ""],
            "NOS_gcFatJets":["NOS_gcFatJets", 10, 0, 10, ""], 
            "NSS_gcFatJets":["NSS_gcFatJets", 10, 0, 10, ""],
            #"gcJet_HT":["gcJet_HT", 30, 0, 5000, "[GeV]"], # transform var
            "gcJet_ST":["gcJet_ST", 30, 0, 5000, "[GeV]"],
            "gcFatJet_nJ":["gcFatJet_nJ", 10, 0, 10, ""],
            "gcFatJet_nT":["gcFatJet_nT", 10, 0, 10, ""],
            "gcFatJet_nW":["gcFatJet_nW", 10, 0, 10, ""],
            "leadingOSFatJet_TorW":["leadingOSFatJet_TorW", 20, 0, 2, ""], # max(gcOSFatJet_pNetT[0], gcOSFatJet_pNetW[0])
            "leadingOSFatJet_tau":["leadingOSFatJet_tau", 20, 0, 1, ""], # min(gcOSFatJet_tau21[0],gcOSFatJet_tau32[0])
            "leadingOSFatJet_TWorQCD":["leadingOSFatJet_TWorQCD", 20, 0, 1, ""], #max(gcOSFatJet_pNetTvsQCD[0], gcOSFatJet_pNetWvsQCD[0])
            #"FatJet_pt_1":["FatJet_pt_1", 50, 0, 1500, "[GeV]"],
            #"FatJet_pt_2":["FatJet_pt_2", 50, 0, 1500, "[GeV]"],
            #"FatJet_sdMass_1":["FatJet_sdMass_1", 50, 0, 500, "[Gev]"],
            #"FatJet_sdMass_2":["FatJet_sdMass_2", 50, 0, 500, "[Gev]"],
#            "tau21_1":["tau21_1", 30, 0, 1, ""], # not as useful as tau21 vs tau 23 for OSFatjets
#            "tau21_2":["tau21_2", 30, 0, 1, ""],
#            "minDR_lep_FatJet":["minDR_lep_FatJet", 50, 0, 5, "[GeV]"],
#            "ptRel_lep_FatJet":["ptRel_lep_FatJet", 50, 0, 500, "[GeV]"],
            "minDR_leadAK8otherAK8":["minDR_leadAK8otherAK8", 20, 0, 10, "[GeV]"],
            "minDR_leadAK4otherAK4":["minDR_leadAK8otherAK8", 20, 0, 10, "[GeV]"],
######            "minDR_AK8s_discrete":["minDR_AK8s_discrete", 14, 0, 14, ""],
######            "minDR_AK4s_discrete":["minDR_AK4s_discrete", 18, 0, 18, ""],
#            "minDR_lep_Jet":["minDR_lep_Jet", 50, 0, 5, "[GeV]"],
#            "ptRel_lep_Jet":["ptRel_lep_Jet", 50, 0, 500, "[GeV]"],
#            "W_pt":["W_pt", 30, 0, 1000, "[GeV]"],
            "W_eta":["W_eta", 40, -4, 4, ""],
            "W_MT":["W_MT", 30, 0, 1500, "[GeV]"],
            "DR_W_lep":["DR_W_lep", 30, 0, 5, ""],
            "minM_lep_Jet":["minM_lep_Jet", 30, 0, 1500, "[GeV]"],
            "minM_lep_Jet_jetID":["minM_lep_Jet_jetID", 15, 0, 15, ""],
            "minM_lep_Jet_TorW":[ "minM_lep_Jet_TorW", 2, 0, 2, ""],
#            "t_pt":["t_pt", 30, 0, 1000, "[GeV]"],
#####            "t_eta":["t_eta", 30, -4, 4, ""],
            "DR_W_b":["DR_W_b", 30, 0, 7, ""],
            "Bprime_chi2":["Bprime_chi2", 30, 0, 1000, ""],
#            "Bdecay_obs":["Bdecay_obs", 4, 0, 4, ""]
}

tags_cases = {"_BdecayCase1":"Bdecay_obs==1",
              "_BdecayCase2":"Bdecay_obs==2",
              "_BdecayCase3":"Bdecay_obs==3",
              "_BdecayCase4":"Bdecay_obs==4",
          }

# 2D phase space choices
var_list = ["NJets_central", "NJets_DeepFlavL", "NJets_forward", "NFatJets", "NOS_gcJets_central", "NSS_gcJets_central", "NOS_gcJets_DeepFlavL", "NSS_gcJets_DeepFlavL", "NOS_gcFatJets", "NSS_gcFatJets", "gcJet_ST", "gcFatJet_nJ", "gcFatJet_nT", "gcFatJet_nW", "leadingOSFatJet_TorW", "leadingOSFatJet_tau", "leadingOSFatJet_TWorQCD", "minDR_leadAK8otherAK8", "minDR_leadAK4otherAK4", "W_eta", "W_MT", "DR_W_lep", "minM_lep_Jet", "minM_lep_Jet_jetID", "minM_lep_Jet_TorW", "DR_W_b", "Bprime_chi2"]

combinations = list(it.combinations(var_list, 2))
print(len(combinations)) # 351

# split into smaller pieces
combinations1 = combinations[:35]
combinations2 = combinations[35:70]
combinations3 = combinations[70:105]
combinations4 = combinations[105:140]
combinations5 = combinations[140:175]
combinations6 = combinations[175:210]
combinations7 = combinations[210:245]
combinations8 = combinations[245:280]
combinations9 = combinations[280:315]
combinations10 = combinations[315:]

# Define functions
gInterpreter.Declare("""     
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}                                                                                                                                        """)

def CreateHistos(Events_tag, combinations, sample):
    histfile.cd()
    for branch1, branch2 in combinations:
        histo_tag = Events_tag.Histo2D(("", "", branches[branch1][1], branches[branch1][2], branches[branch1][3], branches[branch2][1], branches[branch2][2], branches[branch2][3]), branch1, branch2, "weights")
        histo_tag.Write(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + case)

def AddHistos(combinations, sampleList):
    for branch1, branch2 in combinations:
        histo1 = histfile.Get(branch1 + "_vs_" + branch2 + "_" + sampleList[0] + "_weighted" + case)
        for i in range(1, len(sampleList)):
            histo =  histfile.Get(branch1 + "_vs_" + branch2 + "_" + sampleList[i] + "_weighted" + case)
            histo1.Add(histo)
        
        histo1.Write(branch1 + "_vs_" + branch2 + "_bkg_weighted" + case)

def CreateFromSamples(sample):
    print("Processing {}".format(sample.prefix))
    
    samplename = sample.samplename.split('/')[1]
    tfiles = readTreeNominal(samplename,step1dir,"Events")

    Events = RDataFrame(tfiles).Define("weights","weights(genWeight,{},{},{})".format(lumi,sample.xsec,sample.nrun)).Filter(tags_cases[case]).Define(
        "leadingOSFatJet_TorW", "max(gcOSFatJet_pNetT[0], gcOSFatJet_pNetW[0])"
    ).Define(
        "leadingOSFatJet_tau", "min(gcOSFatJet_tau21[0],gcOSFatJet_tau32[0])"
    ).Define("leadingOSFatJet_TWorQCD", "max(gcOSFatJet_pNetTvsQCD[0], gcOSFatJet_pNetWvsQCD[0])")
    
    # Create histogram from all possible 2d phases
    CreateHistos(Events, combinations1, sample.prefix)
    CreateHistos(Events, combinations2, sample.prefix)
    CreateHistos(Events, combinations3, sample.prefix)
    CreateHistos(Events, combinations4, sample.prefix)
    CreateHistos(Events, combinations5, sample.prefix)
    #CreateHistos(Events, combinations6, sample.prefix)
    #CreateHistos(Events, combinations7, sample.prefix)
    #CreateHistos(Events, combinations8, sample.prefix)
    #CreateHistos(Events, combinations9, sample.prefix)
    #CreateHistos(Events, combinations10, sample.prefix)

##################
# Get Histograms #
##################
### Two histograms of interest: _bkg and _Bp1400 ######
bkgList = ["QCDHT3002018UL", "QCDHT5002018UL", "QCDHT7002018UL", "QCDHT10002018UL", "QCDHT15002018UL", "QCDHT20002018UL", 
           "TTToSemiLeptonic2018UL", 
           "WJetsHT2002018UL", "WJetsHT4002018UL", "WJetsHT6002018UL", "WJetsHT8002018UL", "WJetsHT12002018UL", "WJetsHT25002018UL"
]

if(getHistos):
    start_time1 = time.time()

    print("Preparing bkgsig_histos2D.root...")
    histfile = TFile.Open("bkgsig_histos2D.root", "RECREATE")

    # Create histograms from samples
    CreateFromSamples(Bprime_M1400_2018UL)
    CreateFromSamples(QCDHT3002018UL)
    CreateFromSamples(QCDHT5002018UL)
    CreateFromSamples(QCDHT7002018UL)
    CreateFromSamples(QCDHT10002018UL)
    CreateFromSamples(QCDHT15002018UL)
    CreateFromSamples(QCDHT20002018UL)
    CreateFromSamples(TTToSemiLeptonic2018UL)
    CreateFromSamples(WJetsHT2002018UL)
    CreateFromSamples(WJetsHT4002018UL)
    CreateFromSamples(WJetsHT6002018UL)
    CreateFromSamples(WJetsHT8002018UL)
    CreateFromSamples(WJetsHT12002018UL)
    CreateFromSamples(WJetsHT25002018UL)

    # Add background histograms together
    AddHistos(combinations1, bkgList)
    AddHistos(combinations2, bkgList)
    AddHistos(combinations3, bkgList)
    AddHistos(combinations4, bkgList)
    AddHistos(combinations5, bkgList)
    #AddHistos(combinations6, bkgList)
    #AddHistos(combinations7, bkgList)
    #AddHistos(combinations8, bkgList)
    #AddHistos(combinations9, bkgList)
    #AddHistos(combinations10, bkgList)

    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)


########
# Plot #
########
print("plotting...")
start_time2 = time.time()
histfile = TFile.Open("bkgsig_histos2D.root", "READ")    

sig = "Bprime_M1400_2018UL"
def plot2D(combinations):
    for branch1, branch2 in combinations:
        # get base histograms
        histfile.cd()
        hist_sig = histfile.Get(branch1 + "_vs_" + branch2 + "_" + sig + "_weighted" + case)
        hist_bkg = histfile.Get(branch1 + "_vs_" + branch2 + "_" + "bkg" + "_weighted" + case)
        hist_bkg_copy = hist_bkg.Clone()
        hist_purity = hist_sig.Clone()
        hist_sensitivity = hist_sig.Clone()

        # set branch names
        xname = branch1+branches[branch1][-1]
        yname = branch2+branches[branch2][-1]

        # plot signal
        c1.cd(1)
        hist_sig.GetXaxis().SetTitle(xname)
        hist_sig.GetYaxis().SetTitle(yname)
        hist_sig.SetTitle("Signal")
        hist_sig.Draw("COLZ")

        # plot backgrounds
        c1.cd(2)
        hist_bkg.GetXaxis().SetTitle(xname)
        hist_bkg.GetYaxis().SetTitle(yname)
        hist_bkg.SetTitle("Backgrounds")
        hist_bkg.Draw("COLZ")

        # plot purity
        hist_purity.Divide(hist_bkg)
        c1.cd(3)
        hist_purity.GetXaxis().SetTitle(xname)
        hist_purity.GetYaxis().SetTitle(yname)
        hist_purity.SetTitle("Signal Purity S/B")
        hist_purity.Draw("COLZ")

        # plot signal sensitivity
        for bin in range(hist_bkg_copy.GetNcells()):
            BinContent = hist_bkg_copy.GetBinContent(bin) + hist_sig.GetBinContent(bin)
            if(BinContent>=0):
                hist_bkg_copy.SetBinContent(bin, sqrt(BinContent))
            else:
                hist_bkg_copy.SetBinContent(bin, 0)
        hist_sensitivity.Divide(hist_bkg_copy)

        c1.cd(4)
        hist_sensitivity.GetXaxis().SetTitle(xname)
        hist_sensitivity.GetYaxis().SetTitle(yname)
        hist_sensitivity.SetTitle("Signal Sensitivity S/sqrt(S+B)")
        hist_sensitivity.Draw("COLZ")

        c1.Modified()
        c1.Update()

        outname = outdir + subdir + branch1 + "_vs_" +branch2 + ".png"
        c1.SaveAs(outname)

if(getPlots):
    # set up canvas     
    c1 = TCanvas("c1", "c1", 600, 600)
    gStyle.SetOptStat(0)
    c1.Divide(2,2)
    c1_1.SetLeftMargin(0.12)
    c1_2.SetLeftMargin(0.12)
    c1_3.SetLeftMargin(0.12)
    c1_4.SetLeftMargin(0.12)
    c1_1.SetRightMargin(0.15)
    c1_2.SetRightMargin(0.15)
    c1_3.SetRightMargin(0.15)
    c1_4.SetRightMargin(0.15)

    plot2D(combinations1)
    plot2D(combinations2)
    plot2D(combinations3)
    plot2D(combinations4)
    plot2D(combinations5)
    #plot2D(combinations6)
    #plot2D(combinations7)
    #plot2D(combinations8)
    #plot2D(combinations9)
    #plot2D(combinations10)


end_time2 = time.time()
print("time elapsed: ", end_time2 - start_time2)    
