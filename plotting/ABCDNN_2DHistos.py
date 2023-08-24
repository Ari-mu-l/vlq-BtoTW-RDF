import os, sys, time
import itertools as it
import numpy as np
from math import sqrt
from ROOT import *

###############
#   Options   #
###############
case = "_BdecayCase1"
getHistos = False
getPlots = True

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Aug2023_2018/"
outdir = os.getcwd()+'/plots_bkgsig2D/'
subdir = case[1:]+"/"
if not os.path.exists(outdir + subdir): os.system('mkdir -p ' + outdir + subdir)

##################
# Weighting info #
##################
samples = {#'Bp800':'Bprime_M800_20UL18_hadd.root',
           #'Bp1000':'Bprime_M1000_20UL18_hadd.root',
           #'Bp1200':'Bprime_M1200_20UL18_hadd.root',
           #'Bp1300':'Bprime_M1300_20UL18_hadd.root',
           'Bp1400':'Bprime_M1400_20UL18_hadd.root',
           #'Bp1500':'Bprime_M1500_20UL18_hadd.root',
           #'Bp1600':'Bprime_M1600_20UL18_hadd.root',
           #'Bp1700':'Bprime_M1700_20UL18_hadd.root',
           #'Bp1800':'Bprime_M1800_20UL18_hadd.root',
           #'Bp2000':'Bprime_M2000_20UL18_hadd.root',
           #'Bp2200':'Bprime_M2200_20UL18_hadd.root',
           #'DYMHT1200':'DYMHT12002018UL_hadd.root',
           #'DYMHT200':'DYMHT2002018UL_hadd.root',
           #'DYMHT2500':'DYMHT25002018UL_hadd.root',
           #'DYMHT400':'DYMHT4002018UL_hadd.root',
           #'DYMHT600':'DYMHT6002018UL_hadd.root',
           #'DYMHT800':'DYMHT8002018UL_hadd.root',
           'QCD1000':'QCDHT10002018UL_hadd.root',
           'QCD1500':'QCDHT15002018UL_hadd.root',
           'QCD2000':'QCDHT20002018UL_hadd.root',
           #'QCD200':'QCDHT2002018UL_hadd.root', # no contribution
           'QCD300':'QCDHT3002018UL_hadd.root',
           'QCD500':'QCDHT5002018UL_hadd.root',
           'QCD700':'QCDHT7002018UL_hadd.root',
           #'STs':'STs2018UL_hadd.root',
           #'STt':'',
           #'STtW':'STtW2018UL_hadd.root',
           #'STtWb':'STtWb2018UL_hadd.root',
           #'STtb':'STtb2018UL_hadd.root',
           #'SingleElA':'SingleElectronRun2018A2018UL_hadd.root',
           #'SingleElB':'SingleElectronRun2018B2018UL_hadd.root',
           #'SingleElC':'SingleElectronRun2018C2018UL_hadd.root',
           #'SingleElD':'SingleElectronRun2018D2018UL_hadd.root',
           #'SingleMuA':'SingleMuonRun2018A2018UL_hadd.root',
           #'SingleMuB':'SingleMuonRun2018B2018UL_hadd.root',
           #'SingleMuC':'SingleMuonRun2018C2018UL_hadd.root',
           #'SingleMuD':'',
           #'TTHB':'TTHB2018UL_hadd.root',
           #'TTHnonB':'TTHnonB2018UL_hadd.root',
           #'TTMT1000':'TTMT10002018UL_hadd.root',
           #'TTMT700':'TTMT7002018UL_hadd.root',
           #'TTTo2L2Nu':'',
           #'TTTo2L2Nu_0':'',
           #'TTToHadronic':'',
           'TTToSemiLeptonic':'TTToSemiLeptonic2018UL_hadd.root',
           #'TTWl':'TTWl2018UL_hadd.root',
           #'TTWq':'TTWq2018UL_hadd.root',           
           #'TTZM10':'TTZM102018UL_hadd.root',
           #'TTZM1to10':'TTZM1to102018UL_hadd.root',
           'WJets1200':'WJetsHT12002018UL_hadd.root',
           'WJets200':'WJetsHT2002018UL_hadd.root',
           'WJets2500':'WJetsHT25002018UL_hadd.root',
           'WJets400':'WJetsHT4002018UL_hadd.root',
           'WJets600':'WJetsHT6002018UL_hadd.root',
           'WJets800':'WJetsHT8002018UL_hadd.root',
           #'WW':'WW2018UL_hadd.root',
           #'WZ':'WZ2018UL_hadd.root',
           #'ZZ':'ZZ2018UL_hadd.root',
}

lumi = 138000.0

nRun = {#"DY":,
        "Bp800":1000000.,
        "Bp1000":1000000.,
        "Bp1200":1000000.,
        "Bp1300":1000000.,
        "Bp1400":1000000.,
        "Bp1500":1000000.,
        "Bp1600":1000000.,
        "Bp1700":1000000.,
        "Bp1800":1000000.,
        "Bp2000":1000000.,
        "Bp2200":1000000.,
        "QCD1000":14394786.,
        "QCD1500":10411831.,
        "QCD2000":5374711.,
        "QCD200":57336623.,
        "QCD300":61609663.,
        "QCD500":49184771.,
        "QCD700":48506751.,
        "TTToSemiLeptonic":476408000.,
        "WJets1200":6481518.,
        "WJets200":58225632.,
        "WJets2500":2097648.,
        "WJets400":7444030.,
        "WJets600":7718765.,
        "WJets800":7306187.,
    }

xsec = {"Bp800":1.0,
        "Bp1000":1.0,
        "Bp1200":1.0,
        "Bp1300":1.0,
        "Bp1400":1.0,
        "Bp1500":1.0,
        "Bp1600":1.0,
        "Bp1700":1.0,
        "Bp1800":1.0,
        "Bp2000":1.0,
        "Bp2200":1.0,
        'QCD200':1712000,
        'QCD300':347700,
        'QCD500':32100,
        'QCD700':6831,
        'QCD1000':1207,
        'QCD1500':119.9,
        'QCD2000':25.24,
        'TTToSemiLeptonic':831.76*0.438,
        'WJets200':359.7*1.21,
        'WJets400':48.91*1.21,
        'WJets600':12.05*1.21,
        'WJets800':5.501*1.21,
        'WJets1200':1.329*1.21,
        'WJets2500':0.03216*1.21,
}

###### Add branches, categories, tags here ####

branches = {"nSignalIsoMu":["nSignalIsoMu", 10, 0, 10, ""],
            "nSignalIsoEl":["nSignalIsoEl", 10, 0, 10, ""],
            "nVetoIsoLep":["nVetoIsoLep", 10, 0, 10, ""],
            "lepton_pt":["lepton_pt", 50, 0, 1000, "[GeV]"],
            "lepton_eta":["lepton_eta", 40, -4, 4, ""],
            "lepton_miniIso":["lepton_miniIso", 50, 0, 0.2, ""],
            "NJets_central":["NJets_central", 20, 0, 20, ""],
            "NJets_DeepFlavL":["NJets_DeepFlavL", 20, 0, 20, ""],
            "NJets_forward":["NJets_forward", 20, 0, 20, ""],
            "NFatJets":["NFatJets", 10, 0, 10, ""],
            "NOS_gcJets_central":["NOS_gcJets_central", 20, 0, 20, ""],
            "NSS_gcJets_central":["NSS_gcJets_central", 20, 0, 20, ""],
            "NOS_gcJets_DeepFlavL":["NOS_gcJets_DeepFlavL", 15, 0, 15, ""],
            "NSS_gcJets_DeepFlavL":["NSS_gcJets_DeepFlavL", 15, 0, 15, ""],
            "NOS_gcFatJets":["NOS_gcFatJets", 10, 0, 10, ""], 
            # no NSS_gcFatJets yet. FIXME
            #"Jet_HT":["Jet_HT", 50, 0, 5000, "[GeV]"],
            #"Jet_ST":["Jet_ST", 50, 0, 5000, "[GeV]"],
            #"FatJet_pt_1":["FatJet_pt_1", 50, 0, 1500, "[GeV]"],
            #"FatJet_pt_2":["FatJet_pt_2", 50, 0, 1500, "[GeV]"],
            #"FatJet_sdMass_1":["FatJet_sdMass_1", 50, 0, 500, "[Gev]"],
            #"FatJet_sdMass_2":["FatJet_sdMass_2", 50, 0, 500, "[Gev]"],
            "dpak8_J_1":["dpak8_J_1", 50, 0, 1, ""],
            "dpak8_J_2":["dpak8_J_2", 50, 0, 1, ""],
            "dpak8_T_1":["dpak8_T_1", 50, 0, 1, ""],
            "dpak8_T_2":["dpak8_T_2", 50, 0, 1, ""],
            "dpak8_W_1":["dpak8_W_1", 50, 0, 1, ""],
            "dpak8_W_2":["dpak8_W_2", 50, 0, 1, ""],
            "dpak8_tag_1":["dpak8_tag_1", 3, 0, 3, ""],
            "dpak8_tag_2":["dpak8_tag_1", 3, 0, 3, ""],
            "nJ_dpak8":["nJ_dpak8", 20, 0, 20, ""],
            "nT_dpak8":["nJ_dpak8", 15, 0, 15, ""],
            "nW_dpak8":["nW_dpak8", 15, 0, 15, ""],
            "pNet_J_1":["pNet_J_1", 50, 0, 1, ""],
            "pNet_J_2":["pNet_J_2", 50, 0, 1, ""],
            "pNet_T_1":["pNet_T_1", 50, 0, 1, ""],
            "pNet_T_2":["pNet_T_2", 50, 0, 1, ""],
            "pNet_W_1":["pNet_W_1", 50, 0, 1, ""],
            "pNet_W_2":["pNet_W_2", 50, 0, 1, ""],
            "pNet_tag_1":["pNet_tag_1", 3, 0, 3, ""],
            "pNet_tag_2":["pNet_tag_2", 3, 0, 3, ""],
            "nJ_pNet":["nJ_pNet", 20, 0, 20, ""],
            "nT_pNet":["nT_pNet", 15, 0, 15, ""],
            "nW_pNet":["nW_pNet", 15, 0, 15, ""],
            "tau21_1":["tau21_1", 50, 0, 1, ""],
            "tau21_2":["tau21_2", 50, 0, 1, ""],
            "minDR_lep_FatJet":["minDR_lep_FatJet", 50, 0, 5, "[GeV]"],
            "ptRel_lep_FatJet":["ptRel_lep_FatJet", 50, 0, 500, "[GeV]"],
            "minDR_leadAK8otherAK8":["minDR_leadAK8otherAK8", 50, 0, 5, "[GeV]"],
            "minDR_lep_Jet":["minDR_lep_Jet", 50, 0, 5, "[GeV]"],
            "ptRel_lep_Jet":["ptRel_lep_Jet", 50, 0, 500, "[GeV]"],
            "W_pt":["W_pt", 50, 0, 1000, "[GeV]"],
            "W_eta":["W_eta", 40, -4, 4, ""],
            "W_MT":["W_MT", 50, 0, 1500, "[GeV]"],
            "DR_W_lep":["DR_W_lep", 50, 0, 5, ""],
            "minM_lep_Jet":["minM_lep_Jet", 50, 0, 1000, "[GeV]"],
            "t_pt":["t_pt", 50, 0, 1000, "[GeV]"],
            "t_eta":["t_eta", 50, -4, 4, ""],
            "DR_W_b":["DR_W_b", 50, 0, 7, ""],
            "Bprime_chi2":["Bprime_chi2", 50, 0, 1000, ""],
            "Bdecay_obs":["Bdecay_obs", 4, 0, 4, ""]
}

categories = {#"DY": ["DYMHT1200", "DYMHT200", "DYMHT2500"],       
              "QCD": ["QCD1000", "QCD1500", "QCD2000", #"QCD200", 
                      "QCD300", "QCD500", "QCD700"],
              "WJets":["WJets1200", "WJets200", "WJets2500", "WJets400", "WJets600", "WJets800"],
              #"TTbar":["TTToSemiLeptonic"],
          }

tags_cases = {"_BdecayCase1":"Bdecay_obs==1",
        "_BdecayCase2":"Bdecay_obs==2",
        "_BdecayCase3":"Bdecay_obs==3",
        "_BdecayCase4":"Bdecay_obs==4",
    }

##########################
# 2D phase space choices #
##########################
var_list = ["nSignalIsoMu", "nSignalIsoEl", "nVetoIsoLep", 
"lepton_pt", "lepton_eta", "lepton_miniIso", "NJets_central", "NJets_DeepFlavL", "NJets_forward", "NFatJets", "NOS_gcJets_central", "NSS_gcJets_central", "NOS_gcJets_DeepFlavL", "NSS_gcJets_DeepFlavL", "NOS_gcFatJets", "dpak8_J_1", "dpak8_J_2", "dpak8_T_1", "dpak8_T_2", "dpak8_W_1", "dpak8_W_2", "dpak8_tag_1", "dpak8_tag_2", "nJ_dpak8", "nT_dpak8", "nW_dpak8", "pNet_J_1", "pNet_J_2", "pNet_T_1", "pNet_T_2", "pNet_W_1", "pNet_W_2", "pNet_tag_1", "pNet_tag_2", "nJ_pNet", "nT_pNet", "nW_pNet", "tau21_1", "tau21_2", "minDR_lep_FatJet", "ptRel_lep_FatJet", "minDR_leadAK8otherAK8", "minDR_lep_Jet", "ptRel_lep_Jet", "W_pt", "W_eta", "W_MT", "DR_W_lep", "minM_lep_Jet", "t_pt", "t_eta", "DR_W_b", "Bprime_chi2"
]

combinations = list(it.combinations(var_list, 2))

print(len(combinations)) #1378

# split into smaller pieces
combinations1 = combinations[:50]
combinations2 =combinations[50:100]
combinations3 =combinations[100:150]
combinations4 =combinations[150:200]
combinations5 =combinations[200:250]
combinations6 =combinations[250:300]
combinations7 =combinations[300:350]
combinations8 =combinations[350:400]
combinations9 =combinations[400:450]
combinations10 =combinations[450:500]
combinations11 =combinations[500:550]
combinations12 =combinations[550:600]
combinations13 =combinations[600:650]
combinations14 =combinations[650:700]
combinations15 =combinations[700:750]
combinations16 =combinations[750:800]
combinations17 =combinations[800:850]
combinations18 =combinations[850:900]
combinations19 =combinations[900:950]
combinations20 =combinations[950:1000]
combinations21 =combinations[1000:1050]
combinations22 =combinations[1050:1100]
combinations23 =combinations[1100:1150]
combinations24 =combinations[1150:1200]
combinations25 =combinations[1200:1250]
combinations26 =combinations[1250:1300]
combinations27 =combinations[1300:]

####################
# Define functions #
####################
def CreateHistos(Events_tag, combinations, sample):
    histfile.cd()
    for branch1, branch2 in combinations:
        histo_tag = Events_tag.Histo2D(("", "", branches[branch1][1], branches[branch1][2], branches[branch1][3], branches[branch2][1], branches[branch2][2], branches[branch2][3]), branch1, branch2, "weights")
        histo_tag.Write(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + case)

def AddHistos(combinations, sampleList):
    for branch1, branch2 in combinations:
        histo1 = histfile.Get(branch1 + "_vs_" + branch2 + "_" + sampleList[0] + "_weighted" + case)
        histo1.Draw()
        for i in range(1, len(sampleList)):
            histo =  histfile.Get(branch1 + "_vs_" + branch2 + "_" + sampleList[i] + "_weighted" + case)
            histo1.Add(histo)
        
        histo1.Write(branch1 + "_vs_" + branch2 + "_bkg_weighted" + case)

def CreateFromSamples(sample):
    print("Processing" + sample)
    filename = indir + samples[sample]
    Events = RDataFrame("Events", filename).Filter("NJets_forward>0 && Bprime_mass>0").Define("weights","weights(genWeight,{},{},{})".format(lumi,xsec[sample],nRun[sample]))
    Events_tag = Events.Filter(tags_cases[case])

    # Create histogram from all possible 2d phases
    CreateHistos(Events_tag, combinations1, sample)
    CreateHistos(Events_tag, combinations2, sample)
    CreateHistos(Events_tag, combinations3, sample)
    #CreateHistos(Events_tag, combinations4, sample)
    #CreateHistos(Events_tag, combinations5, sample)
    #CreateHistos(Events_tag, combinations6, sample)
    #CreateHistos(Events_tag, combinations7, sample)
    #CreateHistos(Events_tag, combinations8, sample)
    #CreateHistos(Events_tag, combinations9, sample)
    #CreateHistos(Events_tag, combinations10, sample)
    #CreateHistos(Events_tag, combinations11, sample)
    #CreateHistos(Events_tag, combinations12, sample)
    #CreateHistos(Events_tag, combinations13, sample)
    #CreateHistos(Events_tag, combinations14, sample)
    #CreateHistos(Events_tag, combinations15, sample)
    #CreateHistos(Events_tag, combinations16, sample)
    #CreateHistos(Events_tag, combinations17, sample)
    #CreateHistos(Events_tag, combinations18, sample)
    #CreateHistos(Events_tag, combinations19, sample)
    #CreateHistos(Events_tag, combinations20, sample)
    #CreateHistos(Events_tag, combinations21, sample)
    #CreateHistos(Events_tag, combinations22, sample)
    #CreateHistos(Events_tag, combinations23, sample)
    #CreateHistos(Events_tag, combinations24, sample)
    #CreateHistos(Events_tag, combinations25, sample)
    #CreateHistos(Events_tag, combinations26, sample)
    #CreateHistos(Events_tag, combinations27, sample)

##################
# Get Histograms #
##################
### Two histograms of interest: _bkg and _Bp1400 ######
bkgList = ["QCD300", "QCD500", "QCD700", "QCD1000", "QCD1500", "QCD2000", "TTToSemiLeptonic", "WJets200", "WJets400", "WJets600", "WJets800", "WJets1200", "WJets2500"]

if(getHistos):
    start_time1 = time.time()

    gInterpreter.Declare("""
    float weights( float genWeight, float lumi, float xsec, float nRun ){
    return genWeight * lumi * xsec / (nRun * abs(genWeight));
    }
    """)

    print("Preparing bkgsig_histos2D.root...")
    histfile = TFile.Open("bkgsig_histos2D.root", "RECREATE")

    # Create histograms from samples
    CreateFromSamples('Bp1400')
    CreateFromSamples('QCD300')
    CreateFromSamples('QCD500')
    CreateFromSamples('QCD700')
    CreateFromSamples('QCD1000')
    CreateFromSamples('QCD1500')
    CreateFromSamples('QCD2000')
    CreateFromSamples('TTToSemiLeptonic')
    CreateFromSamples('WJets200')
    CreateFromSamples('WJets400')
    CreateFromSamples('WJets600')
    CreateFromSamples('WJets800')
    CreateFromSamples('WJets1200')
    CreateFromSamples('WJets2500')

    # Add background histograms together
    AddHistos(combinations1, bkgList)
    AddHistos(combinations2, bkgList)
    AddHistos(combinations3, bkgList)
    #AddHistos(combinations4, bkgList)
    #AddHistos(combinations5, bkgList)
    #AddHistos(combinations6, bkgList)
    #AddHistos(combinations7, bkgList)
    #AddHistos(combinations8, bkgList)
    #AddHistos(combinations9, bkgList)
    #AddHistos(combinations10, bkgList)
    #AddHistos(combinations11, bkgList)
    #AddHistos(combinations12, bkgList)
    #AddHistos(combinations13, bkgList)
    #AddHistos(combinations14, bkgList)
    #AddHistos(combinations15, bkgList)
    #AddHistos(combinations16, bkgList)
    #AddHistos(combinations17, bkgList)
    #AddHistos(combinations18, bkgList)
    #AddHistos(combinations19, bkgList)
    #AddHistos(combinations20, bkgList)
    #AddHistos(combinations21, bkgList)
    #AddHistos(combinations22, bkgList)
    #AddHistos(combinations23, bkgList)
    #AddHistos(combinations24, bkgList)
    #AddHistos(combinations25, bkgList)
    #AddHistos(combinations26, bkgList)
    #AddHistos(combinations27, bkgList)

    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)


########
# Plot #
########
print("plotting...")
start_time2 = time.time()
histfile = TFile.Open("bkgsig_histos2D.root", "READ")    

sig = "Bp1400"
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
    #plot2D(combinations4)
    #plot2D(combinations5)
    #plot2D(combinations6)
    #plot2D(combinations7)
    #plot2D(combinations8)
    #plot2D(combinations9)
    #plot2D(combinations10)
    #plot2D(combinations11)
    #plot2D(combinations12)
    #plot2D(combinations13)
    #plot2D(combinations14)
    #plot2D(combinations15)
    #plot2D(combinations16)
    #plot2D(combinations17)
    #plot2D(combinations18)
    #plot2D(combinations19)
    #plot2D(combinations20)
    #plot2D(combinations21)
    #plot2D(combinations22)
    #plot2D(combinations23)
    #plot2D(combinations24)
    #plot2D(combinations25)
    #plot2D(combinations26)
    #plot2D(combinations27)


end_time2 = time.time()
print("time elapsed: ", end_time2 - start_time2)
