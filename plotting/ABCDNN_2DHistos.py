import os, sys, time
import itertools as it
import numpy as np
from ROOT import *

###############
#   Options   #
###############
getHistos = True

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Aug2023_2018/"
outdir = os.getcwd()+'/plots_bkgsig2D/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

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
            "NOS_gcFatJets":["NOS_gcFatJets", 10, 0, 10, ""], # no NSS_gcFatJets yet. FIXME
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
#"lepton_pt", "lepton_eta", "lepton_miniIso", "NJets_central", "NJets_DeepFlavL", "NJets_forward", "NFatJets", "NOS_gcJets_central", "NSS_gcJets_central", "NOS_gcJets_DeepFlavL", "NSS_gcJets_DeepFlavL", "NOS_gcFatJets", "dpak8_J_1", "dpak8_J_2", "dpak8_T_1", "dpak8_T_2", "dpak8_W_1", "dpak8_W_2", "dpak8_tag_1", "dpak8_tag_2", "nJ_dpak8", "nT_dpak8", "nW_dpak8", "pNet_J_1", "pNet_J_2", "pNet_T_1", "pNet_T_2", "pNet_W_1", "pNet_W_2", "pNet_tag_1", "pNet_tag_2", "nJ_pNet", "nT_pNet", "nW_pNet", "tau21_1", "tau21_2", "minDR_lep_FatJet", "ptRel_lep_FatJet", "minDR_leadAK8otherAK8", "minDR_lep_Jet", "ptRel_lep_Jet", "W_pt", "W_eta", "W_MT", "DR_W_lep", "minM_lep_Jet", "t_pt", "t_eta", "DR_W_b", "Bprime_chi2"
]

combinations = list(it.combinations(var_list, 2))

####################
# Define functions #
####################
def CreateHistos(Events, tags, branches, sample, tag):
    Events_tag = Events.Filter(tags[tag])

    for branch1, branch2 in combinations:
        histo_tag = Events_tag.Histo2D(("", "", branches[branch1][1], branches[branch1][2], branches[branch1][3], branches[branch2][1], branches[branch2][2], branches[branch2][3]), branch1, branch2, "weights")

        histfile.cd()
        histo_tag.Write(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + tag)
        print(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + tag)
    
def AddHistos(branches, sampleList, tag):
    histo1 = histfile.Get(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + tag)
    for i in range(1, len(sampleList)):
        histo =  histfile.Get(branch1 + "_vs_" + branch2 + "_" + sampleList[i] + "_weighted" + tag)
        histo1.Add(histo)
        
    histo1.Write(branch1 + "_vs_" + branch2 + "_bkg_weighted" + tag)

##################
# Get Histograms #
##################
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

    for sample in samples:
        print("Processing" + sample)
        filename = indir + samples[sample]
        Events = RDataFrame("Events", filename).Filter("NJets_forward>0 && Bprime_mass>0").Define("weights","weights(genWeight,{},{},{})".format(lumi,xsec[sample],nRun[sample]))
        
        CreateHistos(Events, tags_cases, branches, sample, "_BdecayCase1")
        #CreateHistos(Events, tags_cases, branches, sample, "_BdecayCase2")
        #CreateHistos(Events, tags_cases, branches, sample, "_BdecayCase3")
        #CreateHistos(Events, tags_cases, branches, sample, "_BdecayCase4")

    for branch1, branch2 in combinations:
        AddHistos(branches, bkgList, "_BdecayCase1")
        #AddHistos(branches, bkgList, "_BdecayCase2")
        #AddHistos(branches, bkgList, "_BdecayCase3")
        #AddHistos(branches, bkgList, "_BdecayCase4")

    histfile.Close()        
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)


    
