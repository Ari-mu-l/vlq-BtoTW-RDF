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

###############
#   Options   #
###############
case = "_BdecayCase2and3"
getHistos = True
getRegions = True
branch1 = "NJets_forward"
branch2 = "NJets_DeepFlavL"
lumi = 138000.0

sigList = ["Bprime_M800_2018UL", "Bprime_M1400_2018UL", "Bprime_M2200_2018UL"]

step1dir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Sep2023_2018/"

outdir = os.getcwd()+'/plots_ABCDnnRegion/'
subdir = case[1:]+"/"
if not os.path.exists(outdir + subdir): os.system('mkdir -p ' + outdir + subdir)

###### Add branches, categories, tags here ####
branches = {
#            "lepton_pt":["lepton_pt", 50, 0, 1000, "[GeV]"],
#            "lepton_eta":["lepton_eta", 40, -4, 4, ""],
#            "lepton_miniIso":["lepton_miniIso", 50, 0, 0.2, ""],
            "NJets_central":["NJets_central", 10, 0, 10, ""],
            "NJets_DeepFlavL":["NJets_DeepFlavL", 8, 0, 8, ""],
            "NJets_forward":["NJets_forward", 8, 0, 8, ""],
            "NFatJets":["NFatJets", 10, 0, 10, ""],
            "NOS_gcJets_central":["NOS_gcJets_central", 10, 0, 10, ""],
            "NSS_gcJets_central":["NSS_gcJets_central", 10, 0, 10, ""],
            "NOS_gcJets_DeepFlavL":["NOS_gcJets_DeepFlavL", 10, 0, 10, ""],
            "NSS_gcJets_DeepFlavL":["NSS_gcJets_DeepFlavL", 10, 0, 10, ""],
            "NOS_gcFatJets":["NOS_gcFatJets", 10, 0, 10, ""], 
            "NSS_gcFatJets":["NSS_gcFatJets", 10, 0, 10, ""],
            #"gcJet_HT":["gcJet_HT", 30, 0, 5000, "[GeV]"], # transform var
            "gcJet_ST":["gcJet_ST", 25, 0, 5000, "[GeV]"],
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
            "minDR_lep_FatJet":["minDR_lep_FatJet", 30, 0, 5, "[GeV]"],
            "ptRel_lep_FatJet":["ptRel_lep_FatJet", 30, 0, 500, "[GeV]"],
            "minDR_leadAK8otherAK8":["minDR_leadAK8otherAK8", 20, 0, 10, "[GeV]"],
            "minDR_leadAK4otherAK4":["minDR_leadAK8otherAK8", 20, 0, 10, "[GeV]"],
######            "minDR_AK8s_discrete":["minDR_AK8s_discrete", 14, 0, 14, ""],
######            "minDR_AK4s_discrete":["minDR_AK4s_discrete", 18, 0, 18, ""],
#            "minDR_lep_Jet":["minDR_lep_Jet", 50, 0, 5, "[GeV]"],
#            "ptRel_lep_Jet":["ptRel_lep_Jet", 50, 0, 500, "[GeV]"],
#            "W_pt":["W_pt", 30, 0, 1000, "[GeV]"],
            "W_eta":["W_eta", 16, -4, 4, ""],
            "W_MT":["W_MT", 30, 0, 1500, "[GeV]"],
            "DR_W_lep":["DR_W_lep", 20, 0, 5, ""],
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
              "_BdecayCase1and4": "(Bdecay_obs==1) || (Bdecay_obs==4)",
              "_BdecayCase2and3": "(Bdecay_obs==2) || (Bdecay_obs==3)",
          }

####################
# Define functions #
####################
gInterpreter.Declare("""     
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}
""")

def CreateHistos(Events_tag, sample):
    histfile.cd()

    histo_tag = Events_tag.Histo2D(("", "", branches[branch1][1], branches[branch1][2], branches[branch1][3], branches[branch2][1], branches[branch2][2], branches[branch2][3]), branch1, branch2, "weights")
    histo_tag.Write(branch1 + "_vs_" + branch2 + "_" + sample + "_weighted" + case)

def AddHistos(sampleList):
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

    CreateHistos(Events, sample.prefix)

##################
# Get Histograms #
##################
### Two histograms of interest: _bkg and _Bp1400 ######
histfile_name =  "ABCDnn_Regions_{}_vs_{}{}.root".format(branch1, branch2, case)

bkgList = ["QCDHT3002018UL", "QCDHT5002018UL", "QCDHT7002018UL", "QCDHT10002018UL", "QCDHT15002018UL", "QCDHT20002018UL", 
           "TTToSemiLeptonic2018UL", 
           "WJetsHT2002018UL", "WJetsHT4002018UL", "WJetsHT6002018UL", "WJetsHT8002018UL", "WJetsHT12002018UL", "WJetsHT25002018UL"
]

if(getHistos):
    start_time1 = time.time()

    print("Preparing {} ...".format(histfile_name))
    histfile = TFile.Open(histfile_name, "RECREATE")

    CreateFromSamples(Bprime_M800_2018UL)
    CreateFromSamples(Bprime_M1400_2018UL)
    CreateFromSamples(Bprime_M2200_2018UL)
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
    AddHistos(bkgList)

    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)


########
# Plot #
########
print("plotting...")
start_time2 = time.time()
histfile = TFile.Open(histfile_name, "READ")    

def plotRegions(line1, line2, line3, line4, line5):
    line1.Draw()
    line2.Draw()
    #line3.Draw()
    line4.Draw()
    #line5.Draw()
        
def integrateRegions(hist_sig, hist_bkg, region, xmin, xmax, ymin, ymax, text, c1):
    integral_sig = hist_sig.Integral(hist_sig.GetXaxis().FindBin(xmin), hist_sig.GetXaxis().FindBin(xmax), hist_sig.GetYaxis().FindBin(ymin), hist_sig.GetYaxis().FindBin(ymax))
    integral_bkg = hist_bkg.Integral(hist_bkg.GetXaxis().FindBin(xmin), hist_bkg.GetXaxis().FindBin(xmax), hist_bkg.GetYaxis().FindBin(ymin), hist_bkg.GetYaxis().FindBin(ymax))

    c1.cd(1)
    
    content = "#color[2]{"+region+"}"
    text.DrawLatex((xmin+xmax)/2, (ymin+ymax)/2, content)

    c1.cd(2)
    text.DrawLatex((xmin+xmax)/2, (ymin+ymax)/2, content)

    print("Region {}: Signal {}, Background {}".format(region, int(integral_sig), int(integral_bkg)))
    return integral_sig, integral_bkg

def plot2D(case, sig):

    hist_sig = histfile.Get(branch1 + "_vs_" + branch2 + "_" + sig + "_weighted" + case)
    hist_bkg = histfile.Get(branch1 + "_vs_" + branch2 + "_" + "bkg" + "_weighted" + case)
    hist_bkg_copy = hist_bkg.Clone()
    hist_purity = hist_sig.Clone()
    hist_sensitivity = hist_sig.Clone()

    # plot signal
    c1.cd(1)
    hist_sig.GetXaxis().SetTitle(xname)
    hist_sig.GetYaxis().SetTitle(yname)
    hist_sig.SetTitle("Signal")
    hist_sig.Draw("COLZ")
    if(getRegions):
        plotRegions(line1, line2, line3, line4, line5)

    # plot background
    c1.cd(2)
    hist_bkg.GetXaxis().SetTitle(xname)
    hist_bkg.GetYaxis().SetTitle(yname)
    hist_bkg.SetTitle("Backgrounds")
    hist_bkg.Draw("COLZ")
    if(getRegions):
        plotRegions(line1, line2, line3, line4, line5)

    # plot purity
    hist_purity.Divide(hist_bkg)
    
    c1.cd(3)
    hist_purity.GetXaxis().SetTitle(xname)
    hist_purity.GetYaxis().SetTitle(yname)
    hist_purity.SetTitle("Signal Purity S/B")
    hist_purity.Draw("COLZ")
    if(getRegions):
        plotRegions(line1, line2, line3, line4, line5)

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
    if(getRegions):
        plotRegions(line1, line2, line3, line4, line5)

    if(getRegions):
        # define region ranges
        Dx_min, Cx_min, Yx_min = 1, 1, 1
        Dx_max, Cx_max, Yx_max = 8, 8, 8 # probably should be changed to inf
        Bx_min, Ax_min, Xx_min = 0, 0, 0
        Bx_max, Ax_max, Xx_max = 1, 1, 1
    
        Dy_min, By_min = 0, 0
        Dy_max, By_max = 3, 3
        Ay_min, Cy_min = 3, 3
        Ay_max, Cy_max = 4, 4
        Xy_min, Yy_min = 4, 4
        Xy_max, Yy_max = 8, 8

        # integrate histograms over regions
        text = TLatex()

        integrateRegions(hist_sig, hist_bkg, "A", Ax_min+0.01, Ax_max-0.01, Ay_min+0.01, Ay_max-0.01, text, c1)
        integrateRegions(hist_sig, hist_bkg, "B", Bx_min+0.01, Bx_max-0.01, By_min+0.01, By_max-0.01, text, c1)
        integrateRegions(hist_sig, hist_bkg, "C", Cx_min+0.01, Cx_max-0.01, Cy_min+0.01, Cy_max-0.01, text, c1)
        integrateRegions(hist_sig, hist_bkg, "D", Dx_min+0.01, Dx_max-0.01, Dy_min+0.01, Dy_max-0.01, text, c1)
        integrateRegions(hist_sig, hist_bkg, "X", Xx_min+0.01, Xx_max-0.01, Xy_min+0.01, Xy_max-0.01, text, c1)
        integrateRegions(hist_sig, hist_bkg, "Y", Yx_min+0.01, Yx_max-0.01, Yy_min+0.01, Yy_max-0.01, text, c1)

    c1.Modified()
    c1.Update()

    outname = outdir + subdir + branch1 + "_vs_" +branch2 + "_" + sig + ".png"
    c1.SaveAs(outname)


# Set up canvas

c1 = TCanvas("c1", "c1", 600, 600)
gStyle.SetOptStat(0)
c1.Divide(2,2)

c1_1.SetLeftMargin(0.12)
c1_1.SetRightMargin(0.12)
c1_2.SetLeftMargin(0.12)
c1_2.SetRightMargin(0.12)
c1_3.SetLeftMargin(0.12)
c1_3.SetRightMargin(0.12)
c1_4.SetLeftMargin(0.12)
c1_4.SetRightMargin(0.12)

# Define lines
if(getRegions):
    line1 = TLine(0, 3, 8, 3)
    line2 = TLine(1, 0, 1, 8)
    line3 = TLine(4, 0, 4, 8)
    line4 = TLine(0, 4, 8, 4)
    line5 = TLine(2, 1, 2, 8)

    line1.SetLineColor(kRed)
    line2.SetLineColor(kRed)
    line3.SetLineColor(kRed)
    line4.SetLineColor(kRed)
    line5.SetLineColor(kRed)

    line1.SetLineWidth(2)
    line2.SetLineWidth(2)
    line3.SetLineWidth(2)
    line4.SetLineWidth(2)
    line5.SetLineWidth(2)

xname = branch1+branches[branch1][-1]
yname = branch2+branches[branch2][-1]
    
plot2D(case, sigList[0])
plot2D(case, sigList[1])
plot2D(case, sigList[2])

end_time2 = time.time()
print("time elapsed: ", end_time2 - start_time2)


