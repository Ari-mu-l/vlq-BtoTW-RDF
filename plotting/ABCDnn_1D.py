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
yLog = False
transform1 = "Bprime_mass"
transform2 = "gcJet_ST"
lumi = 138000.0

sigList = ["Bprime_M800_2018UL", "Bprime_M1400_2018UL", "Bprime_M2200_2018UL"]

step1dir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Sep2023_2018/"

outdir = os.getcwd()+'/plots_ABCDnnRegion1D/'
subdir = case[1:]+"/"
if not os.path.exists(outdir + subdir): os.system('mkdir -p ' + outdir + subdir)

Regions = {"A":"(NJets_DeepFlavL==3) && (NJets_forward==0)",
           "B":"(NJets_DeepFlavL<3)  && (NJets_forward==0)",
           "C":"(NJets_DeepFlavL==3) && (NJets_forward>=1)",
           "D":"(NJets_DeepFlavL<3)  && (NJets_forward>=1)",
           "X":"(NJets_DeepFlavL>=4) && (NJets_forward==0)",
           "Y":"(NJets_DeepFlavL>=4) && (NJets_forward>=1)",
}

###### Add branches, categories, tags here ####
branches = {"gcJet_ST":["gcJet_ST", 100, 0, 5000, "[GeV]"],
            "Bprime_mass":["Bprime_mass", 100, 0, 6000, "[GeV]"],
}

tags_cases = {"_BdecayCase1":"Bdecay_obs==1",
              "_BdecayCase2":"Bdecay_obs==2",
              "_BdecayCase3":"Bdecay_obs==3",
              "_BdecayCase4":"Bdecay_obs==4",
              "_BdecayCase1and4": "(Bdecay_obs==1) || (Bdecay_obs==4)",
              "_BdecayCase2and3": "(Bdecay_obs==2) || (Bdecay_obs==3)",
          }


# Define functions
gInterpreter.Declare("""     
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}
""")

def CreateHistos(Events_tag, sample, region):
    # get events in the region
    Events_region = Events_tag.Filter(Regions[region])
    #print("Number of events (weighted) in region {}: {}".format(region, int(Events_region.Count())))
    # get histograms in the region
    histfile.cd()
    histo_tag1 = Events_region.Histo1D(("", "", branches[transform1][1], branches[transform1][2], branches[transform1][3]), transform1, "weights")
    histo_tag2 = Events_region.Histo1D(("", "", branches[transform2][1], branches[transform2][2], branches[transform2][3]), transform2, "weights")
    histo_tag1.Write("{}_{}_weighted{}_{}".format(transform1, sample, case, region))
    histo_tag2.Write("{}_{}_weighted{}_{}".format(transform2, sample, case, region))

def CreateFromSamples(sample):
    prefix = sample.prefix
    print("Processing {}".format(prefix))
    # read files
    samplename = sample.samplename.split('/')[1]
    tfiles = readTreeNominal(samplename,step1dir,"Events")
    # load RDF
    Events_tag = RDataFrame(tfiles).Define("weights","weights(genWeight,{},{},{})".format(lumi,sample.xsec,sample.nrun)).Filter(tags_cases[case])
    
    # create histos
    CreateHistos(Events_tag, prefix, "A")
    CreateHistos(Events_tag, prefix, "B")
    CreateHistos(Events_tag, prefix, "C")
    CreateHistos(Events_tag, prefix, "D")
    CreateHistos(Events_tag, prefix, "X")
    CreateHistos(Events_tag, prefix, "Y")

def AddHistos(bkg, region):
    bkgList = bkgDir[bkg]
    histo1 = histfile.Get("{}_{}_weighted{}_{}".format(transform1, bkgList[0], case, region))
    histo2 = histfile.Get("{}_{}_weighted{}_{}".format(transform2, bkgList[0], case, region))
    for i in range(1, len(bkgList)):
        histo1_i =  histfile.Get("{}_{}_weighted{}_{}".format(transform1, bkgList[i], case, region))
        histo2_i =  histfile.Get("{}_{}_weighted{}_{}".format(transform2, bkgList[i], case, region))
        histo1.Add(histo1_i)
        histo2.Add(histo2_i)
    histo1.Write("{}_{}_weighted{}_{}".format(transform1, bkg, case, region))
    histo2.Write("{}_{}_weighted{}_{}".format(transform2, bkg, case, region))

# Get Histograms
histfile_name =  "ABCDnn_1D_{}_{}{}.root".format(transform1, transform2, case)  
bkgDir = {"QCD":["QCDHT3002018UL", "QCDHT5002018UL", "QCDHT7002018UL", "QCDHT10002018UL", "QCDHT15002018UL", "QCDHT20002018UL"], 
          "TTToSemiLeptonic":["TTToSemiLeptonic2018UL"], 
          "WJets": ["WJetsHT2002018UL", "WJetsHT4002018UL", "WJetsHT6002018UL", "WJetsHT8002018UL", "WJetsHT12002018UL", "WJetsHT25002018UL"],
      }

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
    # add histos
    for region in Regions:
        AddHistos("QCD", region)
        AddHistos("TTToSemiLeptonic", region)
        AddHistos("WJets", region)
    # close file
    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)

# Plot
print("plotting...")
start_time2 = time.time()
# Define colors
colors_sig = {"Bprime_M800_2018UL":kRed,
              "Bprime_M1000_2018UL":kMagenta,
              "Bp1200":kBlue,
              "Bp1300":kCyan,
              "Bprime_M1400_2018UL":kViolet,
              "Bp1500":kPink+10,
              "Bp1600":kOrange+1,
              "Bp1700":kRed+3,
              "Bp1800":kPink+10,
              "Bp2000":kMagenta+3,
              "Bprime_M2200_2018UL":kGreen
          }
colors_bkg = {"QCD":40,
              "WJets":41,
              "TTToSemiLeptonic":42
          }

def plot1D_bkgsig(transform, region):
    Legend = TLegend(0.6, 0.7, 0.9, 0.9)
    histo_stack_sig = THStack(transform, transform)
    histo_stack_bkg = THStack(transform, transform)
    if(transform=="Bprime_mass"):
        histo_sum_bkg = TH1D("","",100,0,6000)
    if(transform=="gcJet_ST"):
        histo_sum_bkg = TH1D("","",100,0,5000)

    xname = transform+branches[transform][-1]
    yname = "Events"
    c1 = TCanvas("c", "c", 700,600)

    gStyle.SetOptStat(0)
    if(yLog):
        gPad.SetLogy()
        outname = outdir + subdir + transform + "_" + region + "_logY.png"
    else:
        outname = outdir + subdir + transform + "_" + region + ".png"

    histfile = TFile.Open(histfile_name, "READ")

    for sig in sigList:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, sig, case, region))
        #histo.Draw()
        #time.sleep(10)
        histo.SetLineColor(colors_sig[sig])
        Legend.AddEntry(histo, sig, 'l')
        histo_stack_sig.Add(histo)

    for bkg in bkgDir:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, region))
        histo_sum_bkg.Add(histo)

    histo_stack_sig.Draw("HIST NOSTACK")
    time.sleep(3)
    rightmax = 1.1*histo_sum_bkg.GetMaximum()
    scale = gPad.GetUymax()/rightmax

    for bkg in bkgDir:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, region))
        histo.Scale(scale)
        #histo.Draw("HIST")                                                                                                                              
        #time.sleep(10)
        histo.SetLineColor(colors_bkg[bkg])
        histo.SetFillColor(colors_bkg[bkg])
        Legend.AddEntry(histo, bkg, 'f')
        histo_stack_bkg.Add(histo)

    histo_stack_bkg.Draw("HIST")
    histo_stack_sig.Draw("HIST NOSTACK SAME")
    
    axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 510, "+L")
    axis.Draw()
    Legend.Draw()
    histo_stack_bkg.GetXaxis().SetTitle(xname)
    histo_stack_bkg.GetYaxis().SetTitle(yname)
    histo_stack_bkg.GetXaxis().SetRangeUser(branches[transform][2], branches[transform][3])
    c1.Modified()
    c1.SaveAs(outname)
    c1.Close()

if bkgOnly == True:

else:
    for region in Regions:
        plot1D_bkgsig(transform1, region)
        plot1D_bkgsig(transform2, region)
