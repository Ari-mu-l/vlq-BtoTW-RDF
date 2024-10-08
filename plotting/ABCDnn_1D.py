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
case = "_BdecayCase1and4"
getHistos = True
addNew = True
withSig = False
withData = True
stacked = False
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

def AddHistos(sampleDir, region):
    for sample in sampleDir:
        sampleList = sampleDir[sample]
        histo1 = histfile.Get("{}_{}_weighted{}_{}".format(transform1, sampleList[0], case, region))
        histo2 = histfile.Get("{}_{}_weighted{}_{}".format(transform2, sampleList[0], case, region))

        for i in range(1, len(bkgList)):
            histo1_i =  histfile.Get("{}_{}_weighted{}_{}".format(transform1, sampleList[i], case, region))
            histo2_i =  histfile.Get("{}_{}_weighted{}_{}".format(transform2, sampleList[i], case, region))
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
dataDir = {"data": ["SingleMuon", "EGamma"]}

if(getHistos):
    start_time1 = time.time()

    print("Preparing {} ...".format(histfile_name))
    if(addNew):
        histfile = TFile.Open(histfile_name, "UPDATE")
    else:
        histfile = TFile.Open(histfile_name, "RECREATE")

    #CreateFromSamples(Bprime_M800_2018UL)
    #CreateFromSamples(Bprime_M1400_2018UL)
    #CreateFromSamples(Bprime_M2200_2018UL)
    #CreateFromSamples(QCDHT3002018UL)
    #CreateFromSamples(QCDHT5002018UL)
    #CreateFromSamples(QCDHT7002018UL)
    #CreateFromSamples(QCDHT10002018UL)
    #CreateFromSamples(QCDHT15002018UL)
    #CreateFromSamples(QCDHT20002018UL)
    #CreateFromSamples(TTToSemiLeptonic2018UL)
    #CreateFromSamples(WJetsHT2002018UL)
    #CreateFromSamples(WJetsHT4002018UL)
    #CreateFromSamples(WJetsHT6002018UL)
    #CreateFromSamples(WJetsHT8002018UL)
    #CreateFromSamples(WJetsHT12002018UL)
    #CreateFromSamples(WJetsHT25002018UL)
    CreateFromSamples(SingleMuon)
    CreateFromSamples(EGamma)

    # add histos
    #for region in Regions:
        #AddHistos(bkgDir, region)

    for region in Regions:
        AddHistos(dataDir, region)
        
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

def plot1D_Stack(transform, regionList, label):
    xname = transform+branches[transform][-1]
    yname = "Events"
    c1 = TCanvas("c", "c", 700,600)
    Legend = TLegend(0.6, 0.7, 0.9, 0.9)

    gStyle.SetOptStat(0)
    if(yLog):
        gPad.SetLogy()
        outname = outdir + subdir + transform + "_" + label + "_logY.png"
    else:
        outname = outdir + subdir + transform + "_" + label + ".png"

    histfile = TFile.Open(histfile_name, "READ")
    
    histo_sum_bkg0 = TH1D("", "", branches[transform][1], branches[transform][2], branches[transform][3])
    histo_sum_bkg1 = TH1D("", "", branches[transform][1], branches[transform][2], branches[transform][3])
    histo_sum_bkg2 = TH1D("", "", branches[transform][1], branches[transform][2], branches[transform][3])
    histo_stack_bkg = THStack(transform, transform)
        
    for bkg in bkgDir:
        histo0 = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, regionList[0]))
        histo1 = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, regionList[1]))
        histo2 = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, regionList[2]))
        histo_sum_bkg0.Add(histo0)
        histo_sum_bkg1.Add(histo1)
        histo_sum_bkg2.Add(histo2)

    histo_sum_bkg0.SetLineColor(40)
    histo_sum_bkg0.SetFillColor(40)
    histo_sum_bkg1.SetLineColor(41)
    histo_sum_bkg1.SetFillColor(41)
    histo_sum_bkg2.SetLineColor(42)
    histo_sum_bkg2.SetFillColor(42)
    Legend.AddEntry(histo_sum_bkg0, regionList[0], 'f')
    Legend.AddEntry(histo_sum_bkg1, regionList[1], 'f')
    Legend.AddEntry(histo_sum_bkg2, regionList[2], 'f')
    histo_stack_bkg.Add(histo_sum_bkg0)
    histo_stack_bkg.Add(histo_sum_bkg1)
    histo_stack_bkg.Add(histo_sum_bkg2)
    
    histo_stack_bkg.Draw("HIST")
    Legend.Draw()
    histo_stack_bkg.GetXaxis().SetTitle(xname)
    histo_stack_bkg.GetYaxis().SetTitle(yname)
    histo_stack_bkg.GetXaxis().SetRangeUser(branches[transform][2], branches[transform][3])
    c1.Modified()
    c1.SaveAs(outname)
    c1.Close()

def histoScale(sigList):
    histfile = TFile.Open(histfile_name, "READ")

    for bkg in bkgDir:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, region))
        histo_sum_bkg.Add(histo)

    for sig in sigList:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, sig, case, region))
        if(withSig):
            histo.SetLineColor(colors_sig[sig])
            Legend.AddEntry(histo, sig, 'l')
        else:
            Legend.AddEntry(histo, sig, 'lep')
        histo_stack_sig.Add(histo)

    histo_stack_sig.Draw("HIST NOSTACK")
    time.sleep(3)
    rightmax = 1.1*histo_sum_bkg.GetMaximum()
    scale = gPad.GetUymax()/rightmax

    for bkg in bkgDir:
        histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, region))
        histo.Scale(scale)

        histo.SetLineColor(colors_bkg[bkg])
        histo.SetFillColor(colors_bkg[bkg])
        Legend.AddEntry(histo, bkg, 'f')
        histo_stack_bkg.Add(histo)

def plot1D_noStack(transform, region):
    Legend = TLegend(0.6, 0.7, 0.9, 0.9)
    histo_stack_sig = THStack(transform, transform)
    histo_stack_bkg = THStack(transform, transform)

    histo_sum_bkg = TH1D("", "", branches[transform][1], branches[transform][2], branches[transform][3])

    xname = transform+branches[transform][-1]
    yname = "Events"
    c1 = TCanvas("c", "c", 700,600)

    gStyle.SetOptStat(0)
    outname = "{}{}{}_{}".format(outdir, subdir, transform, region)
    if(withData):
        outname += "data"
    if(yLog):
        gPad.SetLogy()
        outname += "_logY"
    outname += ".png"

    histfile = TFile.Open(histfile_name, "READ")
    
    if(withSig):
        histoScale(sigList)
    elif(withData):
        histoScale(dataDir)
    else:
        for bkg in bkgDir:
            histo = histfile.Get("{}_{}_weighted{}_{}".format(transform, bkg, case, region))
            
            histo.SetLineColor(colors_bkg[bkg])
            histo.SetFillColor(colors_bkg[bkg])
            Legend.AddEntry(histo, bkg, 'f')
            histo_stack_bkg.Add(histo)

    histo_stack_bkg.Draw("HIST")

    if(withSig):
        histo_stack_sig.Draw("HIST NOSTACK SAME")
    
        axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 510, "+L")
        axis.Draw()
    if(withData):
        histo_stack_sig.Draw("P SAME")

        axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 510, "+L")
        axis.Draw()


    Legend.Draw()
    histo_stack_bkg.GetXaxis().SetTitle(xname)
    histo_stack_bkg.GetYaxis().SetTitle(yname)
    histo_stack_bkg.GetXaxis().SetRangeUser(branches[transform][2], branches[transform][3])
    c1.Modified()
    c1.SaveAs(outname)
    c1.Close()

if (stacked):
    plot1D_Stack(transform1, ["A", "B", "X"], "ABX")
    plot1D_Stack(transform1, ["C", "D", "Y"], "CDY")
    plot1D_Stack(transform2, ["A", "B", "X"], "ABX")
    plot1D_Stack(transform2, ["C", "D", "Y"], "CDY")
else:
    for region in Regions:
        plot1D_noStack(transform1, region)
        plot1D_noStack(transform2, region)
