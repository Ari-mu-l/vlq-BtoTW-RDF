import os, sys
from numpy import linspace
from array import array
from ROOT import *

getHistos = True
#getHistos = False

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Jul2023/Run2018_Jul2023/"
outdir = os.getcwd()+'/plots_background/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

samples = {#'DYMHT1200':'DYMHT12002018UL_hadd.root',
           #'DYMHT200':'DYMHT2002018UL_hadd.root',
           #'DYMHT2500':'DYMHT25002018UL_hadd.root',
#           'DYMHT400':'DYMHT4002018UL_hadd.root',
#           'DYMHT600':'DYMHT6002018UL_hadd.root',
#           'DYMHT800':'DYMHT8002018UL_hadd.root',
           'QCD1000':'QCDHT10002018UL_hadd.root',
           'QCD1500':'QCDHT15002018UL_hadd.root',
           'QCD2000':'QCDHT20002018UL_hadd.root',
           'QCD200':'QCDHT2002018UL_hadd.root',
           'QCD300':'QCDHT3002018UL_hadd.root',
           'QCD500':'QCDHT5002018UL_hadd.root',
           'QCD700':'QCDHT7002018UL_hadd.root',           
#           'STs':'STs2018UL_hadd.root',
           #'STt':'',
#           'STtW':'STtW2018UL_hadd.root',
#           'STtWb':'STtWb2018UL_hadd.root',
#           'STtb':'STtb2018UL_hadd.root',
#           'SingleElA':'SingleElectronRun2018A2018UL_hadd.root',
#           'SingleElB':'SingleElectronRun2018B2018UL_hadd.root',
#           'SingleElC':'SingleElectronRun2018C2018UL_hadd.root',
#           'SingleElD':'SingleElectronRun2018D2018UL_hadd.root',
#           'SingleMuA':'SingleMuonRun2018A2018UL_hadd.root',
#           'SingleMuB':'SingleMuonRun2018B2018UL_hadd.root',
#           'SingleMuC':'SingleMuonRun2018C2018UL_hadd.root',
           #'SingleMuD':'',
#           'TTHB':'TTHB2018UL_hadd.root',
#           'TTHnonB':'TTHnonB2018UL_hadd.root',
#           'TTMT1000':'TTMT10002018UL_hadd.root',
#           'TTMT700':'TTMT7002018UL_hadd.root',
           #'TTTo2L2Nu':'',
           #'TTTo2L2Nu_0':'',
           #'TTToHadronic':'',
#           'TTToSemiLeptonic':'TTToSemiLeptonic2018UL_hadd.root',
#           'TTWl':'TTWl2018UL_hadd.root',
#           'TTWq':'TTWq2018UL_hadd.root',
#           'TTZM10':'TTZM102018UL_hadd.root',
#           'TTZM1to10':'TTZM1to102018UL_hadd.root',
           'WJets1200':'WJetsHT12002018UL_hadd.root',
           'WJets200':'WJetsHT2002018UL_hadd.root',
           'WJets2500':'WJetsHT25002018UL_hadd.root',
           'WJets400':'WJetsHT4002018UL_hadd.root',
           'WJets600':'WJetsHT6002018UL_hadd.root',
           'WJets800':'WJetsHT8002018UL_hadd.root',
#           'WW':'WW2018UL_hadd.root',
#           'WZ':'WZ2018UL_hadd.root',
#           'ZZ':'ZZ2018UL_hadd.root',
}

#################
# Get histogram #
#################

###### Add branches here ####

branches = {"M_reco": ["Bprime_mass", 50, 0, 4000, "[GeV]"],
            "pt_reco": ["Bprime_pt", 50, 0, 4000, "[GeV]"],
            "Number of opposite-side jets": ["NOS_gcJets_central", 15, 0, 15, ""],
            "Number of same-side jets": ["NSS_gcJets_central", 15, 0, 15, ""],
            "Number of opposite-side b-tagged jets": ["NOS_gcJets_DeepFlavL", 15, 0, 15, ""],
            "Number of same-side b-tagged jets": ["NSS_gcJets_DeepFlavL", 15, 0, 15, ""],
}

categories = {#"DY": ["DYMHT1200", "DYMHT200", "DYMHT2500"],
              "QCD": ["QCD1000", "QCD1500", "QCD2000", "QCD200", "QCD300", "QCD500", "QCD700"],
              "WJets":["WJets1200","WJets200","WJets2500","WJets400","WJets600","WJets800"]
          }

lumi = 138000.0
nRun = {#"DY":,
        "QCD1000":15230975.,
        "QCD1500":11887406.,
        "QCD2000":5710430.,
        "QCD200":61542214.,
        "QCD300":56214199.,
        "QCD500":61097673.,
        "QCD700":47314826.,
        "WJets1200":6481518.,
        "WJets200":58225632.,
        "WJets2500":2097648.,
        "WJets400":7444030.,
        "WJets600":7718765.,
        "WJets800":7306187.,
    }

xsec = {
        'QCD300':347700,
        'QCD500':32100,
        'QCD700':6831,
        'QCD1000':1207,
        'QCD1500':119.9,
        'QCD2000':25.24,
        'Wjets200':359.7*1.21,
        'Wjets400':48.91*1.21,
        'Wjets600':12.05*1.21,
        'Wjets800':5.501*1.21,
        'Wjets1200':1.329*1.21,
        'Wjets2500':0.03216*1.21,
}

if(getHistos):
    histfile = TFile.Open("background_histos.root", "RECREATE")

    for sample in samples:
        filename = indir+samples[sample]
        Events = RDataFrame("Events", filename)
        
        for branch in branches:
            nbins = branches[branch][1]
            bin_lo = branches[branch][2]
            bin_hi = branches[branch][3]
            
            histo = Events.Histo1D((branch, branch, nbins, bin_lo, bin_hi), branches[branch][0])
            
            histfile.cd()
            
            histo.Write(branch+"_"+sample)
    
    for category in categories:
        sampleList = categories[category]
        if(len(sampleList)!=0):
            for branch in branches:
                histo1 = histfile.Get(branch + "_" + sampleList[0])
                for i in range(1,len(sampleList)):
                    histo =  histfile.Get(branch + "_" + sampleList[i])
                    histo1.Add(histo)
        
                histo1.Write(branch + "_" + category)
        
    histfile.Close()

########
# Plot #
########
histfile = TFile.Open("background_histos.root", "READ")

background_list = ["QCD", "WJets"]

#colors = {"Bp800":kOrange,
#          "Bp1000":kGreen-8,
#          "Bp1200":kCyan-8,
#          "Bp1300":kRed-10,
#          "Bp1400":kOrange-9,
#          "Bp1500":kGreen-6,
#          "Bp1600":kBlue-10,
#          "Bp1700":kRed-9,
#          "Bp1800":kMagenta-10,
#          "Bp2000":kBlue-9,
#          "Bp2200":kRed-7
#}

colors = [kOrange, kRed, kBlue-9, kGreen+2]
#fillstyles = [3002, 3006, 3375, 3357]

for branch in branches:
    c1 = TCanvas("c", "c", 600,400)
    gStyle.SetOptStat(0)

    Legend = TLegend(0.6, 0.7, 0.9, 0.9)
    histo_stack = THStack(branch, branch)

    count=0
    for bkg in background_list:
        histo = histfile.Get(branch + "_" + bkg)
        
        histo.Scale(1/histo.Integral())

        histo.SetLineColor(colors[count])
        histo.SetFillColor(colors[count])
        histo.SetFillStyle(3002)
        Legend.AddEntry(histo, bkg, 'f')

        histo_stack.Add(histo)
        count+=1
            
    xname = branch+branches[branch][-1]
    yname = "Normalized frequency (Unweighted)"
        
    histo_stack.Draw("HIST")
    Legend.Draw()
    c1.Update()
        
    histo_stack.GetXaxis().SetTitle(xname)
    histo_stack.GetYaxis().SetTitle(yname)
    histo_stack.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])
        
    c1.Modified()
    c1.SaveAs(outdir + "histos_" + branch + "_" + ".png")
