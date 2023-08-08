import os, sys
from numpy import linspace
from array import array
from ROOT import *

getHistos = True
#getHistos = False

indir = "root://cmseos.fnal.gov//store/user/kjohnso/BtoTW_Jul2023/Run2018_Jul2023/"
outdir = os.getcwd()+'/plots_signal/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

samples = {'Bp800':'Bprime_M800__20UL18/RDF_BprimeBtoTW_M-800_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1000':'Bprime_M1000_20UL18/RDF_BprimeBtoTW_M-1000_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1200':'Bprime_M1200_20UL18/RDF_BprimeBtoTW_M-1200_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1300':'Bprime_M1300_20UL18/RDF_BprimeBtoTW_M-1300_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1400':'Bprime_M1400_20UL18/RDF_BprimeBtoTW_M-1400_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1500':'Bprime_M1500_20UL18/RDF_BprimeBtoTW_M-1500_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1600':'Bprime_M1600_20UL18/RDF_BprimeBtoTW_M-1600_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1700':'Bprime_M1700_20UL18/RDF_BprimeBtoTW_M-1700_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp1800':'Bprime_M1800_20UL18/RDF_BprimeBtoTW_M-1800_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp2000':'Bprime_M2000_20UL18/RDF_BprimeBtoTW_M-2000_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root',
           'Bp2200':'Bprime_M2200_20UL18/RDF_BprimeBtoTW_M-2200_NWALO_TuneCP5_13TeV-madgraph-pythia8_finalsel_0.root'
}

#################
# Get histogram #
#################

###### Add branches here ####

branches = {"M_reco": ["Bprime_mass", 50, 0, 4000, "[GeV]"],
            "Number of opposite-side jets": ["NOS_gcJets_central", 10, 0, 10, ""],
            "Number of same-side jets": ["NSS_gcJets_central", 7, 0, 7, ""],
            "Number of opposite-side b-tagged jets": ["NOS_gcJets_DeepFlavL", 7, 0, 7, ""],
            "Number of same-side b-tagged jets": ["NSS_gcJets_DeepFlavL", 3, 0, 3, ""],
}

if(getHistos):
    histfile = TFile.Open("signal_histos.root", "RECREATE")

    for sample in samples:
        filename = indir+samples[sample]
        Events = RDataFrame("Events", filename)

        Events_t = Events.Filter("trueLeptonicT==1 && trueLeptonicW==0")
        Events_W = Events.Filter("trueLeptonicT==0 && trueLeptonicW==1")
        Events_noLep = Events.Filter("trueLeptonicT==0 && trueLeptonicW==0")
        Events_diLep = Events.Filter("trueLeptonicT==1 && trueLeptonicW==1")

        for branch in branches:
            nbins = branches[branch][1]
            bin_lo = branches[branch][2]
            bin_hi = branches[branch][3]
            
            histo = Events.Histo1D((branch,             branch,              nbins, bin_lo, bin_hi), branches[branch][0])
            histo_t = Events_t.Histo1D((branch,         branch+"_LeptonicT", nbins, bin_lo, bin_lo), branches[branch][0])
            histo_W = Events_W.Histo1D((branch,         branch+"_LeptonicW", nbins, bin_lo, bin_lo), branches[branch][0])
            histo_noLep = Events_noLep.Histo1D((branch, branch+"_noLep",     nbins, bin_lo, bin_lo), branches[branch][0])
            histo_diLep = Events_diLep.Histo1D((branch, branch+"_diLep",     nbins, bin_lo, bin_hi), branches[branch][0])
            
            histfile.cd()
            
            histo.Write(branch+"_"+sample)
            histo_t.Write(branch+"_LeptonicT_"+sample)
            histo_W.Write(branch+"_LeptonicW_"+sample)
            histo_noLep.Write(branch+"_noLep_"+sample)
            histo_diLep.Write(branch+"_diLep_"+sample)
            
    histfile.Close()

########
# Plot #
########
histfile = TFile.Open("signal_histos.root", "READ")

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

tags = ["", "_LeptonicT", "_LeptonicW", "_noLep", "_diLep"]

isFirst = True
for branch in branches:
    for tag in tags:
        count, nplot = 0, 0
        for sample in samples:
            if(count==0):
                 c1 = TCanvas("c"+tag, "c"+tag, 600,400)
                 gStyle.SetOptStat(0)

                 Legend = TLegend(0.6, 0.7, 0.9, 0.9)
                 histo_stack = THStack(branch+tag, branch+tag)

            if(len(tag)!=0):
                histo = histfile.Get(branch + tag + "_" + sample)
            else:
                histo = histfile.Get(branch + "_" + sample)

            histo.Scale(1/histo.Integral())

            histo.SetLineColor(colors[count])
            histo.SetFillColor(colors[count])
            histo.SetFillStyle(3002)
            Legend.AddEntry(histo, sample, 'f')

            histo_stack.Add(histo)
            count+=1
            
            if(count==4 or (nplot==3 and count==3)):
                xname = branch+branches[branch][-1]
                yname = "Normalized frequency (Unweighted)"

                histo_stack.Draw("HIST NOSTACK")
                Legend.Draw()
                c1.Update()

                histo_stack.GetXaxis().SetTitle(xname)
                histo_stack.GetYaxis().SetTitle(yname)
                histo_stack.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])

                c1.Modified()
                c1.SaveAs(outdir + "histos_" + branch + tag + "_" + str(nplot) + ".png")
                nplot+=1
                count=0
