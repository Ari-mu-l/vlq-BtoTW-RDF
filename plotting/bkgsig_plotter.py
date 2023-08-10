import os, sys
from numpy import linspace
from array import array
from ROOT import *

getHistos = True
#getHistos = False

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Jul2023/Run2018_Jul2023/"
outdir = os.getcwd()+'/plots_bkgsig/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

samples = {'Bp800':'Bprime_M800__20UL18_hadd.root',
           'Bp1000':'Bprime_M1000_20UL18_hadd.root',
           'Bp1200':'Bprime_M1200_20UL18_hadd.root',
           'Bp1300':'Bprime_M1300_20UL18_hadd.root',
           'Bp1400':'Bprime_M1400_20UL18_hadd.root',
           'Bp1500':'Bprime_M1500_20UL18_hadd.root',
           'Bp1600':'Bprime_M1600_20UL18_hadd.root',
           'Bp1700':'Bprime_M1700_20UL18_hadd.root',
           'Bp1800':'Bprime_M1800_20UL18_hadd.root',
           'Bp2000':'Bprime_M2000_20UL18_hadd.root',
           'Bp2200':'Bprime_M2200_20UL18_hadd.root',
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

#################
# Get histogram #
#################
  
###### Add branches here ####                                                                                                                    

branches = {"M_reco": ["Bprime_mass", 50, 0, 4000, "[GeV]"],
            "pt_reco": ["Bprime_pt", 50, 0, 4000, "[GeV]"],
            #"Number of opposite-side jets": ["NOS_gcJets_central", 15, 0, 15, ""],
            #"Number of same-side jets": ["NSS_gcJets_central", 15, 0, 15, ""],
            #"Number of opposite-side b-tagged jets": ["NOS_gcJets_DeepFlavL", 15, 0, 15, ""],
            #"Number of same-side b-tagged jets": ["NSS_gcJets_DeepFlavL", 15, 0, 15, ""],
}

categories = {#"DY": ["DYMHT1200", "DYMHT200", "DYMHT2500"],       
              "QCD": ["QCD1000", "QCD1500", "QCD2000", #"QCD200", 
                      "QCD300", "QCD500", "QCD700"],
              "WJets":["WJets1200", "WJets200", "WJets2500", "WJets400", "WJets600", "WJets800"],
              #"TTbar":["TTToSemiLeptonic"],
          }

if(getHistos):

    gInterpreter.Declare("""
    float weights( float Generator_weight, float lumi, float xsec, float nRun ){
    return Generator_weight * lumi * xsec / (nRun * abs(Generator_weight));
    }
    """)

    print("Preparing bkgsig_histos.root...")
    histfile = TFile.Open("bkgsig_histos.root", "RECREATE")

    for sample in samples:
        filename = indir + samples[sample]
        Events = RDataFrame("Events", filename).Define("Generator_weight","(float) 1.0").Define("weights","weights(Generator_weight,{},{},{})".format(lumi,xsec[sample],nRun[sample]))

        for branch in branches:
            # set histogram bins
            nbins = branches[branch][1]
            bin_lo = branches[branch][2]
            bin_hi = branches[branch][3]
            
            histo = Events.Histo1D((branch, branch, nbins, bin_lo, bin_hi), branches[branch][0], "weights")
            
            histfile.cd()
            histo.Write(branch+"_"+sample+"_weighted")


    for category in categories:
        sampleList = categories[category]
        if(len(sampleList)>1):
            for branch in branches:
                histo1 = histfile.Get(branch + "_" + sampleList[0]+"_weighted")

                for i in range(1,len(sampleList)):
                    branch + "_" + sampleList[i]+"_weighted"
                    histo =  histfile.Get(branch + "_" + sampleList[i]+"_weighted")
                    histo1.Add(histo)
        
                histo1.Write(branch + "_" + category+"_weighted")

    histfile.Close()

########
# Plot #
########
print("plotting...")
histfile = TFile.Open("bkgsig_histos.root", "READ")

sig_list = ["Bp800", "Bp1000", "Bp1200", "Bp1300", "Bp1400", "Bp1500", "Bp1600", "Bp1700", "Bp1800", "Bp2000", "Bp2200"]
bkg_list = ["QCD", "WJets", "TTToSemiLeptonic"]

colors_sig = {"Bp800":kRed,
              "Bp1000":kMagenta,
              "Bp1200":kBlue,
              "Bp1300":kCyan,
              "Bp1400":kViolet,
              "Bp1500":kPink+10,
              "Bp1600":kOrange+1,
              "Bp1700":kRed+3,
              "Bp1800":kPink+10,
              "Bp2000":kMagenta+3,
              "Bp2200":kGreen
}

colors_bkg = {"QCD":40,
              "WJets":41,
              "TTToSemiLeptonic":42
}
#fillstyles = [3002, 3006, 3375, 3357]

for branch in branches:
    c1 = TCanvas("c", "c", 600,400)
    gStyle.SetOptStat(0)

    Legend = TLegend(0.6, 0.7, 0.9, 0.9)
    histo_stack_sig = THStack(branch, branch)
    histo_stack_bkg = THStack(branch, branch)

    for sig in sig_list:
        histo = histfile.Get(branch + "_" + sig + "_weighted")
        
        histo.SetLineColor(colors_sig[sig])
        #histo.SetFillColor(colors_sig[sample])
        #histo.SetFillStyle(3002)                                                                                                        
        Legend.AddEntry(histo, sig, 'l')
        histo_stack_sig.Add(histo)

    for bkg in bkg_list:
        histo = histfile.Get(branch + "_" + bkg + "_weighted")
        
        histo.SetLineColor(colors_bkg[bkg])
        histo.SetFillColor(colors_bkg[bkg])
        #histo.SetFillStyle(3002)

        Legend.AddEntry(histo, bkg, 'f')
        histo_stack_bkg.Add(histo)
            
    xname = branch+branches[branch][-1]
    yname = "Normalized frequency (Unweighted)"
        
    #histo_stack_sig.Draw("HIST NOSTACK SAME")
    histo_stack_bkg.Draw("HIST")
    histo_stack_sig.Draw("HIST NOSTACK SAME")
    Legend.Draw()
    c1.Update()
        
    histo_stack_sig.GetXaxis().SetTitle(xname)
    histo_stack_sig.GetYaxis().SetTitle(yname)
    histo_stack_sig.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])
        
    c1.Modified()
    c1.SaveAs(outdir + "histos_" + branch + ".png")
