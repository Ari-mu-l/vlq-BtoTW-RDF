import os, sys, time
from numpy import linspace
from array import array
from ROOT import *

###############
#   Options   #
###############
getHistos = True
yLog = False
withTags = True
noTag = False
plotSeparate = False

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Aug2023_2018/"
outdir = os.getcwd()+'/plots_bkgsig/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

##################
# Weighting info #
##################
samples = {'Bp800':'Bprime_M800_20UL18_hadd.root',
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

###### Add branches, categories, tags here ####

branches = {"M_reco": ["Bprime_mass", 50, 0, 4000, "[GeV]"],
            #"pt_reco": ["Bprime_pt", 50, 0, 4000, "[GeV]"],
            #"Number of opposite-side jets": ["NOS_gcJets_central", 15, 0, 15, ""],
            #"Number of same-side jets": ["NSS_gcJets_central", 15, 0, 15, ""],
            #"Number of opposite-side b-tagged jets": ["NOS_gcJets_DeepFlavL", 15, 0, 15, ""],
            #"Number of same-side b-tagged jets": ["NSS_gcJets_DeepFlavL", 15, 0, 15, ""],
            #"Number of opposite-side fatjets": ["NOS_gcFatJets", 15, 0, 15, ""],
}

categories = {#"DY": ["DYMHT1200", "DYMHT200", "DYMHT2500"],       
              "QCD": ["QCD1000", "QCD1500", "QCD2000", #"QCD200", 
                      "QCD300", "QCD500", "QCD700"],
              "WJets":["WJets1200", "WJets200", "WJets2500", "WJets400", "WJets600", "WJets800"],
              #"TTbar":["TTToSemiLeptonic"],
          }

tags_general = {"_BdecayCase1":"Bdecay_obs==1",
                "_BdecayCase2":"Bdecay_obs==2",
                "_BdecayCase3":"Bdecay_obs==3",
                "_BdecayCase4":"Bdecay_obs==4",
                "_BdecayAll":"Bdecay_obs>=1",
                "_BdecayOther":"Bdecay_obs<1",
                "_tTag":"NSS_gcJets_DeepFlavL>0",
                "_WTag":"NSS_gcJets_DeepFlavL==0"
            }

tags_signal = {"_trueLepT":"trueLeptonicT==1 && trueLeptonicW==0 && leptonicParticle==1",
               "_trueLepW":"trueLeptonicT==0 && trueLeptonicW==1 && leptonicParticle==0",
               "_falseLepT":"trueLeptonicT==0 && trueLeptonicW==1 && leptonicParticle==1",
               "_falseLepW":"trueLeptonicT==1 && trueLeptonicW==0 && leptonicParticle==0"
           }

####################
# Define functions #
####################
def CreateHistos(Events, tags, branches, sample):
    for tag in tags:
        Events_tag = Events.Filter(tags[tag])
        for branch in branches:
            # set histogram bins
            nbins = branches[branch][1]
            bin_lo = branches[branch][2]
            bin_hi = branches[branch][3]

            histo_tag = Events_tag.Histo1D((branch, branch, nbins, bin_lo, bin_hi), branches[branch][0], "weights")

            histfile.cd()
            histo_tag.Write(branch + "_" + sample + "_weighted" + tag)
    
def AddHistos(branches, sampleList, tags):
    for tag in tags:
        histo1 = histfile.Get(branch + "_" + sample + "_weighted" + tag)
        for i in range(1, len(sampleList)):
            histo =  histfile.Get(branch + "_" + sampleList[i] + "_weighted" + tag)
            histo1.Add(histo)

        histo1.Write(branch + "_" + category+"_weighted"+tag)

##################
# Get Histograms #
##################
if(getHistos):
    start_time1 = time.time()

    gInterpreter.Declare("""
    float weights( float genWeight, float lumi, float xsec, float nRun ){
    return genWeight * lumi * xsec / (nRun * abs(genWeight));
    }
    """)

    print("Preparing bkgsig_histos.root...")
    histfile = TFile.Open("bkgsig_histos.root", "RECREATE")

    for sample in samples:
        filename = indir + samples[sample]
        Events = RDataFrame("Events", filename).Filter("NJets_forward>0 && Bprime_mass>0").Define("weights","weights(genWeight,{},{},{})".format(lumi,xsec[sample],nRun[sample]))

        if(withTags):
            CreateHistos(Events, tags_general, branches, sample)
        
            if("Bp" in sample):
                CreateHistos(Events, tags_signal, branches, sample)

        if(noTag):
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
                if(withTags):
                    AddHistos(branches, sampleList, tags_general)
                if(noTag):
                    AddHistos(branches, sampleList, {"":""})

    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", start_time1 - end_time1)

########
# Plot #
########
print("plotting...")
start_time2 = time.time()
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

def plot(bkgTag, sigTag):
    for branch in branches:
        Legend = TLegend(0.6, 0.7, 0.9, 0.9)
        histo_stack_sig = THStack(branch, branch)
        histo_stack_bkg = THStack(branch, branch)

        for sig in sig_list:
            histo = histfile.Get(branch + "_" + sig + "_weighted" + sigTag) 
        
            histo.SetLineColor(colors_sig[sig])
                 
            Legend.AddEntry(histo, sig, 'l')
            histo_stack_sig.Add(histo)

        for bkg in bkg_list:
            histo = histfile.Get(branch + "_" + bkg + "_weighted" + bkgTag)
        
            histo.SetLineColor(colors_bkg[bkg])
            histo.SetFillColor(colors_bkg[bkg])

            Legend.AddEntry(histo, bkg, 'f')
            histo_stack_bkg.Add(histo)

        xname = branch+branches[branch][-1]
        yname = "Count"

        if(plotSeparate):
            c_bkg = TCanvas("c", "c", 600,400)
            gStyle.SetOptStat(0)

            if(yLog):
                gPad.SetLogy()

            histo_stack_bkg.Draw("HIST")
            Legend.Draw()
            
            histo_stack_bkg.GetXaxis().SetTitle(xname)
            histo_stack_bkg.GetYaxis().SetTitle(yname)
            histo_stack_bkg.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])

            c_bkg.Modified()

            if(yLog):
                outname = outdir + "histos_" + branch + bkgTag + "bkg_logY.png"
            else:
                outname = outdir + "histos_" + branch + bkgTag + "bkg.png"

            c_bkg.SaveAs(outname)
            c_bkg.Close()

            c_sig = TCanvas("c", "c", 600,400)
            gStyle.SetOptStat(0)

            if(yLog):
                gPad.SetLogy()

            histo_stack_sig.Draw("HIST NOSTACK")
            Legend.Draw()

            histo_stack_sig.GetXaxis().SetTitle(xname)
            histo_stack_sig.GetYaxis().SetTitle(yname)
            histo_stack_sig.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])

            c_sig.Modified()

            if(yLog):
                outname = outdir + "histos_" + branch + sigTag + "sig_logY.png"
            else:
                outname = outdir + "histos_" + branch + sigTag + "sig.png"
            c_sig.SaveAs(outname)

        else:
            c1 = TCanvas("c", "c", 600,400)
            gStyle.SetOptStat(0)

            if(yLog):
                gPad.SetLogy()
    
            histo_stack_bkg.Draw("HIST")
            histo_stack_sig.Draw("HIST NOSTACK SAME")
            Legend.Draw()
            c1.Update()
        
            histo_stack_bkg.GetXaxis().SetTitle(xname)
            histo_stack_bkg.GetYaxis().SetTitle(yname)
            histo_stack_bkg.GetXaxis().SetRangeUser(branches[branch][2], branches[branch][3])
        
            c1.Modified()

            if(yLog):
                outname = outdir + "histos_" + branch + sigTag + "_logY.png"
            else:
                outname = outdir + "histos_" + branch + sigTag + ".png"
            c1.SaveAs(outname)

if(withTags):
    for tag in tags_general:
        plot(tag, tag)

    if(~plotSeparate):
        plot("_tTag", "_trueLepT")
        plot("_tTag", "_falseLepT")
        plot("_WTag", "_trueLepW")
        plot("_WTag", "_falseLepW")

if(noTag):
    plot("", "")

end_time2 = time.time()
print("time elapsed: ", start_time2 - end_time2)
