import os, sys
import numpy as np
import matplotlib.pyplot as plt
#from ROOT import gStyle, TH1F, TFile, TTree, TCanvas, THStack, RDataFrame
from ROOT import *

###########################
#   Weight Histograms     #
###########################

indir = "root://cmseos.fnal.gov//store/user/xshen/BtoTW_Jul2023/LeptonChecks/QCDBp_record/"
outdir = os.getcwd()+'/histos_MET/'
if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

samples = {"Bprime800":"Bprime800.root",
           "Bprime1400":"Bprime1400.root", 
           "Bprime2000":"Bprime2000.root",
           "QCD200":"QCD200.root",
           "QCD300":"QCD300.root",
           "QCD500":"QCD500.root",
           "QCD700":"QCD700.root",
           "QCD1500":"QCD1500.root",
           "QCD1000":"QCD1000.root",
           "QCD2000":"QCD2000.root",
}

lumi = 138000.0
nRun = {'Bprime800':99800.,
        'Bprime1400':99600.,
        'Bprime2000':99600.,  # fixme, need a sum-of-gen-weights calculated...                   
        'QCD200':61542214.,
        'QCD300':56214199.,
        'QCD500':61097673.,
        'QCD700':47314826.,
        'QCD1000':15230975.,
        'QCD1500':11887406.,
        'QCD2000':5710430.
}
xsec = {'Bprime800':1.0,'Bprime1400':1.0,'Bprime2000':1.0, #signals all the same at 1pb for now, predictions vary                                                                                          
        'QCD200':1712000,
        'QCD300':347700,
        'QCD500':32100,
        'QCD700':6831,
        'QCD1000':1207,
        'QCD1500':119.9,
        'QCD2000':25.24
}


tfile = TFile.Open("histos_MET.root","RECREATE")

for sample in samples:
    print(sample)
    filename = indir+samples[sample]
    infile = TFile.Open(filename, "READ")
    ftree = infile.Get("NoLepIso_Events")
    ftree.GetEntry(0)
    Generator_weight = ftree.Generator_weight

    weight = Generator_weight*lumi*xsec[sample]/(nRun[sample]*abs(Generator_weight))

    MET_pt = infile.Get("MET_pt")
    MET_pt.Scale(weight)

    tfile.cd()
    MET_pt.Write(sample)

#######################
#   Make MET plots    #
#######################

#histoFile = TFile.Open(outfile, "READ")

METbkg = THStack("MET_bkg", "stacked bkg")
SUMbkg = TH1D("","",300,0,1500)

c1 = TCanvas("c1", "c1", 700,600)

Legend = TLegend(600,500,600,500)

Bprime800 = tfile.Get('Bprime800')
Bprime1400 = tfile.Get('Bprime1400')
Bprime2000 = tfile.Get('Bprime2000')
QCD200 = tfile.Get('QCD200')
QCD300 = tfile.Get('QCD300')
QCD500 = tfile.Get('QCD500')
QCD700 = tfile.Get('QCD700')
QCD1000 = tfile.Get('QCD1000')
QCD1500 = tfile.Get('QCD1500')
QCD2000 = tfile.Get('QCD2000')

Bprime800.SetLineColor(kGreen+1)                         
Bprime1400.SetLineColor(kRed+2)
Bprime2000.SetLineColor(kOrange+2)

QCD200.SetLineColor(41)
QCD300.SetLineColor(42)
QCD500.SetLineColor(43)
QCD700.SetLineColor(44)
QCD1000.SetLineColor(45)
QCD1500.SetLineColor(46)
QCD2000.SetLineColor(47)

QCD200.SetFillColor(41)
QCD300.SetFillColor(42)
QCD500.SetFillColor(43)
QCD700.SetFillColor(44)
QCD1000.SetFillColor(45)
QCD1500.SetFillColor(46)
QCD2000.SetFillColor(47)

Bprime800.SetLineWidth(2)
Bprime1400.SetLineWidth(2)
Bprime2000.SetLineWidth(2)

Bprime800.Draw("HIST")
Bprime1400.Draw("HIST SAME")
Bprime2000.Draw("HIST SAME")

c1.Update()

SUMbkg.Add(QCD200)
SUMbkg.Add(QCD300)
SUMbkg.Add(QCD500)
SUMbkg.Add(QCD700)
SUMbkg.Add(QCD1000)
SUMbkg.Add(QCD1500)
SUMbkg.Add(QCD2000)

rightmax = 1.1*SUMbkg.GetMaximum()
scale = gPad.GetUymax()/rightmax

QCD200.Scale(scale)
QCD300.Scale(scale)
QCD500.Scale(scale)
QCD700.Scale(scale)
QCD1000.Scale(scale)
QCD1500.Scale(scale)
QCD2000.Scale(scale)

METbkg.Add(QCD200)
METbkg.Add(QCD300)
METbkg.Add(QCD500)
METbkg.Add(QCD700)
METbkg.Add(QCD1000)
METbkg.Add(QCD1500)
METbkg.Add(QCD2000)

METbkg.Draw("HIST")

Bprime800.Draw("HIST SAME")
Bprime1400.Draw("HIST SAME")
Bprime2000.Draw("HIST SAME")

axis = TGaxis(gPad.GetUxmax(),gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0, rightmax, 510, "+L")

axis.Draw()

Legend.AddEntry(Bprime800, 'Bprime800', 'l')
Legend.AddEntry(Bprime1400, 'Bprime1400', 'l')
Legend.AddEntry(Bprime2000, 'Bprime2000', 'l')
Legend.AddEntry(QCD200, 'QCD200', 'f')
Legend.AddEntry(QCD300, 'QCD300', 'f')
Legend.AddEntry(QCD500, 'QCD500', 'f')
Legend.AddEntry(QCD700, 'QCD700', 'f')
Legend.AddEntry(QCD1000, 'QCD1000', 'f')
Legend.AddEntry(QCD1500, 'QCD1500', 'f')
Legend.AddEntry(QCD2000, 'QCD2000', 'f')

Legend.Draw()

c1.Update()
c1.SaveAs("histos_MET.png")

############################
#     Efficiency plot      #
############################

Nbkg = np.zeros(Ncuts)
cutList = np.linspace(50,550,Ncuts)

tfile.cd()

for sample in samples:
    MET_pt = tfile.Get(sample)
    if('Bprime' in sample):
        for i in range(Ncuts):
            Nsigs[sample][i] = MET_pt.Integral(MET_pt.FindFixBin(cutList[i]), MET_pt.FindFixBin(1500))
    else:
        for i in range(Ncuts):
            Nbkg[i] += MET_pt.Integral(MET_pt.FindFixBin(cutList[i]), MET_pt.FindFixBin(1500))

plt.figure()

for sample in Nsigs:
    plt.plot(cutList[:250], Nsigs[sample][:250]/(Nbkg[:250]**0.5), 'o', markersize=1, label=sample) # check if they are in diff colors                                                                                  

plt.legend()
plt.show()

