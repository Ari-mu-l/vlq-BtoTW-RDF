import os, time
import numpy as np
from math import sqrt
from ROOT import *
from samples import *
from utils import *
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

cutBranch = "W_MT"
nbin = 100
bin_min = 0
bin_max = 2000
getHisto = False
bkgAll = True
lumi = 138000.0

step1dir = "root://cmseos.fnal.gov//store/user/jmanagan/BtoTW_Sep2023_2018/"

outdir = os.getcwd()+'/optimizationScan/W_MT/'
if not os.path.exists(outdir): os.system('mkdir -p ' + outdir)

histfile_name =  "optimizationScan_{}.root".format(cutBranch)

# Define functions
gInterpreter.Declare("""     
float weights( float genWeight, float lumi, float xsec, float nRun ){
return genWeight * lumi * xsec / (nRun * abs(genWeight));
}
""")

def CreateFromSamples(sample):
    prefix = sample.prefix
    print("Processing {}".format(prefix))
    # read files
    samplename = sample.samplename.split('/')[1]
    tfiles = readTreeNominal(samplename,step1dir,"Events")
    # load RDF
    Events = RDataFrame(tfiles).Define("weights","weights(genWeight,{},{},{})".format(lumi,sample.xsec,sample.nrun))
    histfile.cd()
    histo = Events.Histo1D(("", "", nbin, bin_min, bin_max), cutBranch, "weights")
    histo.Write("{}_{}".format(cutBranch, prefix))
    
def AddHistos(process):
    sampleList = bkgDir[process]
    for sample in sampleList:
        histo0 = histfile.Get("{}_{}".format(cutBranch, sample))
        for i in range(1, len(sampleList)):
            histo =  histfile.Get("{}_{}".format(cutBranch, sample))
            histo0.Add(histo)
        
        histo0.Write("{}_{}".format(cutBranch, process))

sigList = ["Bprime_M800_2018UL", "Bprime_M1400_2018UL", "Bprime_M2200_2018UL"]

bkgDir = {"QCD":["QCDHT3002018UL", "QCDHT5002018UL", "QCDHT7002018UL", "QCDHT10002018UL", "QCDHT15002018UL", "QCDHT20002018UL"], 
          "TTToSemiLeptonic":["TTToSemiLeptonic2018UL"], 
          "WJets": ["WJetsHT2002018UL", "WJetsHT4002018UL", "WJetsHT6002018UL", "WJetsHT8002018UL", "WJetsHT12002018UL", "WJetsHT25002018UL"],
      }

if(getHisto):
    print("Preparing {} ...".format(histfile_name))
    start_time1 = time.time()

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

    for bkg in bkgDir:
        AddHistos("QCD")
        AddHistos("TTToSemiLeptonic")
        AddHistos("WJets")

    histfile.Close()
    end_time1 = time.time()
    print("time elapsed: ", end_time1 - start_time1)

# Plot
print("plotting...")
start_time2 = time.time()

histfile = TFile.Open(histfile_name, "READ")

Nbkg_before = np.ones(nbin)
Nsig_before = np.ones((3, nbin))
Nbkg_after = np.zeros(nbin)
Nsig_after = np.zeros((3, nbin))

bin_del = int((bin_max-bin_min)/nbin)
cutList = np.arange(bin_min, bin_max, bin_del) + bin_del

histfile.cd()

for i in range(3):
    hist_sig = histfile.Get("{}_{}".format(cutBranch, sigList[i]))
    Nsig_before[i,:]*=(hist_sig.Integral())
    for j in range(nbin):
        Nsig_after[i][j] = hist_sig.Integral(hist_sig.FindFixBin(bin_min+1), hist_sig.FindFixBin(cutList[j]-1))

hist_bkg = histfile.Get("{}_{}".format(cutBranch, "QCD"))
if(bkgAll):
    hist_bkg.Add(histfile.Get("{}_{}".format(cutBranch, "TTToSemiLeptonic")))
    hist_bkg.Add(histfile.Get("{}_{}".format(cutBranch, "WJets")))

for j in range(nbin):
    Nbkg_after[j] = hist_bkg.Integral(hist_bkg.FindFixBin(bin_min+1), hist_bkg.FindFixBin(cutList[j]-1))

Nbkg_before*=(hist_bkg.Integral())

sig800  = Nsig_after[0]/np.sqrt(Nbkg_after)
sig1400 = Nsig_after[1]/np.sqrt(Nbkg_after)
sig2200 = Nsig_after[2]/np.sqrt(Nbkg_after)

eff800  = Nsig_after[0]/Nsig_before[0]
eff1400 = Nsig_after[1]/Nsig_before[1]
eff2200 = Nsig_after[2]/Nsig_before[2]

fig, (ax1, ax2) = plt.subplots(1,2)
fig.set_size_inches(12, 6)
if(bkgAll):
    fig.suptitle('Optimization scan against major backgrounds'.format(cutBranch))
else:
    fig.suptitle('Optimization scan against QCD background'.format(cutBranch))

ax1.plot(cutList, sig800,  label = 'Bp800')
ax1.plot(cutList, sig1400, label = 'Bp1400')
ax1.plot(cutList, sig2200, label = 'Bp2200')
ax2.plot(cutList, eff800,  label = 'Bp800')
ax2.plot(cutList, eff1400, label = 'Bp1400')
ax2.plot(cutList, eff2200, label = 'Bp2200')

ax1.set_xlabel('W_MT Cut [GeV]')
ax2.set_xlabel('W_MT Cut [GeV]')

ax1.set_ylabel('S/$\sqrt{B}$')
ax2.set_ylabel('Signal Efficiency')

ax1.set(xlim=(bin_min, bin_max), ylim=(0,40))
ax2.set(xlim=(bin_min, bin_max), ylim=(0,1.2))

ax1.xaxis.set_major_locator(MultipleLocator(200))
ax2.xaxis.set_major_locator(MultipleLocator(200))

ax1.grid(which='major', alpha=0.5)
ax2.grid(which='major', alpha=0.5)

#ax1.plot(np.ones(10)*120, np.linspace(0, max(np.max(sig800), np.max(sig1400), np.max(sig2200)), 10), color='r')
#ax2.plot(np.ones(10)*120, np.linspace(0, 1, 10), color='r')

ax1.legend()
ax2.legend()
plt.show()

selection_min = 100/bin_del - 1
selection_max = 400/bin_del - 1

print("In the interested region: ")
for i in range(selection_max-selection_min):
    print("W_MT<{} GeV".format(cutList[selection_min+i]))
    print("Bp800: Significance = {}, Signal Efficiency = {}".format(sig800[selection_min+i], eff800[selection_min+i]))
    print("Bp1400: Significance = {}, Signal Efficiency = {}".format(sig1400[selection_min+i], eff1400[selection_min+i]))
    print("Bp2200: Significance = {}, Signal Efficiency = {}".format(sig2200[selection_min+i], eff2200[selection_min+i]))
    print(" ")
