#Arguements: iPlot, region, isCategorized, directory, blind, yLog, rebinning

#python -u plotTemplatesPaper.py DnnTprime SR True templatesSRCR_Nov2021_TT False True 0p2
#python -u plotTemplatesPaper.py HTdnnL CR True templatesCR_Nov2021_TT False False 0p2

#plotList='lepPt lepEta lepPhi lepIso JetEta JetPt NJetsCentral NJetsForward NBJets MET HT ST FatJetEta FatJetPt Tau21 SoftDrop probj probt probw deeptag NFatJets nT nW minDR_twoAK8s tmass tpt Wdrlep tdrWb isLepW minMlj PtRel PtRelAK8 minDR minDRAK8 BpMass BpPt BpEta BpPhi BpDeltaR BpPtBal BpChi2 MLP_HT500_Bprime MLP_HT500_TTbar MLP_HT500_WJets MLP_HT250_Bprime MLP_HT250_TTbar MLP_HT250_WJets'
plotList='mlp_bprime mlp_ttbar mlp_wjets'
for iPlot in $plotList; do
    echo $iPlot
    python -u -b histoDrawTrees.py $iPlot
done

