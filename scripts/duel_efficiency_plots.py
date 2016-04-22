#!/usr/bin/python

# The first line allows for somebody to execute this script using a command
# like '> ./run_analysis.py'
# without it you need to use '> python run_analysis.py'
#
# This is a simple script allowing you to run a simple analysis of multiple trees.
# It requires three arguments ... the base file name, the first run number, and the last run number
#
# For example to analyses 100 events from the MIND simulation execute
# > ./run_analysis.py "/data/neutrino05/rbayes/Rec_out/mubarCC/mubarCC_" 1000 1099

# For more information on python programming in  general refer to python.org
# For more information on using root with python see the documentation at root.cern.ch

import sys
from ROOT import TFile, TTree, TH1D, TH2D, TChain, TCanvas, gDirectory

def analysis(basename1, basename2, first, last, maxE, qsel1, qsel2):
    # this initializes a TChain object to read a tree file
    tc = [TChain('tree'), TChain('tree')]

    # here is a simple for loop
    for i in range(int(first), int(last)+1):
        # define the filename
        filename = basename1
        filename += str(i)
        filename += '.root'
        tc[0].Add(filename)
        filename = basename2
        filename += str(i)
        filename += '.root'
        tc[1].Add(filename)

    #cuts = 'TrajectoryNo > 0   \
    #        && abs(evVertex[2]) < 50000 \
    #        && Fitted[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj]==1 \
    #        && TrajVertex[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj] < 49000 \
    #        && initrangqP[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj]/RecMom[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj] > 0.25 \
    #        && fittedNode[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj]/InMu[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj] > 0.8 \
    #        && abs(ErrMom[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj]/RecMom[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj]) < 20.0'
    cuts = 'classif_NumTrajectory > 0 '
    
    # Cuts that do not make sens for a mini-mind analysis
    # 
    # && InMu[Fitted[LongMuTraj] ? LongMuTraj : NoMuTraj] > 20'

    
    # Now we can do whatever analysis we want on it
    # Consider first making a simple histogram of the momenta
    hEtru = [TH1D('hEtru0', ';True Energy (GeV);Entries per 200 keV', 50, 0., float(maxE)), TH1D('hEtru1', ';True Energy (GeV);Entries per 200 keV', 50, 0., float(maxE))]
    hE = [TH1D('hE0', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE)),TH1D('hE1', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))]
    hErec = [TH1D('hErec0', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE)),TH1D('hErec1', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))]
    hEqID = [TH1D('hEqID0', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE)),TH1D('hEqID1', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))]
    hEff  = [TH1D('hEff0', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE)),TH1D('hEff1', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))]
    hback = [TH1D('hback0', ';True Energy (GeV);Fractional Background', 50, 0., float(maxE)),TH1D('hback1', ';True Energy (GeV);Fractional Background', 50, 0., float(maxE))]
    hppull = [TH2D('hppull0', ';True Momentum (MeV/c); (p_{rec} - p_{true}) / #sigma_{p}', 50, 0., float(maxE), 50, -0.5, 0.5),TH2D('hppull1', ';True Momentum (MeV/c); (p_{rec} - p_{true}) / #sigma_{p}', 50, 0., float(maxE), 50, -0.5, 0.5)]
    hkpull   = [TH2D('hkpull0', ';True Momentum (MeV/c);(q/p_{rec} - q/p_{true}) / #sigma_{q/p}', 50, 0., float(maxE), 50, -0.5, 0.5),TH2D('hkpull1', ';True Momentum (MeV/c);(q/p_{rec} - q/p_{true}) / #sigma_{q/p}', 50, 0., float(maxE), 50, -0.5, 0.5)]
    # Create a canvas to catch the output of the tree analysis
    can_E = TCanvas('can_E','can_E')
    tc[0].Draw('MC_PrimaryEng*0.001>>hE0', cuts)
    tc[1].Draw('MC_PrimaryEng*0.001>>hE1', cuts)
    hE[0].Sumw2()
    hE[1].Sumw2()
    tc[0].Draw('MC_PrimaryEng*0.001>>hEtru0')
    tc[1].Draw('MC_PrimaryEng*0.001>>hEtru1')
    hEtru[0].Sumw2()
    hEtru[0].Sumw2()
    # The first option selects the variable and redirects the results to a histogram
    # The second option is a list of conditions on the accepted events. These are the cuts.

    can_Erec = TCanvas('can_Erec','can_Erec')
    tc[0].Draw('MC_PrimaryEng*0.001>>hErec0', cuts + '&& traj_Fitted[classif_LongMuTraj]')    
    tc[1].Draw('MC_PrimaryEng*0.001>>hErec1', cuts + '&& traj_Fitted[classif_LongMuTraj]')
    hErec[0].Sumw2()
    hErec[1].Sumw2()
    
    can_EqID = TCanvas('can_EqID','can_EqID')
    tc[0].Draw('MC_PrimaryEng*0.001>>hEqID0', cuts + '&& traj_Fitted[classif_LongMuTraj]==1 && traj_Charge[classif_LongMuTraj]=='+str(qsel1))
    hEqID[0].Sumw2()
    tc[1].Draw('MC_PrimaryEng*0.001>>hEqID1', cuts + '&& traj_Fitted[classif_LongMuTraj]==1 && traj_Charge[classif_LongMuTraj]=='+str(qsel2))
    hEqID[1].Sumw2()

    can_EqID.cd()
    for i in [0,1]:
        hEqID[i].Divide(hErec[i])
        hEqID[i].SetMinimum(0.0)
        hEqID[i].SetMaximum(1.0)
        hEqID[i].SetMarkerStyle(20+i)
        hEqID[i].SetMarkerColor(i+1)
    hEqID[0].Draw('p')
    hEqID[1].Draw('psame')

    can_Erec.cd()
    
    for i in [0,1]:
        hErec[i].Divide(hE[i])
        hErec[i].SetMinimum(0.0)
        hErec[i].SetMaximum(1.0)
        hErec[i].SetMarkerStyle(20+i)
        hErec[i].SetMarkerColor(i+1)
    hErec[0].Draw('p')
    hErec[1].Draw('psame')

    can_E.cd()
    for i in [0,1]:
        hE[i].Divide(hEtru[i])
        hE[i].SetMinimum(0.0)
        hE[i].SetMaximum(1.0)
        hE[i].SetMarkerStyle(20+i)
        hE[i].SetMarkerColor(i+1)
    hE[0].Draw('p')
    hE[1].Draw('psame')


    can_eff = TCanvas('can_eff','can_eff')
    tc[0].Draw('MC_PrimaryEng*0.001>>hEff0',cuts + ' && traj_Charge[classif_LongMuTraj]=='+str(qsel1))
    hEff[0].Sumw2()
    print hEff[0].Integral(), "\t", hEtru[0].Integral(), "\t", hEff[0].Integral()/hEtru[0].Integral()
    
    hEff[0].Divide(hEtru[0])
    hEff[0].SetStats(0)
    hEff[0].Draw()

    can_back = TCanvas('can_back','can_back')
    tc[0].Draw('MC_PrimaryEng*0.001>>hback0',cuts + ' && traj_Charge[classif_LongMuTraj]!='+str(qsel1))
    hback[0].Sumw2()
    hback[0].Divide(hEtru[0])
    hback[0].SetStats(0)
    hback[0].Draw()
    print hback[0].Integral(), "\t", hEtru[0].Integral(), "\t", hback[0].Integral()/hEtru[0].Integral()

    can_Erec.Print('tree_Erec.eps')
    can_EqID.Print('tree_EqID.eps')
    can_E.Print('tree_Eeff.eps')

    can_ppull = TCanvas("can_ppull","can_ppull")
    can_ppull.SetLogz()
    tc[0].Draw('(traj_Mom[0] - MCtr_Mom[0])/traj_ErrMom[0]:abs(MCtr_Mom[0])*0.001>>hppull1','traj_Fitted[classif_LongMuTraj]==1','colz')
    can_ppull.Print('ppull.eps')
    

    can_kpull = TCanvas("can_kpull","can_kpull")
    
    tc[0].Draw('(traj_Charge[0]/traj_Mom[0] - MCtr_Charge[0]/MCtr_Mom[0])/(traj_ErrMom[0]*traj_Charge[0]/traj_Mom[0]/traj_Mom[0]):abs(MCtr_Mom[0])*0.001>>hkpull0',
            'traj_Fitted[classif_LongMuTraj]==1','colz')
    tc[1].Draw('(traj_Charge[0]/traj_Mom[0] - MCtr_Charge[0]/MCtr_Mom[0])/(traj_ErrMom[0]*traj_Charge[0]/traj_Mom[0]/traj_Mom[0]):abs(MCtr_Mom[0])*0.001>>hkpull1',
            'traj_Fitted[classif_LongMuTraj]==1','colz')
    hkpull[0].FitSlicesY()
    hkpull[1].FitSlicesY()
    can_kpull.SetLogz()
    can_kpull.Print('kpull.eps')
    gDirectory.Get('hkpull0_1').SetMarkerColor(1)
    gDirectory.Get('hkpull0_1').SetMarkerStyle(20)
    gDirectory.Get('hkpull1_1').SetMarkerColor(2)
    gDirectory.Get('hkpull1_1').SetMarkerStyle(21)
    gDirectory.Get('hkpull0_1').Draw('p')
    gDirectory.Get('hkpull1_1').Draw('psame')
    gDirectory.Get('hkpull0_1').SetMinimum(-0.1)
    gDirectory.Get('hkpull0_1').SetMaximum(0.1)
    gDirectory.Get('hkpull0_1').SetStats(0)
    gDirectory.Get('hkpull1_1').SetStats(0)
    can_kpull.Print('kpull_mean.eps')
    gDirectory.Get('hkpull0_2').SetMarkerColor(1)
    gDirectory.Get('hkpull0_2').SetMarkerStyle(20)
    gDirectory.Get('hkpull1_2').SetMarkerColor(2)
    gDirectory.Get('hkpull1_2').SetMarkerStyle(21)
    gDirectory.Get('hkpull0_2').Draw('p')
    gDirectory.Get('hkpull1_2').Draw('psame')
    gDirectory.Get('hkpull0_2').SetMinimum(0.0)
    gDirectory.Get('hkpull0_2').SetMaximum(0.2)
    gDirectory.Get('hkpull0_2').SetStats(0)
    gDirectory.Get('hkpull1_2').SetStats(0)
    can_kpull.Print('kpull_sigma.eps')
    au = raw_input('>')
    
if __name__=='__main__':
    if len(sys.argv) == 8:
        analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], 
                 sys.argv[5], sys.argv[6], sys.argv[7])
    else:
        '''
        An incorrect number of arguments have been used

        use:
          ./run_analysis <absolute path of tree files> <first run number> <last run number>

        '''

        
