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

def analysis(basename, first, last, maxE, qsel):
    # this initializes a TChain object to read a tree file
    tc = TChain('tree')

    # here is a simple for loop
    for i in range(int(first), int(last)+1):
        # define the filename
        filename = basename
        filename += str(i)
        filename += '.root'
        tc.Add(filename)

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
    hEtru = TH1D('hEtru', ';True Energy (GeV);Entries per 200 keV', 50, 0., float(maxE))
    hE = TH1D('hE', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))
    hErec = TH1D('hErec', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))
    hEqID = TH1D('hEqID', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))
    hEff  = TH1D('hEff', ';True Energy (GeV);Fractional Efficiency', 50, 0., float(maxE))
    hback = TH1D('hback', ';True Energy (GeV);Fractional Background', 50, 0., float(maxE))
    hppull = TH2D('hppull', ';True Momentum (MeV/c); (p_{rec} - p_{true}) / #sigma_{p}', 50, 0., float(maxE), 50, -0.5, 0.5)
    hkpull   = TH2D('hkpull', ';True Momentum (MeV/c);(q/p_{rec} - q/p_{true}) / #sigma_{q/p}', 50, 0., float(maxE), 50, -0.5, 0.5)
    # Create a canvas to catch the output of the tree analysis
    can_E = TCanvas('can_E','can_E')
    tc.Draw('MC_PrimaryEng*0.001>>hE', cuts)
    hE.Sumw2()
    tc.Draw('MC_PrimaryEng*0.001>>hEtru')
    hEtru.Sumw2()
    # The first option selects the variable and redirects the results to a histogram
    # The second option is a list of conditions on the accepted events. These are the cuts.

    can_Erec = TCanvas('can_Erec','can_Erec')
    tc.Draw('MC_PrimaryEng*0.001>>hErec', cuts + '&& traj_Fitted[classif_LongMuTraj]')
    hErec.Sumw2()
    
    can_EqID = TCanvas('can_EqID','can_EqID')
    tc.Draw('MC_PrimaryEng*0.001>>hEqID', cuts + '&& traj_Fitted[classif_LongMuTraj]==1 && traj_Charge[classif_LongMuTraj]=='+str(qsel))
    hEqID.Sumw2()

    can_EqID.cd()
    hEqID.Divide(hErec)
    hEqID.SetMinimum(0.0)
    hEqID.SetMaximum(1.0)
    hEqID.Draw()

    can_Erec.cd()

    hErec.Divide(hE)
    hErec.SetMinimum(0.0)
    hErec.SetMaximum(1.0)
    hErec.Draw()

    can_E.cd()
    hE.Divide(hEtru)
    hE.SetMinimum(0.0)
    hE.SetMaximum(1.0)
    hE.Draw()


    can_eff = TCanvas('can_eff','can_eff')
    tc.Draw('MC_PrimaryEng*0.001>>hEff',cuts + ' && traj_Charge[classif_LongMuTraj]=='+str(qsel))
    hEff.Sumw2()
    print hEff.Integral(), "\t", hE.Integral(), "\t", hEff.Integral()/hE.Integral()
    
    hEff.Divide(hE)
    hEff.SetStats(0)
    hEff.Draw()

    can_back = TCanvas('can_back','can_back')
    tc.Draw('MC_PrimaryEng*0.001>>hback',cuts + ' && traj_Charge[classif_LongMuTraj]!='+str(qsel))
    hback.Sumw2()
    hback.Divide(hE)
    hback.SetStats(0)
    hback.Draw()
    print hback.Integral(), "\t", hE.Integral(), "\t", hback.Integral()/hE.Integral()

    can_Erec.Print('tree_Erec.eps')
    can_EqID.Print('tree_EqID.eps')
    can_E.Print('tree_Eeff.eps')

    can_ppull = TCanvas("can_ppull","can_ppull")
    can_ppull.SetLogz()
    tc.Draw('(traj_Mom[0] - MCtr_Mom[0])/traj_ErrMom[0]:abs(MCtr_Mom[0])*0.001>>hppull','traj_Fitted[classif_LongMuTraj]==1','colz')
    can_ppull.Print('ppull.eps')
    

    can_kpull = TCanvas("can_kpull","can_kpull")
    
    tc.Draw('(traj_Charge[0]/traj_Mom[0] - MCtr_Charge[0]/MCtr_Mom[0])/(traj_ErrMom[0]*traj_Charge[0]/traj_Mom[0]/traj_Mom[0]):abs(MCtr_Mom[0])*0.001>>hkpull',
            'traj_Fitted[classif_LongMuTraj]==1','colz')
    hkpull.FitSlicesY()
    can_kpull.SetLogz()
    can_kpull.Print('kpull.eps')
    gDirectory.Get('hkpull_1').Draw()
    gDirectory.Get('hkpull_1').SetMinimum(-0.1)
    gDirectory.Get('hkpull_1').SetMaximum(0.1)
    can_kpull.Print('kpull_mean.eps')
    gDirectory.Get('hkpull_2').Draw()
    gDirectory.Get('hkpull_2').SetMinimum(0.0)
    gDirectory.Get('hkpull_2').SetMaximum(0.2)
    can_kpull.Print('kpull_sigma.eps')
    au = raw_input('>')
    
if __name__=='__main__':
    if len(sys.argv) == 6:
        analysis(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        '''
        An incorrect number of arguments have been used

        use:
          ./run_analysis <absolute path of tree files> <first run number> <last run number>

        '''

        
