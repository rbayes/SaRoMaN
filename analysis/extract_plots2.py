#!/usr/bin/env python

import sys, os
from ROOT import TFile, TH3D, TH2D, TH1D, TCanvas, TLegend, TStyle, gROOT

def nuSTORMCC1list():

    return [["nu_muCC_1/muCC_app.root", "#nu_{#mu} CC Signal"],
            ["nu_muCC_1/mubarCC_app.root", "#bar{#nu}_{#mu} CC mis-ID Background"],
            ["nu_muCC_1/mubarNC_app.root", "#bar{#nu}_{#mu} NC Background"],
            ["nu_muCC_1/eCC_app.root", "#nu_{e} CC Background"]]
#            ["nuSTORM_CC/eNC.root", "#nu_{e} NC Background"]]

def nuSTORMCC2list():

    return [["nu_muCC_2/muCC_app.root", "#nu_{#mu} CC Signal"],
            ["nu_muCC_2/mubarCC_app.root", "#bar{#nu}_{#mu} CC mis-ID Background"],
            ["nu_muCC_2/mubarNC_app.root", "#bar{#nu}_{#mu} NC Background"],
            ["nu_muCC_2/eCC_app.root", "#nu_{e} CC Background"]]
#            ["nuSTORM_CC/eNC.root", "#nu_{e} NC Background"]]

def nuSTORMCC3list():

    return [["nu_muCC_3/muCC_app.root", "#nu_{#mu} CC Signal"],
            ["nu_muCC_3/mubarCC_app.root", "#bar{#nu}_{#mu} CC mis-ID Background"],
            ["nu_muCC_3/mubarNC_app.root", "#bar{#nu}_{#mu} NC Background"],
            ["nu_muCC_3/eCC_app.root", "#nu_{e} CC Background"]]
#            ["nuSTORM_CC/eNC.root", "#nu_{e} NC Background"]]
def nuSTORMCC4list():

    return [["nu_muCC_4/muCC_app.root", "#nu_{#mu} CC Signal"],
            ["nu_muCC_4/mubarCC_app.root", "#bar{#nu}_{#mu} CC mis-ID Bkgd."],
            ["nu_muCC_4/mubarNC_app.root", "#bar{#nu}_{#mu} NC Background"],
            ["nu_muCC_4/eCC_app.root", "#nu_{e} CC Background"]]
#            ["nuSTORM_CC/eNC.root", "#nu_{e} NC Background"]]


def MakeStyle():
    style = TStyle('Display','Display')
    style.SetTitleSize(0.06,'X')
    style.SetTitleSize(0.06,'Y')
    style.SetLabelSize(0.06,'X')
    style.SetLabelSize(0.06,'Y')
    style.SetLegendFillColor(10)
    style.SetPadBottomMargin(0.15)
    style.SetPadTopMargin(0.05)
    style.SetPadRightMargin(0.10)
    style.SetPadLeftMargin(0.15)
    style.SetPalette(1)
    style.SetFrameFillColor(10)
    style.SetHistFillColor(10)
    style.SetOptStat(0)
    return style

def plotEffbyMethod(filespec, dir, outsuffix, methodlist, islog=0):
    
    f = TFile(filespec[0], 'READ')
    l = TLegend(0.7,0.69,0.99,0.99)
    l.SetFillColor(10)
    hm = []
    he = []
    color = 1
    for m in methodlist:
        print filespec[0], 'recEff' + m
        h = f.Get('recEff' + m)
        heff = h.ProjectionY()
        heff.SetLineColor(color)
        heff.SetLineWidth(2)
        heff.SetLineStyle(color)
        heff.SetMarkerStyle(19 + color)
        heff.SetMarkerColor(color)
        heff.SetStats(0)
        h.SetStats(0)
        h.GetXaxis().SetLabelSize(0.06)
        h.GetXaxis().SetTitleSize(0.06)
        h.GetYaxis().SetLabelSize(0.06)
        h.GetYaxis().SetTitleSize(0.06)
        heff.GetXaxis().SetLabelSize(0.06)
        heff.GetXaxis().SetTitleSize(0.06)
        heff.GetYaxis().SetLabelSize(0.06)
        heff.GetYaxis().SetTitleSize(0.06)
        hm.append(h)
        he.append(heff)
        l.AddEntry(he[-1], m + ' Method', 'lp')
        if(color!=2 and color!=4 and color!=9):
            color += 1
        else:
            color += 2
    
    c = TCanvas()
    c.SetBottomMargin(0.15)
    c.SetLeftMargin(0.15)
    c.SetTopMargin(0.05)
    if islog==1: c.SetLogy()
    if islog==1: he[0].GetYaxis().SetRangeUser(1e-8,0.01)
    else : he[0].GetYaxis().SetRangeUser(1e-6,1.0)
    he[0].Draw()
    for h in he[1:]:
        h.Draw('same')
    l.Draw('same')
    c.Print(dir+'/Eff_'+outsuffix+'.eps')

    im = 0
    for h in he:
        h.Draw()
        c.Print(dir+'/'+methodlist[im]+'Eff_'+outsuffix+'.eps')
        im += 1
    
    if islog==1: c.SetLogy(0)
    
    im = 0
    for h in hm:
        h.Draw('colz')
        c.Print(dir+'/'+methodlist[im]+'recEff_'+outsuffix+'.eps')
        im += 1



def plotBackgroundsForMethod(filespeclist, histname, dir, outsuffix, max=0, proj=1):
    flist = []
    hlist = []
    hpylist = []
    hHitlist = []
    i = 0
        
    c1 = TCanvas()
    c1.SetBottomMargin(0.15)
    c1.SetTopMargin(0.05)
    c1.SetLeftMargin(0.15)
    c1.SetRightMargin(0.05)
    l1 = TLegend(0.2,0.7,0.45,0.925)
    l1.SetFillColor(10)
    print filespeclist[0][0], histname
    f = TFile(filespeclist[0][0], 'READ')
    flist.append(f)
    h = f.Get(histname)
    h.SetName(histname + '_' + str(i))
    norm = h.GetEntries()
    if norm == 0: norm == 1
    if proj==0 and len(filespeclist) > 1: h.Scale(1./norm)
    hlist.append(h)
    if proj == 1:
        hy = h.ProjectionY()
    else:
        hy = h
    hpylist.append(hy)
    hpylist[-1].GetXaxis().SetLabelSize(0.06)
    hpylist[-1].GetXaxis().SetTitleSize(0.06)
    hpylist[-1].GetYaxis().SetLabelSize(0.06)
    hpylist[-1].GetYaxis().SetTitleSize(0.06)
    hpylist[-1].GetYaxis().SetTitle("Fractional Efficiency")
    # hpylist[-1].GetXaxis().SetTitle("Fractional Variation in #Delta E")
    hpylist[-1].SetStats(0)
    hpylist[-1].SetLineColor(i + 1)
    hpylist[-1].SetMarkerStyle(i + 20)
    hpylist[-1].SetMarkerColor(i + 1)
    hpylist[-1].SetLineWidth(2)
    hpylist[-1].SetLineStyle(i + 1)
    l1.AddEntry(hpylist[-1], filespeclist[0][1], 'lp')
    if proj == 1: hpylist[-1].GetYaxis().SetRangeUser(5e-7, 50)
    elif proj == 0 and not max == 0: hpylist[-1].GetYaxis().SetRangeUser(1e-5, max)
    c1.cd()
    hpylist[-1].Draw()

    i+=1
    for filespec in filespeclist[1:]:
        print filespec
        f = TFile(filespec[0], 'READ')
        flist.append(f)
        h = f.Get(histname)
        h.SetName(histname + '_' + str(i))
        norm = h.GetEntries()
        if norm == 0: norm = 1
        if proj==0: h.Scale(1./norm)
        hlist.append(h)
        if proj == 1:
            hy = hlist[-1].ProjectionY()
        else:
            hy = hlist[-1]
        h.SetName(histname + '_py' + str(i))
        hpylist.append(hy)
        hpylist[-1].SetStats(0)
        hpylist[-1].SetLineColor(i+1)
        hpylist[-1].SetMarkerStyle(i + 20)
        hpylist[-1].SetMarkerColor(i+1)
        hpylist[-1].SetLineWidth(2)
        hpylist[-1].SetLineStyle(i + 1)
        l1.AddEntry(hpylist[-1], filespec[1], 'lp')
        hpylist[-1].Draw('same')
        i+=1
        if i==4: i+=1
    
    if len(filespeclist) > 1: l1.Draw()
    # if proj == 1:
    c1.SetLogy()
    c1.Print(dir+'/'+outsuffix+'.eps')
    
    
def extractMigrationMatrix(filename, histname, outfile):
    
    o = open(outfile, 'w')
    f = TFile(filename)
    h = f.Get(histname)
    Ny = h.GetNbinsY()
    Nx = h.GetNbinsX()
    nye = Ny-1
    for i in range(Nx):
        s = '{0, %d' % nye
        o.write(s)
        for j in range(Ny):
            s = ', %e' % h.GetBinContent(i+1,j+1)
            o.write(s)
        if i == Nx-1:
            o.write('};\n')
        else:
            o.write('}:\n')
    
def writeTEXMigrationMatrix(filename, histname, outfile, caption):
    
    o = open(outfile, 'w')
    f = TFile(filename)
    h = f.Get(histname)
    Ny = h.GetNbinsY()
    Nx = h.GetNbinsX()
    nye = Ny-1
    o.write('\\begin{table}\n')
    o.write('\centering\n')
    o.write('\\caption{Probability of observing ' + caption +' in units of $10^{-3}$} \n')
    s = '\\begin{tabular}{c|'
    for i in range(Ny): s += 'c'
    s += '}\n'
    o.write(s)
    o.write('\\hline\n')
    s = ' '
    for j in range(Ny): s += ' & %.1f-%.1f GeV' % (h.GetYaxis().GetBinLowEdge(j+1),h.GetYaxis().GetBinUpEdge(j+1)) 
    s += '\\\\\n\\hline\n'
    o.write(s)
    for i in range(Nx):
        elow = h.GetXaxis().GetBinLowEdge(i+1)
        ehigh = h.GetXaxis().GetBinUpEdge(i+1)
        s = '%.1f-%.1f GeV ' % (elow, ehigh)
        o.write(s)
        for j in range(Ny):
            cont = h.GetBinContent(i+1,j+1) * 1000.0
            s = '& %.3f ' % cont
            o.write(s)
        o.write('\\\\\n')
        
    o.write('\\hline\n')
    o.write('\\end{tabular}\n')
    o.write('\\end{table}\n')
    

def main(case):
    ll = []
    knnMMdir = ''
    dir = ''
    print case
    if int(case)==0:
        ll = nuSTORMCC1list()
        dir = '/afs/phas.gla.ac.uk/user/r/rbayes/nuSTORM/analysis/nu_muCC_1/'
        type = ['MuApp','MuBarCCBack','MuBarNCBack','ECCBack']
        methodlist = ['KNN','BDT','MLPBNN','BDTG']
    if int(case)==1:
        ll = nuSTORMCC2list()
        dir = '/afs/phas.gla.ac.uk/user/r/rbayes/nuSTORM/analysis/nu_muCC_2/'
        type = ['MuApp','MuBarCCBack','MuBarNCBack','ECCBack']
        methodlist = ['KNN','BDT','MLPBNN','BDTG']
    if int(case)==2:
        ll = nuSTORMCC3list()
        dir = '/afs/phas.gla.ac.uk/user/r/rbayes/nuSTORM/analysis/nu_muCC_3/'
        type = ['MuApp','MuBarCCBack','MuBarNCBack','ECCBack']
        methodlist = ['KNN','BDT','MLPBNN','BDTG']
    if int(case)==3:
        ll = nuSTORMCC4list()
        dir = '/afs/phas.gla.ac.uk/user/r/rbayes/nuSTORM/analysis/nu_muCC_4/'
        type = ['MuApp','MuBarCCBack','MuBarNCBack','ECCBack']
        methodlist = ['KNN','BDT','MLPBNN','BDTG']
    else:
        print "Do not recognize case."
        return 1
        
    d = os.path.dirname(dir)
    if not os.path.exists(d):
        os.makedirs(d)
    plotEffbyMethod(ll[0], dir + 'aplots', type[0], methodlist)
    plotEffbyMethod(ll[1], dir + 'aplots', type[1], methodlist, 1)
    plotEffbyMethod(ll[2], dir + 'aplots', type[2], methodlist, 1)
    plotEffbyMethod(ll[3], dir + 'aplots', type[3], methodlist, 1)

    for m in methodlist:
        plotBackgroundsForMethod(ll[0:2], 'Eff' + m, dir + 'aplots', 'Eff_'+m, 1.0, 0.0)
        plotBackgroundsForMethod(ll, 'recEff' + m, dir + 'aplots', 'allEff_'+m)
        plotBackgroundsForMethod(ll[1:], 'recEff' + m, dir + 'aplots', 'Bkgd_'+m)
        if int(case) != 3:
            plotBackgroundsForMethod(ll, 'Nhits_' + m, dir + 'aplots', 'Nhits_'+m, 0.2, 0)
            plotBackgroundsForMethod(ll, 'Rp_' + m, dir + 'aplots', 'Rp_'+m, 0, 0)
            plotBackgroundsForMethod(ll, 'ErrqP_' + m, dir + 'aplots', 'ErrqP_'+m, 0, 0)
            plotBackgroundsForMethod(ll, 'MeanEDep_' + m, dir + 'aplots', 'MeanEDep_'+m, 0, 0)
            plotBackgroundsForMethod(ll, 'EngVar_' + m, dir + 'aplots', 'EngVar_'+m, 0.08, 0)
            plotBackgroundsForMethod(ll, 'Qt_' + m, dir + 'aplots', 'Qt_'+m, 0, 0)
        for i in range(len(ll)): extractMigrationMatrix(ll[i][0], 'recEff' + m, dir + m + '_MM/'+type[i]+'.txt')
        # for i in range(len(ll)): writeTEXMigrationMatrix(ll[i][0], 'recEff' + m, dir + m + '_TM/'+type[i]+'.tex', ll[i][1])
    
    
    return 0

if __name__ == "__main__":
    print sys.argv[1]
    main(sys.argv[1])
