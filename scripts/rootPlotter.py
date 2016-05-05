import getopt, sys
from ROOT import TFile, TTree, TH1D, TH2D, TGraph, TMarker, TLine, TBox, TCanvas

class eventPlotter:
    '''
    Class takes reconstructed output tree from the saroman package
    '''
    def __init__(self):
        '''
        
        '''
        self.initevent = 0
        self.mode      = 0
        self.infile    = ''
        self.outfile   = ''

    def runPlotter(self, argv):
        
        try:
            opts, args = getopt.getopt(argv, 'fs', ['infile=','outfile=','e='])
        except getopt.GetoptError:
            print "rootPlotter.py --infile="
            print "               --outfile="
            print "               --e= event to view the indicated event"
            print "               -f to scan failed events"
            print "               -s to scan all events"
        print opts, args
        for opt, arg in opts:
            print opt, arg
            if opt == '--infile':
                self.infile = arg
                print opt, arg
            if opt == '--outfile':
                self.outfile = arg
                print opt, arg
            if opt == '--e':
                # view a single event
                self.initevent = int(arg)
                self.mode  = 1 # view a single event 
                print opt, arg
            if opt == '-f':
                # scan through failed events
                self.mode  = 2 # scan failed events only
                print opt, arg
            if opt == '-s':
                # scan through all events (from initevent)
                self.mode = 0
                print opt, arg
        if self.infile == "":
            print "No input file provided."
            print "rootPlotter.py -F input file"
            print "               -O output file"
            print "               -e event to view the indicated event"
            print "               -f to scan failed events"
            print "               -s to scan all events"
            sys.exit()
        if self.outfile == "":
            print "No output file provided. Will use default output name."
            self.outfile = "Event_"+str(self.event)+".eps"
        if self.mode == 1:
            self.generateSingleEvent()
        elif self.mode == 2:
            self.loopThroughFailedEvents()
        elif self.mode == 0:
            self.loopThroughAllEvents()

    def generateSingleEvent(self):
        '''
        Plot single event from input file, print the event content and exit
        '''
        InFile = TFile(self.infile)
        tree = InFile.Get("tree")
        
        self.generateEventPlot(tree, self.initevent)
    
    def loopThroughFailedEvents(self):
        '''
        Plot loop through events from input file, print the event content and continue
        '''
        InFile = TFile(self.infile)
        tree = InFile.Get("tree")
        
        entries = tree.GetEntriesFast()
        for ievt in xrange(self.initevent, entries):
            ientry = tree.LoadTree( ievt )
            if ientry < 0:
                break
            nb = tree.GetEntry( ievt )
            # ntraj = event.classif_NumTrajectory
            if nb < 0:
                continue
            failed = 0
            for traj_fit in tree.traj_Fitted:
                if traj_fit == 0:
                    failed = 1
            if failed == 1:
                self.generateEventPlot(tree, ievt)
 
    def loopThroughFailedEvents(self):
        '''
        Loop through events from input file, plot the event content and continue
        '''
        InFile = TFile(self.infile)
        tree = InFile.Get("tree")
        
        entries = tree.GetEntriesFast()
        for ievt in xrange(self.initevent, entries):
            ientry = tree.LoadTree( ievt )
            if ientry < 0:
                break
            nb = tree.GetEntry( ievt )
            # ntraj = event.classif_NumTrajectory
            if nb < 0:
                continue
            self.generateEventPlot(tree, ievt)

    def generateEventPlot(self, tree, evt):
        '''
        Plot the trajectories from the given event.
        '''
        c = TCanvas()
        fitted2d = TH2D("fitted2d","fitted2d",1000, -1800,1800, 1000, -1000,1000)
        raw_mu2d = TH2D("raw_mu2d","raw_mu2d",1000, -1800,1800, 1000, -1000,1000)
        raw_other2d = TH2D("raw_other2d",
                           "; Z Position (mm); Y Position (mm)",
                           1000, -1800,1800, 1000, -1000,1000)
        
        fittedx= TH2D("fittedx","fittedx",1000, -1800,1800, 1000, -1000,1000)
        raw_mux = TH2D("raw_mux","raw_mux",1000, -1800,1800, 1000, -1000,1000)
        raw_otherx = TH2D("raw_otherx","; Z Position (mm); Y Position (mm)",
                          1000, -1800,1800, 1000, -1000,1000)

        evtcut = "MC_Evt=="+str(evt)
        tree.Draw("trajNode_YPos:trajNode_ZPos>>fitted2d",
                  evtcut + " && traj_Fitted==1","COLZ")
        tree.Draw("raw_Ymeas:raw_Zmeas>>raw_mu2d",
                  evtcut + " && raw_MotherProp>=0.5","COLZ")
        tree.Draw("raw_Ymeas:raw_Zmeas>>raw_other2d",
                  evtcut + " && raw_MotherProp<0.5","COLZ")
  
        tree.Draw("trajNode_XPos:trajNode_ZPos>>fittedx","traj_Fitted==1","COLZ")
        tree.Draw("raw_Xmeas:raw_Zmeas>>raw_mux","raw_MotherProp>=0.5","COLZ")
        tree.Draw("raw_Xmeas:raw_Zmeas>>raw_otherx","raw_MotherProp<0.5","COLZ")

        
        fitted2d.SetMarkerColor(4)
        fitted2d.SetMarkerStyle(24)
        
        raw_mu2d.SetMarkerColor(42)
        raw_mu2d.SetMarkerStyle(21)
        
        raw_other2d.SetMarkerColor(3)
        raw_other2d.SetMarkerStyle(21)
        
        fittedx.SetMarkerColor(4)
        fittedx.SetMarkerStyle(24)
        
        raw_mux.SetMarkerColor(42)
        raw_mux.SetMarkerStyle(21)
        
        raw_otherx.SetMarkerColor(2)
        raw_otherx.SetMarkerStyle(21)
        
  
        #  Draw the detector
  
        ymax = 1000
        ymin =-1000
  
        box1 =  TBox(-1630.5,ymin,-1600.5,ymax)
        box1.SetFillColor(1)
        box1.Draw()
        box2 =  TBox(-1100,ymin,-1070,ymax)
        box2.SetFillColor(1)
        box2.Draw("same")
        box3 =  TBox(-960,ymin,-930,ymax)
        box3.SetFillColor(1)
        box3.Draw("same")
        box4 =  TBox(-430,ymin,-400,ymax)
        box4.SetFillColor(1)
        box4.Draw("same")
        box5 =  TBox(-360,ymin,-330,ymax)
        box5.SetFillColor(1)
        box5.Draw("same")
        box6 =  TBox(-295,ymin,-265,ymax)
        box6.SetFillColor(1)
        box6.Draw("same")
        box7 =  TBox(-230,ymin,-200,ymax)
        box7.SetFillColor(1)
        box7.Draw("same")
        box8 =  TBox(-125,ymin,-95,ymax)
        box8.SetFillColor(1)
        box8.Draw("same")
        box9 =  TBox(-20,ymin,10,ymax)
        box9.SetFillColor(1)
        box9.Draw("same")
        box10 =  TBox(85,ymin,115,ymax)
        box10.SetFillColor(1)
        box10.Draw("same")
        box11 =  TBox(225,ymin,255,ymax)
        box11.SetFillColor(1)
        box11.Draw("same")
        box12 =  TBox(365,ymin,395,ymax)
        box12.SetFillColor(1)
        box12.Draw("same")
        box13 =  TBox(505,ymin,535,ymax)
        box13.SetFillColor(1)
        box13.Draw("same")
        box14 =  TBox(645,ymin,675,ymax)
        box14.SetFillColor(1)
        box14.Draw("same")
        box15 =  TBox(790,ymin,820,ymax)
        box15.SetFillColor(1)
        box15.Draw("same")
        box16 =  TBox(930,ymin,960,ymax)
        box16.SetFillColor(1)
        box16.Draw("same")
        box17 =  TBox(1070,ymin,1100,ymax)
        box17.SetFillColor(1)
        box17.Draw("same")
        box18 =  TBox(1600,ymin,1630,ymax)
        box18.SetFillColor(1)
        box18.Draw("same")
        
        raw_mu2d.Draw("same");
        
        raw_other2d.Draw("same");
        fitted2d.Draw("same");
        
        au = raw_input('\'print\' to save event, \'exit\' to end, press enter key to continue >>')
        if au == 'print':
            self.outfile += "_evt"+str(evt)+".eps";
            c.Print(self.outfile);
        elif au == 'exit':
            sys.exit()
        
        
if __name__=='__main__':
    
    plotter = eventPlotter()
    plotter.runPlotter(sys.argv[1:])
