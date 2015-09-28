
import os
import math
import subprocess
import shutil

import print_config

class saroman:

    def __init__(self):
        #Set up paths #
        self.home = os.getenv("HOME") 
        self.exec_base = os.path.join(self.home, 'SaRoMaN')
        self.out_base  = os.path.join(self.home, 'out')
        self.scripts_dir = os.path.join(self.exec_base, 'saroman')
        self.third_party_support = self.home + "/nuSTORM/third_party"

        self.print_config=print_config.print_config()

        #Should be implemented as input values#
        self.train_sample = 0
        self.part = 'mu'#'14'
        self.pid = 14
        self.seed = 10
        self.Nevts = 10
        self.inttype = 'All'
        self.Bfield = 1.0

        #General class variables#
        self.ASeed = str(self.seed + 1000)

    def check_make_dir(self,dirname):
        '''
        Check if the target directories dirname exist, if not create it.
        '''
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    def print_file(self,filename,data):
        '''
        Print data to file filename.
        '''
        outfile = open(filename,'w+')
        outfile.write(data)
        outfile.close()

    
    def print_outdata_file(self,filename,data):
        '''
        Print call outdata to file filename.
        '''
        outfile = open(filename,'w+')
        subprocess.call(data, stdout=outfile)
        outfile.close()

    def shell_source(self, script):
        '''
        Run script as source and update environment accordingly
        '''
        pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE,shell=True)
        output = pipe.communicate()[0]
        env = dict((line.split("=",1) for line in output.splitlines()))
        os.environ.update(env)

    def set_environment(self):
        '''
        Set up the environment required to run.
        '''
        genie_support_ext = self.third_party_support + "/genie_source/src/scripts/build/ext"
        
        os.environ['GENIE_SUPPORT_EXT'] = genie_support_ext
        os.environ['THIRD_PARTY_SUPPORT'] = self.third_party_support 

        # Start by sourcing the root installation
        self.shell_source(genie_support_ext + "/v5_34_10/download/root/bin/thisroot.sh")
        
        # Now add the GENIE libraries to path and the ld library
        os.environ['GENIE'] = self.home + "/nuSTORM/third_party/genie_source"
        os.environ['GENIE_INCDIR']=self.third_party_support + "/include/GENIE"
        os.environ['GENIE_LIBDIR']=self.third_party_support + "/lib"
        os.environ['PATH']+= os.pathsep + self.third_party_support + "/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/lib"
        # Now add the supporting libraries
        # PYTHIA
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + genie_support_ext + "/v6_424/lib"
        # log4cpp
        os.environ['PATH']+= os.pathsep + self.third_party_support+ "/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep+self.third_party_support + "/libs"
        # LHAPDF
        os.environ['LHAPATH']=self.home + "/nuSTORM/third_party/genie_source/data/evgen/pdfs"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + genie_support_ext + "/v5_7_0/stage/lib"
        # CLHEP
        clhep_base_dir = self.home + "/clhep"
        os.environ['CLHEP_BASE_DIR'] = clhep_base_dir 
        os.environ['PATH']+= os.pathsep + clhep_base_dir + "/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + clhep_base_dir + "/lib"
        # GEANT4
        g4install = self.third_party_support + "/geant4.10.00-install"
        os.environ['G4INSTALL']=g4install
        os.environ['G4WORKDIR']=g4install
        os.environ['G4SYSTEM']="Linux-g++"
        os.environ['G4INCLUDE']=g4install + "/include/Geant4"
        os.environ['G4LIB']=g4install + "/lib64"
        os.environ['G4LEVELGAMMADATA']=g4install + "/share/Geant4-10.0.0/data/PhotonEvaporation3.0"
        os.environ['G4RADIOACTIVEDATA']=g4install + "/share/Geant4-10.0.0/data/RadiativeDecay4.0"
        os.environ['G4LEDATA']=g4install + "/share/Geant4-10.0.0/data/G4EMLOW6.35"
        os.environ['G4NEUTRONHPDATA']=g4install + "/share/Geant4-10.0.0/data/G4NDL4.4"
        os.environ['G4ABLADATA']=g4install + "/share/Geant4-10.0.0/data/G4ABLA3.0"
        os.environ['G4SAIDXSDATA']=g4install + "/share/Geant4-10.0.0/data/G4SAIDDATA1.1"
        os.environ['G4ENSDFSTATE']=g4install + "/share/Geant4-10.0.0/data/G4ENSDFSTATE1.0"
        os.environ['G4NEUTRONXSDATA']=g4install + "/share/Geant4-10.0.0/data/G4NEUTRONXS1.4"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + g4install + "/lib64"
        # Support packages
        # CRY simulation information 
        os.environ['CRYPATH']=self.third_party_support + "/cry_v1.7"
        os.environ['CRY_INCDIR']=self.third_party_support + "/cry_v1.7/include"
        os.environ['CRY_LIBDIR']=self.third_party_support + "/cry_v1.7/lib"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/cry_v1.7/lib"
        # BHEP
        os.environ['BHEP_LIB']=self.third_party_support + "/bhep/lib"
        os.environ['BHEP_PATH']=self.third_party_support + "/bhep/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/bhep/lib"
        os.environ['PATH']+= os.pathsep + "$BHEP_PATH"
        
#
#    def set_MIND_parameters(self, Type, xdim, ydim, zdim, \
#                           FeThickness, SciThickness, AlThickness,Gap):

    def submit_run(self, seed, Nevts, pid, inttype, BField):
        self.seed = seed
        self.Nevts = Nevts
        self.BField = BField
        # particle id
        pdgid = {14:'mu',-14:'mubar',12:'e',-12:'ebar',16:'tau',-16:'taubar'}
        part = ''
        if pid % 2 == 0 and \
               pid / 2. > 5. and pid / 2. < 9.:
            # the particle is a neutrino
            self.part = pdgid[pid]
        else:
            print "Invalid particle specified"
            return 1
        # special interaction type for training purposes only
        if inttype.find('Train') > 0:
            self.train_sample = 1
            if inttype.find('CC') > 0:
                inttype = 'CC'
            if inttype.find('NC') > 0:
                inttype = 'NC'

        self.run_genie()
        self.run_simulation()
        self.run_digitization()
        self.run_reconstruction()

    def run_genie(self):
        genieOutBase = os.path.join(self.out_base, 'genie_samples')
        genieOutDir = os.path.join(genieOutBase,'nd_'+self.part+self.inttype)
        genieOutLog = os.path.join(genieOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
        FeTargetCode="1000260560"
        C12HTargetCode='1000060120[0.922582],1000010010[0.077418]'
        equation ='1/(0.0069*x) +1'
        vE ='0.01,4.0'

        # Check if the target directories exist if not create it
        self.check_make_dir(genieOutBase)
        self.check_make_dir(genieOutDir)

        command = ['gevgen','-r',str(self.seed),'-n',str(self.Nevts),'-p',str(self.pid),'-t',FeTargetCode,
                             '-e',vE,'-f',equation]
        command2= ['gevgen','-r',self.ASeed,'-n',str(self.Nevts),'-p',str(self.pid),'-t',C12HTargetCode,
                             '-e',vE,'-f',equation]

        if self.inttype != 'All':
            command += ['--event-generator-list',self.inttype]
            command2 += ['--event-generator-list',self.inttype]

        command += ['--seed',str(self.seed),'--cross-sections',self.scripts_dir+'/xsec_Fe56_splines.xml']
        command2 += ['--seed',self.ASeed,'--cross-sections',self.scripts_dir+'/xsec_C12+H1_splines.xml']

        print subprocess.list2cmdline(command)
        self.print_outdata_file(genieOutLog,command)

        geniefile='gntp.'+str(self.seed)+'.ghep.root'
        geniedest=genieOutDir+'/ev0_'+str(self.seed)+'_'+str(self.pid)+'_'+str(FeTargetCode)+'_'+str(self.Nevts)+'.root'
        shutil.move(geniefile, geniedest)

        print subprocess.list2cmdline(command2)
        self.print_outdata_file(genieOutLog,command2)
        
        geniefile='gntp.'+self.ASeed+'.ghep.root'
        geniedest=genieOutDir+'/ev0_'+self.ASeed+'_'+str(self.pid)+'_'+str(C12HTargetCode)+'_'+str(self.Nevts)+'.root'
        shutil.move(geniefile, geniedest)
            

    def run_simulation(self):
        mindG4OutBase = os.path.join(self.out_base, 'G4_out')
        mindG4OutConfig = os.path.join(self.out_base, 'G4_config')
        mindG4OutDir = os.path.join(mindG4OutBase,'nd_'+self.part+self.inttype)
        mindG4ConfigDir = os.path.join(mindG4OutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.check_make_dir(mindG4OutBase)
        self.check_make_dir(mindG4OutConfig)
        self.check_make_dir(mindG4OutDir)
        self.check_make_dir(mindG4ConfigDir)

        mindG4config = os.path.join(mindG4ConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.config')
        self.print_config.print_mindG4_config(vars(self),mindG4config)

        mindG4OutLog = os.path.join(mindG4OutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')

        self.shell_source(self.third_party_support+'/geant4.10.00-install/bin/geant4.sh')

        command = [self.exec_base+'/sciNDG4/mindG4',mindG4config]
        print subprocess.list2cmdline(command)
        self.print_outdata_file(mindG4OutLog,command)

    def run_digitization(self):
        digiOutBase = os.path.join(self.out_base, 'digi_out')
        digiOutConfig = os.path.join(self.out_base, 'digi_param')
        digiOutDir = os.path.join(digiOutBase,'nd_'+self.part+self.inttype)
        digiConfigDir = os.path.join(digiOutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.check_make_dir(digiOutBase)
        self.check_make_dir(digiOutConfig)
        self.check_make_dir(digiOutDir)
        self.check_make_dir(digiConfigDir)

        digiConfig = os.path.join(digiConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.digi.param')
        self.print_config.print_digi_config(vars(self),digiConfig)

        digiOutLog = os.path.join(digiOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
        command = [self.exec_base+"/digi_ND/examples/simple_smear",digiConfig]
        print subprocess.list2cmdline(command)
        self.print_outdata_file(digiOutLog,command)

    def run_reconstruction(self):
        recOutBase = os.path.join(self.out_base, 'rec_out')
        recOutConfig = os.path.join(self.out_base, 'rec_param')
        recOutDir = os.path.join(recOutBase,'nd_'+self.part+self.inttype)
        recConfigDir = os.path.join(recOutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.check_make_dir(recOutBase)
        self.check_make_dir(recOutConfig)
        self.check_make_dir(recOutDir)
        self.check_make_dir(recConfigDir)

        recConfig = os.path.join(recConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.rec.param')
        self.print_config.print_rec_config(vars(self),recConfig)

        recOutLog = os.path.join(recOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
        command = [self.exec_base+"/mind_rec/examples/fit_tracks",recConfig, str(self.Nevts)]
        print subprocess.list2cmdline(command)
        self.print_outdata_file(recOutLog,command)

if __name__ == "__main__":

    s=saroman()
    s.set_environment()
    s.run_genie()
    s.run_simulation()
    s.run_digitization()
    s.run_reconstruction()


