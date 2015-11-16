#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 4/11-2015
#Created on 25/9-2015
#######################################################################################################################
#General python import
#######################################################################################################################
import os
import math
import subprocess
import shutil
import sys
import getopt
import time
#sys.path.append('pythonlib')
#######################################################################################################################
#Importing own python files
#######################################################################################################################
from pythonlib import print_config
from pythonlib import field_map_generator
from pythonlib import handle_third_party
from pythonlib import xml_parser
#import print_config
#import field_map_generator
#import handle_third_party
#######################################################################################################################
#Class generation
#######################################################################################################################
class saroman:
    '''
    Class saroman handles calls to the simulation software.
    Setup:
    Make sure the paths are correctly set. self.exec_base, self.out_base, self.scripts_dir and self.third_party_support.
    Also make sure all the values set in __init__ are correct for the geometry you want to simulate.
    If you do not need to install the third party software -> set self.need_third_party_install = False
    If you do not need to build our software -> set self.need_own_install = False

    How to run:
    At the moment, set up everything in saroman.py then simply run with python saroman.py
    Command line commands will be implemented in the future.
    '''

    '''
    Input same as submit script
    perhaps make it possible to set up variables in another py file and import saroman.
    '''

    def __init__(self):
        #Set up paths #
        self.home = os.getenv("HOME") 
        self.exec_base = os.path.join(self.home, 'SaRoMaN')
        self.out_base  = os.path.join(self.home, 'out')
        self.scripts_dir = os.path.join(self.exec_base, 'saroman')
        self.third_party_support = '/data/neutrino05/phallsjo/third_party'
        self.xml_file_path = os.path.join(self.exec_base,'MIND.gdml')
        self.parsed_file_path  = os.path.join(self.exec_base,'parsedGdml.log')

        #General flags
        self.need_third_party_install = False
        self.need_own_install = False
        self.generate_field_map = True # If false remember to change self.field_map_name to point to your field map!
        self.parse_gdml = True

        #Should be implemented as input values#
        self.train_sample = 0
        self.part = 'mu-'#'14'
        self.pid = -14
        self.seed = 10
        self.Nevts = 10
        self.inttype = 'CC'
        self.Bfield = 1.5

        #Mind geometry
        #Different types of geometry, 3 represents a rectangular detector.
        self.MIND_type = 3#0   # Cylinder
        self.MIND_xdim = 0.96#7.0 # m
        self.MIND_ydim = 0.96#6.0 # m
        self.MIND_zdim = 2.0#13.0 # m
        #Not used for rectangular detector
        self.MIND_vertex_xdim = 0#2.0 # m
        self.MIND_vertex_ydim = 0#2.0 # m
        self.MIND_vertex_zdim = 0#2.0 # m
        #Used as the size of the steal plates on each side of the detector.
        self.MIND_ear_xdim = 2.54#3.5-0.96#0.4393 # m
        self.MIND_ear_ydim = 1.04#2.0-0.96#2.8994 # m
        self.MIND_bore_diameter = 0.2 # m
        #Mind internal dimensions
        #de_dx found at pdg.lbl.gov/2015/AtomicNuclearProperties
        self.MIND_active_mat = 'G4_POLYSTYRENE'
        self.MIND_active_de_dx = 0.2052 #MeV/mm
        self.MIND_thickness_active = 1.5 # cm
        self.MIND_thickness_sigma = self.MIND_thickness_active / math.sqrt(12)
        self.MIND_width_active = 1.5 #cm
        self.MIND_width_sigma = self.MIND_width_active / math.sqrt(12)
        self.MIND_rad_length_active = 413.1 #mm
        self.MIND_npanels = 2 #Describe howmany 'parts' the magnetic field has.
        self.MIND_active_layers = 1 #1
        self.MIND_passive_mat = 'G4_Fe'
        self.MIND_passive_de_dx = 1.143 #MeV/mm
        self.MIND_thickness_passive = 3.0#1.5 # cm
        self.MIND_rad_length_passive = 17.58 #mm
        self.MIND_bracing_mat = 'G4_Al'
        self.MIND_bracing_de_dx = 0.4358 #MeV/mm
        self.MIND_thickness_bracing = 0.1 # cm
        self.MIND_thickness_air = 0.5 # cm
        self.MIND_thickness_air_mm = 10*self.MIND_thickness_air
        self.MIND_rad_length_air = 303.9 #mm
        self.MIND_min_eng_at_plane = 0#0.000016 #MeV
        self.MIND_module_length = (self.MIND_thickness_active*self.MIND_active_layers +
                                   self.MIND_thickness_passive+self.MIND_thickness_bracing +
                                   self.MIND_thickness_air)

        self.MIND_module_de_dx = (self.MIND_thickness_active*self.MIND_active_layers*self.MIND_active_de_dx +
                                  self.MIND_thickness_passive*self.MIND_passive_de_dx +
                                  self.MIND_thickness_bracing*self.MIND_bracing_de_dx)/self.MIND_module_length

        #Print config object, used to generate config files correctly
        #Set to either single_particle generation or generation through genie.
        self.GenerationMode = 'SINGLE_PARTICLE' #GENIE
        self.print_config=print_config(self.GenerationMode)

        #Setup for field_map_generator.py
        self.CreateFieldMap = True
        self.field_map_generator = field_map_generator(self.Bfield,self.MIND_ydim+self.MIND_ear_ydim,
            self.MIND_xdim+self.MIND_ear_xdim, self.MIND_npanels)
        self.field_map_name = 'field_map_test.res'
        self.field_map_folder = self.out_base
        self.field_map_full_name =os.path.join(self.field_map_folder,self.field_map_name)

        #Setup for handle_third_party.py
        self.handle_third_party = handle_third_party(self.exec_base,self.third_party_support)

        #Setup for xml_parser.py
        self.xml_parser = xml_parser(self.xml_file_path,self.parsed_file_path)
        self.useGDML = 0
        if self.parse_gdml:
            self.useGDML = 1

        #General class variables
        self.ASeed = str(self.seed + 1000)

        #Config file variables not generalized
        self.config_digi_seed = 107311191
        self.config_digi_eng_res = 0.11 #(%%)
        self.config_digi_gaus_sigma = 1.0 #(cm)

        self.config_rec_likelihood = 1 #1 if we require output of likelihood info
        self.config_rec_max_multOcc_plane = 10
        self.config_rec_plane_occupancy = 10
        self.config_rec_min_iso_prop = 0.8
        self.config_rec_min_check_nodes = 3
        self.config_rec_max_seed_hits = 20
        self.config_rec_min_seed_hits = 7
        self.config_rec_pat_rec = 1 # Do Pattern recognition (0=false, 1=true)
        self.config_rec_fac_refit = 10000
        self.config_rec_refit = 1 # fit data twice (0=false, 1=true)
        self.config_rec_z_cut = 300 #cm
        self.config_rec_x_cut = 50 #cm
        self.config_rec_y_cut = 50 #cm
        self.config_rec_low_fit_cut2 = 0.6
        self.config_rec_low_fit_cut0 = 0.8
        self.config_rec_high_pass_hits = 500
        self.config_rec_low_pass_hits = 6
        self.config_rec_max_n_trajs = 4
        self.config_rec_max_coincidence = 0.5
        self.config_rec_accept_chi = 20
        self.config_rec_max_traj = 40
        self.config_rec_max_sep = 7 #cm
        self.config_rec_min_blobOcc = 2
        self.config_rec_max_blobSkip = 0.2
        self.config_rec_min_use_prop = 0.7
        self.config_rec_max_consec_missed_planes = 3
        self.config_rec_pat_rec_max_outliners = 10
        self.config_rec_pat_rec_max_chi = 20
        self.config_rec_chi2node_max = 20
        self.config_rec_max_outliners = 5
        self.config_rec_chi2fit_max = 50
        # verbosities for recpack services
        self.config_rec_vfit = 0
        self.config_rec_vmat = 0
        self.config_rec_vnav = 0
        self.config_rec_vmod = 0
        self.config_rec_vsim = 0
        self.config_rec_model = 'particle/helix'
        self.config_rec_kfitter = 'kalman'
        self.config_rec_gen_seed = 373940592
        self.config_rec_boxX = 1.5 #cm
        self.config_rec_do_clust = 1
        self.config_rec_detect = 'tracking'
        self.config_rec_step_size = 1 #cm
        self.config_rec_pos_res = 1.5 #cm
        self.config_rec_meas_type = 'xyz'
        self.config_rec_WLSatten = 5000
        # relative density, Sc/Fe, AIR/Sc.
        self.config_rec_rel_denSI = 0.135
        self.config_rec_rel_denAS = 0.001

        self.config_mindG4_min_kilEng = 300 #MeV
        self.config_mindG4_prod_cut = 30 #mm
        self.config_mindG4_vertex_location = 'GAUSS'
        self.config_mindG4_isUniform = 0
        
        
        

        

################################################
    def Check_make_dir(self,dirname):
        '''
        Check if the target directories dirname exist, if not create it.
        '''
        if not os.path.exists(dirname):
            os.makedirs(dirname)

    def Generate_field_map(self):
        if(self.CreateFieldMap):
            self.field_map_generator.Print_field_to_file(self.field_map_full_name)

    def Clean_up_own(self):
        '''
        Clean up our own software, use before building and before committing to git.
        '''
        #sciNDG4
        command = [self.third_party_support+'/bin/scons','-c']
        subprocess.call(command, cwd = self.exec_base+'/sciNDG4')

        #digi_ND
        subprocess.call('make clean', shell=True, cwd = self.exec_base+'/digi_ND') 
        command = self.exec_base+'/digi_ND/cleanup.sh'
        subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/digi_ND')
        #mind_rec
        subprocess.call('make clean', shell=True, cwd = self.exec_base+'/mind_rec') 

        command = self.exec_base+'/mind_rec/cleanup.sh'
        subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/mind_rec')

    def Config_and_build_own(self):
        '''
        Build all the own libraries that are used for SaRoMaN
        '''
        #digi_ND
        #run configure and autogen in that context.
        #command = self.exec_base+'/digi_ND/autogen.sh'
        #print command
        #subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/digi_ND')
        #subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/digi_ND')
        #command = self.exec_base+'/digi_ND/configure'
        #print command
        #subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/digi_ND')
        #subprocess.call('make', shell=True, cwd = self.exec_base+'/digi_ND')

        #mind_rec
        #run configure and autogen in that context.
        command = self.exec_base+'/mind_rec/autogen.sh'
        print command
        subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/mind_rec')
        subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/mind_rec')

        command = self.exec_base+'/mind_rec/configure'
        print command
        subprocess.call('bash %s' %command, shell=True, cwd = self.exec_base+'/mind_rec')
        subprocess.call('make', shell=True, cwd = self.exec_base+'/mind_rec')    
        
        #sciNDG4
        #command = [self.third_party_support+'/bin/scons']
        #print subprocess.list2cmdline(command)
        #subprocess.call(command, cwd = self.exec_base+'/sciNDG4', env=os.environ)
    '''        
    def Create_folder_structure(self,name,ending):
        OutBase = os.path.join(self.out_base, name+'_out')
        OutConfig = os.path.join(self.out_base, name+ending)
        OutDir = os.path.join(OutBase,'nd_'+self.part+self.inttype)
        ConfigDir = os.path.join(OutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.Check_make_dir(OutBase)
        self.Check_make_dir(OutConfig)
        self.Check_make_dir(OutDir)
        self.Check_make_dir(ConfigDir)
    '''

    def Download_config_and_build_third_party(self):
        '''
        Handles the third party software
        '''
        self.handle_third_party.Download_and_install_genie()
        self.handle_third_party.Download_and_install_geant()
        self.handle_third_party.Download_and_install_depencencies_digi()
        self.handle_third_party.Download_and_install_depencencies_rec()
        self.handle_third_party.Download_and_install_scons()

    def Handle_commandline_input(self,argv):
        '''
        Handles commandline flags, CIO are implemented, if there are no flags the code is run depending on how the variables are setup in __init__
        '''
        if argv==[]:
            if self.need_third_party_install:
                self.handle_third_party.Download_and_install_genie_depencencies()
            self.Set_environment()
            if self.need_third_party_install:
                self.Download_config_and_build_third_party()
            if self.need_own_install:
                self.Clean_up_own()
                self.Config_and_build_own()
            if self.generate_field_map:
                self.Generate_field_map()
            if self.parse_gdml:
                self.xml_parser.Parse_file()

            self.Run_genie()
            self.Run_simulation()
            self.Run_digitization()
            self.Run_reconstruction()

        else:
            try:
                opts, args = getopt.getopt(argv,"CIO")
            except getopt.GetoptError:
                print 'saronman.py -C to clean, -I to install, and -O to build own, to run when setups is ok have no flags'
                sys.exit(2)
            for opt, arg in opts:
                if opt == '-C':
                    self.Set_environment()
                    self.Clean_up_own()
                if opt == '-I':
                    self.handle_third_party.Download_and_install_genie_depencencies()
                    self.Set_environment()
                    self.Download_config_and_build_third_party()
                    self.Clean_up_own()
                    self.Config_and_build_own()
                if opt== '-O':
                    self.Set_environment()
                    #self.Clean_up_own()
                    self.Config_and_build_own()

    def Print_file(self,filename,data):
        '''
        Print data to file filename.
        '''
        outfile = open(filename,'w+')
        outfile.write(data)
        outfile.close()

    def Print_outdata_file(self,filename,command):
        '''
        Print the command to stdout then print call output to file filename.
        '''
        start = time.time()
        print subprocess.list2cmdline(command)
        outfile = open(filename,'w+')
        subprocess.call(command, stdout=outfile)
        outfile.close()
        elapsed = (time.time()-start)
        print 'Time to run process: %s seconds' % elapsed

    def Shell_source(self, script):
        '''
        Run script as source and update environment accordingly
        '''
        pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE,shell=True)
        output = pipe.communicate()[0]
        #print output
        env = dict((line.split("=",1) for line in output.splitlines()))
        os.environ.update(env)

    def Set_environment(self):
        '''
        Set up the environment required to run.
        '''
        genie_support_ext = self.third_party_support + "/genie_source/src/scripts/build/ext"
        
        #os.environ['GENIE_SUPPORT_EXT'] = genie_support_ext
        os.environ['THIRD_PARTY_SUPPORT'] = self.third_party_support

        #ROOT
        self.Shell_source(self.third_party_support + "/root/bin/thisroot.sh")
        #print os.environ

        #GENIE
        os.environ['GENIE'] = self.third_party_support + "/genie2.8.6"
        os.environ['GENIE_INCDIR']=self.third_party_support + "/genie-install/include/GENIE"
        os.environ['GENIE_LIBDIR']=self.third_party_support + "/genie-install/lib"
        os.environ['PATH']+= os.pathsep + self.third_party_support + "/genie-install/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/genie-install/lib"
        # Now add the supporting libraries
        # PYTHIA
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/pythia/v6_428/lib"
        # log4cpp
        #os.environ['PATH']+= os.pathsep + self.third_party_support+ "/bin"
        #os.environ['LD_LIBRARY_PATH']+= os.pathsep+self.third_party_support + "/libs"
        os.environ['PATH']+= os.pathsep + self.third_party_support+ "/log4cpp-install/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep+self.third_party_support + "/log4cpp-install/lib"
        # LHAPDF
        #os.environ['LHAPATH']=self.third_party_support +"/genie_source/data/evgen/pdfs"
        os.environ['LHAPATH']=self.third_party_support +"/genie2.8.6/data/evgen/pdfs"
        #os.environ['LD_LIBRARY_PATH']+= os.pathsep + genie_support_ext + "/v5_7_0/stage/lib"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/lhapdf-5.9.1-install/lib"
        # CLHEP
        clhep_base_dir = self.third_party_support + "/install/clhep-2.1.4.1"
        os.environ['PATH']+= os.pathsep + clhep_base_dir + "/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + clhep_base_dir + "/lib"
        # GEANT4
        version = "Geant4-10.0.1"
        g4install = self.third_party_support + "/install"
        os.environ['G4INSTALL']=g4install
        os.environ['G4WORKDIR']=g4install
        os.environ['G4SYSTEM']="Linux-g++"
        os.environ['G4INCLUDE']=g4install + "/include/Geant4"
        os.environ['G4LIB']=g4install + "/lib64"
        os.environ['G4LEVELGAMMADATA']=g4install + "/share/"+version+"/data/PhotonEvaporation3.0"
        os.environ['G4RADIOACTIVEDATA']=g4install + "/share/"+version+"/data/RadiativeDecay4.0"
        os.environ['G4LEDATA']=g4install + "/share/"+version+"/data/G4EMLOW6.35"
        os.environ['G4NEUTRONHPDATA']=g4install + "/share/"+version+"/data/G4NDL4.4"
        os.environ['G4ABLADATA']=g4install + "/share/"+version+"/data/G4ABLA3.0"
        os.environ['G4SAIDXSDATA']=g4install + "/share/"+version+"/data/G4SAIDDATA1.1"
        os.environ['G4ENSDFSTATE']=g4install + "/share/"+version+"/data/G4ENSDFSTATE1.0"
        os.environ['G4NEUTRONXSDATA']=g4install + "/share/"+version+"/data/G4NEUTRONXS1.4"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + g4install + "/lib64"
        # Support packages
        # CRY simulation information 
        os.environ['CRYPATH']=self.third_party_support + "/cry_v1.7/src"
        os.environ['CRY_INCDIR']=self.third_party_support + "/cry_v1.7/src"
        os.environ['CRY_LIBDIR']=self.third_party_support + "/cry_v1.7/lib"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/cry_v1.7/lib"
        # BHEP
        os.environ['BHEP_LIB']=self.third_party_support + "/bhep-install/lib"
        os.environ['BHEP_PATH']=self.third_party_support + "/bhep-install/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/bhep-install/lib"
        os.environ['PATH']+= os.pathsep + self.third_party_support + "/bhep-install/bin"

        #recpack
        os.environ['PATH']+= os.pathsep + self.third_party_support + "/recpack-install/bin"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/recpack-install/lib"

        #xerces
        os.environ['XERCES_LIBDIR']= self.third_party_support + "/install/lib"
        os.environ['XERCES_INCDIR']= self.third_party_support + "/install/include"
        os.environ['LD_LIBRARY_PATH']+= os.pathsep + self.third_party_support + "/install/lib"

        #Expat
        os.environ['EXPAT_LIBDIR']= self.third_party_support + "/install/lib"
        os.environ['EXPAT_INCDIR']= self.third_party_support + "/install/include"

        #print os.environ
        
#
#    def set_MIND_parameters(self, Type, xdim, ydim, zdim, \
#                           FeThickness, SciThickness, AlThickness,Gap):

    def Submit_run(self, seed, Nevts, pid, inttype, BField):
        '''
        Unused at the moment
        '''
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

        self.Run_genie()
        self.Run_simulation()
        self.Run_digitization()
        self.Run_reconstruction()

    def Run_genie(self):

        if(self.GenerationMode != 'SINGLE_PARTICLE'):
            
            genieOutBase = os.path.join(self.out_base, 'genie_samples')
            genieOutDir = os.path.join(genieOutBase,'nd_'+self.part+self.inttype)
            genieOutLog = os.path.join(genieOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
            FeTargetCode="1000260560"
            C12HTargetCode='1000060120[0.922582],1000010010[0.077418]'
            equation ='1/(0.0069*x) +1'
            vE ='0.01,4.0'

            # Check if the target directories exist if not create it
            self.Check_make_dir(genieOutBase)
            self.Check_make_dir(genieOutDir)

            command = ['gevgen','-r',str(self.seed),'-n',str(self.Nevts),'-p',str(self.pid),'-t',FeTargetCode,
                             '-e',vE,'-f',equation]
            command2= ['gevgen','-r',self.ASeed,'-n',str(self.Nevts),'-p',str(self.pid),'-t',C12HTargetCode,
                             '-e',vE,'-f',equation]

            if self.inttype != 'All':
                command += ['--event-generator-list',self.inttype]
                command2 += ['--event-generator-list',self.inttype]

                command += ['--seed',str(self.seed),'--cross-sections',self.scripts_dir+'/xsec_Fe56_splines.xml']
                command2 += ['--seed',self.ASeed,'--cross-sections',self.scripts_dir+'/xsec_C12+H1_splines.xml']

                self.Print_outdata_file(genieOutLog,command)

                geniefile='gntp.'+str(self.seed)+'.ghep.root'
                geniedest=genieOutDir+'/ev0_'+str(self.seed)+'_'+str(self.pid)+'_'+str(FeTargetCode)+'_'+str(self.Nevts)+'.root'
                shutil.move(geniefile, geniedest)

                self.Print_outdata_file(genieOutLog,command2)
        
                geniefile='gntp.'+self.ASeed+'.ghep.root'
                geniedest=genieOutDir+'/ev0_'+self.ASeed+'_'+str(self.pid)+'_'+str(C12HTargetCode)+'_'+str(self.Nevts)+'.root'
                shutil.move(geniefile, geniedest)
            
    def Run_simulation(self):
        #Used to remove the .gdml file which is created by Genie each run.
        gdmlFile = os.path.join(self.exec_base,'MIND_detector.gdml')

        if os.path.exists(gdmlFile):
            os.remove(gdmlFile)

        mindG4OutBase = os.path.join(self.out_base, 'G4_out')
        mindG4OutConfig = os.path.join(self.out_base, 'G4_config')
        mindG4OutDir = os.path.join(mindG4OutBase,'nd_'+self.part+self.inttype)
        mindG4ConfigDir = os.path.join(mindG4OutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.Check_make_dir(mindG4OutBase)
        self.Check_make_dir(mindG4OutConfig)
        self.Check_make_dir(mindG4OutDir)
        self.Check_make_dir(mindG4ConfigDir)

        mindG4config = os.path.join(mindG4ConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.config')
        self.print_config.print_mindG4_config(vars(self),mindG4config)

        mindG4OutLog = os.path.join(mindG4OutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')

        self.Shell_source(self.third_party_support+"/install/bin/geant4.sh")

        #Uncomment to run visually
        #command = [self.exec_base+'/sciNDG4/mindG4','-v',mindG4config]
        #subprocess.call(command)
        
        command = [self.exec_base+'/sciNDG4/mindG4',mindG4config]
        self.Print_outdata_file(mindG4OutLog,command)

    def Run_digitization(self):
        digiOutBase = os.path.join(self.out_base, 'digi_out')
        digiOutConfig = os.path.join(self.out_base, 'digi_param')
        digiOutDir = os.path.join(digiOutBase,'nd_'+self.part+self.inttype)
        digiConfigDir = os.path.join(digiOutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.Check_make_dir(digiOutBase)
        self.Check_make_dir(digiOutConfig)
        self.Check_make_dir(digiOutDir)
        self.Check_make_dir(digiConfigDir)

        digiConfig = os.path.join(digiConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.digi.param')
        self.print_config.print_digi_config(vars(self),digiConfig)

        digiOutLog = os.path.join(digiOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
        command = [self.exec_base+"/digi_ND/examples/simple_smear",digiConfig]
        self.Print_outdata_file(digiOutLog,command)

    def Run_reconstruction(self):
        recOutBase = os.path.join(self.out_base, 'rec_out')
        recOutConfig = os.path.join(self.out_base, 'rec_param')
        recLikelihoods = os.path.join(self.out_base, 'likelihoods')
        recOutDir = os.path.join(recOutBase,'nd_'+self.part+self.inttype)
        recConfigDir = os.path.join(recOutConfig,'nd_'+self.part+self.inttype)

        # Check if the target directories exist if not create it
        self.Check_make_dir(recOutBase)
        self.Check_make_dir(recOutConfig)
        self.Check_make_dir(recLikelihoods)
        self.Check_make_dir(recOutDir)
        self.Check_make_dir(recConfigDir)

        recConfig = os.path.join(recConfigDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.rec.param')
        self.print_config.print_rec_config(vars(self),recConfig)

        recOutLog = os.path.join(recOutDir,'nd_'+self.part+self.inttype+'_'+str(self.seed)+'.log')
        command = [self.exec_base+"/mind_rec/examples/fit_tracks",recConfig, str(self.Nevts)]
        self.Print_outdata_file(recOutLog,command)
        print 'Completed reconstruction'

#######################################################################################################################
#File specific functions
#######################################################################################################################

if __name__ == "__main__":
    s=saroman()
    s.Handle_commandline_input(sys.argv[1:])

#######################################################################################################################
