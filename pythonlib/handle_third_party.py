#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 27/10-2015
#Created on 13/10-2015
#######################################################################################################################
#General python import
#######################################################################################################################
import os
import subprocess
#######################################################################################################################
#Importing own python files
#######################################################################################################################

#######################################################################################################################
#Class generation
#######################################################################################################################
class handle_third_party:

	def __init__(self,exec_base,third_party_support):
		self.FNULL = open(os.devnull, 'w') # So that output is not printed to console
		self.exec_base = exec_base
		self.third_party_support = third_party_support

################################################
	  	#def Clean_up_install_directory(self):
		#Delete all tar and extra files which are no longer needed.

	def Download_and_install_genie_depencencies(self):
		self.Check_make_dir(self.third_party_support)            
		#Log4cpp
		print 'Installing Log4cpp...'
		command = ['cvs','-d',':pserver:anonymous@log4cpp.cvs.sourceforge.net:/cvsroot/log4cpp','login']
		print 'Please press enter when the password request comes\n'
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		command = ['cvs','-d',':pserver:anonymous@log4cpp.cvs.sourceforge.net:/cvsroot/log4cpp','-z3','co','log4cpp']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		command = self.third_party_support+'/log4cpp/autogen.sh'
		subprocess.call('bash %s' %command, shell=True, cwd = self.third_party_support+'/log4cpp',stdout=self.FNULL)
		command = self.third_party_support+'/log4cpp/configure --prefix='+self.third_party_support+'/log4cpp-install'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/log4cpp',stdout=self.FNULL)
		subprocess.call('gmake', shell=True, cwd = self.third_party_support+'/log4cpp',stdout=self.FNULL)
		subprocess.call('gmake install', shell=True, cwd = self.third_party_support+'/log4cpp',stdout=self.FNULL)

		print 'Log4cpp was installed successfully'

		#LHAPDF
		print 'Installing LHAPDF...'
		command = ['wget','http://www.hepforge.org/archive/lhapdf/lhapdf-5.9.1.tar.gz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		command = ['tar','xzvf','lhapdf-5.9.1.tar.gz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		command = self.third_party_support+'/lhapdf-5.9.1/configure --prefix='+self.third_party_support+'/lhapdf-5.9.1-install'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/lhapdf-5.9.1',stdout=self.FNULL)
		subprocess.call('gmake', shell=True, cwd = self.third_party_support+'/lhapdf-5.9.1',stdout=self.FNULL)
		subprocess.call('gmake install', shell=True, cwd = self.third_party_support+'/lhapdf-5.9.1',stdout=self.FNULL)

		print 'LHAPDF was installed successfully'

		#Pythia6
		print 'Installing Pythia6...'
		self.Check_make_dir(self.third_party_support + '/pythia')
		command = ['wget','http://home.fnal.gov/~rhatcher/build_pythia6.sh']
		#command = ['scp',self.exec_base+'/build_pythia6.sh', self.third_party_support + '/pythia']
		subprocess.call(command,cwd = self.third_party_support,stdout=self.FNULL)
		self.Shell_source_no_environ(self.third_party_support + '/pythia'+'/build_pythia6.sh',self.third_party_support + '/pythia')
		
		print 'Pythia6 was installed successfully'

		#ROOT
		print 'Installing ROOT...'
		command = ['git','clone','http://root.cern.ch/git/root.git']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)       
		#command = ['git','tag','-l'] #Is this just visual?
		#subprocess.call(command, cwd = self.third_party_support + '/root')
		command = ['git','checkout','-b','v5-34-34','v5-34-34']
		subprocess.call(command, cwd = self.third_party_support + '/root',stdout=self.FNULL)
		command = self.third_party_support+'/root/configure linuxx8664gcc --enable-pythia6 --with-pythia6-libdir='+self.third_party_support+'/pythia/v6_428/lib --enable-python'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/root',stdout=self.FNULL)
		subprocess.call('make', shell=True, cwd = self.third_party_support+'/root',stdout=self.FNULL)
		
		print 'ROOT was installed successfully'

		#self.Shell_source(self.third_party_support + '/root/bin'+'/thisroot.sh',self.third_party_support + '/root')

	def Download_and_install_genie(self):
		#Genie
		print 'Installing Genie...'

		command = ['svn','co','http://genie.hepforge.org/svn/generator/branches/R-2_8_6','genie2.8.6']
		subprocess.call(command, cwd = self.third_party_support, stdout=self.FNULL)

		#os.environ['GENIE'] = self.third_party_support + "/genie2.8.6" #Is also done in saroman.py
		command = [self.third_party_support+'/genie2.8.6/configure','--prefix='+self.third_party_support+'/genie-install',
		'--enable-flux-drivers','--enable-geom-drivers','--enable-test','--enable-mueloss','--with-log4cpp-inc='+self.third_party_support+'/log4cpp-install/include',
		' --with-log4cpp-lib=' +self.third_party_support+'/log4cpp-install/lib','--with-pythia6-lib=' +self.third_party_support+'/pythia/v6_428/lib', 
		'--with-lhapdf-inc=' +self.third_party_support+'/lhapdf-5.9.1-install/include','--with-lhapdf-lib=' +self.third_party_support+'/lhapdf-5.9.1-install/lib',
		'--enable-numi','--enable-atmo','--enable-eventserver']

		subprocess.call(command, cwd = self.third_party_support+'/genie2.8.6',stdout=self.FNULL)
		#subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/genie2.8.6')

		subprocess.call('make', shell=True, cwd = self.third_party_support+'/genie2.8.6',stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support+'/genie2.8.6',stdout=self.FNULL)

		print 'Genie was installed successfully'

	def Download_and_install_geant(self):

		self.Check_make_dir(self.third_party_support + '/source')
		self.Check_make_dir(self.third_party_support + '/build')
		self.Check_make_dir(self.third_party_support + '/install')

		# Bash 28 xerces-c-3.1.2               
		self.Mice_script_install_emulator('xerces-c-3.1.2','xerces-c-3.1.2.tar.gz','http://apache.mirror.anlx.net//xerces/c/3/sources/xerces-c-3.1.2.tar.gz')
		# Bash 29 expat-2.1.0   
		self.Mice_script_install_emulator('expat-2.1.0','expat-2.1.0.tar.gz','http://downloads.sourceforge.net/expat/expat-2.1.0.tar.gz')

		#Bhep v3r0p0
		print 'Installing bhep...'
		command = ['svn','export','svn://next.ific.uv.es/svn/bhep/tags/v3r0p0','bhep']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)

		#self.Shell_source(self.third_party_support + '/root/bin'+'/thisroot.sh',self.third_party_support + '/root') #Make sure rootdir is pointing in the right place

		command = self.third_party_support+'/bhep/bhep3/autogen.sh'
		subprocess.call('bash %s' %command, shell=True, cwd = self.third_party_support+'/bhep/bhep3',stdout=self.FNULL)

		command = self.third_party_support+'/bhep/bhep3/configure --prefix='+self.third_party_support+'/bhep-install'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/bhep/bhep3',stdout=self.FNULL)
		subprocess.call('make', shell=True, cwd = self.third_party_support+'/bhep/bhep3',stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support+'/bhep/bhep3',stdout=self.FNULL)

		print 'Bhep was installed successfully'
		#Cry 1.7
		print 'Installing cry...'
		command = ['wget','http://nuclear.llnl.gov/simulation/cry_v1.7.tar.gz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		command = ['tar','xzvf','cry_v1.7.tar.gz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)
		subprocess.call('make', shell=True, cwd = self.third_party_support+'/cry_v1.7',stdout=self.FNULL)

		print 'Cry was installed successfully'

		#Geant4
		print 'Installing Geant4...'
		directory = 'geant4.10.00.p01'
		filename = directory + '.tar.gz'
		url = 'http://www.geant4.org/geant4/support/source/' + filename
		
		self.Check_make_dir(self.third_party_support + '/source/'+directory)
		self.Check_make_dir(self.third_party_support + '/build/'+directory)

		command = ['wget','--directory-prefix='+self.third_party_support+'/source',url]
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL) 

		command = ['tar','xzvf',self.third_party_support+'/source/'+filename,'-C',self.third_party_support+'/source']
		subprocess.call(command, cwd = self.third_party_support +'/source',stdout=self.FNULL) 

		command = ['tar','xzvf',self.third_party_support+'/source/'+filename]
		subprocess.call(command, cwd = self.third_party_support +'/build/'+directory,stdout=self.FNULL) 
		
		command = ['cmake','-DCMAKE_INSTALL_PREFIX='+self.third_party_support+'/install','-DGEANT4_INSTALL_DATA=ON','-DGEANT4_INSTALL_DATADIR'
		'-DGEANT4_USE_GDML=ON',
		'-DXERCESC_ROOT_DIR='+self.third_party_support+'/install','-DGEANT4_USE_OPENGL_X11=ON','-DEXPAT_INCLUDE_DIR='+self.third_party_support+'/install/include',
		'-DEXPAT_LIBRARY='+self.third_party_support+'/install/lib/libexpat.so',self.third_party_support+'/source/'+directory]
		subprocess.call(command, cwd = self.third_party_support +'/build/'+directory,stdout=self.FNULL)

		subprocess.call('make', shell=True, cwd = self.third_party_support +'/build/'+directory,stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support +'/build/'+directory,stdout=self.FNULL)       

		print 'Geant4 was installed successfully'    

	def Download_and_install_depencencies_digi(self):
		#CLHEP-2.1.4.1
		print 'Installing CLHEP...'
		command = ['wget','--directory-prefix='+self.third_party_support+'/source','http://proj-clhep.web.cern.ch/proj-clhep/DISTRIBUTION/tarFiles/clhep-2.1.4.1.tgz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)

		command = ['tar','xvfz',self.third_party_support+'/source/clhep-2.1.4.1.tgz','-C',self.third_party_support+'/build']
		subprocess.call(command, cwd = self.third_party_support +'/source',stdout=self.FNULL)

		command = ['cmake','-DCMAKE_INSTALL_PREFIX='+self.third_party_support+'/install/'+'clhep-2.1.4.1',self.third_party_support+'/build/2.1.4.1/CLHEP']
		subprocess.call(command, cwd = self.third_party_support +'/install',stdout=self.FNULL)

		subprocess.call('make', shell=True, cwd = self.third_party_support+'/install',stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support+'/install',stdout=self.FNULL)

		print 'CLHEP was installed successfully'

	def Download_and_install_depencencies_rec(self):
		#RECPACK-v1r2p0
		print 'Installing RECPACK...'
		#Set up dependencies correctly, just needed for testing during writing!  
		#clhep_base_dir = self.third_party_support + "/install/clhep-2.1.4.1"
		#os.environ['PATH']+= os.pathsep + clhep_base_dir + "/bin"
		#os.environ['LD_LIBRARY_PATH']= os.pathsep + clhep_base_dir + "/lib" 
		#self.Shell_source(self.third_party_support + '/root/bin'+'/thisroot.sh',self.third_party_support + '/root')

		command = ['svn','export','svn://next.ific.uv.es/svn/recpack/alpha/trunk/recpack','recpack-source']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)

		command = self.third_party_support+'/recpack-source/autogen.sh'
		subprocess.call('bash %s' %command, shell=True, cwd = self.third_party_support+'/recpack-source',stdout=self.FNULL)
		command = self.third_party_support+'/recpack-source/autogen.sh'
		subprocess.call('bash %s' %command, shell=True, cwd = self.third_party_support+'/recpack-source',stdout=self.FNULL)

		command = self.third_party_support+'/recpack-source/configure --prefix='+self.third_party_support+'/recpack-install'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/recpack-source',stdout=self.FNULL)
		subprocess.call('make', shell=True, cwd = self.third_party_support+'/recpack-source',stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support+'/recpack-source',stdout=self.FNULL)

		print 'RECPACK was installed successfully'

	def Download_and_install_scons(self):
		print 'Installing SCONS 1.2.0...'
		command = ['wget','http://downloads.sourceforge.net/project/scons/scons/1.2.0/scons-1.2.0.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fscons%2Ffiles%2Fscons%2F1.2.0%2F']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)

		command = ['tar','xvfz',self.third_party_support+'/scons-1.2.0.tar.gz']
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL)

		command = ['python','setup.py','install','--prefix='+self.third_party_support]
		subprocess.call(command, cwd = self.third_party_support +'/scons-1.2.0',stdout=self.FNULL)
		  
	def Mice_script_install_emulator(self,directory,filename,url):
		print ('Installing ' + directory + '...')

		command = ['wget','--directory-prefix='+self.third_party_support+'/source',url]
		subprocess.call(command, cwd = self.third_party_support,stdout=self.FNULL) 

		command = ['tar','xzvf',self.third_party_support+'/source/'+filename,'-C',self.third_party_support+'/build']
		subprocess.call(command, cwd = self.third_party_support +'/source',stdout=self.FNULL) 

		command = self.third_party_support+'/build/'+directory+'/configure --prefix='+self.third_party_support+'/install'
		subprocess.call('bash %s'%command, shell=True, cwd = self.third_party_support+'/build/'+directory,stdout=self.FNULL)

		subprocess.call('make', shell=True, cwd = self.third_party_support+'/build/'+directory,stdout=self.FNULL)
		subprocess.call('make install', shell=True, cwd = self.third_party_support+'/build/'+directory,stdout=self.FNULL)

		print (directory +' Genie was installed successfully')

	def Shell_source(self, script,where):
		'''
		Run script as source and update environment accordingly
		'''
		pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE,shell=True,cwd = where)
		output = pipe.communicate()[0]
		env = dict((line.split("=",1) for line in output.splitlines()))
		os.environ.update(env)
	
	def Shell_source_no_environ(self, script,where):
		'''
		Run script as source and update environment accordingly
		'''
		pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE,shell=True,cwd = where)
		output = pipe.communicate()[0]
		
	def Check_make_dir(self,dirname):
		'''
		Check if the target directories dirname exist, if not create it.
		'''
		if not os.path.exists(dirname):
			os.makedirs(dirname)
		

#######################################################################################################################
#File specific functions
#######################################################################################################################
#if __name__ == "__main__":
#    home = os.getenv("HOME")
#    test = '/data/neutrino05/phallsjo'
#    s=handle_third_party(home+'/SaRoMaN',test+'/test2')
#    s.Download_and_install_genie_dependencies()
#    s.Download_and_install_genie()
#    s.Download_and_install_geant()
#    s.Download_and_install_depencencies_digi()
#    s.Download_and_install_depencencies_rec()
