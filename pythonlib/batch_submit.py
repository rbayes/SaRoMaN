import os
import subprocess
from time import sleep

class batch_submit:

	def __init__(self):
		self.seed = 1
		self.home = os.getenv("HOME")
		self.working_dir = os.path.join(self.home,'SaRoMaN')
		self.out_dir = os.path.join(self.home,'batch2')
		self.instances = 20
		#self.configfile =
		self.batchfile = self.working_dir + "/run_"+str(self.seed)+".job"


	def print_file(self, filename, data):
		outfile = open(filename, 'w+')
		outfile.write(data)
		outfile.close()

	def print_batch_submission(self):
		batch_data = '''
#PBS -N saroman.py_%(seed)d
#PBS -q medium6
#PBS -l walltime=5:59:00
#PBS -e %(out_dir)s/saroman.py_%(seed)d.err
#PBS -o %(out_dir)s/saroman.py_%(seed)d.out

python %(working_dir)s/saroman.py -B %(seed)d
'''% vars(self)

		#--configuration_file %(configfile)s
        
		self.print_file(self.batchfile, batch_data)
        
	def generate_submission(self):
        
		for self.seed in range(0,self.instances):
			#self.configfile = self.working_dir + "/maus_config_"+str(self.seed)
			#self.print_config()
			self.print_batch_submission()
			cmd = ['qsub', self.batchfile]
			q = subprocess.Popen(cmd)
			sleep(3) #Wait 3 sec to make sure that all processes are submitted correctly.


if __name__ == "__main__":
	b = batch_submit()
	b.generate_submission()

	#REMOVE JOBS FROM QSTAT
	#for x in range(1795878,1796077):
		#string = str(x) + '.offler'
		#cmd = ['qdel', string]
		#subprocess.Popen(cmd)
		#sleep(0.5)
	
