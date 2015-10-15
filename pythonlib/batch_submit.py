import os
import subprocess

class batch_submit:

	def __init__(self):
		self.seed = 1
		self.home = os.getenv("HOME")
		self.working_dir = os.path.join(self.home,'SaRoMaN')
		self.instances = 100
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
#PBS -e %(working_dir)s/saroman.py_%(seed)d.err
#PBS -o %(working_dir)s/saroman.py_%(seed)d.out

python %(working_dir)s/saroman.py 
'''% vars(self)

		#--configuration_file %(configfile)s
        
		self.print_file(self.batchfile, batch_data)
        
	def generate_submission(self):
        
		for self.seed in range(1,self.instances):
			#self.configfile = self.working_dir + "/maus_config_"+str(self.seed)
			#self.print_config()
			self.print_batch_submission()
			cmd = ['qsub', self.batchfile]
			q = subprocess.Popen(cmd)


if __name__ == "__main__":

    b = batch_submit()
    b.generate_submit()
