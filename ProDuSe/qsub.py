import os
import subprocess
import re
import random

class Qsub:

	previous_job_ids=[]
	previous_job_names=[]
	previous_job_commands={}

	def __init__(self, job_name, virtual_memory='4G', wait_on_previous_job=False, source_directory=""):

		# Job Attributes		
		self.job_name=job_name
		self.virtual_memory=virtual_memory
		self.wait_on_previous_job=wait_on_previous_job
                self.source_directory=source_directory


	def change(self, job_name='', virtual_memory='', wait_on_previous_job=''):
		
		# If you would like to, change qsub specific attributes
		# May want to add more to this list
		if not job_name == '':
			self.job_name=job_name
		if not virtual_memory=='':
			self.virtual_memory=virtual_memory
		if not wait_on_previous_job=='':
			self.wait_on_previous_job=wait_on_previous_job	


	def run(self, commands):

		# Store commands in a temporary .sh file:
		if not os.path.exists('./tmp'):
			os.makedirs('./tmp')

		tmp_file = '/'.join(['.', 'tmp', self.job_name])
		tmp_file = '.'.join([tmp_file, 'sh'])

		# Make sure the tmp file does not exist, if it does, generate a random number and append it to the end
		prev_tmp_file = tmp_file
		while os.path.isfile(tmp_file):
			tmp_file = '.'.join([prev_tmp_file[:-3], str(random.randint(1000000000, 9999999999)), 'sh'])

		# Write the commands to this new tmp file
		tmp_handler = open(tmp_file, 'w')
                if not self.source_directory == "":
                    tmp_handler.write(' '.join(["source", self.source_directory, ";"]))
		tmp_handler.write(' '.join(commands))
		tmp_handler.close()

		# Change this file to be an executable
		os.chmod(tmp_file, 777)

		# Create tmp command
		tmp_command = '/'.join(['.', tmp_file])

                subprocess.call(['chmod', '777', tmp_command])

		hold_jid = ""
		if self.wait_on_previous_job:
			hold_jid = ','.join(self.previous_job_ids)

		qsub_command=[
			'qsub',			'-cwd',
			'-b',			'y',
			'-N',			self.job_name
                        ]
                if not hold_jid == "":
                        qsub_command = qsub_command + ['-hold_jid', hold_jid]
                qsub_command = qsub_command + [
			'-l',			''.join(['h_vmem=',self.virtual_memory]),
			tmp_command
			]


		# Call the Qsub while storing the stdout
		qsub_output = subprocess.Popen(qsub_command, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
                err, std = qsub_output.communicate()


		# Grab the id from the stdout
		id = re.search("[0-9]+", err).group()


		# Store information about this job
		self.previous_job_ids.append(id)
		self.previous_job_names.append(self.job_name)
		self.previous_job_commands[id] = commands
