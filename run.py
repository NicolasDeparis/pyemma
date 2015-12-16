import os

import step
import param

class Run:
	def __init__(self,folder):
		self.folder=folder
		self.data_folder=folder+"data/"
		
		for folder in os.listdir(self.data_folder):
			try:				
				stepnum=int(folder)							
			except ValueError:
				continue
				
			key="step_%05d"%stepnum
			val= step.Step(stepnum, self.data_folder)
			setattr(self,key,val)
		
		self.param=param.Param(self.folder)
