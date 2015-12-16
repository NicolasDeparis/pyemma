import os
import numpy as np
import field

class Part : 
	def __init__(self, number,folder, isStar):	
		"""			
			Create a part object from EMMA output
		"""
		
		self.number=number		
		self.isStar	= isStar
				
		if isStar:
			self.type="star_"
		else:
			self.type="part_"
			
		path = "%s%05d/"%(folder,number)
		for cur_folder in  os.listdir(path):
			if "part_" in cur_folder:
				key=cur_folder.strip("part_").replace(".","_")
				val= field.Field(folder,number,cur_folder)
				setattr(self,key,val)
