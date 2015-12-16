import os
import hop
import field

class Step:
	def __init__(self,number,folder):
		"""			
		Create a step object			
		"""
	
		self.number=number
		self.folder=folder
				
		self.part=Fields(number,folder,0)
		self.star=Fields(number,folder,1)
		self.grid=Fields(number,folder,2)
		
		self.hop=hop.Hop(number,folder)


class Fields:
	def __init__(self, number,folder, sets_type):
		"""			
			Create a set of fields object
		"""
		
		self.number=number		
		self.sets_type	= sets_type
				
		if sets_type==0:
			self.type="star_"
		elif sets_type==1:
			self.type="part_"
		elif sets_type==2:
			self.type="grid_"
			
		path = "%s%05d/"%(folder,number)
		for cur_folder in  os.listdir(path):
			if self.type in cur_folder:
				key=cur_folder[5:].replace(".","_")
				val= field.Field(folder,number,cur_folder)
				setattr(self,key,val)
