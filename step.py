import grid
import part

class Step:
	def __init__(self,number,folder):
		"""			
		Create a step object from EMMA output 					
		"""
	
		self.number=number
		self.folder=folder
		
		self.grid=grid.Grid(number,folder)		
		self.part=part.Part(number,folder,0)
		self.star=part.Part(number,folder,1)			
