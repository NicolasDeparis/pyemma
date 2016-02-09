
# coding: utf-8

# In[1]:

import copy
import pickle


# In[2]:

class Run():
    def __init__(self):
        pass

class Runset():
    def __init__(self):
        self.runs = []
        
    def add(self,run):
        self.runs.append(run)
        
    def get_description(self):
        for i in range(len(self.runs)):
            print( "%02d %s"%(i,self.runs[i].description) )
    def get_folder(self):
        for i in range(len(self.runs)):
            print( "%02d %s"%(i,self.runs[i].folder) )
        
    def save(self):
        name="db"
        with open(name, 'wb') as output:
            pickle.dump(self, output,-1)
            
    def load(self):
        name="db"
        with open(name, 'rb') as input:
            tmp=pickle.load(input)
        self.runs = tmp.runs


# In[3]:

runset=Runset()
r=Run()


# In[4]:

r.description="run de base -> m2"
r.folder="/home/deparis/curie_data/data/8_8_gather_6/"
r.colors="Red"
r.labels="m2 SN1"
runset.add(copy.copy(r))

r.description="masse d'etoile *8 -> m1"
r.folder="/home/deparis/curie_data/data/8_8_gather_7/"
r.colors="Blue"
r.labels="m1 SN1"
runset.add(copy.copy(r))

r.description="même masse d'étoile mais pas de SN"
r.folder="/home/deparis/curie_data/data/8_8_gather_8/"
r.colors="Green"
r.labels="m2 SN0"
runset.add(copy.copy(r))

r.description="masse d'etoile *8 mais pas de SN"
r.folder="/home/deparis/curie_data/data/8_8_gather_9/"
r.colors="Cyan"
r.labels="m1 SN0"
runset.add(copy.copy(r))

r.description="masse d'étoile /8 -> m3"
r.folder="/home/deparis/curie_data/data/8_8_gather_10/"
r.colors="#cc9900"
r.labels="m3 SN1"
runset.add(copy.copy(r))

r.description="même masse d'étoile mais SN thermique"
r.folder="/home/deparis/curie_data/data/8_8_gather_11/"
r.colors="#cc9900"
r.labels="m2 SN th"
runset.add(copy.copy(r))

r.description="même masse d'étoile mais SN kin simple sans eject"
r.folder="/home/deparis/curie_data/data/8_8_gather_12/"
r.colors="#cc9900"
r.labels="m2 SN th"
runset.add(copy.copy(r))

r.description="même masse d'étoile mais SN kin simple avec eject"
r.folder="/home/deparis/curie_data/data/8_8_gather_13/"
r.colors="Blue"
r.labels="SN kin simple"
runset.add(copy.copy(r))

r.description="même masse d'étoile mais sans tirage de Poisson __ DOUTE"
r.folder="/home/deparis/curie_data/data/8_8_gather_14/"
r.colors="Blue"
r.labels="no random"
runset.add(copy.copy(r))

r.description="tout pareil mais avec flux __ PROBLEME"
r.folder="/home/deparis/curie_data/data/8_8_gather_15/"
r.colors="Blue"
r.labels="no random"
runset.add(copy.copy(r))

r.description="RERUN tout pareil mais avec flux"
r.folder="/home/deparis/curie_data/data/8_8_gather_16/"
r.colors="Blue"
r.labels="no random"
runset.add(copy.copy(r))


# In[5]:

runset.save()


# In[6]:

# runset2=Runset()
# runset2.load()
# runset2.get_desc()

