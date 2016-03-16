import pickle

class Run():
    """
    create a run database
    """
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

    def save(self, name="db"):
        with open(name, 'wb') as output:
            pickle.dump(self, output,-1)

    def load(self, name="db"):
        with open(name, 'rb') as input:
            tmp=pickle.load(input)
        self.runs = tmp.runs
