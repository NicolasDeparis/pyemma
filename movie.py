# coding: utf-8

import os
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

class MovieField:
    """
    movie physical field
    """

    def __init__(self,path):
        self._path=path

    def _read(self,file):
        """
        read EMMA movie snap
        """
        with open(file, "rb") as f:
            x = np.fromfile(f, dtype=np.int32  ,count=1)[0]
            y = np.fromfile(f, dtype=np.int32  ,count=1)[0]
            a = np.fromfile(f, dtype=np.float32  ,count=1)[0]
            data=np.fromfile(f, dtype=np.float32  ,count=x*y).reshape(y,x)
            return x,y,a,data

#     def draw_one(self,file_path):
#         plt.figure(figsize=(10,10))
#         x,y,a,data =readMovie(file_path)
#         fig=plt.imshow(np.log10(data), interpolation='nearest',extent=(0,1,0,1),origin='lower')
#         plt.axis("off")
#         fig.axes.get_xaxis().set_visible(False)
#         fig.axes.get_yaxis().set_visible(False)

    def draw_all(self):
        """
        draw all movie snap into img/ folder
        """

        #make img dir
        try:
            import os
            os.mkdir("%simg/"%(self._path))
        except FileExistsError:
            pass

        #getting the file list
        files = np.sort(os.listdir(self._path))

        #draw each file
        for i,file in enumerate(files):
                if file != "img":
                    plt.clf()
                    file_path = "%s%s"%(self._path,file)
                    print(file_path)

                    x,y,a,data = self._read(file_path)

                    data = np.log10(data)
                    data = (data-np.min(data)) / (np.max(data)-np.min(data))

                    img = Image.fromarray(np.uint8(plt.cm.hot(data)*255))
                    img.save("%simg/%07d.png"%(self._path,i))

class Movie:
    """
    Global movie class : scan the movie folder
    """
    def __init__(self,path):        
        for field in os.listdir(path):
            setattr(self,field,MovieField("%s%s/"%(path,field)))
