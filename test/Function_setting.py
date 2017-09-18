import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy.spatial.distance import cdist
from matplotlib import pyplot as plt

class Parameters:
    def __init__(self, molecule_number, atomic_number, method, mode ):
        self.molecule_number = molecule_number
        self.atomic_number = atomic_number
        self.method = method
        self.mode = mode
        self.path = None 
        self.dimension = None
        self.component = None
        self.path = None
        self.testingnumber = None
        self.point_interval = None
        self.startneighbor = None
        self.endneighbor = None
        self.neighbor_interval =None

    def set_output_path(self, path = '/home/sky8/structure/'):
        self.path = path

    def getdata(self, filepath):
        if self.mode == "Internal":
            coordi=np.loadtxt(filepath)
            co = coordi[0 : self.molecule_number]
        elif self.mode == 'Cartisian':
            coordi = np.loadtxt(filepath)
            co = coordi[0 : self.atomic_number*self.molecule_number]
            co = co.reshape(self.molecule_number, 3*self.atomic_number)
        self.dimension = co.shape[1]
        return co

    def set_testParameter(self, component, testingnumber, point_interval, startneighbor, endneighbor, neighbor_interval=1):
        self.component = component
        self.testingnumber = testingnumber
        self.point_interval = point_interval
        self.startneighbor = startneighbor
        self.endneighbor = endneighbor
        self.neighbor_interval = neighbor_interval


def Manifold_convergence_test(func, co, par):
    fig = plt.figure()
    for j in range(par.startneighbor, par.endneighbor, par.neighbor_interval):
        output = func(n_components=par.component, n_neighbors=j).fit_transform(co)
        print("standrdok")
        output = output[0:100]
        standard = cdist(output, output, 'euclidean')
        convergence = []
        F_lis_1_15 = []
        F_lis_2_25 = []
        F_lis_16_39 = []
        F_lis_51_92 = []
        F_lis_61_86 = []
        F_lis = []
        F_lis = [F_lis_1_15, F_lis_2_25, F_lis_16_39, F_lis_51_92, F_lis_61_86]
        F_name = ['1_15', '2_25', '16_39', '51_92', '61_86']
        for i in range((par.testingnumber-100)//par.point_interval):
            coordinate = co[0:10 ** 2 + i * par.point_interval ]
            output = func(n_components=par.component, n_neighbors=j).fit_transform(coordinate)
            output = output[0:100]
            F = cdist(output, output, 'euclidean')
            for k in range(5):
                s = F_name[k].split('_')
                F_lis[k].append(F[int(s[0])][int(s[1])])
            di = np.sum(np.absolute(F - standard))
            convergence.append(di)
            if i == 10:
                print(10)
            if i == 200:
                print(200)
            if i == 400:
                print(400)
            if i == 800:
                print(800)
            if i == 900:
                print(900)
        print("finish part")


        for i in range(5):
            plt.scatter( list(range(100,par.testingnumber,par.point_interval)),F_lis[i], s=1, alpha=0.5,figure = fig)
            plt.savefig("%s%s_%s_%s_neighbor%d_dimension%d_%d.png" % (par.path,par.mode,par.method,F_name[i],j,par.dimension,par.component),dpi=600)
            plt.clf()

        plt.scatter( list(range(100,par.testingnumber,par.point_interval)),convergence, s=1, alpha=0.5,figure = fig)
        plt.savefig("%s%s_%s_convergence_neighbor%d_dimension%d_%d.png" % (par.path,par.mode,par.method,j,par.dimension,par.component),dpi=600)
        plt.clf()
