import numpy as np
from scipy import stats
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

    def set_neighbor(self, component, testingnumber, point_interval, startneighbor=None, endneighbor=None, neighbor_interval=1 ):
        self.component = component
        self.testingnumber = testingnumber
        self.point_interval = point_interval
        self.startneighbor = startneighbor
        self.endneighbor = endneighbor
        self.neighbor_interval = neighbor_interval

    def set_component(self, neighbor, component=None, smallest=None, biggest=None, component_interval=1):
        self.neighbor = neighbor
        self.component = component
        self.smallest = smallest
        self.biggest = biggest
        self.component_interval = component_interval




def Manifold_convergence_test(func, co, par, **kwargs):
    fig = plt.figure()
    for j in range(par.startneighbor, par.endneighbor, par.neighbor_interval):
        output = func(n_components=par.component, n_neighbors=j, **kwargs).fit_transform(co)
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
            output = func(n_components=par.component, n_neighbors=j, **kwargs).fit_transform(coordinate)
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

def manifold_convergence_simple_test(func, co, par, **kwargs):
    fig = plt.figure()
    output = func(n_components=par.component, **kwargs).fit_transform(co)
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
        output = func(n_components=par.component, **kwargs).fit_transform(coordinate)
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
        plt.savefig("%s%s_%s_%s_dimension%d_%d.png" % (par.path,par.mode,par.method,F_name[i],par.dimension,par.component),dpi=600)
        plt.clf()
    plt.scatter( list(range(100,par.testingnumber,par.point_interval)),convergence, s=1, alpha=0.5,figure = fig)
    plt.savefig("%s%s_%s_convergence_dimension%d_%d.png" % (par.path,par.mode,par.method,par.dimension,par.component),dpi=600)
    plt.clf()

def manifold_residual_variance(func, coordinate, par ,**kwargs):
    residual_variance = []
    fig = plt.figure()
    for i in range(par.smallest,par.biggest+1):
        iso = func(n_components=i , n_neighbors=par.neighbor, **kwargs)
        iso.fit(coordinate)
        embedding = iso.embedding_
        Dy = cdist(embedding, embedding , 'euclidean')
        Dm = iso.dist_matrix_
        Dy = Dy.ravel()
        Dm = Dm.ravel()
        slope, intercept, r_value, p_value, std_err = stats.linregress(Dm,Dy)
        R_squared = r_value**2
        residual = 1-R_squared
        residual_variance.append(residual)
    plt.scatter(list(range(par.smallest,par.biggest+1,par.component_interval)),residual_variance, s=10, alpha=0.5,figure=fig)
    plt.plot()
    plt.xlabel('Dimension')
    plt.ylabel('Residual Variance')
    plt.savefig('%s%s_%s_residual_variance_neighbor%d' % (par.path,par.mode,par.method,par.neighbor))
    plt.clf()

