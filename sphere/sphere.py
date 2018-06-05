"""
Plotting/sampling points on a sphere.
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

from vMF import *


__all__ = ['Sphere', 'vMFSample', 'Coord3D']


class Coord3D():
    """
    Store 3D coordinates.
    """

    def __init__(self, x, y, z):
        """
        Store 3D coordinates.
        """
        self.x = x
        self.y = y
        self.z = z


class Sphere():
    """
    Plotting/sampling points on a sphere.
    """

    def __init__(self):
        """
        Plotting points on a sphere.
        """
        self.sphere = []
        self.samples = []
        self.distribution = ['uniform', 'vMF']
        
        
    def _make_sphere(self, radius):
        """
        Get mesh for a sphere of a given radius.
        """
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0 * pi:100j]
        x = radius * sin(phi) * cos(theta)
        y = radius * sin(phi) * sin(theta)
        z = radius * cos(phi)

        self.sphere = Coord3D(x, y, z)
        return self.sphere


    def _draw_sphere(self, ax, radius):
        """
        Draw a sphere on a given axis.
        """
        sphere = self._make_sphere(radius - 0.01) # subtract a little so points show on sphere
        ax.plot_surface(sphere.x, sphere.y, sphere.z,  
                        rstride = 1, cstride = 1, 
                        color = sns.xkcd_rgb["light grey"],
                        alpha = 0.5,
                        linewidth = 0)
    

    def plot(self, radius = 1, data = None):
        """
        Plot a sphere with samples superposed, if supplied.
        
        :param radius: radius of the base sphere
        :param data: list of sample objects 
        """

        sns.set_style('dark')
        fig = plt.figure(figsize = [10, 10])
        ax = fig.add_subplot(111, projection = '3d')

        # plot the sphere
        self._draw_sphere(ax, radius)

        if data == None:
            data = self.samples

        is_list = True
        try:
            N = len(data)
        except:
            N = 1
            is_list = False

        palette = sns.color_palette("GnBu_d", N)
        colors = [palette[i] for i in reversed(range(N))]
           
        # plot data, if supplied
        if is_list:
            data_check = data[0]
        else:
            data_check = data
            
        if type(data_check) is vMFSample:
            i = 0
            if is_list:
                for d in data:
                    ax.scatter(d.x, d.y, d.z, s = 50, alpha = 0.7,
                               label = '$\kappa = $' + str(d.kappa), color = colors[i])
                    i += 1
            else:
                ax.scatter(data.x, data.y, data.z, s = 50, alpha = 0.7,
                           label = '$\kappa = $' + str(data.kappa), color = colors[i])
                    
                    
        elif type(data_check) is Coord3D:
            i = 0
            if is_list:
                for d in data:
                    ax.scatter(d.x, d.y, d.z, s = 50, alpha = 0.7,
                               label = 'uniform samples', color = palette[i])
                    i += 1
            else:
                ax.scatter(data.x, data.y, data.z, s = 50, alpha = 0.7,
                           label = 'uniform samples', color = palette[i])
                    
                    
        else:
            print('Error: data type not recognised')

        ax.set_axis_off()
        ax.legend(bbox_to_anchor = [0.65, 0.75])


    def sample(self, n_samples, radius = 1, distribution = "uniform", mu = None, kappa = None):
        """
        Sample points on a spherical surface.
        """

        # uniform
        if distribution == self.distribution[0]:
            u = np.random.uniform(0, 1, n_samples)
            v = np.random.uniform(0, 1, n_samples)

            theta = 2 * np.pi * u
            phi = np.arccos(2 * v - 1)

            # convert to cartesian
            x = radius * np.cos(theta) * np.sin(phi)
            y = radius * np.sin(theta) * np.sin(phi) 
            z = radius * np.cos(phi) 
            self.samples = Coord3D(x, y, z)

        # vMF
        elif distribution == self.distribution[1]:
            try:
                s = sample_vMF(mu, kappa, n_samples)
            except:
                print('Error: mu and kappa must be defined when sampling from vMF')
                return
            self.samples = vMFSample(s, kappa)
                
        else:
            print('Error: sampling distribution not recognised (try \'uniform\' or \'vMF\')')

        return self.samples
        
        
class vMFSample():
    """
    Store 3D coordinates from the sample_vMF function
    and kappa value.
    """

    def __init__(self, sample, kappa):
        """
        Store 3D coordinates from the sample_vMF function
        and kappa value.
        """
        tp = np.transpose(sample)
        self.x = tp[0]
        self.y = tp[1]
        self.z = tp[2]
        self.kappa = kappa
        self.sample = sample
