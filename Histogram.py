import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.stats import laplace as lap


class histogram:
    def __init__(self, data, bars=10, fit='NORM'):
        self.histogram, self.abscissa, self.fit_parameter, self.data_fit= [], [], [], []
        self.data, self.bars = data, bars
        self.mean, self.std = 0, 0
        self.abscissa, self.data_hist = [], []
        self.getHist(data, bars, fit)


    def getHist(self, data, bars, fit):
        """
        Calcul les carac de l'histograme
        Met en place l'abscisse en fonction du nombre de bâtons
        Fait le fitting des 'donnees' où donnees_param [0] = moyenne et donnees_param[1] = STD
        trace le fit suivant une loi normale de paramêtre moyenne et STD
        """
        self.histogram, bin_edges = np.histogram(data, 10)
        self.abscissa = np.linspace(min(bin_edges), max(bin_edges), bars)

        if (fit == 'NORM'):
            self.fit_parameter = norm.fit(data)
            self.data_fit = norm.pdf(self.abscissa, self.fit_parameter[0], self.fit_parameter[1])
        elif (fit == 'LAP'):
            self.fit_parameter = lap.fit(data)
        else :
            self.data_fit = None


    def plotHistogram(self, title='Title', drawfit=True):
        """
        Déssine l'histogramme des 'donnees' avec un nombre de 'batons' de bâtons
        """
        self.data_hist, bins, patches = plt.hist(self.data, self.bars, normed=0, facecolor='b', alpha=0.20)
        if self.data_fit is not None:
            self.data_fit = (self.data_fit / max(self.data_fit)) * max(self.data_hist)

        if (drawfit == True):
            plt.plot(self.abscissa, self.data_fit, 'r')
            plt.title(title)
            plt.xlabel(title)
            plt.ylabel('Occurence')
            plt.ylim([0, 1.05 * max(self.data_fit)])
            plt.show()
        else:
            plt.title(title)
            plt.xlabel(title)
            plt.ylabel('Occurence')
            #plt.ylim([0, 1.05 * max(self.data_fit)])
            plt.show()

    def fitParameters (self):
        print('Mean =', self.fit_parameter[0])
        print('STD =', self.fit_parameter[1])