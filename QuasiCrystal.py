import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from NanoCavity import nanoCavity
#from Histogram import histogram
#from progressBar import printProgressBar

PI = np.pi

class criticalPt:
    def __init__(self, xc, yc, type):
        self.xc = xc
        self.yc = yc
        self.type = type


class quasiCrystal():
    def __init__(self, listGrating, stretching=-1, resolution=0.5, dimensionX=3, dimensionY=3, xOffset=0, yOffset=0, percentageOfDiffForDoppel=0.5,
                 filtering=True):
        self.dimensionX, self.dimensionY = dimensionX, dimensionY
        self.arrayMin, self.arrayMax, self.arraySad = [], [], []
        self.distMinMin, self.distMinMax, self.distMinSad = [], [], []
        self.xOffset, self.yOffset = xOffset, yOffset
        self.criticalPt = []
        self.listWell = []
        self.stretching = stretching

        self.Grating = listGrating
        self.resolution = resolution
        self.percentageOfDiffForDoppel = percentageOfDiffForDoppel

        print("Searching for critical point")
        self.findCriticalPoint()
        print("Found %d critical point" % len(self.criticalPt))
        print("Testing extremum type")
        self.assignExtremumTypeToCriticalPt()
        print("Found %d min" % len(self.arrayMin))
        print("Found %d max" % len(self.arrayMax))
        print("Found %d saddle" % len(self.arraySad))

        print("Creating well object")
        self.populateNanoCavities()
        if (filtering):
            print('Nbr of wells before filtering =', len(self.listWell))
            self.filterWellDoppelganger()
            print('Nbr of well types =', len(self.listWell))


    def findCriticalPoint(self):

        def distance(xc, yc, pt):
            return math.sqrt((xc - pt[0]) ** 2 + (yc - pt[1]) ** 2)

        def findZeroDerivative(p):
            x, y = p

            def partialDx(x, y):
                sol = 0
                for g in self.Grating:
                    sol += g.partialDx(x, y)
                return sol

            def partialDy(x, y):
                sol = 0
                for g in self.Grating:
                    sol += g.partialDy(x, y)
                return sol

            return (partialDx(x, y), partialDy(x, y))

        def jacobian(p):
            x, y = p
            Dxx = 0
            for g in self.Grating:
                Dxx += g.partialDxx(x, y)
            Dxy = 0
            for g in self.Grating:
                Dxy += g.partialDxy(x, y)
            Dyx = 0
            for g in self.Grating:
                Dyx += g.partialDyx(x, y)
            Dyy = 0
            for g in self.Grating:
                Dyy += g.partialDyy(x, y)
            return [[Dxx, Dxy], [Dyx, Dyy]]

        nbPointX = int(self.dimensionX / self.resolution)
        nbPointY = int(self.dimensionY / self.resolution)

        self.criticalPt = []
        x0 = self.xOffset
        for idxX in range(nbPointX):
            x0 += self.resolution
            y0 = self.yOffset
            for idxY in range(nbPointY):
                y0 += self.resolution
                xc, yc = fsolve(func=findZeroDerivative, x0=(x0, y0), args=(), fprime=jacobian)
                isOK = False
                if (xc < self.dimensionX + self.xOffset) and (xc >  self.xOffset):
                    if (yc < self.dimensionY + self.yOffset) and (yc >  self.yOffset) :
                        # Test if the critical point is not on the edges of the sample
                        # if (xc > (self.dimensionX + + self.xOffset) * 0.02) and (xc < 0.98 * (self.dimensionX + + self.xOffset)):
                        #     if (yc > (self.dimensionY + self.yOffset) * 0.02) and (yc < 0.98 * (self.dimensionY + self.yOffset)):
                        isOK = True
                        # Test us the critical point is not already in the self.criticalPt list +/- a fraction of the resolution
                        for pt in self.criticalPt:
                            if (distance(xc, yc, pt) < self.resolution / 5.0):
                                isOK = False
                if (isOK):
                    self.criticalPt.append([xc, yc, "n"])

    def assignExtremumTypeToCriticalPt(self):
        """
        Pour chaque couple de self.listCriticalPt, trouve qui est Max, min ou saddle et les place dans un array
        avec D = dxx * dyy - (dxy * dxy)
        si D>0 et dxx<0 alors (xc,yc) est un maximum
        si D>0 et dxx>0 alors (xc,yc) est un minimum
        si D<0 alors (xc,yc) est un minimum en scelle de cheval (saddle)
        si D=0 alors le couple (xc,yc) est quelconque
        """
        listMin = []
        listMax = []
        listSad = []

        for point in self.criticalPt:
            D = 0
            sDxx = sDyy = sDxy = 0
            z=0
            for g in self.Grating:
                sDxx += g.partialDxx(point[0], point[1])
                sDyy += g.partialDyy(point[0], point[1])
                sDyy += g.partialDxy(point[0], point[1])
                #D = sDxx * sDyy - (sDxy) ** 2
                z += g.getGrattingValue(point[0], point[1])

            D = sDxx * sDyy - (sDxy) ** 2
            if (D > 0) and (sDxx > 0):
                #Minimum
                # Test if the minimum is not on the edges of the sample
                if (point[0] > (self.dimensionX + + self.xOffset) * 0.02) and (point[0] < 0.98 * (self.dimensionX + + self.xOffset)):
                    if (point[1] > (self.dimensionY + self.yOffset) * 0.02) and (point[1] < 0.98 * (self.dimensionY + self.yOffset)):
                        listMin.append((point[0], point[1], z))
                        point[2] = 'm'
            elif (D > 0) and (sDxx < 0):
                #Maximum
                listMax.append((point[0], point[1], z))
                point[2] = 'M'
            elif (D < 0):
                #Saddle
                listSad.append((point[0], point[1], z))
                point[2] = 's'

        self.arrayMin = np.asarray(listMin)
        self.arrayMax = np.asarray(listMax)
        self.arraySad = np.asarray(listSad)

        """
        Uniformisation des hauteur de l'echantillon
        Toute les valeurs sont comprise entre 0 et 1
        Si on les multiplie par la hauteur calculée par les courbes Bessel ou AFM
        alors l'axe des z n'est plus arbitraire (mais n'est pas absolu)
        """
        # self.minmin = min(self.arrayMin[:, 2])
        # self.maxmax = max(self.arrayMax[:, 2])
        # self.arrayMin[:, 2] = (self.arrayMin[:, 2] + np.abs(self.minmin) / self.maxmax)
        # self.arrayMax[:, 2] = (self.arrayMax[:, 2] + np.abs(self.minmin) / self.maxmax)
        # self.arraySad[:, 2] = (self.arraySad[:, 2] + np.abs(self.minmin) / self.maxmax)

    def populateNanoCavities(self):
        """
        Fait appel à la class Well : calcule et trouve les carac des puits de l'échantillon
        """
        for idxMin in range(np.size(self.arrayMin[:, 0])):
            #def __init__(self, listGrating, minLoc, arrayMin, arrayMax, arraySad, drawResolution=0.05):
            self.listWell.append(
                nanoCavity(self.Grating, idxMin, self.arrayMin, self.arrayMax, self.arraySad, self.resolution / 10))

    def filterWellDoppelganger(self):
        """
        Pour chaque "if", on regarde si "i" et "j" est identique ; à chaque "if" on compare "i" et "j"
        Si "non" on passe au "j" suivant
        Si "oui" on passe au "if" suivant et si tout les "if" sont vrai alors les puits sont identiques
        Si les puit sont identiques on conserve le premier puits observé en guise de référence
        """

        #First : eliminate all the well that has no neighbors and that have failed the MVEE analysis
        nbOfwell = len(self.listWell)
        i = 0
        while i < nbOfwell:
            if (self.listWell[i].ellipsoid == -1) \
                    or (self.listWell[i].nbOfNeighbors == 0):
                self.listWell.pop(i)
                # nbOfwell -= 1 ?
                nbOfwell = len(self.listWell)
            else:
                i += 1

        #Secondly two wells are declared the same when  they have:
        #   -  the same number of neighbors
        #   -  the same number of Max
        #   -  the same number of Saddle
        #
        #
        #   -   For each neighbor :
        #           Is the normalized distance difference


        # Je prends une distance  pour le puit i. Puis-je trouver une distance dans j qui est égale à self.percentageOfDiffForDoppel
        # La distance est un critère pas assez contraignant, à cause des symétries. Le mieux est d'utiliser le resultat de la MVEE.

        nbOfwell = len(self.listWell)
        i = 0   #currentWell
        while (i < nbOfwell):
            j = i + 1   #Closest Well, since it is in principle the more similar to well i
            while (j < nbOfwell):
                doFilter = False

                #   -  the same number of neighbors
                #   -  the same number of Max
                #   -  the same number of Saddle
                if (self.listWell[j].nbOfNeighbors == self.listWell[i].nbOfNeighbors) \
                        and (self.listWell[j].nbOfNeighborsMax == self.listWell[i].nbOfNeighborsMax) \
                        and (self.listWell[j].nbOfNeighborsSad == self.listWell[i].nbOfNeighborsSad):

                    i_ellipsoid = self.listWell[i].ellipsoid
                    j_ellipsoid = self.listWell[j].ellipsoid

                    #TODO cela ne marche pas si on a une sphère à cause des approx numérique
                    relativeEccentricityDiff = abs(i_ellipsoid.eccentricity - j_ellipsoid.eccentricity) / max(i_ellipsoid.eccentricity, j_ellipsoid.eccentricity)
                    relativeEllipticityDiff = abs(i_ellipsoid.ellipticity - j_ellipsoid.ellipticity) / max(i_ellipsoid.ellipticity, j_ellipsoid.ellipticity)
                    relativeVolumeDiff = abs(i_ellipsoid.volume - j_ellipsoid.volume) / max(i_ellipsoid.volume, j_ellipsoid.volume)

                    if (relativeEccentricityDiff < self.percentageOfDiffForDoppel) and (relativeEllipticityDiff < self.percentageOfDiffForDoppel) and (relativeVolumeDiff < self.percentageOfDiffForDoppel) :
                        doFilter = True

                    # # Testing if
                    # distDiff = np.zeros(len(self.listWell[i].dist2Min))
                    # for k in range(len(self.listWell[i].dist2Min)):
                    #     # normalized distance difference for each neighbors
                    #     distDiff[k] = (abs((self.listWell[i].dist2Min[k] - self.listWell[j].dist2Min[k])) / max(
                    #         self.listWell[i].dist2Min[k], self.listWell[j].dist2Min[k]))
                    #
                    # # Find the point where the normalized distance difference is inferior to our threshold self.percentageOfDiffForDoppel
                    # idxMax = np.where(distDiff[:] < self.percentageOfDiffForDoppel)
                    # if (len(idxMax[0]) == self.listWell[i].nbOfNeighbors):
                    #     doFilter = True

                if (doFilter):
                    self.listWell.pop(j)
                    nbOfwell = len(self.listWell)

                #NB : if we have filtered th j-th well, no need to do j +=1
                j += 1

            i += 1

    # def getNanoCavitiesStats(self, basicinfo=True, MVEE=True):
    #     """
    #     Affiche les dimension des puits après les avoir filtrer si demander
    #     """
    #
    #     for i in range(len(self.listWell)):
    #         if (self.listWell[i].nbOfNeighbors == 0):
    #             print('>> Well n°', i, '(', self.listWell[i].minLoc, ')')
    #             print('* No informations available *')
    #         else:
    #             print('>> Well n°', i, '(', self.listWell[i].minLoc, ')')
    #             print('( Xc = ', self.listWell[i].xc, '; Yc = ', self.listWell[i].yc, '; Zc = ', self.listWell[i].zc,
    #                   ')')
    #             print(' Nb of neighbors :', self.listWell[i].nbOfNeighbors, '( Max :',
    #                   self.listWell[i].nbOfNeighborsMax, '; Sad :',
    #                   self.listWell[i].nbOfNeighborsSad, ')')
    #             if (basicinfo):
    #                 print('- Basic informations :')
    #                 print(' min height =', self.listWell[i].heightMin)
    #                 print(' max height =', self.listWell[i].heightMax)
    #                 print(' min width  =', self.listWell[i].widthMin)
    #                 print(' max width  =', self.listWell[i].widthMax)
    #
    #             if (MVEE):
    #                 if (self.listWell[i].ellipsoid != -1):
    #                     print('- Ellipsoid informations :')
    #                     print(' eccentricity =', self.listWell[i].ellipsoid.eccentricity)
    #                     print(' ellipticity  =', self.listWell[i].ellipsoid.ellipticity)
    #                     print(' a diameter   =', self.listWell[i].ellipsoid.a)
    #                     print(' b diameter   =', self.listWell[i].ellipsoid.b)
    #                     print(' c diameter   =', self.listWell[i].ellipsoid.c)
    #                 else:
    #                     print('* No ellispoid informations available *')
    #
    #         print('\n')


    def drawContour(self, nbOfLevel=50, fname=None):
        X = np.arange(self.xOffset, self.dimensionX + self.xOffset, self.resolution)
        Y = np.arange(self.yOffset, self.dimensionY + self.yOffset, self.resolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.Grating[0].getGrattingValue(X, Y)
        for i in range(1, len(self.Grating)):
            Z += self.Grating[i].getGrattingValue(X, Y)



        plt.figure()
        cp = plt.contourf(X, Y, Z, nbOfLevel, cmap=cm.jet)
        plt.colorbar(cp)
        plt.title('Filled Contours Plot')
        plt.xlabel('x ')
        plt.ylabel('y ')

        if fname is not None:
            plt.savefig(fname, dpi=600)
        plt.show()

    def drawAllSample(self, dots=True, coat='full', coatalpha=0.5, ellipses=False, wellidx=[], ellipsesalpha=0.5):
        fig = plt.figure()
        ax = Axes3D(fig)
        plt.hold(True)

        X = np.arange(self.xOffset, self.dimensionX + self.xOffset, self.resolution)
        Y = np.arange(self.yOffset, self.dimensionY + self.yOffset, self.resolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.Grating[0].getGrattingValue(X, Y)
        for i in range(1, len(self.Grating)):
            Z += self.Grating[i].getGrattingValue(X, Y)


        #TODO le scaling, ici on considere l'amplitude max à partir de la somme des réseaux s'ils sont tous en phase
        ax.set_zlim(-len(self.Grating) - 0.01, len(self.Grating) + 0.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        if (dots):
            for crPt in self.criticalPt:
                zc = self.Grating[0].getGrattingValue(crPt[0], crPt[1])
                for i in range(1, len(self.Grating)):
                    zc += self.Grating[i].getGrattingValue(crPt[0], crPt[1])
                if crPt[2] == 'm':
                    ax.scatter(crPt[0], crPt[1], zc, color='b')
                elif crPt[2] == 's':
                    ax.scatter(crPt[0], crPt[1], zc, color='g')
                elif crPt[2] == 'M':
                    ax.scatter(crPt[0], crPt[1], zc, color='r')

        if (ellipses):
            sin, cos = np.sin, np.cos
            for idx in wellidx:
                if (self.listWell[idx].ellipsoid) and (self.listWell[idx].ellipsoid != -1):
                    # rx, ry, rz = 1. / np.sqrt(self.listWell[idx].ellipsoid.D)
                    u, v = np.mgrid[0:2 * PI:20j, -PI / 2:PI / 2:10j]

                    def ellipse(RX, RY, RZ, u, v):
                        x = RX * cos(u) * cos(v)
                        y = RY * sin(u) * cos(v)
                        z = RZ * sin(v)
                        return x, y, z

                    RX = self.listWell[idx].ellipsoid.a / 2
                    RY = self.listWell[idx].ellipsoid.b / 2
                    RZ = self.listWell[idx].ellipsoid.c / 2
                    E = np.dstack(ellipse(RX, RY, RZ, u, v))
                    E = np.dot(E, self.listWell[idx].ellipsoid.V) + self.listWell[idx].ellipsoid.center
                    x, y, z = np.rollaxis(E, axis=-1)
                    ax.plot_surface(x, y, z, cstride=1, rstride=1, color='y', alpha=ellipsesalpha)

        if (coat == 'full'):
            surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.coolwarm, alpha=coatalpha, linewidth=0,
                                   antialiased=False)
            fig.colorbar(surf, shrink=0.5, aspect=5)
        elif (coat == 'wire'):
            ax.plot_surface(X, Y, Z, rstride=2, cstride=2, alpha=coatalpha)

        plt.show()

    # def saveFile(self, name='G?', basicinfo=True, MVEE=True):
    #     name1 = name + '_(Basic_stretch' + str(self.stretching * 2) + 'mm_dimX' + str(self.dimensionX) + '_dimY' + str(
    #         self.dimensionY) + '_res' + str(self.resolution) + '_filtering' + str(self.percentageOfDiffForDoppel) + ')'
    #     file = open('{0}.txt'.format(name1), "w")
    #
    #     for i in range(len(self.listWell)):
    #         if (self.listWell[i].nbOfNeighbors != 0) and (self.listWell[i].heightMin != -1):
    #             file.write('X = ; ' + str(self.listWell[i].xc) + '; Y = ;' + str(self.listWell[i].yc) + '; Z = ;' + str(
    #                 self.listWell[i].zc))
    #             file.write('; Neigh = ;' + str(self.listWell[i].nbOfNeighbors))
    #         if (basicinfo):
    #             file.write('; Hmin = ;' + str(self.listWell[i].heightMin))
    #             file.write('; Hmax = ;' + str(self.listWell[i].heightMax))
    #             file.write('; Wmin = ;' + str(self.listWell[i].widthMin))
    #             file.write('; Wmax = ;' + str(self.listWell[i].widthMax))
    #         file.write('\n')
    #     file.close()
    #
    #     name2 = name + '_(Ellipse_stretch' + str(self.stretching * 2) + 'mm_dimX' + str(
    #         self.dimensionX) + '_dimY' + str(
    #         self.dimensionY) + '_resolution' + str(self.resolution) + '_filtering' + str(self.percentageOfDiffForDoppel) + ')'
    #     file = open('{0}.txt'.format(name2), "w")
    #
    #     for i in range(len(self.listWell)):
    #         if (MVEE):
    #             if (self.listWell[i].ellipsoid != -1):
    #                 file.write(
    #                     'X = ; ' + str(self.listWell[i].xc) + '; Y = ;' + str(self.listWell[i].yc) + '; Z = ;' + str(
    #                         self.listWell[i].zc))
    #                 file.write('; Neigh = ;' + str(self.listWell[i].nbOfNeighbors))
    #                 file.write('; ellipticity = ;' + str(self.listWell[i].ellipsoid.ellipticity))
    #                 file.write('; Aaxis = ;' + str(self.listWell[i].ellipsoid.a))
    #                 file.write('; Baxis = ;' + str(self.listWell[i].ellipsoid.b))
    #                 file.write('; Caxis = ;' + str(self.listWell[i].ellipsoid.c))
    #         file.write('\n')
    #     file.close()
    #
    #     print('\n')

