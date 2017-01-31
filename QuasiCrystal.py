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
    """
    Main class of the program representing the quasiCrystal as a sum of sinusoidal Gratings.
    Parameters :
    - listSinusoidalGrating : list of sinusodialGraing object
    - minPtSearchResolution :
    - sizeX : Size in µm in the X direction of the quasiCrystal to reconstruct (and analyse) from a sum of sinusoidal Gratings
    - sizeY : same as previous parameter, but for Y direction
    - xOffset : shift in µm for the X direction from the center
    - yOffset : shift in µm for the Y direction from the center
    -
    """
    def __init__(self, listSinusoidalGrating, minPtSearchResolution=0.5, sizeX=3, sizeY=3, xOffset=0, yOffset=0, maxRelativeDiffForEquivalentFilter=0.5,
                 isFilterEquivalentCavities=True):

        self.sizeX, self.sizeY = sizeX, sizeY
        self.arrayMin, self.arrayMax, self.arraySad = [], [], []
        self.distMinMin, self.distMinMax, self.distMinSad = [], [], []
        self.xOffset, self.yOffset = xOffset, yOffset
        self.listCriticalPt = []
        self.listNanoCavities = []

        self.gratings = listSinusoidalGrating
        self.minPtSearchResolution = minPtSearchResolution
        self.maxRelativeDiffForEquivalentFilter = maxRelativeDiffForEquivalentFilter

        print("Searching for critical point")
        self.findCriticalPoint()
        print("Found %d critical point" % len(self.listCriticalPt))
        print("Testing extremum type")
        self.assignExtremumTypeToCriticalPt()
        print("Found %d min" % len(self.arrayMin))
        print("Found %d max" % len(self.arrayMax))
        print("Found %d saddle" % len(self.arraySad))

        print("Creating well object")
        self.populateNanoCavities()
        if (isFilterEquivalentCavities):
            print('Nbr of wells before isFilterEquivalentCavities =', len(self.listNanoCavities))
            self.filterEquivalentCavities()
            print('Nbr of well types =', len(self.listNanoCavities))


    def findCriticalPoint(self):
        """
        Scan the quasiCrystal surface with a resolution of self.minPtSearchResolution in order to find local critical point,
        i.e. points where the first 2D derivative is zero.
        Put this points in self.listCriticalPt, after being sure that this critical point is not already there
        """

        def distance(xc, yc, pt):
            return math.sqrt((xc - pt[0]) ** 2 + (yc - pt[1]) ** 2)

        def findZeroDerivative(p):
            x, y = p

            def partialDx(x, y):
                sol = 0
                for g in self.gratings:
                    sol += g.partialDx(x, y)
                return sol

            def partialDy(x, y):
                sol = 0
                for g in self.gratings:
                    sol += g.partialDy(x, y)
                return sol

            return (partialDx(x, y), partialDy(x, y))

        def jacobian(p):
            x, y = p
            Dxx = 0
            for g in self.gratings:
                Dxx += g.partialDxx(x, y)
            Dxy = 0
            for g in self.gratings:
                Dxy += g.partialDxy(x, y)
            Dyx = 0
            for g in self.gratings:
                Dyx += g.partialDyx(x, y)
            Dyy = 0
            for g in self.gratings:
                Dyy += g.partialDyy(x, y)
            return [[Dxx, Dxy], [Dyx, Dyy]]

        nbPointX = int(self.sizeX / self.minPtSearchResolution)
        nbPointY = int(self.sizeY / self.minPtSearchResolution)

        self.listCriticalPt = []
        x0 = self.xOffset
        for idxX in range(nbPointX):
            x0 += self.minPtSearchResolution
            y0 = self.yOffset
            for idxY in range(nbPointY):
                y0 += self.minPtSearchResolution
                xc, yc = fsolve(func=findZeroDerivative, x0=(x0, y0), args=(), fprime=jacobian)
                isOK = False
                if (xc < self.sizeX + self.xOffset) and (xc >  self.xOffset):
                    if (yc < self.sizeY + self.yOffset) and (yc >  self.yOffset) :
                        isOK = True
                        # Test if the critical point is not already in the self.listCriticalPt list +/- a fraction of the minPtSearchResolution
                        for pt in self.listCriticalPt:
                            if (distance(xc, yc, pt) < self.minPtSearchResolution / 5.0):
                                isOK = False
                if (isOK):
                    self.listCriticalPt.append([xc, yc, "n"])

    def assignExtremumTypeToCriticalPt(self):
        """
        For each critical Point inside self.listCriticalPt we perform a second derivative test in order to know if the
        critical point is a minimum, a maximum or a saddle point

        With D = dxx * dyy - (dxy * dxy)
        if D>0 and dxx<0 then (xc,yc) is a maximum
        if D>0 and dxx>0 then (xc,yc) is a minimum
        if D<0 then (xc,yc) is a saddle point
        if D=0 then (xc,yc) is any kind...
        """
        listMin = []
        listMax = []
        listSad = []

        for criticalPt in self.listCriticalPt:
            xc = criticalPt[0]
            yc = criticalPt[1]
            sDxx = sDyy = sDxy = 0
            z=0
            for g in self.gratings:
                sDxx += g.partialDxx(xc, yc)
                sDyy += g.partialDyy(xc, yc)
                sDyy += g.partialDxy(xc, yc)
                z += g.getGrattingValue(xc, yc)

            D = sDxx * sDyy - (sDxy) ** 2
            if (D > 0) and (sDxx > 0):
                #Minimum
                # Test if the minimum is not on the edges of the sample, if it is the case, discard this point.
                if (xc > (self.sizeX + + self.xOffset) * 0.02) and (xc < 0.98 * (self.sizeX + + self.xOffset)):
                    if (yc > (self.sizeY + self.yOffset) * 0.02) and (yc < 0.98 * (self.sizeY + self.yOffset)):
                        listMin.append((xc, yc, z))
                        criticalPt[2] = 'm'
            elif (D > 0) and (sDxx < 0):
                #Maximum
                listMax.append((xc, yc, z))
                criticalPt[2] = 'M'
            elif (D < 0):
                #Saddle
                listSad.append((xc, yc, z))
                criticalPt[2] = 's'

        self.arrayMin = np.asarray(listMin)
        self.arrayMax = np.asarray(listMax)
        self.arraySad = np.asarray(listSad)

    def populateNanoCavities(self):
        """
        Explore the nanocavity around a minimum on the quasiCrystal Surface.
        """

        nbOfMin = np.size(self.arrayMin[:, 0])
        for idxMin in range(nbOfMin):
            self.listNanoCavities.append(
                nanoCavity(self.gratings, idxMin, self.arrayMin, self.arrayMax, self.arraySad, self.minPtSearchResolution / 10.0))

    def filterEquivalentCavities(self):
        """
        Search and remove among the minimum list, the cavity that are similar.
        Due to the symmetry, if we increase the size of the quasiCrystal (sizeX and sizeX parameters) we will artificially increase the number of cavities.
        """

        #First : eliminate all the cavities that has no neighbors and that have failed the MVEE analysis
        nbOfCavities = len(self.listNanoCavities)
        i = 0
        while i < nbOfCavities:
            if (self.listNanoCavities[i].ellipsoid == -1) \
                    or (self.listNanoCavities[i].nbOfNeighbors == 0):
                self.listNanoCavities.pop(i)
                # nbOfwell -= 1 ?
                nbOfCavities = len(self.listNanoCavities)
            else:
                i += 1

        #Secondly two cavities are declared equivalent when they have:
        #   -  the same number of neighbors
        #   -  the same number of Max
        #   -  the same number of Saddle
        #   -  The relative difference in ellipsoid volume, eccentricity and ellipticiy are smaller than self.maxRelativeDiffForEquivalentFilter


        nbOfCavities = len(self.listNanoCavities)
        i = 0   #currentCavity
        while (i < nbOfCavities):
            j = i + 1   #Closest cavity, since it is in principle the more similar to cavity i
            while (j < nbOfCavities):
                doFilter = False

                #   -  the same number of neighbors
                #   -  the same number of Max
                #   -  the same number of Saddle
                if (self.listNanoCavities[j].nbOfNeighbors == self.listNanoCavities[i].nbOfNeighbors) \
                        and (self.listNanoCavities[j].nbOfNeighborsMax == self.listNanoCavities[i].nbOfNeighborsMax) \
                        and (self.listNanoCavities[j].nbOfNeighborsSad == self.listNanoCavities[i].nbOfNeighborsSad):

                    i_ellipsoid = self.listNanoCavities[i].ellipsoid
                    j_ellipsoid = self.listNanoCavities[j].ellipsoid

                    #TODO cela ne marche pas si on a une sphère à cause des approx numérique
                    relativeEccentricityDiff = abs(i_ellipsoid.eccentricity - j_ellipsoid.eccentricity) / max(i_ellipsoid.eccentricity, j_ellipsoid.eccentricity)
                    relativeEllipticityDiff = abs(i_ellipsoid.ellipticity - j_ellipsoid.ellipticity) / max(i_ellipsoid.ellipticity, j_ellipsoid.ellipticity)
                    relativeVolumeDiff = abs(i_ellipsoid.volume - j_ellipsoid.volume) / max(i_ellipsoid.volume, j_ellipsoid.volume)

                    if (relativeEccentricityDiff < self.maxRelativeDiffForEquivalentFilter) and (relativeEllipticityDiff < self.maxRelativeDiffForEquivalentFilter) and (relativeVolumeDiff < self.maxRelativeDiffForEquivalentFilter) :
                        doFilter = True

                if (doFilter):
                    self.listNanoCavities.pop(j)
                    nbOfCavities = len(self.listNanoCavities)

                #NB : if we have filtered th j-th well, no need to do j +=1
                j += 1

            i += 1

    #########################

    # def getNanoCavitiesStats(self, basicinfo=True, MVEE=True):
    #     """
    #     Affiche les dimension des puits après les avoir filtrer si demander
    #     """
    #
    #     for i in range(len(self.listNanoCavities)):
    #         if (self.listNanoCavities[i].nbOfNeighbors == 0):
    #             print('>> Well n°', i, '(', self.listNanoCavities[i].minIdx, ')')
    #             print('* No informations available *')
    #         else:
    #             print('>> Well n°', i, '(', self.listNanoCavities[i].minIdx, ')')
    #             print('( Xc = ', self.listNanoCavities[i].xc, '; Yc = ', self.listNanoCavities[i].yc, '; Zc = ', self.listNanoCavities[i].zc,
    #                   ')')
    #             print(' Nb of neighbors :', self.listNanoCavities[i].nbOfNeighbors, '( Max :',
    #                   self.listNanoCavities[i].nbOfNeighborsMax, '; Sad :',
    #                   self.listNanoCavities[i].nbOfNeighborsSad, ')')
    #             if (basicinfo):
    #                 print('- Basic informations :')
    #                 print(' min height =', self.listNanoCavities[i].heightMin)
    #                 print(' max height =', self.listNanoCavities[i].heightMax)
    #                 print(' min width  =', self.listNanoCavities[i].widthMin)
    #                 print(' max width  =', self.listNanoCavities[i].widthMax)
    #
    #             if (MVEE):
    #                 if (self.listNanoCavities[i].ellipsoid != -1):
    #                     print('- Ellipsoid informations :')
    #                     print(' eccentricity =', self.listNanoCavities[i].ellipsoid.eccentricity)
    #                     print(' ellipticity  =', self.listNanoCavities[i].ellipsoid.ellipticity)
    #                     print(' a diameter   =', self.listNanoCavities[i].ellipsoid.a)
    #                     print(' b diameter   =', self.listNanoCavities[i].ellipsoid.b)
    #                     print(' c diameter   =', self.listNanoCavities[i].ellipsoid.c)
    #                 else:
    #                     print('* No ellispoid informations available *')
    #
    #         print('\n')


    def drawContour(self, nbOfLevel=50, fname=None):
        X = np.arange(self.xOffset, self.sizeX + self.xOffset, self.minPtSearchResolution)
        Y = np.arange(self.yOffset, self.sizeY + self.yOffset, self.minPtSearchResolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.gratings[0].getGrattingValue(X, Y)
        for i in range(1, len(self.gratings)):
            Z += self.gratings[i].getGrattingValue(X, Y)



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

        X = np.arange(self.xOffset, self.sizeX + self.xOffset, self.minPtSearchResolution)
        Y = np.arange(self.yOffset, self.sizeY + self.yOffset, self.minPtSearchResolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.gratings[0].getGrattingValue(X, Y)
        for i in range(1, len(self.gratings)):
            Z += self.gratings[i].getGrattingValue(X, Y)


        #TODO le scaling, ici on considere l'amplitude max à partir de la somme des réseaux s'ils sont tous en phase
        ax.set_zlim(-len(self.gratings) - 0.01, len(self.gratings) + 0.01)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        if (dots):
            for crPt in self.listCriticalPt:
                zc = self.gratings[0].getGrattingValue(crPt[0], crPt[1])
                for i in range(1, len(self.gratings)):
                    zc += self.gratings[i].getGrattingValue(crPt[0], crPt[1])
                if crPt[2] == 'm':
                    ax.scatter(crPt[0], crPt[1], zc, color='b')
                elif crPt[2] == 's':
                    ax.scatter(crPt[0], crPt[1], zc, color='g')
                elif crPt[2] == 'M':
                    ax.scatter(crPt[0], crPt[1], zc, color='r')

        if (ellipses):
            sin, cos = np.sin, np.cos
            for idx in wellidx:
                if (self.listNanoCavities[idx].ellipsoid) and (self.listNanoCavities[idx].ellipsoid != -1):
                    # rx, ry, rz = 1. / np.sqrt(self.listNanoCavities[idx].ellipsoid.D)
                    u, v = np.mgrid[0:2 * PI:20j, -PI / 2:PI / 2:10j]

                    def ellipse(RX, RY, RZ, u, v):
                        x = RX * cos(u) * cos(v)
                        y = RY * sin(u) * cos(v)
                        z = RZ * sin(v)
                        return x, y, z

                    RX = self.listNanoCavities[idx].ellipsoid.a / 2
                    RY = self.listNanoCavities[idx].ellipsoid.b / 2
                    RZ = self.listNanoCavities[idx].ellipsoid.c / 2
                    E = np.dstack(ellipse(RX, RY, RZ, u, v))
                    E = np.dot(E, self.listNanoCavities[idx].ellipsoid.V) + self.listNanoCavities[idx].ellipsoid.center
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
    #     name1 = name + '_(Basic_stretch' + str(self.stretching * 2) + 'mm_dimX' + str(self.sizeX) + '_dimY' + str(
    #         self.sizeY) + '_res' + str(self.minPtSearchResolution) + '_filtering' + str(self.maxRelativeDiffForEquivalentFilter) + ')'
    #     file = open('{0}.txt'.format(name1), "w")
    #
    #     for i in range(len(self.listNanoCavities)):
    #         if (self.listNanoCavities[i].nbOfNeighbors != 0) and (self.listNanoCavities[i].heightMin != -1):
    #             file.write('X = ; ' + str(self.listNanoCavities[i].xc) + '; Y = ;' + str(self.listNanoCavities[i].yc) + '; Z = ;' + str(
    #                 self.listNanoCavities[i].zc))
    #             file.write('; Neigh = ;' + str(self.listNanoCavities[i].nbOfNeighbors))
    #         if (basicinfo):
    #             file.write('; Hmin = ;' + str(self.listNanoCavities[i].heightMin))
    #             file.write('; Hmax = ;' + str(self.listNanoCavities[i].heightMax))
    #             file.write('; Wmin = ;' + str(self.listNanoCavities[i].widthMin))
    #             file.write('; Wmax = ;' + str(self.listNanoCavities[i].widthMax))
    #         file.write('\n')
    #     file.close()
    #
    #     name2 = name + '_(Ellipse_stretch' + str(self.stretching * 2) + 'mm_dimX' + str(
    #         self.sizeX) + '_dimY' + str(
    #         self.sizeY) + '_resolution' + str(self.minPtSearchResolution) + '_filtering' + str(self.maxRelativeDiffForEquivalentFilter) + ')'
    #     file = open('{0}.txt'.format(name2), "w")
    #
    #     for i in range(len(self.listNanoCavities)):
    #         if (MVEE):
    #             if (self.listNanoCavities[i].ellipsoid != -1):
    #                 file.write(
    #                     'X = ; ' + str(self.listNanoCavities[i].xc) + '; Y = ;' + str(self.listNanoCavities[i].yc) + '; Z = ;' + str(
    #                         self.listNanoCavities[i].zc))
    #                 file.write('; Neigh = ;' + str(self.listNanoCavities[i].nbOfNeighbors))
    #                 file.write('; ellipticity = ;' + str(self.listNanoCavities[i].ellipsoid.ellipticity))
    #                 file.write('; Aaxis = ;' + str(self.listNanoCavities[i].ellipsoid.a))
    #                 file.write('; Baxis = ;' + str(self.listNanoCavities[i].ellipsoid.b))
    #                 file.write('; Caxis = ;' + str(self.listNanoCavities[i].ellipsoid.c))
    #         file.write('\n')
    #     file.close()
    #
    #     print('\n')

