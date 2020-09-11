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
        """
        Critical point where the derivative of the profile equals zero
        :param xc: x position
        :param yc: y position
        :param type: minimum, maximum, saddle
        """
        self.xc = xc
        self.yc = yc
        self.type = type


class quasiCrystal():
    """
    Main class of the program representing the quasiCrystal as a sum of sinusoidal Gratings.
    """
    def __init__(self, list_sinusoidal_grating, min_pt_search_resolution=0.5, sizeX=3, sizeY=3, xOffset=0, yOffset=0, max_relative_diff_for_equivalent_filter=0.1,
                 is_filter_equivalent_cavities=True, name=None):
        """
        :param list_sinusoidal_grating: list of sinusodialGrating object
        :param min_pt_search_resolution: Resolution in µm used during the systematic scan of critical point (derivative equals zero)
        :param sizeX: Size in µm in the X direction of the quasiCrystal to reconstruct (and analyse) from a sum of sinusoidal Gratings
        :param sizeY: same as previous parameter, but for Y direction
        :param xOffset: shift in µm for the X direction from the center
        :param yOffset: shift in µm for the Y direction from the center
        :param max_relative_diff_for_equivalent_filter:
        :param is_filter_equivalent_cavities:
        """

        self.sizeX, self.sizeY = sizeX, sizeY
        self.arrayMin, self.arrayMax, self.arraySad = [], [], []
        self.distMinMin, self.distMinMax, self.distMinSad = [], [], []
        self.xOffset, self.yOffset = xOffset, yOffset
        self.list_critical_pt = []
        self.list_nano_cavities = []

        self.name = name

        self.gratings = list_sinusoidal_grating
        self.min_pt_search_resolution = min_pt_search_resolution
        self.max_relative_diff_for_equivalent_filter = max_relative_diff_for_equivalent_filter

        print("Searching for critical point")
        self.findCriticalPoint()
        print("Found %d critical point" % len(self.list_critical_pt))
        print("Testing extremum type")
        self.assignExtremumTypeToCriticalPt()
        print("Found %d min" % len(self.arrayMin))
        print("Found %d max" % len(self.arrayMax))
        print("Found %d saddle" % len(self.arraySad))

        print("Creating nanoCavities objects")
        self.populateNanoCavities()
        print('Nbr of nanoCavities found =', len(self.list_nano_cavities))
        if (is_filter_equivalent_cavities):
            print('Nbr of nanoCavities before filtering =', len(self.list_nano_cavities))
            self.filterEquivalentCavities()
            print('Nbr of nanoCavities types =', len(self.list_nano_cavities))
            # self.print_info_cavities()

    def createPerfectSRG(self, nbOfGrating, depth=1, groove=1, negative=False):
        #FIXME
        if negative:
            depth = - depth

        s = np.ones(nbOfGrating) * depth
        p = np.ones(nbOfGrating) * groove
        phase = PI / nbOfGrating
        list_gratings = []
        for i in range(nbOfGrating):
            g = sinusoidalGrating(s[i], i * phase, p[i])
            list_gratings.append(g)

        return quasiCrystal(list_gratings, res, dimX, dimY, offset, offset, tol, is_filter_equivalent_cavities=True)

    def findCriticalPoint(self):
        """
        Scan the quasiCrystal surface with a resolution of self.min_pt_search_resolution in order to find local critical point,
        i.e. points where the first 2D derivative is zero.
        Put this points in self.list_critical_pt, after being sure that this critical point is not already there
        """

        def distance(xc, yc, pt):
            return math.sqrt((xc - pt[0])**2 + (yc - pt[1])**2)

        def findZeroDerivative(p):
            x, y = p

            # def partialDx(x, y):
            #     sol = 0
            #     for g in self.gratings:
            #         sol += g.partialDx(x, y)
            #     return sol
            #
            # def partialDy(x, y):
            #     sol = 0
            #     for g in self.gratings:
            #         sol += g.partialDy(x, y)
            #     return sol

            Dx = 0
            for g in self.gratings:
                Dx += g.partialDx(x, y)

            Dy = 0
            for g in self.gratings:
                Dy += g.partialDy(x, y)

            return np.array([Dx, Dy])

        def jacobian(p):
            x, y = p
            # x, y = p[0], p[1]
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

        nbPointX = int(self.sizeX / self.min_pt_search_resolution)
        nbPointY = int(self.sizeY / self.min_pt_search_resolution)

        self.list_critical_pt = []
        x0 = self.xOffset
        for idxX in range(nbPointX):
            x0 += self.min_pt_search_resolution
            y0 = self.yOffset
            for idxY in range(nbPointY):
                y0 += self.min_pt_search_resolution
                # res = fsolve(func=findZeroDerivative, x0=np.array([x0, y0]), fprime=jacobian)
                # res = fsolve(func=findZeroDerivative, x0=np.array([x0, y0]))
                # print(res)
                xc, yc = fsolve(func=findZeroDerivative, x0=(x0, y0), args=(), fprime=jacobian)
                isOK = False
                if (xc < self.sizeX + self.xOffset) and (xc > self.xOffset):
                    if (yc < self.sizeY + self.yOffset) and (yc > self.yOffset):
                        isOK = True
                        # Test if the critical point is not already in the self.listCriticalPt list +/- a fraction of the minPtSearchResolution
                        for pt in self.list_critical_pt:
                            if (distance(xc, yc, pt) < self.min_pt_search_resolution / 5.0):
                                isOK = False
                if (isOK):
                    self.list_critical_pt.append([xc, yc, "n"])

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

        for criticalPt in self.list_critical_pt:
            xc = criticalPt[0]
            yc = criticalPt[1]
            sDxx = sDyy = sDxy = 0
            z = 0
            for g in self.gratings:
                sDxx += g.partialDxx(xc, yc)
                sDyy += g.partialDyy(xc, yc)
                sDyy += g.partialDxy(xc, yc)
                z += g.getGrattingValue(xc, yc)

            D = sDxx * sDyy - (sDxy) ** 2
            if (D > 0) and (sDxx > 0):
                # Minimum
                # Test if the minimum is not on the edges of the sample, if it is the case, discard this point.
                if (xc > self.sizeX*0.02 + self.xOffset) and (xc < 0.98 * self.sizeX - self.xOffset):
                    if (yc > self.sizeY * 0.02 + self.yOffset ) and (yc < 0.98 * self.sizeY - self.yOffset):
                        listMin.append((xc, yc, z))
                        criticalPt[2] = 'm'
            elif (D > 0) and (sDxx < 0):
                # Maximum
                listMax.append((xc, yc, z))
                criticalPt[2] = 'M'
            elif D < 0:
                # Saddle
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
            self.list_nano_cavities.append(
                nanoCavity(self.gratings, idxMin, self.arrayMin, self.arrayMax, self.arraySad, self.min_pt_search_resolution / 10.0))

    def filterEquivalentCavities(self):
        """
        Search and remove among the minimum list, the cavity that are similar.
        Due to the symmetry, if we increase the size of the quasiCrystal (sizeX and sizeX parameters) we will artificially increase the number of cavities.
        """

        # First : eliminate all the cavities that has no neighbors and that have failed the MVEE analysis
        nb_of_cavities = len(self.list_nano_cavities)
        i = 0
        while i < nb_of_cavities:
            if (self.list_nano_cavities[i].ellipsoid == -1) \
                    or (self.list_nano_cavities[i].nbOfNeighbors == 0):
                self.list_nano_cavities.pop(i)
                # nbOfwell -= 1 ?
                nb_of_cavities = len(self.list_nano_cavities)
            else:
                i += 1


        # Secondly remove cavities that are on the edge of the reconstructed surface
        # We search if the distance to the edge is inferior to the distance to the closest minimum

        nb_of_cavities = len(self.list_nano_cavities)
        i = 0
        while i < nb_of_cavities:
            d = np.sqrt(self.list_nano_cavities[i].distNextMinSquared)
            xc, yc = self.list_nano_cavities[i].xc, self.list_nano_cavities[i].yc
            if (xc + d > self.sizeX + self.xOffset) or (xc - d < -self.sizeX - self.xOffset) or (yc + d > self.sizeY + self.yOffset) or (yc - d < -self.sizeY - self.yOffset):
                self.list_nano_cavities.pop(i)
                nb_of_cavities = len(self.list_nano_cavities)
            else:
                i += 1

        # Thirdly two cavities are declared equivalent when they have:
        #   -  the same number of neighbors
        #   -  the same number of Max
        #   -  the same number of Saddle
        #   -  The relative difference in ellipsoid volume, eccentricity and ellipticiy are smaller than self.maxRelativeDiffForEquivalentFilter


        nb_of_cavities = len(self.list_nano_cavities)
        i = 0   #currentCavity
        while (i < nb_of_cavities):
            j = i + 1
            while (j < nb_of_cavities):
                do_filter = False

                #   -  the same number of neighbors
                #   -  the same number of Max
                #   -  the same number of Saddle
                if self.list_nano_cavities[j].nbOfNeighbors == self.list_nano_cavities[i].nbOfNeighbors:
                    if self.list_nano_cavities[j].nbOfNeighborsMax == self.list_nano_cavities[i].nbOfNeighborsMax:
                        if self.list_nano_cavities[j].nbOfNeighborsSad == self.list_nano_cavities[i].nbOfNeighborsSad:
                            i_ellipsoid = self.list_nano_cavities[i].ellipsoid
                            j_ellipsoid = self.list_nano_cavities[j].ellipsoid

                            relative_eccentricity_diff = abs(i_ellipsoid.eccentricity - j_ellipsoid.eccentricity) / max(i_ellipsoid.eccentricity, j_ellipsoid.eccentricity)
                            # relative_ellipticity_diff = abs(i_ellipsoid.ellipticity - j_ellipsoid.ellipticity) / max(i_ellipsoid.ellipticity, j_ellipsoid.ellipticity)
                            relative_volume_diff = abs(i_ellipsoid.volume - j_ellipsoid.volume) / max(i_ellipsoid.volume, j_ellipsoid.volume)

                            threshold = self.max_relative_diff_for_equivalent_filter
                            if relative_eccentricity_diff < threshold:
                                if relative_volume_diff < threshold:
                                    do_filter = True

                if do_filter:
                    self.list_nano_cavities.pop(j)
                    nb_of_cavities = len(self.list_nano_cavities)
                else:
                    #NB : if we have filtered th j-th well, no need to do j +=1
                    j += 1

            i += 1

    def get_cavities_statistics(self):
        eccentricity_list = []
        volume_list = []

        listOrientationX = []
        listOrientationY = []
        listOrientationZ = []

        for cavity in self.list_nano_cavities:
            ellipsoid = cavity.ellipsoid
            eccentricity_list.append(ellipsoid.eccentricity)
            volume_list.append(ellipsoid.volume)

            o_X, o_Y, o_Z = ellipsoid.assessOrientationWithStretchConstraint()

            listOrientationX.append(o_X)
            listOrientationY.append(o_Y)
            listOrientationZ.append(o_Z)

        return np.mean(eccentricity_list), np.std(eccentricity_list), np.mean(volume_list), np.std(volume_list), np.mean(listOrientationX), np.std(listOrientationX), np.mean(listOrientationY), np.std(listOrientationY), np.mean(listOrientationZ), np.std(listOrientationZ)

    def drawContour(self, nb_of_level=50, fname=None, title=None, dots_critical=False, dots_cavity=False):
        """
        Display function for getting a countour profile od the surface
        :param nb_of_level: Nb of Level for the countour
        :param fname: file name
        :param title: Title for the graph
        :return:
        """
        X = np.arange(self.xOffset, self.sizeX + self.xOffset, self.min_pt_search_resolution)
        Y = np.arange(self.yOffset, self.sizeY + self.yOffset, self.min_pt_search_resolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.gratings[0].getGrattingValue(X, Y)
        for i in range(1, len(self.gratings)):
            Z += self.gratings[i].getGrattingValue(X, Y)


        plt.figure()
        # cp = plt.contourf(X, Y, Z, nb_of_level, cmap=cm.jet)
        cp = plt.contourf(X, Y, Z, nb_of_level, cmap=cm.Greys)
        #Only one level
        #cp = plt.contourf(X, Y, Z, [0, 0.1])
        plt.colorbar(cp)
        #plt.title('Filled Contours Plot')
        plt.xlabel('x (µm)', size=20)
        plt.ylabel('y (µm)', size=20)
        if title is not None:
            plt.title(title, size=20)

        if dots_critical:
            for crPt in self.list_critical_pt:
                if crPt[2] == 'm':
                    plt.scatter(crPt[0], crPt[1], color='b')
                elif crPt[2] == 's':
                    plt.scatter(crPt[0], crPt[1], color='g')
                elif crPt[2] == 'M':
                    plt.scatter(crPt[0], crPt[1], color='r')

        if dots_cavity:
            for cavity in self.list_nano_cavities:
                plt.scatter(cavity.xc, cavity.yc, color='b')

        if fname is not None:
            plt.savefig(fname, dpi=600)
        plt.show()

        return cp

    def drawAllQuasiCrystal(self, dots=True, coat='full', coatalpha=0.5, ellipsoid=False, wellidx=[], ellipsesalpha=0.5, file_save_Path=None, title=None):
        """
        Draw the qausi-crystal surface with the critical points and the ellipsoids obtain from MVEE

        :param dots: Draw the criticals point
        :param coat: How to draw the surface, full or wire
        :param coatalpha:
        :param ellipsoid: Draw the ellipsoid obtained via MVEE in the nanocavities
        :param wellidx:
        :param ellipsesalpha: alpha blending for the ellipse
        :param file_save_Path: File path for saving the graph
        :return:
        """
        fig = plt.figure()
        ax = Axes3D(fig)
        # plt.hold(True)

        X = np.arange(self.xOffset, self.sizeX + self.xOffset, self.min_pt_search_resolution)
        Y = np.arange(self.yOffset, self.sizeY + self.yOffset, self.min_pt_search_resolution)
        X, Y = np.meshgrid(X, Y)

        Z = self.gratings[0].getGrattingValue(X, Y)
        for i in range(1, len(self.gratings)):
            Z += self.gratings[i].getGrattingValue(X, Y)


        #TODO le scaling, ici on considere l'amplitude max à partir de la somme des réseaux s'ils sont tous en phase
        #ax.set_zlim(-len(self.gratings) - 0.01, len(self.gratings) + 0.01)
        # ax.zaxis.set_major_locator(LinearLocator(10))
        # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        ax.set_zticks([])

        ax.set_xlabel("x /µm")
        ax.set_ylabel("y /µm")

        ax.view_init(60, 120)

        if dots:
            for crPt in self.list_critical_pt:
                zc = self.gratings[0].getGrattingValue(crPt[0], crPt[1])
                for i in range(1, len(self.gratings)):
                    zc += self.gratings[i].getGrattingValue(crPt[0], crPt[1])
                if crPt[2] == 'm':
                    ax.scatter(crPt[0], crPt[1], zc, color='b')
                elif crPt[2] == 's':
                    ax.scatter(crPt[0], crPt[1], zc, color='g')
                elif crPt[2] == 'M':
                    ax.scatter(crPt[0], crPt[1], zc, color='r')

        if ellipsoid:
            sin, cos = np.sin, np.cos
            for nano_cavity in self.list_nano_cavities:
                ellipsoid = nano_cavity.ellipsoid
                if ellipsoid is not None and ellipsoid != -1:
                    # rx, ry, rz = 1. / np.sqrt(self.listNanoCavities[idx].ellipsoid.D)
                    u, v = np.mgrid[0:2 * PI:20j, -PI / 2:PI / 2:10j]

                    # cf parametric description of an ellipsoid in wikipedia
                    def ellipse(RX, RY, RZ, u, v):
                        x = RX * cos(u) * cos(v)
                        y = RY * sin(u) * cos(v)
                        z = RZ * sin(v)
                        return x, y, z

                    RX = ellipsoid.rx
                    RY = ellipsoid.ry
                    RZ = ellipsoid.rz

                    # cosmetic
                    # RX /=2
                    # RY /= 2
                    # RZ /= 2

                    # Stack arrays in sequence depth wise (along third axis). Takes a sequence of
                    # arrays and stack them along the third axis to make a single array
                    E = np.dstack(ellipse(RX, RY, RZ, u, v))
                    # self.ellipsoid.V comes from the SVD and contains the normalized eigen vectors
                    E = np.dot(E, ellipsoid.V) + ellipsoid.center
                    x, y, z = np.rollaxis(E, axis=-1)
                    ax.plot_surface(x, y, z, cstride=1, rstride=1, color='y', alpha=ellipsesalpha)

        if (coat == 'full'):
            surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.coolwarm, alpha=coatalpha, linewidth=0,
                                   antialiased=True)
            # fig.colorbar(surf, shrink=0.5, aspect=5)

        elif (coat == 'wire'):
            ax.plot_surface(X, Y, Z, rstride=2, cstride=2, alpha=coatalpha)

        if title is not None:
            plt.title(title, size=20)

        if file_save_Path is not None:
            plt.savefig(file_save_Path, dpi=300)
        plt.show()

    def draw_cavities(self):
        i = 0
        if self.name is not None:
            base_name = self.name
        else:
            base_name = ""
        for cavity in self.list_nano_cavities:
            fname = base_name + "_cavity_num" + str(i) +".png"
            cavity.draw(dots=False, coat='full', coatAlpha=0.1, drawEllipsoid=True, ellipseAlpha=0.2, fname=fname, is_draw_princ_axes=True)
            i += 1


    def print_info_cavities(self):
        for i in range(len(self.list_nano_cavities)):
            print("cavity nb :", i)
            self.list_nano_cavities[i].print_info()



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

