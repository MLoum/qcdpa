import matplotlib.pyplot as plt
import numpy as np
from cavityEllipsoid import cavityEllipsoid
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

PI, sin, cos = np.pi, np.sin, np.cos


class nanoCavity():
    """
    Class representing the neighborhood of a local minimum in the quasiCrystal.
    This neighborhood is composed of :
        - a list of the closest MaximumPoint
        - a list of the closest Saddle Point
        - a ellipsoid that is enclosing all the previous listed points

    "Closest" means a distance smaller than the distance from the local minimum to the closest other minimum on the quasiCrystal
    (their positions are inside the arrayMin list)
    """
    def __init__(self, listGrating, minIdx, arrayMin, arrayMax, arraySad, drawResolution=0.05):
        """
        :param listGrating:
        :param minIdx:
        :param arrayMin:
        :param arrayMax:
        :param arraySad:
        :param drawResolution:
        """
        self.minIdx = minIdx
        self.Grating = listGrating
        self.drawResolution = drawResolution
        self.xc, self.yc, self.zc = arrayMin[minIdx, :]
        self.nbOfNeighborsMax, self.nbOfNeighborsSad, self.nbOfNeighbors = 0, 0, 0
        self.distNextMinSquared, self.heightMin, self.heightMax, self.widthMax, self.widthMin = 0, 0, 0, 0, 0
        self.neighbors, self.dist2Min,  = [], []
        self.ellipsoid = None

        self.findDistLim(arrayMin)
        self.findNeighbors(arrayMin, arrayMax, arraySad)
        self.fitEllipsoidInCavity()



    def findDistLim(self, arrayMin):
        """
        Find the closest other minimum from the position of the nanocavities
        This will define the upper limit for the distance between the max and saddle neighbors

        On calcul le carré pour eviter d'avoir à calculer la racine. Comme on cherche juste le plus proche,le carré est suffisant.
        """
        idxMin = self.minIdx


        #TODO à tester
        #    D²   =  (      x_c²        -            x_i² )²     +   (      y_c²        -            y_i² )²
        # distLim = (arrayMin[idxMin, 0]  - arrayMin[:, 0]) ** 2 +  (arrayMin[idxMin, 1]  - arrayMin[:, 1]) ** 2

        x_c, y_c = self.xc, self.yc
        x_i = arrayMin[:, 0]
        y_i = arrayMin[:, 1]

        DsquareArray = (x_c - x_i) ** 2 +  (y_c - y_i) ** 2


        # distLim = arrayMin[idxMin, 0] ** 2 + arrayMin[idxMin, 1] ** 2 + arrayMin[:, 1] ** 2 + arrayMin[:, 0] ** 2 \
        #           - 2 * (arrayMin[idxMin, 0] * arrayMin[:, 0] + arrayMin[idxMin, 1] * arrayMin[:, 1])

        DsquareArray[:].sort()


        #NB : distLim[0] is the distance between the cavity and itself and is of course D²=0
        self.distNextMinSquared = DsquareArray[1]

    def findNeighbors(self, arrayMin, arrayMax, arraySad):
        """
        Pour trouver les plus proche voisins il suffit de calculer la distance entre le point iMin et les points Max
        et aussi le point iMin et les Saddle en enmployant la même formule que dans 'findDistLim'.
        Enduite on cherche les indices des distances Max ou Sad inférieur à la distance limite
        Enfin on crée un tableau 'neighbors' contenant les x et y des points au indices concernés
        :param arrayMin:
        :param arrayMax:
        :param arraySad:
        :return:
        """

        x_c, y_c = self.xc, self.yc

        x_max_i = arrayMax[:, 0]
        y_max_i = arrayMax[:, 1]
        dist2AllMax = (x_c - x_max_i) ** 2 +  (y_c - y_max_i) ** 2


        # dist2AllMax = arrayMin[self.minIdx, 0] ** 2 + arrayMin[self.minIdx, 1] ** 2 + arrayMax[:, 0] ** 2 \
        #               + arrayMax[:, 1] ** 2 - 2 * (
        #     arrayMin[self.minIdx, 0] * arrayMax[:, 0] + arrayMin[self.minIdx, 1] * arrayMax[:, 1])
        idxMax = np.where(dist2AllMax < self.distNextMinSquared)
        neighborsMax = arrayMax[idxMax[0], :]

        x_sad_i = arraySad[:, 0]
        y_sad_i = arraySad[:, 1]
        dist2AllSad = (x_c - x_sad_i) ** 2 +  (y_c - y_sad_i) ** 2


        #dist2AllSad = arrayMin[self.minIdx, 0] ** 2 + arrayMin[self.minIdx, 1] ** 2 + arraySad[:, 0] ** 2 + arraySad[:, 1] ** 2 - 2 * (arrayMin[self.minIdx, 0] * arraySad[:, 0] + arrayMin[self.minIdx, 1] * arraySad[:, 1])
        idxSad = np.where(dist2AllSad < self.distNextMinSquared)
        neighborsSad = arraySad[idxSad[0], :]

        self.nbOfNeighborsMax = np.size(idxMax[0])
        self.nbOfNeighborsSad = np.size(idxSad[0])
        self.nbOfNeighbors = self.nbOfNeighborsMax + self.nbOfNeighborsSad

        self.neighbors = []

        for j in range(len(neighborsMax)):
            #TODO A 2D tester si point plus près et plus haut qui "masque" ce voisin.
            # toAppend = False
            self.neighbors.append((neighborsMax[j, 0], neighborsMax[j, 1], neighborsMax[j, 2]))
        for i in range(len(neighborsSad)):
            self.neighbors.append((neighborsSad[i, 0], neighborsSad[i, 1], neighborsSad[i, 2]))
        self.neighbors = np.asarray(self.neighbors)

        self.nbOfNeighbors = len(self.neighbors)

        if (self.nbOfNeighbors != 0):
            self.dist2Min = (self.xc - self.neighbors[:, 0]) ** 2 + (self.yc - self.neighbors[:, 1]) ** 2

    def fitEllipsoidInCavity(self):
        """
        Using the Minimum Volume Enclosure Ellipsoid (MVEE) algorithm, we estimate the volume of the current cavity by
        fitting an ellipsoid enclosing all the neighbors points (max or saddle) and also the minimum point of the cavity.
        """

        # 3 points is the limit for the Minimum Volume Enclosure Ellipsoid (3 points give a plane).
        if (len(self.neighbors) >= 3):
            #add the minimum of the cavity  to the MVEE points
            center = np.array([[self.xc, self.yc, self.zc]])
            mveePoints = np.vstack((self.neighbors, center))

            # mveePoints = self.neighbors
            self.ellipsoid = cavityEllipsoid(mveePoints)
        else:
            self.ellipsoid = -1



    def draw(self, dots=True, coat='full', coatAlpha=0.5, drawEllipsoid=False, ellipseAlpha=0.5, title=None, fname=None, is_draw_princ_axes=False):
        """
        :param dots:
        :param coat:
        :param coatAlpha:
        :param drawEllipsoid:
        :param ellipseAlpha:
        :return:
        """
        if (self.nbOfNeighbors != 0):
            fig = plt.figure()
            ax = Axes3D(fig)
            # plt.hold(True)

            distNextMin = np.sqrt(self.distNextMinSquared)
            X = np.arange(self.xc - distNextMin, self.xc + distNextMin, distNextMin / 20.)
            Y = np.arange(self.yc - distNextMin, self.yc + distNextMin, distNextMin / 20.)
            X, Y = np.meshgrid(X, Y)

            Z = self.Grating[0].getGrattingValue(X, Y)
            for i in range(1, len(self.Grating)):
                Z += self.Grating[i].getGrattingValue(X, Y)

            ax.zaxis.set_major_locator(LinearLocator(10))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


            if (dots):
                for pt in self.neighbors:
                    ax.scatter(pt[0], pt[1], pt[2], color='r')

                ax.scatter(self.xc, self.yc, self.zc, color='k')



            # Draw the ellipsoid :
            if (drawEllipsoid) and (self.ellipsoid != -1):
                u, v = np.mgrid[0:2 * PI:20j, -PI / 2:PI / 2:10j]


                # cf parametric description of an ellipsoid in wikipedia
                def ellipse(RX, RY, RZ, u, v):
                    x = RX * cos(u) * cos(v)
                    y = RY * sin(u) * cos(v)
                    z = RZ * sin(v)
                    return x, y, z


                RX = self.ellipsoid.rx
                RY = self.ellipsoid.ry
                RZ = self.ellipsoid.rz

                #cosmetic
                # RX /=2
                # RY /= 2
                # RZ /= 2

                #Stack arrays in sequence depth wise (along third axis). Takes a sequence of
                # arrays and stack them along the third axis to make a single array
                E = np.dstack(ellipse(RX, RY, RZ, u, v))
                # self.ellipsoid.V comes from the SVD
                E = np.dot(E, self.ellipsoid.V) + self.ellipsoid.center
                x, y, z = np.rollaxis(E, axis=-1)
                ax.plot_surface(x, y, z, cstride=1, rstride=1, alpha=ellipseAlpha, color='y')

                center = self.ellipsoid.center
                #draw the center for the ellipsoid :
                ax.scatter(center[0], center[1], center[2], color='k')

                if is_draw_princ_axes:
                    # Draw the principal axes
                    print("pp axes")
                    # for v in self.ellipsoid.principalAxes:
                    #     print (v)
                    #     # ax.plot3D(xs=[center[0] - v[0], v[0] + center[0]], ys=[center[1] - v[1], v[1] + center[1]], zs=[center[2] - v[2], v[2] + center[2]], color='red', alpha=0.8, lw=5)
                    #     first_vector = v * self.ellipsoid.rx
                    #     ax.scatter(first_vector[0], first_vector[1], first_vector[2], color='b')
                    #     first_vector = v * self.ellipsoid.ry
                    #     ax.scatter(first_vector[0], first_vector[1], first_vector[2], color='b')
                    #     ax.scatter(first_vector[0], first_vector[1], first_vector[2], color='b')
                    x_vector = center + self.ellipsoid.principalAxes[0]*self.ellipsoid.rx
                    ax.scatter(x_vector[0], x_vector[1], x_vector[2], color='r')
                    y_vector = center + self.ellipsoid.principalAxes[1]*self.ellipsoid.ry
                    ax.scatter(y_vector[0], y_vector[1], y_vector[2], color='g')
                    z_vector = center + self.ellipsoid.principalAxes[2]*self.ellipsoid.rz
                    ax.scatter(z_vector[0], z_vector[1], z_vector[2], color='b')

            # Draw the topography of the nanoCavity:
            if (coat == 'full'):
                surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=coatAlpha, cmap=cm.coolwarm, linewidth=0,
                                       antialiased=False)
                fig.colorbar(surf, shrink=0.5, aspect=5)

                # Draw countour at the bottom
                # cset = ax.contourf(X, Y, Z, cmap=cm.Greys)

            elif (coat == 'wire'):
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=coatAlpha)



        ax.set_xlabel('X axis (µm)')
        ax.set_ylabel('Y axis (µm)')
        ax.set_zlabel('Z axis (nm)')

        if title is not None:
            plt.title(title)

        if fname is not None:
            plt.savefig(fname, dpi=300)

        plt.show()


    def print_info(self):
        print("Center : ", self.xc, self.yc, self.zc)
        print("nb Neighbors : ", len(self.neighbors))
        print("nb Neighbors max : ", self.nbOfNeighborsMax)
        print("nb Neighbors sad : ", self.nbOfNeighborsSad)
        if (self.ellipsoid is not None) or self.ellipsoid != -1:
            # print("ellipticity :", self.ellipsoid.ellipticity)
            print("volume :", self.ellipsoid.volume)
            print("eccentricity :", self.ellipsoid.eccentricity)