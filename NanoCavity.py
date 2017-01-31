import matplotlib.pyplot as plt
import numpy as np
from cavityEllipsoid import cavityEllipsoid
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

PI, sin, cos = np.pi, np.sin, np.cos


class nanoCavity():
    def __init__(self, listGrating, minLoc, arrayMin, arrayMax, arraySad, drawResolution=0.05):
        """
        Classe "puit" qui correspond  au puit centré en un minimum minLoc de coordonées
        xc = arrayMin[minLoc][0] et yc = arrayMin[minLoc][1] et qui a pour attribut le nombre de ses voisin max
        ou saddles, la disatnce les séparant du min local, la largeur et la hauteur du puit
        ainsi que la distance au prochain min.
        """
        self.minLoc = minLoc
        self.Grating = listGrating
        self.drawResolution = drawResolution
        self.xc = arrayMin[minLoc, 0]
        self.yc = arrayMin[minLoc, 1]
        self.zc = arrayMin[minLoc, 2]
        self.nbOfNeighborsMax, self.nbOfNeighborsSad, self.nbOfNeighbors = 0, 0, 0
        self.distNextMinSquared, self.heightMin, self.heightMax, self.widthMax, self.widthMin = 0, 0, 0, 0, 0
        self.neighbors, self.dist2Min, self.ellipsoid = [], [], []

        self.findDistLim(arrayMin)
        self.findNeighbors(arrayMin, arrayMax, arraySad)
        self.wellDimensions()



    def findDistLim(self, arrayMin):
        """
        Calcule la distance au carré  entre ce min et tout les autre min puis récupère la plus petite valeure
        D² = x_c² - 2.x_i.x_c + x_i²     +     y_c² - 2.y_i.y_c + y_i²

        On calcul le carré pour eviter d'avoir à calculer la racine. Comme on cherche juste le plus proche,le carré est suffisant.
        """
        idxMin = self.minLoc


        #TODO à tester
        #    D²   =  (      x_c²        -            x_i² )²     +   (      y_c²        -            y_i² )²
        # distLim = (arrayMin[idxMin, 0]  - arrayMin[:, 0]) ** 2 +  (arrayMin[idxMin, 1]  - arrayMin[:, 1]) ** 2

        #x_c = arrayMin[idxMin, 0]
        #y_c = arrayMin[idxMin, 1]
        #x_i = arrayMin[:, 1]
        #y_i = arrayMin[:, 1]

        #DsquareArray = (x_c - x_i) ** 2 +  (y_c - y_i) ** 2


        distLim = arrayMin[idxMin, 0] ** 2 + arrayMin[idxMin, 1] ** 2 + arrayMin[:, 1] ** 2 + arrayMin[:, 0] ** 2 \
                  - 2 * (arrayMin[idxMin, 0] * arrayMin[:, 0] + arrayMin[idxMin, 1] * arrayMin[:, 1])

        distLim[:].sort()

        # A 2D tester si point plus près et plus haut qui "masque" ce voisin.

        # Pourquoi le 2/3 ? et  distLim[1], parce que le 0 c'est lui même avec D²=0
        #self.distLim = (2 / 3) * distLim[1]
        #parce que le 0 c'est lui même avec D²=0
        self.distNextMinSquared = distLim[1]

    def findNeighbors(self, arrayMin, arrayMax, arraySad):
        """
        Pour trouver les plus proche voisins il suffit de calculer la distance entre le point iMin et les points Max
        et aussi le point iMin et les Saddle en enmployant la même formule que dans 'findDistLim'.
        Enduite on cherche les indices des distances Max ou Sad inférieur à la distance limite
        Enfin on crée un tableau 'neighbors' contenant les x et y des points au indices concernés
        """
        #x_c = arrayMin[idxMin, 0]
        #y_c = arrayMin[idxMin, 1]


        #x_max_i = arrayMax[:, 0]
        #y_max_i = arrayMax[:, 1]
        #dist2AllMax = (x_c - x_max_i) ** 2 +  (y_c - y_max_i) ** 2



        dist2AllMax = arrayMin[self.minLoc, 0] ** 2 + arrayMin[self.minLoc, 1] ** 2 + arrayMax[:, 0] ** 2 \
                      + arrayMax[:, 1] ** 2 - 2 * (
            arrayMin[self.minLoc, 0] * arrayMax[:, 0] + arrayMin[self.minLoc, 1] * arrayMax[:, 1])
        idxMax = np.where(dist2AllMax < self.distNextMinSquared)
        neighborsMax = arrayMax[idxMax[0], :]

        #x_sad_i = arrayMax[:, 0]
        #y_sad_i = arrayMax[:, 1]
        #dist2AllSad = (x_c - x_sad_i) ** 2 +  (y_c - y_sad_i) ** 2


        dist2AllSad = arrayMin[self.minLoc, 0] ** 2 + arrayMin[self.minLoc, 1] ** 2 + arraySad[:, 0] ** 2 \
                      + arraySad[:, 1] ** 2 - 2 * (
            arrayMin[self.minLoc, 0] * arraySad[:, 0] + arrayMin[self.minLoc, 1] * arraySad[:, 1])
        idxSad = np.where(dist2AllSad < self.distNextMinSquared)
        neighborsSad = arraySad[idxSad[0], :]

        self.nbOfNeighborsMax = np.size(idxMax[0])
        self.nbOfNeighborsSad = np.size(idxSad[0])
        self.nbOfNeighbors = self.nbOfNeighborsMax + self.nbOfNeighborsSad

        self.neighbors = []

        for j in range(len(neighborsMax)):
            #TODO
            # A 2D tester si point plus près et plus haut qui "masque" ce voisin.
            toAppend = False
            self.neighbors.append((neighborsMax[j, 0], neighborsMax[j, 1], neighborsMax[j, 2]))
        for i in range(len(neighborsSad)):
            self.neighbors.append((neighborsSad[i, 0], neighborsSad[i, 1], neighborsSad[i, 2]))
        self.neighbors = np.asarray(self.neighbors)

        self.nbOfNeighbors = len(self.neighbors)

        if (self.nbOfNeighbors != 0):
            self.dist2Min = (self.xc - self.neighbors[:, 0]) ** 2 + (self.yc - self.neighbors[:, 1]) ** 2

    def wellDimensions(self):
        """
        On choisit la taille comme suit :
        La hauteur est le plus petit G(x,y) entre les Saddle et les Max
        les dimensions 2D sont données par l'ellipse fittée sur tout les voisins (Sad et Max)
        le codes est expliquer et inspirer de :
        http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html

        On utilise plutôt maintenant la MVEE
        """

        # Pourquoi 5, y-a-til une limite à la MVEE ?
        #if (len(self.neighbors) >= 5):
        # A 2 on obtient Singular matrix ac 3 points on obtient un plan.
        if (len(self.neighbors) >= 3):
            # La MVEE donne aussi les dimension rx, ry, rz de l'ellipse englobante
            #def __init__(self, points, tol=1 / 1000, size='half'):
            #mveePoint = self.neighbors + [self.xc, self.yc, self.zc]

            #add the minimum point to the MVEE points
            #print(np.shape(self.neighbors))
            center = np.array([[self.xc, self.yc, self.zc]])
            #print(np.shape(center))
            #print(np.shape(self.neighbors))
            #mveePoint = np.concatenate(self.neighbors, center)
            mveePoints = np.vstack((self.neighbors, center))

            #print(mveePoints)

            self.ellipsoid = cavityEllipsoid(mveePoints)
        else:
            self.ellipsoid = -1

        #Obsolete -> replaced by MVEE
        # if (len(self.neighbors) != 0):
        #     # La largeur min est la distance minimale non nulle entre le min et le Max/Sad (dans le plans x0y)
        #     neighborsDist = self.dist2Min
        #     neighborsDist[:].sort()
        #     self.widthMin = np.sqrt(neighborsDist[0])
        #     self.widthMax = np.sqrt(neighborsDist[-1])
        #     # La hauteur min est la distance minimale non nulle entre le min et le Max/Sad (dans le plans x0y)
        #     neighborsHeight = self.neighbors[:, 2]
        #     neighborsHeight[:].sort()
        #     self.heightMin = abs(neighborsHeight[0] - self.zc)
        #     self.heightMax = abs(neighborsHeight[-1] - self.zc)
        #
        # elif (self.nbOfNeighbors == 0):
        #     self.heightMax, self.heightMin = -1, -1
        #     self.widthMax, self.widthMin = -1, -1

    def drawWell(self, dots=True, coat='full', coatAlpha=0.5, ellipse=False, ellipseAlpha=0.5):
        if (self.nbOfNeighbors != 0):
            fig = plt.figure()
            ax = Axes3D(fig)
            plt.hold(True)

            #TODO surface pas en accord avec les points.

            distNextMin = np.sqrt(self.distNextMinSquared)
            X = np.arange(self.xc - distNextMin, self.xc + distNextMin, distNextMin / 20.)
            Y = np.arange(self.yc - distNextMin, self.yc + distNextMin, distNextMin / 20.)
            X, Y = np.meshgrid(X, Y)

            Z = self.Grating[0].getGrattingValue(X, Y)
            for i in range(1, len(self.Grating)):
                Z += self.Grating[i].getGrattingValue(X, Y)

            #ax.set_zlim(-len(self.Grating) - 0.01, len(self.Grating) + 0.01)
            ax.zaxis.set_major_locator(LinearLocator(10))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

            # # Dessine les Points :
            # #Oops we don't know anymore if we have a sad or a max (can't be a min because we have chosen the max and sad with a distance inferior to the closest min
            # if (dots):
            #     tab = []
            #     tab.append((self.xc, self.yc, 'm'))
            #     idx = 0
            #     while idx <= self.nbOfNeighborsMax:
            #         tab.append((self.neighbors[idx, 0], self.neighbors[idx, 1], 'M'))
            #         idx += 1
            #     # while idx <= self.nbOfNeighbors:
            #     #     tab.append((self.neighbors[idx, 0], self.neighbors[idx, 1], 's'))
            #     #     idx += 1
            #     #tab.append((self.xc, self.yc, 'm'))
            #
            #     for Pt in tab:
            #         zc = self.Grating[0].getGrattingValue(Pt[0], Pt[1])
            #         for i in range(1, len(self.Grating)):
            #             zc += self.Grating[i].getGrattingValue(Pt[0], Pt[1])
            #         if Pt[2] == 'm':
            #             ax.scatter(Pt[0], Pt[1], zc, color='b')
            #         elif Pt[2] == 's':
            #             ax.scatter(Pt[0], Pt[1], zc, color='g')
            #         elif Pt[2] == 'M':
            #             ax.scatter(Pt[0], Pt[1], zc, color='r')

            if (dots):
                for pt in self.neighbors:
                    ax.scatter(pt[0], pt[1], pt[2], color='r')

                ax.scatter(self.xc, self.yc, self.zc, color='k')



            # Dessine l'ellipse :
            if (ellipse) and (self.ellipsoid != -1):
                u, v = np.mgrid[0:2 * PI:20j, -PI / 2:PI / 2:10j]


                #cf parametric description in wiki
                def ellipse(RX, RY, RZ, u, v):
                    x = RX * cos(u) * cos(v)
                    y = RY * sin(u) * cos(v)
                    z = RZ * sin(v)
                    return x, y, z

                #Est-ce que cela ne donne pas les demi demi Axe ?
                RX = self.ellipsoid.rx
                RY = self.ellipsoid.ry
                RZ = self.ellipsoid.rz

                #Stack arrays in sequence depth wise (along third axis). Takes a sequence of arrays and stack them along the third axis to make a single array
                E = np.dstack(ellipse(RX, RY, RZ, u, v))
                # self.ellipsoid.V provient de la SVD
                E = np.dot(E, self.ellipsoid.V) + self.ellipsoid.center
                x, y, z = np.rollaxis(E, axis=-1)
                ax.plot_surface(x, y, z, cstride=1, rstride=1, alpha=ellipseAlpha, color='y')

                center = self.ellipsoid.center
                #dessine le centre de l'ellipse :
                ax.scatter(center[0], center[1], center[2], color='k')

                #Dessine les axes principaux

                for v in self.ellipsoid.principalAxes:
                   ax.plot3D([center[0],v[0] + center[0]], [center[1],v[1] + center[1]], [center[2],v[2] + center[2]], color='red', alpha=0.8, lw=3)




            # Dessine la topographie :
            if (coat == 'full'):
                surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=coatAlpha, cmap=cm.coolwarm, linewidth=0,
                                       antialiased=False)
                fig.colorbar(surf, shrink=0.5, aspect=5)
            elif (coat == 'wire'):
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=coatAlpha)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show()


    def printInfo(self):
        print("Center : ", self.xc, self.yc, self.zc)
        print("Neighbors : ", self.neighbors)