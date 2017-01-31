import numpy.linalg as la
import numpy as np

PI = np.pi

#TODO cite Source !
#Greatly inspired from :
#http://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python

class cavityEllipsoid():
    def __init__(self, points, tol=1 / 1000.0, size='half'):
        self.rx, self.ry, self.rz = 0, 0, 0
        self.ellipseX, self.ellipseY, self.ellipseZ = 0, 0, 0
        self.A, self.center = [], []
        self.deltaX, self.deltaY, self.deltaZ = [], [], []
        self.a, self.b, self.c, self.volume = 0, 0, 0, 0
        self.eccentricity, self.ellipticity = 0, 0
        self.size = size

        self.MVEE(points, tol)
        #self.dist2centre(points)
        self.parameters()

    def MVEE(self, points, tol):
        """
        Finds the ellipse equation in "center form"
        (x-c).T * A * (x-c) = 1
        """
        N, d = points.shape
        Q = np.column_stack((points, np.ones(N))).T
        u = np.ones(N) / N
        err = 1.0 + tol

        while err > tol:
            # assert u.sum() == 1 # invariant
            X = np.dot(np.dot(Q, np.diag(u)), Q.T)
            M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
            jdx = np.argmax(M)
            step_size = (M[jdx] - d - 1.0) / ((d + 1) * (M[jdx] - 1.0))
            new_u = (1 - step_size) * u
            new_u[jdx] += step_size
            err = la.norm(new_u - u)
            u = new_u
        c = np.dot(u, points)
        A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points) - np.multiply.outer(c, c)) / d
        #return A, c

        self.center = c
        self.A = A

    def dist2centre(self, points):
        """
        Obtain distance to center coordinates for the entire x,y,z array passed.
        """
        self.delta_x, self.delta_y, self.delta_z = (points[:, 0] - self.center[0]), (points[:, 1] - self.center[1]), (
            points[:, 2] - self.center[2])
        # Distance entre l'ellipse et l'origine du repère :
        self.dist = np.sqrt(self.delta_x ** 2 + self.delta_y ** 2 + self.delta_z ** 2)

    def parameters(self):
        """
        From wikipédia
        (x-V)T A (x-V) = 1
        De plus, les vecteurs propres de A définissent les axes de l'ellipsoïde et les valeurs propres de A sont égales à l'inverse du carré des demi-axes
        :return:
        """
        # Données extirpé de la matrice solution de l'ellipse (où D est la diagonal de A et composé des axes a, b, c

        #alculating the SVD consists of finding the eigenvalues and eigenvectors of A At and At A.
        # The eigenvectors of At A make up the columns of V , the eigenvectors of A At  make up the columns of U.
        # Also, the singular values in S are square roots of eigenvalues from A At or At A.
        # The singular values are the diagonal entries of the S matrix and are arranged in descending order.
        # The singular values are always real numbers. If the matrix A is a real matrix, then U and V are also real.

        U, self.D, self.V = la.svd(self.A)

        self.eigenValue, self.principalAxes = la.eig(self.A)

        #w, v = la.eig(A)

        #inverse du carré des ->demi<--axes à partir des elements diagonaux de la matrice A
        self.rx, self.ry, self.rz = 1. / np.sqrt(self.D)






        #self.a, self.b, self.c = max(self.rx, self.ry), min(self.rx, self.ry), self.rz
        self.a, self.b, self.c = self.rx, self.ry, self.rz

        #Attention les valeurs propres données par la svd ou le eig sont-ils dans le même ordre (et surtout les mêmes).
        if self.rx > self.ry:
            princAxeX = self.principalAxes[0]
            princAxeY = self.principalAxes[1]
        else :
            princAxeX = self.principalAxes[1]
            princAxeY = self.principalAxes[0]


        #self.VecZ = self.principalAxes[2]

        # Eccentricité et Ellipticité dans le plan :

        #L'ellipticité est une mesure de l'aplatissement d'une ellipse. Elle est comprise entre les valeurs 0 et 1,
        # le premier cas correspondant à un cercle et le second à une ellipse infiniment allongée, c'est - à - dire un segment.


        if self.a < self.b :
            self.eccentricity = np.sqrt(self.b ** 2 - self.a ** 2) / self.b
            self.ellipticity = 1 - self.a / self.b
        else :
            self.eccentricity = np.sqrt(self.a ** 2 - self.b ** 2) / self.a
            self.ellipticity = 1 - self.b / self.a
        self.volume = (4/3)*PI*self.a*self.b*self.c

    def printData(self):
        print("pos : ",self.center)
        print("a : ", self.a)
        print("b : ", self.b)
        print("c : ", self.c)
        print("dir a : ", self.principalAxes[0])
        print("dir b : ", self.principalAxes[1])
        print("dir c : ", self.principalAxes[2])

    def testOrientation(self, v):
        if(self.a > self.b):
            alignementFactor = abs(np.dot(self.principalAxes[0], v))
        else:
            alignementFactor = abs(np.dot(self.principalAxes[1], v))
        #pondération par l'ellipticité
        return (self.ellipticity  * alignementFactor)


    def testOrientation_old(self, v):
        """
        Ne donne pas les résultats attendus, pb de normalisation ?
        """
        axeA = abs(self.a * np.dot(self.principalAxes[0], v))
        axeB = abs(self.b * np.dot(self.principalAxes[1], v))
        axeC = abs(self.c * np.dot(self.principalAxes[2], v))
        return (axeA + axeB + axeC)

