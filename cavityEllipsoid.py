import numpy.linalg as la
import numpy as np

PI = np.pi



class cavityEllipsoid():
    """
    Fits an ellipsoid using MVEE inside a list of point
    Greatly inspired from :
    http://stackoverflow.com/questions/14016898/port-matlab-bounding-ellipsoid-code-to-python
    """
    def __init__(self, points, tol=1 / 1000.0, size='half'):
        self.rx, self.ry, self.rz = 0, 0, 0
        self.ellipseX, self.ellipseY, self.ellipseZ = 0, 0, 0
        self.A, self.center = [], []
        self.deltaX, self.deltaY, self.deltaZ = [], [], []
        self.a, self.b, self.c, self.volume = 0, 0, 0, 0
        self.eccentricity, self.ellipticity = 0, 0
        self.size = size

        self.MVEE(points, tol)
        self.getEllipsoidParameters()

    def MVEE(self, points, tol):
        """
        Moshtagh, N. (2005). Minimum volume enclosing ellipsoids. Convex Optimization, (January), 1–9. http://doi.org/10.1.1.116.7691
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


    def getEllipsoidParameters(self):
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

    def printParameters(self):

        print("pos : ",self.center)
        print("a : ", self.a)
        print("b : ", self.b)
        print("c : ", self.c)
        print("dir a : ", self.principalAxes[0])
        print("dir b : ", self.principalAxes[1])
        print("dir c : ", self.principalAxes[2])

