import numpy as np

PI = np.pi


class sinusoidalGrating:
    """
    Implement the sinusoidal 2D phase grating with :
    - s : peak to peak amplitude
    - alpha : angle of the grating with respect to the x axis
    - p : pitch of the grating in µm
    """
    def __init__(self, s, alpha, p):
        self.s = s
        self.alpha = alpha
        self.p = p  #in µm

    def getPos(self, x, y):
        return 2 * PI / self.p * (x * np.cos(self.alpha) + y * np.sin(self.alpha))

    def getGrattingValue(self, X, Y):
        return self.s * np.cos(self.getPos(X, Y))

    def partialDx(self, x, y):
        return - 2 * PI * self.s / self.p * np.cos(self.alpha) * np.sin(self.getPos(x, y))

    def partialDy(self, x, y):
        return - 2 * PI * self.s / self.p * np.sin(self.alpha) * np.sin(self.getPos(x, y))

    def partialDxx(self, x, y):
        return -(2 * PI * np.cos(self.alpha) / self.p) ** 2 * self.getGrattingValue(x, y)

    def partialDyy(self, x, y):
        return -(2 * PI * np.sin(self.alpha) / self.p) ** 2 * self.getGrattingValue(x, y)

    def partialDxy(self, x, y):
        return -(2 * PI / self.p) ** 2 * np.sin(self.alpha) * np.cos(self.alpha) * self.getGrattingValue(x, y)

    def partialDyx(self, x, y):
        return self.partialDxy(x, y)
