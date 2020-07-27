import numpy as np
import matplotlib.pyplot as plt
from sinusoidalGrating import sinusoidalGrating
from QuasiCrystal import quasiCrystal
from matplotlib import cm
import matplotlib.colors as colors
import shelve

PI = np.pi

class StretchSeries:

    def __init__(self, stretch_dists=None, exp_data=None, name=None, min_pt_search_resolution=0.5, sizeX=3, sizeY=3, xOffset=0, yOffset=0, max_relative_diff_for_equivalent_filter=0.1,
                 is_filter_equivalent_cavities=True):
        """

        :param stretch_dists:
        :param exp_data: if exp_data is a string, then it implicitly means that this is the path to a shelve file
        :param name:
        :param min_pt_search_resolution: in µm
        :param sizeX: in µm
        :param sizeY: in µm
        :param xOffset: in µm
        :param yOffset: in µm
        :param max_relative_diff_for_equivalent_filter:
        :param is_filter_equivalent_cavities:
        """
        if isinstance(exp_data, str):
            self.load_state(exp_data)
            return
        self.exp_data = exp_data
        self.name = name
        self.stretch_dists = stretch_dists
        self.min_pt_search_resolution = min_pt_search_resolution
        self.sizeX, self.sizeY = sizeX, sizeY
        self.xOffset, self.yOffset = xOffset, yOffset
        self.max_relative_diff_for_equivalent_filter = max_relative_diff_for_equivalent_filter
        self.is_filter_equivalent_cavities = is_filter_equivalent_cavities
        self.nb_stretch = self.exp_data.nb_stretch
        self.nb_gratings = self.exp_data.nb_gratings
        self.SRG_dict = {}
        self.create_SRG_from_exp_data()

    def create_SRG_from_exp_data(self):
        for i_s in range(self.nb_stretch):
            print("Stretch nb %d :" % i_s)
            list_gratings = []
            amps = self.exp_data.amps[i_s]
            angles = self.exp_data.angles[i_s]
            pitches = self.exp_data.pitches[i_s]
            for j_g in range(self.nb_gratings):
                list_gratings.append(sinusoidalGrating(amps[j_g], angles[j_g], pitches[j_g]))

            self.SRG_dict[i_s] = quasiCrystal(list_gratings, self.min_pt_search_resolution, self.sizeX, self.sizeY, self.xOffset, self.yOffset, self.max_relative_diff_for_equivalent_filter, self.is_filter_equivalent_cavities)

    def save_state(self, savefile_path):
        self.shelf = shelve.open(savefile_path, 'n')  # n for new

        self.shelf['exp_data'] = self.exp_data
        self.shelf['name'] = self.name
        self.shelf['min_pt_search_resolution'] = self.min_pt_search_resolution
        self.shelf['sizeX'], self.shelf['sizeY'] = self.sizeX, self.sizeY
        self.shelf['xOffset'], self.shelf['yOffset'] = self.xOffset, self.yOffset
        self.shelf['max_relative_diff_for_equivalent_filter'] = self.max_relative_diff_for_equivalent_filter
        self.shelf['is_filter_equivalent_cavities'] = self.is_filter_equivalent_cavities
        self.shelf['nb_stretch'] = self.nb_stretch
        self.shelf['nb_gratings'] = self.nb_gratings
        self.shelf['SRG_dict'] = self.SRG_dict

        print(list(self.shelf.keys()))

        self.shelf.close()

    def load_state(self, load_file_path):
        self.shelf = shelve.open(load_file_path)

        print (list(self.shelf.keys()))

        self.exp_data = self.shelf['exp_data']
        self.name = self.shelf['name']
        self.min_pt_search_resolution = self.shelf['min_pt_search_resolution']
        self.sizeX, self.sizeY = self.shelf['sizeX'], self.shelf['sizeY']
        self.xOffset, self.yOffset = self.shelf['xOffset'], self.shelf['yOffset']
        self.max_relative_diff_for_equivalent_filter = self.shelf['max_relative_diff_for_equivalent_filter']
        self.is_filter_equivalent_cavities = self.shelf['is_filter_equivalent_cavities']
        self.nb_stretch = self.shelf['nb_stretch']
        self.nb_gratings = self.shelf['nb_gratings']
        self.SRG_dict = self.shelf['SRG_dict']

        self.shelf.close()

    def drawContourFromStretchSerie(self, dots_critical=False, dots_cavity=False):
        for i_s in range(self.nb_stretch):
            qc = self.SRG_dict[i_s]
            title = "Stretch = " + str(i_s) +  "mm"
            fname = self.name + "_contour_s" + str(i_s) + ".png"
            qc.drawContour(75, fname, title, dots_critical, dots_cavity)

    def drawCavitiesFromStretchSerie(self):
        for i_s in range(self.nb_stretch):
            qc = self.SRG_dict[i_s]
            title = "Stretch = " + str(i_s) +  "mm"
            fname = self.name + "_cavities_s" + str(i_s) + ".png"
            qc.drawAllQuasiCrystal(dots=False, file_save_Path=fname, ellipses=True, title=title)

    def draw_3D_FromStretchSerie(self, nbStetch, amp, angle, pas):

        for i_s in range(nbStetch):
            qc = self.SRG_dict[i_s]
            title = "Stretch = " + str(i_s) +  "mm"
            fname = self.name + str(i_s) + ".png"
            qc.drawAllQuasiCrystal(dots=False, coat='full', coatalpha=0.5, ellipses=False, wellidx=[], ellipsesalpha=0.5, file_save_Path=fname)

    def histogram_3D_OfCavityEllipticy(self):
        for i_s in range(self.nb_stretch):
            listVolume = []
            listEllipticity = []
            qc = self.SRG_dict[i_s]
            for cavity in qc.list_nano_cavities:
                ellipsoid = cavity.ellipsoid
                listVolume.append(ellipsoid.volume)
                listEllipticity.append(ellipsoid.ellipticity)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            hist, xedges, yedges = np.histogram2d(listVolume, listEllipticity, bins=(4, 4))
            xpos, ypos = np.meshgrid(xedges[:-1] + xedges[1:], yedges[:-1] + yedges[1:])

            xpos = xpos.flatten() / 2.
            ypos = ypos.flatten() / 2.
            zpos = np.zeros_like(xpos)

            dx = xedges[1] - xedges[0]
            dy = yedges[1] - yedges[0]
            dz = hist.flatten()

        # For the6 faces
        # When coloring the faces of the boxes specifically, this is the order of the coloring:
        #
        # -Z (bottom of box)
        # +Z (top of box)
        # -Y
        # +Y
        # -X
        # +X


            offset = dz + np.abs(dz.min())
            fracs = offset.astype(float) / offset.max()
            # norm = colors.Normalize(fracs.min(), fracs.max())
            # colors = cm.jet(norm(fracs))

            # ax.bar3d(xpos,ypos,zpos,1,1,dz, color=colors)
            col = ['b', 'r', 'b', 'b', 'b', 'b']
            col = 'r'
            # ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, zsort='average')
            ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
            plt.xlabel("Volume")
            plt.ylabel("Ellipticity")
            ax.set_zlabel("Occurence")
            plt.show()

    def get_evolution_of_cavity_statistics(self):
        #FIXME
        stretch_value = [0, 1, 2, 3, 4, 5, 6]
        list_volume_mean, list_volume_std  = [], []
        list_ellipticity_mean, list_ellipticity_std = [], []
        for i_s in range(self.nb_stretch):
            qc = self.SRG_dict[i_s]
            e_mean, e_std, v_mean, v_std, oX_mean, oX_std, oY_mean, oY_std, oZ_mean, oZ_std = qc.get_cavities_statistics()
            list_volume_mean.append(v_mean)
            list_volume_std.append(v_std)
            list_ellipticity_mean.append(e_mean)
            list_ellipticity_std.append(e_std)

        plt.errorbar(stretch_value, list_volume_mean, yerr=list_volume_std)
        plt.title("Evolution Volume")
        name = self.name + "evol_volume" + ".png"
        plt.savefig(name)
        plt.show()

        plt.errorbar(stretch_value, list_ellipticity_mean, yerr=list_ellipticity_std)
        plt.title("Evolution Ellipticity")
        name = self.name + "evol_ell" + ".png"
        plt.savefig(name)
        plt.show()



    def statisticsOnAlignementWithStretchConstraint(self, amp, pas) :
        arrayStretch = np.arange(0, nbStetch, 1)

        listMeanOrientationX = []
        listMeanOrientationY = []
        listMeanOrientationZ = []

        listStdDevOrientationX = []
        listStdDevOrientationY = []
        listStdDevOrientationZ = []

        vecX = [1, 0, 0]
        vecY = [0, 1, 0]
        vecZ = [0, 0, 1]

        for i_s in range(nbStetch):
            print("Stretch nb :", i_s, "over :", nbStetch)
            #create Sample
            listGratings = []
            phase = PI / nb_gratings
            for j_g in range(nb_gratings):
                listGratings.append(sinusoidalGrating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
            qc = quasiCrystal(listGratings, res, dimX, dimY, offset, offset, tol, filtering=True)

            listOrientationX = []
            listOrientationY = []
            listOrientationZ = []

            for cavity in qc.list_nano_cavities:
                ellipsoid = cavity.ellipsoid
                listOrientationX.append(ellipsoid.assessOrientationWithStretchConstraint(vecX))
                listOrientationY.append(ellipsoid.assessOrientationWithStretchConstraint(vecY))
                listOrientationZ.append(ellipsoid.assessOrientationWithStretchConstraint(vecZ))

            listMeanOrientationX.append(np.mean(np.asarray(listOrientationX)))
            listMeanOrientationY.append(np.mean(np.asarray(listOrientationY)))
            listMeanOrientationZ.append(np.mean(np.asarray(listOrientationZ)))

            listStdDevOrientationX.append(np.std(np.asarray(listOrientationX)))
            listStdDevOrientationY.append(np.std(np.asarray(listOrientationY)))
            listStdDevOrientationZ.append(np.std(np.asarray(listOrientationZ)))

        # print(arrayStretch)
        # print(listMeanOrientationX)
        plt.errorbar(arrayStretch, listMeanOrientationX, yerr=listStdDevOrientationX, marker='o')
        plt.show()
        plt.errorbar(arrayStretch, listMeanOrientationY, yerr=listStdDevOrientationY, marker='o')
        plt.show()
        plt.errorbar(arrayStretch, listMeanOrientationZ, yerr=listStdDevOrientationZ, marker='o')
        plt.show()



    def getCenterCavityWithEllips(self, amp, pas, angle):
        dim = 5
        offsetXY = -2.5
        dimX, dimY, offset, res, tol = dim, dim, offsetXY, 0.1, 0.1

        #firstGrating
        n = 0
        listGratings = []
        nb_gratings = 3
        for j_g in range(nb_gratings):
            listGratings.append(sinusoidalGrating(amp[n, j_g], angle[n, j_g], pas[n, j_g]))
        qc = quasiCrystal(listGratings, res, dimX, dimY, offset, offset, tol, is_filter_equivalent_cavities=False)

        for cavity in qc.list_nano_cavities:
            if (-0.2 < cavity.xc < 0.2) and (-0.2 < cavity.yc < 0.2):
                cavity.draw(dots=False, coat='full', coatAlpha=0.5, drawEllipsoid=False, ellipseAlpha=0.2)

#getCenterCavityWithEllips(amp3G, p3G, angle3G)


    def evolutionOverStretching(self, amp, pas, angle, nbStetch):
        arrayStretch = np.arange(0, nbStetch, 1)

        dim = 5
        offsetXY = -2.5
        dimX, dimY, offset, res, tol = dim, dim, offsetXY, 0.1, 0.1

        listMeanEllicticity = []
        listMeanVolume = []
        listMeanEccentricity = []

        listStdDevEllicticity = []
        listStdDevVolume = []
        listStdDevEccentricity = []

        for i_s in range(nbStetch):
            print("Stretch nb :", i_s, "over :", nbStetch)
            # create Sample
            listGratings = []
            for j_g in range(nb_gratings):
                # listGratings.append(sinusoidalGrating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
                listGratings.append(sinusoidalGrating(amp[i_s, j_g], angle3G[i_s, j_g], pas[i_s, j_g]))

            qc = quasiCrystal(listGratings, res, dimX, dimY, offset, offset, tol, is_filter_equivalent_cavities=True)

            listEllicticity = []
            listVolume = []
            listEccentricity = []

            for cavity in qc.list_nano_cavities:
                cavity.draw(dots=True, coat='full', coatAlpha=0.5, drawEllipsoid=True, ellipseAlpha=0.1)
                ellipsoid = cavity.ellipsoid
                listEllicticity.append(ellipsoid.ellipticity)
                listVolume.append(ellipsoid.volume)
                listEccentricity.append(ellipsoid.eccentricity)

            listMeanEllicticity.append(np.mean(np.asarray(listEllicticity)))
            listMeanVolume.append(np.mean(np.asarray(listVolume)))
            listMeanEccentricity.append(np.mean(np.asarray(listEccentricity)))

            listStdDevEllicticity.append(np.std(np.asarray(listEllicticity)))
            listStdDevVolume.append(np.std(np.asarray(listVolume)))
            listStdDevEccentricity.append(np.std(np.asarray(listEccentricity)))


        # print(arrayStretch)
        # print(listMeanOrientationX)
        plt.errorbar(arrayStretch, listMeanEllicticity, yerr=listStdDevEllicticity, marker='o')
        plt.show()
        plt.errorbar(arrayStretch, listMeanVolume, yerr=listStdDevVolume, marker='o')
        plt.show()
        plt.errorbar(arrayStretch, listMeanEccentricity, yerr=listStdDevEccentricity, marker='o')
        plt.show()

#evolutionOverStretching(amp3G, p3G, angle3G, 6)




    def evolutionOverStretchingOneCavity(self, amp, pas, angle, nbStetch):


        dim = 2
        offsetXY = -1
        dimX, dimY, offset, res, tol = dim, dim, offsetXY, 0.1, 0.1

        listEllicticity = []
        listVolume = []
        listEccentricity = []

        for i_s in range(nbStetch):
            print("Stretch nb :", i_s, "over :", nbStetch)
            # create Sample
            listGratings = []
            for j_g in range(nb_gratings):
                # listGratings.append(sinusoidalGrating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
                listGratings.append(sinusoidalGrating(amp[i_s, j_g], angle3G[i_s, j_g], pas[i_s, j_g]))

            qc = quasiCrystal(listGratings, res, dimX, dimY, offset, offset, tol, is_filter_equivalent_cavities=True)
            cs = qc.drawContour(75)

            #Extract level curve
            #print(cs.collections[0])
            # print(cs.collections[0].get_paths())
            #p = cs.collections[0].get_paths()[0]
            maxLenPath = 0
            for path in cs.collections[0].get_paths():
                extent = path.get_extents()

                if extent.contains(0, 0):
                    p = path

            v = p.vertices
            x = v[:, 0]
            y = v[:, 1]

            plt.plot(x,y)
            plt.show()

            def fitEllipse(x, y):
                x = x[:, np.newaxis]
                y = y[:, np.newaxis]
                D = np.hstack((x * x, x * y, y * y, x, y, np.ones_like(x)))
                S = np.dot(D.T, D)
                C = np.zeros([6, 6])
                C[0, 2] = C[2, 0] = 2;
                C[1, 1] = -1
                E, V = eig(np.dot(inv(S), C))
                n = np.argmax(np.abs(E))
                a = V[:, n]
                return a

            def ellipse_center(a):
                b, c, d, f, g, a = a[1] / 2, a[2], a[3] / 2, a[4] / 2, a[5], a[0]
                num = b * b - a * c
                x0 = (c * d - b * f) / num
                y0 = (a * f - b * d) / num
                return np.array([x0, y0])

            def ellipse_angle_of_rotation(a):
                b, c, d, f, g, a = a[1] / 2, a[2], a[3] / 2, a[4] / 2, a[5], a[0]
                return 0.5 * np.arctan(2 * b / (a - c))

            def ellipse_axis_length(a):
                b, c, d, f, g, a = a[1] / 2, a[2], a[3] / 2, a[4] / 2, a[5], a[0]
                up = 2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g)
                down1 = (b * b - a * c) * ((c - a) * np.sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a))
                down2 = (b * b - a * c) * ((a - c) * np.sqrt(1 + 4 * b * b / ((a - c) * (a - c))) - (c + a))
                res1 = np.sqrt(up / down1)
                res2 = np.sqrt(up / down2)
                return np.array([res1, res2])

            a = fitEllipse(x, y)

            axisLength = ellipse_axis_length(a)
            print(axisLength)
            longAxe = max(axisLength[0], axisLength[1])
            minAxe = min(axisLength[0], axisLength[1])
            #print(longAxe, minAxe)
            listEllicticity.append(1 - minAxe/longAxe)


        arrayStretch = np.arange(0, len(listEllicticity), 1)

        plt.plot(arrayStretch, listEllicticity, marker='o')
        plt.show()
        print(listEllicticity)

        # qc.drawAllQuasiCrystal()

        #     for cavity in qc.listNanoCavities:
        #         if (-0.2 < cavity.xc < 0.2) and (-0.2 < cavity.yc < 0.2):
        #             cavity.draw(dots=True, coat='full', coatAlpha=0.5, drawEllipsoid=True, ellipseAlpha=0.1)
        #             ellipsoid = cavity.ellipsoid
        #             listEllicticity.append(ellipsoid.ellipticity)
        #             listVolume.append(ellipsoid.volume)
        #             listEccentricity.append(ellipsoid.eccentricity)
        #
        # # print(arrayStretch)
        # # print(listMeanOrientationX)
        # arrayStretch = np.arange(0, len(listEllicticity), 1)
        # plt.plot(arrayStretch, listEllicticity, marker='o')
        # plt.show()
        # plt.plot(arrayStretch, listVolume, marker='o')
        # plt.show()
        # plt.plot(arrayStretch, listEccentricity, marker='o')
        # plt.show()
