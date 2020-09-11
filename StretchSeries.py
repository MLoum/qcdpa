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
            qc.drawAllQuasiCrystal(dots=False, file_save_Path=fname, ellipsoid=True, title=title)

    def draw_3D_FromStretchSerie(self, nbStetch, amp, angle, pas):

        for i_s in range(nbStetch):
            qc = self.SRG_dict[i_s]
            title = "Stretch = " + str(i_s) +  "mm"
            fname = self.name + str(i_s) + ".png"
            qc.drawAllQuasiCrystal(dots=False, coat='full', coatalpha=0.5, ellipsoid=False, wellidx=[], ellipsesalpha=0.5, file_save_Path=fname)

    def histogram_3D_volume_eccentricity(self):
        for i_s in range(self.nb_stretch):
            listVolume = []
            list_eccentricity = []
            qc = self.SRG_dict[i_s]
            for cavity in qc.list_nano_cavities:
                ellipsoid = cavity.ellipsoid
                listVolume.append(ellipsoid.volume)
                list_eccentricity.append(ellipsoid.eccentricity)

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            hist, xedges, yedges = np.histogram2d(listVolume, list_eccentricity, bins=(4, 4))
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
            plt.ylabel("Eccentricity")
            ax.set_zlabel("Occurence")
            plt.title("Stretch = " + str(i_s) + "mm")
            fname = self.name + "_histo3D_stretch_" + str(i_s) + ".png"
            plt.savefig(fname, dpi=300)
            plt.show()

    def get_evolution_of_cavity_statistics(self):
        #FIXME
        stretch_value = [0, 1, 2, 3, 4, 5, 6]
        list_volume_mean, list_volume_std = [], []
        list_eccentricity_mean, list_eccentricity_std = [], []
        for i_s in range(self.nb_stretch):
            qc = self.SRG_dict[i_s]
            e_mean, e_std, v_mean, v_std, oX_mean, oX_std, oY_mean, oY_std, oZ_mean, oZ_std = qc.get_cavities_statistics()
            list_volume_mean.append(v_mean)
            list_volume_std.append(v_std)
            list_eccentricity_mean.append(e_mean)
            list_eccentricity_std.append(e_std)

        plt.errorbar(stretch_value, list_volume_mean, yerr=list_volume_std, fmt="ro")
        plt.title("Evolution of volume with stretch")
        plt.xlabel("stretch / mm")
        plt.ylabel("Volume / $µm^3$ x10-3")
        name = self.name + "evol_volume" + ".png"
        plt.savefig(name, dpi=300)
        plt.show()

        plt.errorbar(stretch_value, list_eccentricity_mean, yerr=list_eccentricity_std, fmt="ro")
        plt.title("Evolution of Eccentricity with stretch")
        plt.xlabel("stretch / mm")
        plt.ylabel("Eccentricity")
        name = self.name + "evol_ell" + ".png"
        plt.savefig(name, dpi=300)
        plt.show()


    def statistics_on_alignement_with_stretch_constraint(self) :
        stretch_value = [0, 1, 2, 3, 4, 5, 6]

        listMeanOrientationX = []
        listMeanOrientationY = []
        listMeanOrientationZ = []

        listStdDevOrientationX = []
        listStdDevOrientationY = []
        listStdDevOrientationZ = []


        for i_s in range(self.nb_stretch):
            # print("Stretch nb :", i_s, "over :", self.nb_stretch)
            qc = self.SRG_dict[i_s]

            listOrientationX = []
            listOrientationY = []
            listOrientationZ = []

            for cavity in qc.list_nano_cavities:
                ellipsoid = cavity.ellipsoid
                o_X, o_Y, o_Z =  ellipsoid.assessOrientationWithStretchConstraint()
                listOrientationX.append(o_X)
                listOrientationY.append(o_Y)
                listOrientationZ.append(o_Z)

            listMeanOrientationX.append(np.mean(np.asarray(listOrientationX)))
            listMeanOrientationY.append(np.mean(np.asarray(listOrientationY)))
            listMeanOrientationZ.append(np.mean(np.asarray(listOrientationZ)))

            listStdDevOrientationX.append(np.std(np.asarray(listOrientationX)))
            listStdDevOrientationY.append(np.std(np.asarray(listOrientationY)))
            listStdDevOrientationZ.append(np.std(np.asarray(listOrientationZ)))

        plt.errorbar(stretch_value, listMeanOrientationX, yerr=listStdDevOrientationX, fmt='ro')
        plt.xlabel("stretch / mm")
        plt.ylabel("Orientation X")
        name = self.name + "orienation_X" + ".png"
        plt.savefig(name, dpi=300)
        plt.show()
        plt.errorbar(stretch_value, listMeanOrientationY, yerr=listStdDevOrientationY, fmt='ro')
        plt.xlabel("stretch / mm")
        plt.ylabel("Orientation Y")
        name = self.name + "orientation_Y" + ".png"
        plt.savefig(name, dpi=300)
        plt.show()
        plt.errorbar(stretch_value, listMeanOrientationZ, yerr=listStdDevOrientationZ, fmt='ro')
        plt.xlabel("stretch / mm")
        plt.ylabel("Orientation Z")
        name = self.name + "orientation_Z" + ".png"
        plt.savefig(name, dpi=300)
        plt.show()



