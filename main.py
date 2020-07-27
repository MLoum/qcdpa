import warnings
import numpy as np
from sinusoidalGrating import sinusoidalGrating
from QuasiCrystal import quasiCrystal
from  ExperimentalData import exp_data_2G, exp_data_3G, exp_data_4G, exp_data_8G
from StretchSeries import StretchSeries



#TODO tester les échantillons avec une très haute symétrie pour voir si le nombre de types de cavités sature

#TODO filger par rapport aux bords que pour les minimums

#TODO chercher la maille elementaire.

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
PI = np.pi

stretch = 1
dim = 10
offsetXY = -dim/2.0
resolution = 0.1
dimX, dimY, offset, res, tol, bars = dim, dim, offsetXY, resolution, 0.1, 30

print("Ici")
stretch_dists_mm = [0, 1, 2, 3, 4, 5, 6]
# _2G_serie = StretchSeries(stretch_dists_mm, exp_data_2G, "2G", min_pt_search_resolution=resolution, sizeX=dim, sizeY=dim, xOffset=offsetXY, yOffset=offsetXY, max_relative_diff_for_equivalent_filter=0.1, is_filter_equivalent_cavities=True)
# _2G_serie.drawContourFromStretchSerie()

# _4G_serie = StretchSeries(stretch_dists_mm, exp_data_4G, "4G", min_pt_search_resolution=resolution, sizeX=dim, sizeY=dim, xOffset=offsetXY, yOffset=offsetXY, max_relative_diff_for_equivalent_filter=0.1, is_filter_equivalent_cavities=True)
# _4G_serie.drawContourFromStretchSerie(dots_cavity=True)
# _4G_serie.save_state("4G_serie_10microns_filter0p1")
# _4G_serie.drawCavitiesFromStretchSerie()
# _4G_serie.get_evolution_of_cavity_statistics()
# _4G_serie.histogram_3D_OfCavityEllipticy()


amps = exp_data_4G.amps[0]
angles = exp_data_4G.angles[0]
pitches = exp_data_4G.pitches[0]
list_gratings = []
for j_g in range(exp_data_4G.nb_gratings):
    list_gratings.append(sinusoidalGrating(amps[j_g], angles[j_g], pitches[j_g]))

qc_4G = quasiCrystal(list_gratings, min_pt_search_resolution=resolution, sizeX=dim, sizeY=dim, xOffset=offsetXY, yOffset=offsetXY, max_relative_diff_for_equivalent_filter=0.2, is_filter_equivalent_cavities=True)

qc_4G.drawContour(dots_cavity=True)
qc_4G.draw_cavities()


# _2G_serie = StretchSeries(stretch_dists_mm, exp_data_2G, "2G", min_pt_search_resolution=resolution, sizeX=dim, sizeY=dim, xOffset=offsetXY, yOffset=offsetXY, max_relative_diff_for_equivalent_filter=0.2, is_filter_equivalent_cavities=True)
# _2G_serie.drawContourFromStretchSerie(dots_cavity=True)
# _2G_serie.save_state("2G_serie_10microns_filter0p1")
# _2G_serie.drawCavitiesFromStretchSerie()
# # _2G_serie.get_evolution_of_cavity_statistics()
# # _2G_serie.histogram_3D_OfCavityEllipticy()



# _4G_serie = StretchSeries(exp_data="G_serie_10microns_filter0p1.dat")
# # _4G_serie.drawContourFromStretchSerie()
# _4G_serie.histogram_3D_OfCavityEllipticy()


# evolutionOverStretchingOneCavity(amp3G, p3G, angle3G, 6)

# nbofGrating = 3
# sample = createPerfectSRG(nbofGrating, depth=1, groove=0.5, negative=True)
# fname = "perfect_SRG_large300_" + str(nbofGrating) + ".png"
# sample.drawContour(75,fname)
#
# for well in sample.listNanoCavities :
#     well.printInfo()
#     well.ellipsoid.printData()
#     well.draw(dots=True, coat='full', coatAlpha=0.2, ellipse=False, ellipseAlpha=0.5)
#     well.draw(dots=True, coat='wire', coatAlpha=0.2, ellipse=True, ellipseAlpha=0.5)



# list_gratings = ([grating0, grating1, grating2, grating3, grating4, grating5, grating6, grating7])
#
# PDMS_sample = sample(list_gratings, stretch, res, dimX, dimY, offset, offset, tol, filtering=True)



# PDMS_sample = quasiCrystal(listGratings, res, dimX, dimY, offset, offset, tol, isFilterEquivalentCavities=True)
# #Test MVEE et dessin
# for well in PDMS_sample.listNanoCavities :
#     well.printInfo()
#     well.ellipsoid.printData()
#     well.draw(dots=True, coat='full', coatAlpha=0.2, drawEllipsoid=False, ellipseAlpha=0.5)
#     well.draw(dots=True, coat='wire', coatAlpha=0.2, drawEllipsoid=True, ellipseAlpha=0.5)

# listVolume = []
# for well in PDMS_sample.listNanoCavities :
#     ellipsoid = well.ellipsoid
#     listVolume.append(ellipsoid.volume)
#
# nbOfBar = 10
# histVol = histogram(listVolume, nbOfBar, fit='None')
# histVol.plotHistogram('hist Volume', False)






#Alpha variable
# xpos = np.arange(0,4,1)
# ypos = np.arange(0,4,1)
# xpos, ypos = np.meshgrid(xpos, ypos)
# xpos = xpos.flatten()
# ypos = ypos.flatten()
# zpos = np.zeros(4*4)
# rho = np.random.random((4,4))
# dx = 0.5 * np.ones_like(zpos)
# dy = dx.copy()
# dz = rho.flatten()
# nrm=mpl.colors.Normalize(-1,1)
# colors=cm.RdBu(nrm(-dz))
# alpha = np.linspace(0.2, 0.95, len(xpos), endpoint=True)
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# for i in range(len(xpos)):
#     ax.bar3d(xpos[i],ypos[i],zpos[i], dx[i], dy[i], dz[i], alpha=alpha[i], color=colors[i], linewidth=0)
# plt.show()




#Test stability of number of Cavities kind
# nbTest = 30
# inc = 2
#
# dimX, dimY, res, tol, bars = 5.0, 5.0, 0.1, 0.1, 30
# nbPuit = []
# dim = []
#
#
# for i in range(nbTest):
#     print ("test : ", i, "over ", nbTest)
#     dimX = dimY = 5 + inc*i
#     dim.append(dimX)
#     print("dim : ", dimX)
#     PDMS_sample = sample(list_gratings, stretch, res, dimX, dimY, tol, filtering=True)
#     nbPuit.append(len(PDMS_sample.listNanoCavities ))
#
# plt.plot(dim, nbPuit)
# plt.show()



