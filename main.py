"""               Modelisation des réseaux et affichage des données sur les puits et particules               """

import warnings
import numpy as np
from Grating import grating
from QuasiCrystal import quasiCrystal
from Histogram import histogram
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt





import os


#TODO tester les échantillons avec une très haute symétrie pour voir si le nombre de types de cavités sature

#TODO filger par rapport aux bords que pour les minimums

#TODO chercher la maille elementaire.

warnings.filterwarnings('ignore', 'The iteration is not making good progress')
PI = np.pi

stretch = 1
dim = 10
offsetXY = -5
dimX, dimY, offset, res, tol, bars = dim, dim, offsetXY, 0.1, 0.1, 30



def createPerfectSRG(nbOfGrating, depth=1, groove=1, negative=False):
    if negative:
        depth = - depth

    s = np.ones(nbOfGrating) * depth
    p = np.ones(nbOfGrating) * groove
    phase = PI / nbOfGrating
    list_gratings = []
    for i in range(nbOfGrating):
        g = grating(s[i], i * phase, p[i])
        list_gratings.append(g)

    return quasiCrystal(list_gratings, stretch, res, dimX, dimY, offset, offset, tol, filtering=True)


# -------------------------------------------------- 8 Gratings ---------------------------------------------------


print('Réseau 8G :')

nb_gratings = 8
nbStetch = 7
phase = PI / nb_gratings

amp = np.array([[0.427947, 0.350812, 0.505739, 0.305864, 0.679745, 0.578136, 0.554883, 1],
                [0.414172, 0.460252, 0.513327, 0.495240, 0.771715, 0.511313, 0.618458, 1],
                [0.513523, 0.434903, 0.828801, 0.390097, 0.892500, 0.487479, 0.690805, 1],
                [0.549956, 0.425569, 0.509538, 0.435393, 1, 0.419731, 0.822757, 0.901456],
                [0.568348, 0.457213, 0.394343, 0.376277, 1, 0.327728, 0.915975, 0.718308],
                [0.696373, 0.607797, 0.667442, 0.549605, 1, 0.567812, 0.893165, 0.980499],
                [0.596690, 0.771279, 0.624219, 0.767799, 1.086060, 0.555316, 1, 0.950397]])

pas = np.array([[0.957, 0.961, 0.961, 0.961, 0.966, 0.977, 0.972, 0.967],
                [0.995, 0.999, 0.978, 0.954, 0.940, 0.948, 0.964, 0.985],
                [1.039, 1.042, 0.995, 0.945, 0.913, 0.921, 0.953, 1.004],
                [1.090, 1.091, 1.010, 0.936, 0.890, 0.898, 0.942, 1.025],
                [1.142, 1.141, 1.028, 0.926, 0.867, 0.873, 0.930, 1.039],
                [1.192, 1.188, 1.037, 0.916, 0.847, 0.852, 0.918, 1.054],
                [1.234, 1.226, 1.039, 0.903, 0.827, 0.833, 0.905, 1.064]])


# # # #-------------------- 4 gratings
#
# print('Réseau 4G :')
#
# nb_gratings = 4
# nbStetch = 7
# phase = PI / nb_gratings
#
#
# amp = np.array([[0.786441, 0.885087, 0.708722, 1],  # 0mm
#                 [0.707390, 0.903390, 0.716984, 1],  # 1mm
#                 [0.782744, 0.935068, 0.748444, 1],  # 2mm
#                 [0.828168, 0.952814, 0.747646, 1],  # 3mm
#                 [0.839347, 0.943083, 0.752176, 1],  # 4mm
#                 [0.836209, 0.898063, 0.742515, 1],  # 5mm
#                 [0.674446, 0.930514, 0.821491, 1]])  # 6mm
#
# pas = np.array([[1.568, 1.702, 1.841, 1.701],  # 0mm
#                 [1.622, 1.727, 1.811, 1.712],  # 1mm
#                 [1.777, 1.775, 1.728, 1.744],  # 2mm
#                 [1.921, 1.803, 1.659, 1.763],  # 3mm
#                 [2.052, 1.821, 1.599, 1.767],  # 4mm
#                 [2.177, 1.832, 1.544, 1.765],  # 5mm
#                 [2.189, 1.812, 1.532, 1.770]])  # 6m



# # #-------------------- 2 gratings FAKE
# print('Réseau 4G :')
#
# nb_gratings = 2
# nbStetch = 7
# phase = PI / nb_gratings
#
#
# amp = np.array([[0.95, 1],  # 0mm
#                 [0.97, 1],  # 1mm
#                 [0.95, 1],  # 2mm
#                 [0.94, 1],  # 3mm
#                 [0.1, 0.97],  # 4mm
#                 [0.98, 1],  # 5mm
#                 [0.95, 1]])  # 6mm
#
# pas = np.array([[1, 1],  # 0mm
#                 [0.9, 1.1],  # 1mm
#                 [0.8, 1.2],  # 2mm
#                 [0.7, 1.3],  # 3mm
#                 [0.6, 1.4],  # 4mm
#                 [0.5, 1.5],  # 5mm
#                 [0.4, 1.6]])  # 6mm
#
#
#
# # #-------------------- 3 gratings FAKE
# print('Réseau 3G  Fake:')
#
# nb_gratings = 3
# nbStetch = 7
# phase = PI / nb_gratings
#
#
# amp = np.array([[0.95, 1, 0.99],  # 0mm
#                 [0.97, 1, 0.99],  # 1mm
#                 [0.95, 0.99, 1],  # 2mm
#                 [0.99, 0.94, 1],  # 3mm
#                 [0.1, 0.99, 0.97],  # 4mm
#                 [0.98, 1, 0.99],  # 5mm
#                 [0.95, 0.99, 1]])  # 6mm
#
# amp = - amp
#
# pas = np.array([[0.99, 1, 0.99],
#                 [1, 1.05, 0.95],   # 1mm
#                 [1.1, 1.1, 0.9],  # 2mm
#                 [1.2, 1.15, 0.85],  # 3mm
#                 [1.3, 1.2, 0.8],  # 4mm
#                 [1.4, 1.25, 0.75],  # 5mm
#                 [1.5, 1.3, 0.7]])  # 6mm


# plt.plot(pas[:,0], 'k')
# plt.plot(pas[:,1], 'b')
# plt.plot(pas[:,2], 'r')
# plt.plot(pas[:,3], 'k')
# plt.plot(pas[:,4], 'b')
# plt.plot(pas[:,5], 'r')
# plt.plot(pas[:,6], 'y')
# plt.plot(pas[:,7], 'c')

# plt.plot(amp[:,0], 'k')
# plt.plot(amp[:,4], 'b')
# plt.plot(amp[:,7], 'r')

# plt.show()



# for i_s in range(nbStetch):
#     print("Stretch nb :", i_s, "over :", nbStetch)
#     #create Sample
#     listGratings = []
#     phase = PI / nb_gratings
#     for j_g in range(nb_gratings):
#         listGratings.append(grating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
#     PDMS_sample = sample(listGratings, i_s, res, dimX, dimY, offset, offset, tol, filtering=True)
#     PDMS_sample.drawContour()
#     fname = "level_s" + str(i_s) + ".jpg"
#     plt.savefig(fname)


# nbofGrating = 4
# sample = createPerfectSRG(nbofGrating)
# fname = "perfect_SRG_jet_" + str(nbofGrating) + ".png"
# sample.drawContour(75,fname)




# for i_s in range(nbStetch):
#     listGratings = []
#     phase = PI / nb_gratings
#     for j_g in range(nb_gratings):
#         listGratings.append(grating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
#     PDMS_sample = sample(listGratings, i_s, res, dimX, dimY, offset, offset, tol, filtering=True)
#     fname = "3_grating_fake_level_30_jet_s" + str(i_s) + ".png"
#     PDMS_sample.drawContour(75, fname)


# arrayStretch = np.arange(0, nbStetch, 1)
#
# listMeanOrientationX = []
# listMeanOrientationY = []
# listMeanOrientationZ = []
#
# listStdDevOrientationX = []
# listStdDevOrientationY = []
# listStdDevOrientationZ = []
#
# vecX = [1, 0, 0]
# vecY = [0, 1, 0]
# vecZ = [0, 0, 1]
#
# for i_s in range(nbStetch):
#     print("Stretch nb :", i_s, "over :", nbStetch)
#     #create Sample
#     listGratings = []
#     phase = PI / nb_gratings
#     for j_g in range(nb_gratings):
#         listGratings.append(grating(amp[i_s, j_g], j_g * phase, pas[i_s, j_g]))
#     PDMS_sample = sample(listGratings, i_s, res, dimX, dimY, offset, offset, tol, filtering=True)
#
#     listOrientationX = []
#     listOrientationY = []
#     listOrientationZ = []
#
#     for well in PDMS_sample.listWell:
#         ellipsoid = well.ellipsoid
#         listOrientationX.append(ellipsoid.testOrientation(vecX))
#         listOrientationY.append(ellipsoid.testOrientation(vecY))
#         listOrientationZ.append(ellipsoid.testOrientation(vecZ))
#
#     listMeanOrientationX.append(np.mean(np.asarray(listOrientationX)))
#     listMeanOrientationY.append(np.mean(np.asarray(listOrientationY)))
#     listMeanOrientationZ.append(np.mean(np.asarray(listOrientationZ)))
#
#     listStdDevOrientationX.append(np.std(np.asarray(listOrientationX)))
#     listStdDevOrientationY.append(np.std(np.asarray(listOrientationY)))
#     listStdDevOrientationZ.append(np.std(np.asarray(listOrientationZ)))
#
# print(arrayStretch)
# print(listMeanOrientationX)
# plt.errorbar(arrayStretch, listMeanOrientationX, yerr=listStdDevOrientationX, marker='o')
# plt.show()
# plt.errorbar(arrayStretch, listMeanOrientationY, yerr=listStdDevOrientationY, marker='o')
# plt.show()
# plt.errorbar(arrayStretch, listMeanOrientationZ, yerr=listStdDevOrientationZ, marker='o')
# plt.show()


# nbofGrating = 3
# sample = createPerfectSRG(nbofGrating, depth=1, groove=0.5, negative=True)
# fname = "perfect_SRG_large300_" + str(nbofGrating) + ".png"
# sample.drawContour(75,fname)
#
# for well in sample.listWell :
#     well.printInfo()
#     well.ellipsoid.printData()
#     well.drawWell(dots=True, coat='full', coatAlpha=0.2, ellipse=False, ellipseAlpha=0.5)
#     well.drawWell(dots=True, coat='wire', coatAlpha=0.2, ellipse=True, ellipseAlpha=0.5)



# list_gratings = ([grating0, grating1, grating2, grating3, grating4, grating5, grating6, grating7])
#
# PDMS_sample = sample(list_gratings, stretch, res, dimX, dimY, offset, offset, tol, filtering=True)



PDMS_sample = quasiCrystal(listGratings, i_s, res, dimX, dimY, offset, offset, tol, filtering=True)
#Test MVEE et dessin
for well in PDMS_sample.listWell :
    well.printInfo()
    well.ellipsoid.printData()
    well.drawWell(dots=True, coat='full', coatAlpha=0.2, ellipse=False, ellipseAlpha=0.5)
    well.drawWell(dots=True, coat='wire', coatAlpha=0.2, ellipse=True, ellipseAlpha=0.5)

# listVolume = []
# for well in PDMS_sample.listWell :
#     ellipsoid = well.ellipsoid
#     listVolume.append(ellipsoid.volume)
#
# nbOfBar = 10
# histVol = histogram(listVolume, nbOfBar, fit='None')
# histVol.plotHistogram('hist Volume', False)



# listVolume = []
# listEllipticity = []
# for well in PDMS_sample.listWell :
#     ellipsoid = well.ellipsoid
#     listVolume.append(ellipsoid.volume)
#     listEllipticity.append(ellipsoid.ellipticity)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
#
# hist, xedges, yedges = np.histogram2d(listVolume, listEllipticity, bins=(4,4))
# xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:])
#
# xpos = xpos.flatten()/2.
# ypos = ypos.flatten()/2.
# zpos = np.zeros_like (xpos)
#
# dx = xedges [1] - xedges [0]
# dy = yedges [1] - yedges [0]
# dz = hist.flatten()
#
# #pour les 6 faces
# # When coloring the faces of the boxes specifically, this is the order of the coloring:
# #
# # -Z (bottom of box)
# # +Z (top of box)
# # -Y
# # +Y
# # -X
# # +X
#
#
# offset = dz + np.abs(dz.min())
# fracs = offset.astype(float)/offset.max()
# norm = colors.Normalize(fracs.min(), fracs.max())
# colors = cm.jet(norm(fracs))
#
# # ax.bar3d(xpos,ypos,zpos,1,1,dz, color=colors)
#
#
# col = ['b', 'r', 'b', 'b', 'b', 'b']
# col = 'r'
# ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors, zsort='average')
# plt.xlabel ("Volume")
# plt.ylabel ("Ellipticity")
# ax.set_zlabel("Occurence")
#
# plt.show()
#



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
#     nbPuit.append(len(PDMS_sample.listWell ))
#
# plt.plot(dim, nbPuit)
# plt.show()





# ---- Sample Analyses :
#PDMS_sample.getWellsStats(basicinfo=True, MVEE=True)

# ---- Draw All Sample :
#PDMS_sample.drawAllSample(dots=False, coat='full', coatalpha=0.5, ellipses=True, wellidx=range(len(PDMS_sample.listWell)), ellipsesalpha=1)
#PDMS_sample.drawAllSample(dots=True, ellipses=True, wellidx=range(len(PDMS_sample.listWell)), ellipsesalpha=1)
#PDMS_sample.drawAllSample(dots=True, ellipses=True, ellipsesalpha=1)

# ---- Draw Well :
#PDMS_sample.listWell[0].drawWell(dots=False, coat='full', coatalpha=0.5, ellipse=True, ellipsealpha=1)

# ---- Save & Plot Results :
#PDMS_sample.saveFile(name='8G', basicinfo=True, MVEE=True)
#title = str(nb_gratings) + 'G stretch' + str(2 * stretch) + 'mm' + ' '
#PDMS_sample.resultsHistogram(title, bars, info='VOLUME')
#PDMS_sample.resultsHistogram(title, bars, info='ELLIPTICITY')

# listHeight, listWidth, listAxis, listEllipticity, listVolume = [], [], [], [], []

# for i in range(len(PDMS_sample.listWell)):
#     if (PDMS_sample.listWell[i].nbOfNeighbors != 0) and (PDMS_sample.listWell[i].heightMin != -1):
#         if (PDMS_sample.listWell[i].ellipsoid != -1):
#             listAxis.append([PDMS_sample.listWell[i].ellipsoid.a, PDMS_sample.listWell[i].ellipsoid.b,
#                              PDMS_sample.listWell[i].ellipsoid.c])
#             listEllipticity.append([PDMS_sample.listWell[i].ellipsoid.ellipticity])
#             listVolume.append([PDMS_sample.listWell[i].ellipsoid.volume])
#
# arrayVolume = np.asarray(listVolume)
# arrayEllipticity = np.asarray(listEllipticity)
#
# volume = np.linspace(0.1, 2.0, bars)
# ellipt = np.linspace(0.0, 1.0, bars)
# z = np.zeros([len(volume), len(ellipt)])
# for i in range(len(volume)):
#     for j in range(len(ellipt)):
#         for k in range(len(arrayEllipticity)):
#             if arrayEllipticity[k] > ellipt[i] and arrayEllipticity[k] < ellipt[i + 1]:
#                 if arrayVolume[k] > volume[j] and arrayVolume[k] > volume[j + 1]:
#                     z[i,j] += 1
# print ('v', volume)
# print ('e', ellipt)
# print ('occ', z)
# -------------------------------------------------- 4 Gratings --------------------------------------------------- #

# print('Réseau 4G :')
#
# nb_gratings = 4
# phase = PI / nb_gratings
#
# s = np.array([[0.336, 0.573, 0.208, 1], [0.304, 0.700, 0.281, 1],
#               [0.394, 0.672, 0.292, 1], [0.201, 0.806, 0.363, 1]])
#
# p = np.array([[1.568, 1.702, 1.841, 1.701], [1.777, 1.775, 1.728, 1.744],
#               [2.052, 1.821, 1.599, 1.767], [2.189, 1.812, 1.532, 1.770]])
#
# grating0 = grating(s[stretch, 0], 0 * phase, p[stretch, 0])
# grating1 = grating(s[stretch, 1], 1 * phase, p[stretch, 1])
# grating2 = grating(s[stretch, 2], 2 * phase, p[stretch, 2])
# grating3 = grating(s[stretch, 3], 3 * phase, p[stretch, 3])
#
# list_gratings = ([grating0, grating1, grating2, grating3])
#
# PDMS_sample = sample(list_gratings, stretch, res, dimX, dimY, tol, filtering=False)
#
# # ---- Sample Analyses :
# #PDMS_sample.getWellsStats(basicinfo=True, MVEE=True)
#
# # ---- Draw All Sample :
# #PDMS_sample.drawAllSample(dots=False, coat='full', coatalpha=0.5, ellipses=True, wellidx=range(len(PDMS_sample.listWell)), ellipsesalpha=1)
# #PDMS_sample.drawAllSample(dots=True, ellipses=True, wellidx=range(len(PDMS_sample.listWell)), ellipsesalpha=1)
# #PDMS_sample.drawAllSample(dots=True, coat='full', ellipses=True, ellipsesalpha=1)
#
# # ---- Draw Well :
# #PDMS_sample.listWell[0].drawWell(dots=True, coat=False, coatalpha=0.5, ellipse=True, ellipsealpha=1)
#
# # ---- Save & Plot Results :
# #PDMS_sample.saveFile(name='4G', basicinfo=True, MVEE=True)
# title = str(nb_gratings) + 'G stretch ' + str(2*stretch) + 'mm' + ' '
# PDMS_sample.resultsHistogram(title, bars, info='VOLUME')
# PDMS_sample.resultsHistogram(title, bars, info='ELLIPTICITY')
