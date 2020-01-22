#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import yaml

################################## Hypothèses d'étude ###################################
# La fusée est infiniment rigide
# La masse est répartie de manière uniforme
# Dans le repère terrestre, le vent est horizontal
# Le vol s'effectue dans une atmosphère standard, où la densité rho est connye pour une altitude donnée

################################## Données de la fusée ##################################
    
# longueur_total = 1000       # longueur totale de la fusée (ogive + corps + ailerons) en mm
# diametre_ext = 80 	      # Diamètre extérieure de la fusée en mm
# longueur_corps = 750        # longueur du corps de la fusée en mm

# masse_fusee_min = 1500      # masse totale de la fusée en g sans le propulseur
# masse_fusee_max = 2400      # masse maximale de la fusée pour stabilité
# masse_corps = 836 g         # masse du corps de la fusée en g
# masse_ailerons = 400 g      # masse des 4 ailerons, 4*100, en g
# masse_bagues_medium = 38 g  # masse des bagues en medium, 7 en tout

# m = 150                     # Emplanture des ailerons en mm
# n = 100                     # saumon des ailerons en mm
# p = 100                     # flèche des ailerons en mm
# e = 130                     # Envergure des ailerons en mm

# epaisseur_ailerons = 2      # épaisseur des ailerons en mm
# nbre_ailerons = 4           # nombre d'ailerons

########################### valeurs connues de la fusee #################################

with open('caracteristique.yml', 'r') as f:
    c = yaml.load(f.read())

l = c['l'] # longueur de l'ogive
dref = c['D'] # diamètre extérieur de la fusée 
dail = c['D'] # diamètre de référence pour les ailerons autour du corps
m = c['m'] # emplanture des ailerons
n = c['n'] # saumon des ailerons
e = c['e'] # envergure des ailerons
p = c['p'] # flèche des ailerons
Q = c['nb_ailerons'] # nombre d'ailerons
Xail = c['Xail'] # longueur du corps (de la pointe de l'ogive à l'encastrure supérieur des ailerons)
taille = 1 # longueur totale de la fusée (en m)

Xcg_sans = c['Xcg_sans'] # distance au centre de masse à partir du haut de l'ogive sans le propulseur
Xcg_avec = c['Xcg_avec'] # distance au centre de masse à partir du haut de l'ogive avec le propulseur


########################### valeurs obtenues avec StabTraj ##############################
# ST = StabTraj : correspond aux données fournies par la feuille de calcul StabTraj
 
Cna_ST = 22.4 # coefficient de portance
Xcpa_ST = 0.808 # position de centre de poussée aérodynamique par rapport à la pointe de l'ogive, en m

CnaAil_ST = 20.4 # coefficient de portance des ailerons
CnaCoiffe_ST = 2.0 # coefficient de portance pour la coiffe ( = (diametre_ogive / diametre_ailerons)**2 )

XcpaAil_ST = 0.878 # position de centre de poussée pour les ailerons
XcpaCoiffe_ST = 0.093 # position de centre de poussée pour la coiffe

MSCna_ST = 2.85 # marge de stabilité

################################ calculs de stabilite ###################################

finesse = taille/dref # rapport de la taille sur le diamètre de la fusée

f = np.sqrt(e**2+(p+(n-m)/2)**2) # dépend de la géométrie des ailerons
CnaAil = (1+dail/(2*e+dail))*(4*Q*(e/dref)**2)/(1+np.sqrt(1+(2*f/(m+n))**2))
XcpaAil = Xail + p*(m+2*n)/(3*(m+n))+1/6*(m+n-m*n/(m+n)) # centre de portance aérodynamique de l'empennage (dépend des dimensions des ailerons de la fusée)
XcpaCoiffe = 7/15*l

# calculs de X CPA et du Cn_alpha
CnaCoiffe = 2
Cna = CnaCoiffe + CnaAil
Xcpa = (XcpaCoiffe*CnaCoiffe + XcpaAil*CnaAil)/(CnaCoiffe + CnaAil)


def MS(Xcg): # Marge statique
    return abs(Xcg_sans-Xcpa)/dref


MSCna = MS(Xcg_sans)*Cna # Couple de rappel de la portance

# tests stabilite à partir de la donnée de la finesse (doit être compris entre 10 et 20 pour stabilité)

if finesse < 10 or finesse > 20:
    if taille < 10:
        print("Fusée trop petite")
    else:
        print("fusée trop grande")
elif Cna < 15 or Cna > 30:
    if Cna < 15:
        print("criteres de stabilite non respecte (Cna<15)")
    elif Cna > 30:
        print("criteres de stabilite non respecte (Cna>30)")
elif MS(Xcg_sans) < 1.5 or MS(Xcg_sans) > 6:
    print(f"criteres de stabilite non respecte : MS = {MS(Xcg_sans)}")
elif MSCna < 30 or MSCna > 100:
    print(f"criteres de stabilite non respecte : MSCna = {MSCna}")
else:
    print("fusee stable")
    print(f"finesse : {finesse} \nf : {f} \nCna : {Cna}")
    print(f"MS : {MS(Xcg_sans)} \nMS.Cna : {MSCna}")
    print(f"CnaAilerons : {CnaAil} \nXcpaAilerons : {XcpaAil}")
    print(f"Xcpa : {Xcpa}")
    

tMax = 15 # taille maximale du vecteur temps utilisée dans l'affichage des courbes
mTotal = 1.9 # masse de la fusée
mCombustible = 0.0843 # masse du combustible
tCombustion = 1.0 # temps de combustion du propulseur


def P(t): # fonction qui calcule la poussée en fonction du temps, à partir des caractéristiques du propulseur (Cesaroni Pro-24 6G 143G-150 "Pandora")
    if t < 0:
        return 0
    elif t < .04:
        return 250-50/.04*t
    elif t < .68:
        return 200-50/(.68-.4)*(t-.4)
    elif t < 0.84:
        return 150-100/(.84-.68)*(t-.68)
    elif t < 1.04:
        return 50-50/(1.04-.84)*(t-.84)
    else:
        return 0


g = 9.81 # accélération de la pesanteur
rho0 = 1.013  # densité de l'air
Cx = 0.6 # coefficient de trainée de la fusée 
Cyb = Cna  # Cy_beta
Cza = Cna  # Cz_alpha
LongueurTube = 0.750 # longueur du tube en mètres
D = 0.080  # Diamètre de la fusée
Dogive = D  # Diamètre de l'ogive
L = 0.130  # Envergure des ailerons
EPaileron = 0.002  # Épaisseur des ailerons
Strainee = np.pi*D**2/4+4*L*EPaileron  # Surface de trainée (maitre couple)
Sreference = 0.006067  # Surface de référence : (np.pi * D**2 / 4)
LongueurRampe = 2.5 # longueur de la rampe de lancement

Yf, Zf, Lf, Mf, Nf = [0]*5  # Force et moment appliqué à la fusée : poussée 
Vvent = 15 # Vitesse du vent
epsilon = 0 # Direction du vent


def F(Vect, t):
	# 6 degrés de liberté au total : 3 translations et 3 rotations
	# A chaque degré de liberté on associe 2 composantes du vecteur d'état : on a donc 12 composantes au total pour le vecteur d'état
	# x, y, z : repère de la fusée
	# u, v, w : vitesse dans le repère lié à la fusée
	# phi, theta, psi : angles liant le repère fusée au repère terrestre
	# p, q, r : vitesse de rotation autour de chaque axe de la fusée
	
    """Fonction de dérivation du vecteur d'état Vect en fonction de t : en intégrant u, v et w dans le référentiel terrestre, on en déduit l'évolution du vecteur position de la fusée dans le repère terrestre"""
    
    x, y, z, u, v, w, phi, theta, psi, p, q, r = Vect
    if t < tCombustion:
        m = mTotal - t*mCombustible/tCombustion
    else:
        m = mTotal - mCombustible
    # coefficients de la matrice d'inertie (diagonale de la matrice)
    A = m*D**2/2 # pondéré par la masse
    B = m*LongueurTube**2/12 # idem
    C = m*LongueurTube**2/12 # idem
    
    ly = MS(Xcg_sans)*D # lacet
    lz = MS(Xcg_sans)*D # tangage

    Vvent_x = Vvent*np.cos(epsilon)*np.cos(psi)*np.cos(theta) - Vvent*np.sin(epsilon)*np.sin(psi) # composante suivant x de la vitesse du vent
    Vvent_y = Vvent*np.cos(epsilon)*(np.sin(psi)*np.cos(phi) + np.sin(phi)*np.sin(theta))+ Vvent*np.sin(epsilon)*np.cos(psi)*np.cos(phi) # composante suivant y de la vitesse du vent
    Vvent_z = Vvent*np.cos(epsilon)*(-np.sin(psi)*np.sin(phi) + np.sin(theta)*np.cos(phi)) - Vvent*np.sin(epsilon)*np.cos(psi)*np.sin(phi) # composante suivant z de la vitesse du vent
    Vrx = -u+Vvent_x # expression de la vitesse relative suivant x
    Vry = -v+Vvent_y # idem suivant y
    Vrz = -w+Vvent_z # idem suivant z

    rho = rho0*(20000+z)/(20000-z) # variation de la masse volumique de l'air
    Xa = -rho*Strainee*Vrx**2*Cx/2 # traînée
    # Xf = 146.7 if t < 0.97 else 0
    Xf = P(t)

    du = 1/m*(Xa+Xf-m*g*np.sin(theta))+r*v-q*w
    if np.sqrt(x**2+y**2+z**2) < LongueurRampe:  # Test si la fusée est toujours dans la rampe de lancement
        dv, dw, dp, dq, dr = [0]*5
    else:
        alpha = -np.arctan(w/u)
        beta = np.arctan(v/u*np.cos(np.arctan(w/u))) 
        Ya = -rho*Sreference*Vry**2*Cyb*beta/2 # force latérale
        Za = -rho*Sreference*Vrz**2*Cza*alpha/2 # portance
        
        Ma = rho*Sreference*Vrz**2*Cza*alpha*lz/2 # moment de tangage
        Na = rho*Sreference*Vry**2*Cyb*beta*ly/2 # moment de lacet
        

        dv = 1/m*(Ya+Yf+m*g*np.cos(theta)*np.sin(phi))+p*w-r*u
        dw = 1/m*(Za+Zf+m*g*np.cos(theta)*np.cos(phi))+q*u-p*v
        
        dp = 1/A*(Lf+(B-C)*q*r)
        dq = 1/B*(Ma+Mf+(C-A)*p*r)
        dr = 1/C*(Na+Nf+(A-B)*p*q)

    return np.array([*RtoR0(Vect), du, dv, dw, p, q, r, dp, dq, dr])

def RtoR0(Vect):
    """Renvoie les vitesses de la fusée dans le repère terrestre"""
    cphi, ctheta, cpsi = np.cos(Vect[6:9])
    sphi, stheta, spsi = np.sin(Vect[6:9])
    T = np.array([[cpsi*ctheta, -spsi*cphi+cpsi*stheta*sphi, spsi*sphi+cpsi*stheta*cphi],
               [spsi*ctheta, cpsi*cphi+spsi*stheta*sphi, -cpsi*sphi+spsi*stheta*cphi],
               [-stheta, ctheta*sphi, ctheta*cphi]])
    return np.dot(T, Vect[3:6])

############ Simulation du vol de la fusée avec différentes direction de vent ###########

if __name__ == "__main__":

    ################
    # Affichage 2D #
    ################
    
    t = np.linspace(0, tMax, 200)
    for epsilon in np.linspace(0, np.pi, 5): # nombre de trajectoire en fonction de la vitesse du vent
        res = odeint(F, np.array([0, 0, 0, 0, 0, 0, 0, 1.396, 0, 0, 0, 0]), t)
        plt.plot(res[:, 0], -res[:, 2])
    plt.title("Trajectoire de la fusée en fonction du temps, sans ouverture de parachute")
    plt.xlabel("Distance parcourue (en mètres)")
    plt.ylabel("Altitude (en mètres)")
    plt.show()

    ################
    # Affichage 3D #
    ################
    
    # mpl.rcParams['legend.fontsize'] = 10
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # t = np.linspace(0, tMax, 100)
    # for epsilon in np.linspace(0, np.pi, 1):
    #     res = odeint(F, np.array([0, 0, 0, 0, 0, 0, 0, np.pi*80/180, 0, 0, 0, 0]), t)
    #     ax.plot(res[:, 0], res[:, 1], -res[:, 2])
    # ax.set_xlabel("Distance parcourue (en mètres)")
    # ax.set_ylabel("Dérive (en mètres)")
    # ax.set_zlabel(r"Altitude (en mètres)")
    # plt.title("Trajectoire de la fusée, sans ouverture de parachute")
    # plt.show()
