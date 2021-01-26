import random
from cmath import log, sqrt

import pandas

# Variance en basse fréquence
def calculateVarLOW(x):
    return log(((x.iloc[len(x) - 1, 2] + x.iloc[len(x) - 1, 1]) / 2) / ((x.iloc[0, 2] + x.iloc[0, 1]) / 2)) ** 2


# Sur 134 jours
def calculateVarHIGH(x):
    sum = 0
    for i in range(134):
        # 96 records avant de retourner a même moment de la journée
        sum = sum + log(((x.iloc[96 * (i + 1), 1] + x.iloc[96 * (i + 1), 2]) / 2) / (
                (x.iloc[96 * i, 1] + x.iloc[96 * i, 2]) / 2)) ** 2
    sum = sum / (134 - 1)
    return sum.real

# Mk2 avec la fomule du cours
def Mk2High(x):
    sum = 0
    for i in range(134):
        sum = sum + abs(((x.iloc[96 * (i + 1), 1] + x.iloc[96 * (i + 1), 2]) / 2) - (
                (x.iloc[96 * i, 1] + x.iloc[96 * i, 2]) / 2)) ** 2
    sum = sum / (134 - 1)
    return sum.real

# Mk prime avec la formule du cours
def MkPRIME2High(x):
    sum = 0
    for i in range(int(134 / 2)):
        sum = sum + abs(
            ((x.iloc[2 * 96 * (i + 1), 1] + x.iloc[2 * 96 * (i + 1), 2]) / 2) - (
                    (x.iloc[2 * 96 * i, 1] + x.iloc[2 * 96 * i, 2]) / 2)) ** 2
    sum = sum / (134 - 1) * 2
    return sum.real

# En basse fréquence
def Mk2Low(x):
    return abs(((x.iloc[len(x) - 1, 2] + x.iloc[len(x) - 1, 1]) / 2) - ((x.iloc[0, 2] + x.iloc[0, 1]) / 2)) ** 2

# basse fréquence
def MkPRIME2Low(x):
    return abs(((x.iloc[len(x) - 1, 2] + x.iloc[len(x) - 1, 1]) / 2) - ((x.iloc[0, 2] + x.iloc[0, 1]) / 2)) ** 2 * 2

# calcul des ck par ondelettes
def calculateCk(x, data, i):
    sum = abs(
        ((data.iloc[int(96 * x * (i + 1)), 1] + data.iloc[int(96 * x * (i + 1)), 2]) / 2) / (
                (data.iloc[int(96 * i * x), 1] + data.iloc[int(96 * x * i), 2]) / 2))
    return sum

# lecture sur fichiers
GBPEUR = pandas.read_csv("Dataset TD5.csv", usecols=[0, 1, 2], header=1)
print(GBPEUR)

SEKEUR = pandas.read_csv("Dataset TD5.csv", usecols=[4, 5, 6], header=1)
print(SEKEUR)

CADEUR = pandas.read_csv("Dataset TD5.csv", usecols=[8, 9, 10], header=1)
print(CADEUR)

# Calcul des exposants de hurst
print("\n\nHurst Exponent")
hurstHighGBPEUR = 1 / 2 * log(MkPRIME2High(GBPEUR) / Mk2High(GBPEUR)) / log(2)
hurstHighSEKEUR = 1 / 2 * log(MkPRIME2High(SEKEUR) / Mk2High(SEKEUR)) / log(2)
hurstHighCADEUR = 1 / 2 * log(MkPRIME2High(CADEUR) / Mk2High(CADEUR)) / log(2)

hurstLowGBPEUR = 1 / 2 * log(MkPRIME2Low(GBPEUR) / Mk2Low(GBPEUR)) / log(2)
hurstLowSEKEUR = 1 / 2 * log(MkPRIME2Low(SEKEUR) / Mk2Low(SEKEUR)) / log(2)
hurstLowCADEUR = 1 / 2 * log(MkPRIME2Low(CADEUR) / Mk2Low(CADEUR)) / log(2)

# Affichages
print(hurstHighGBPEUR.real)
print(hurstHighSEKEUR.real)
print(hurstHighCADEUR.real)
print(hurstLowGBPEUR.real)
print(hurstLowSEKEUR.real)
print(hurstLowCADEUR.real)

# Affichage des varainces
print("\nGBP - EUR")
dailyGBPEURLOW = sqrt(calculateVarLOW(GBPEUR)).real
dailyGBPEURHIGH = sqrt(calculateVarHIGH(GBPEUR)).real
annualGBPEURHIGH = (dailyGBPEURHIGH * sqrt(135) ** hurstHighGBPEUR).real

print(dailyGBPEURLOW)
print(dailyGBPEURHIGH)
print(annualGBPEURHIGH)

print("\n\nSKE - EUR")
dailySEKEURLOW = sqrt(calculateVarLOW(SEKEUR)).real
dailySEKEURHIGH = sqrt(calculateVarHIGH(SEKEUR)).real
annualSEKEURHIGH = (dailySEKEURHIGH * sqrt(135) ** hurstHighSEKEUR).real

print(dailySEKEURLOW)
print(dailySEKEURHIGH)
print(annualSEKEURHIGH)

print("\n\nCAD - EUR")
dailyCADEURLOW = sqrt(calculateVarLOW(CADEUR)).real
dailyCADEURHIGH = sqrt(calculateVarHIGH(CADEUR)).real
annualCADEURHIGH = (dailyCADEURHIGH * sqrt(135) ** hurstHighCADEUR).real

print(dailyCADEURLOW)
print(dailyCADEURHIGH)
print(annualCADEURHIGH)

# Liste stockant les covariances et variances
listCovGBPCAD = []
listCovSEKCAD = []
listCovSEKGBP = []
listVarGBP = []
listVarSEK = []
listVarCAD = []
listT = []

# 200 points aléatoires pour calculer l'intégrale
for i in range(200):
    cjGBP = 0
    cjSEK = 0
    cjCAD = 0
    t = random.uniform(0, 1)
    for j in range(int(134)):
        cjGBP = cjGBP + calculateCk(t ** hurstHighGBPEUR.real, GBPEUR, j)
        cjSEK = cjSEK + calculateCk(t ** hurstHighSEKEUR.real, SEKEUR, j)
        cjCAD = cjCAD + calculateCk(t ** hurstHighCADEUR.real, CADEUR, j)
    cjlGBP = cjGBP / 134
    cjlCAD = cjCAD / 134
    cjlSEK = cjSEK / 134
    totGBPCAD = 0
    totSEKCAD = 0
    totSEKGBP = 0
    totGBP2 = 0
    totSEK2 = 0
    totCAD2 = 0
    for j in range(int(134)):
        # Calcul des variances covariances
        totGBPCAD = totGBPCAD + (calculateCk(t ** hurstHighGBPEUR.real, GBPEUR, j) - cjlGBP) * (
                calculateCk(t ** hurstHighCADEUR.real, CADEUR, j) - cjlCAD)
        totSEKCAD = totSEKCAD + (calculateCk(t ** hurstHighSEKEUR.real, SEKEUR, j) - cjlSEK) * (
                calculateCk(t ** hurstHighCADEUR.real, CADEUR, j) - cjlCAD)
        totSEKGBP = totSEKGBP + (calculateCk(t ** hurstHighGBPEUR.real, GBPEUR, j) - cjlGBP) * (
                calculateCk(t ** hurstHighSEKEUR.real, SEKEUR, j) ** 2 - cjlSEK)
        totGBP2 = totGBP2 + (calculateCk(t ** hurstHighGBPEUR.real, GBPEUR, j) - cjlGBP) ** 2
        totCAD2 = totCAD2 + (calculateCk(t ** hurstHighCADEUR.real, CADEUR, j) - cjlCAD) ** 2
        totSEK2 = totSEK2 + (calculateCk(t ** hurstHighSEKEUR.real, SEKEUR, j) - cjlSEK) ** 2
    totGBPCAD = totGBPCAD / 134
    totSEKCAD = totSEKCAD / 134
    totSEKGBP = totSEKGBP / 134
    totSEK2 = totSEK2 / 134
    totGBP2 = totGBP2 / 134
    totCAD2 = totCAD2 / 134
    # Ajout des variances et covariances simulées
    listVarCAD.append(totCAD2)
    listVarSEK.append(totSEK2)
    listVarGBP.append(totGBP2)
    listCovSEKCAD.append(totSEKCAD)
    listCovSEKGBP.append(totSEKGBP)
    listCovGBPCAD.append(totGBPCAD)
    listT.append(t)

# Fichier avec les simulatons
df = pandas.DataFrame(list(zip(listT, listVarGBP, listVarSEK, listVarCAD, listCovSEKGBP, listCovSEKCAD, listCovGBPCAD)),
                      columns=['Window', 'Variance GBP', 'Variance SEK', 'Variance CAD', 'Covariance SEK - GBP',
                               'Covariance SEK - CAD', 'Covariance GBP - CAD'])
df.to_csv('TD5 Variance - Covariance.csv')

# valeurs minimales de variances et de covariances
minGBP = min(df['Variance GBP'].tolist())
minSEK = min(df['Variance SEK'].tolist())
minCAD = min(df['Variance CAD'].tolist())
minSEKGBP = min(df['Covariance SEK - GBP'].tolist())
minSEKCAD = min(df['Covariance SEK - CAD'].tolist())
minGBPCAD = min(df['Covariance GBP - CAD'].tolist())

# Les temps d'investissement pris en cosidération pour variance minimale
minTGBP = df[df['Variance GBP'] == minGBP].iloc[0,0]
minTSEK = df[df['Variance SEK'] == minSEK].iloc[0,0]
minTCAD = df[df['Variance CAD'] == minCAD].iloc[0,0]
minTSEKGBP = df[df['Covariance SEK - GBP'] == minSEKGBP].iloc[0,0]
minTSEKCAD = df[df['Covariance SEK - CAD'] == minSEKCAD].iloc[0,0]
minTGBPCAD = df[df['Covariance GBP - CAD'] == minGBPCAD].iloc[0,0]

# Stockage dans une dataframe et output avec nos résultats
df2 = pandas.DataFrame(list(zip(['Min VAR GBP', 'Min VAR SEK', 'Min VAR CAD', 'Min Cov SEK GBP', 'Min Cov SEK CAD', 'Min Cov GBP CAD'], [minGBP, minSEK, minCAD, minSEKGBP, minSEKCAD, minGBPCAD],
                                [minTGBP, minTSEK, minTCAD, minTSEKGBP, minTSEKCAD, minTGBPCAD])),
                      columns = ['Name', 'Value', 'Window'])
df2.to_csv('TD5 minimum values.csv')

