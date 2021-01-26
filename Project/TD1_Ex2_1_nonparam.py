import csv
from cmath import exp, pi, sqrt
import random

import pandas

# fonction de densité kernel
def kernelDensity(x):
    return (1.0 / sqrt(2 * pi) * exp(-x * x / 2)).real

# fonction intégrale de kernel
def kernelCDF(yi, y0):
    sum_total = 0.0
    for i in range(300):
        # valeur aléatoires dans l'interval pour pouvoir calculer l'intégrale
        random_x = random.uniform(y0, yi)
        # Intégrale par MonteCarlo avec 300 points
        sum_total = sum_total + kernelDensity(random_x)
    # Résultat de l'intégralepar Monte Carlo
    return ((yi - y0) / 300 * sum_total).real

# Fonction de répartition
def FApproximation(data_list, x, h):
    sum_integral = 0.0
    for i in range(len(data_list)):
        # x de - 3 à (x - xi)/h car -3 3 contient 99% des valeurs de la fonction gaussienne (kernel gaussien ici)
        sum_integral = sum_integral + kernelCDF((x - data_list[i]) / h, -3)
    # Intégrale par monte carlo
    result = sum_integral * 1 / len(data_list)
    return result.real

# liste des quantiles pour lesquels nous voulons estimer une var
listQuantile = [0.999, 0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7]
# lecture sur fichier
dataframe = pandas.read_csv("MR_Ex2_MonteCarlo.csv")
print(dataframe)
# Liste des prix
listPrice = dataframe['Price return'].tolist()
print(listPrice)
list1 = []
list2 = []
list3 = []
list4 = []

# Listes des VaT
VaRlisth1 = []
VaRlisth2 =[]
VaRlisth3 = []

# Pour chacun des h faire une approximation de F
for i in range(len(listPrice)):
    print(i)
    # Avec la liste des prix (historique)
    x = listPrice[i]
    list1.append(x.real)
    list2.append(FApproximation(listPrice, x, 0.001))
    list3.append(FApproximation(listPrice, x, 0.003))
    list4.append(FApproximation(listPrice, x, 0.01))
# Créer un fichier output avec toutes les valeurs générées
df = pandas.DataFrame(list(zip(list1, list2,list3, list4)), columns=['returns', 'CDF for h = 0.001', 'CDF for h = 0.003', 'CDF for h = 0.01'])
# Trier par return
df = df.sort_values(by=['returns'])
print(df)
df.to_csv('TD1 outputPrices.csv')

# Création des VaR avec la liste des quantiles
for i in range (len(listQuantile)) :
    # VaR = F^-1(1 - alpha)
    VaRlisth1.append(-1 * df[df['CDF for h = 0.001'] >= 1 - listQuantile[i]].iloc[0,0])
    VaRlisth2.append(-1 * df[df['CDF for h = 0.003'] >= 1 - listQuantile[i]].iloc[0,0])
    VaRlisth3.append(-1 * df[df['CDF for h = 0.01'] >= 1 - listQuantile[i]].iloc[0,0])

# Création du fichier output avec les VaR
dfVar = pandas.DataFrame(list(zip(listQuantile, VaRlisth1, VaRlisth2, VaRlisth3)),
                         columns=['Quantile', 'VaR for h = 0.001',
                                    'VaR for h = 0.003', 'VaR for h = 0.01'])
dfVar.to_csv('TD1 outputVaR.csv')

# Réinitialisation des listes
list1 = []
list2 = []
list3 = []
list4 = []
VaRlisth1 = []
VaRlisth2 =[]
VaRlisth3 = []

# Création de 2500 pour avoir un modèle plus précis
for i in range(2500):
    print(i)
    # returns entre -1 et 1
    x = random.uniform(-1, 1)
    list1.append(x.real)
    list2.append(FApproximation(listPrice, x, 0.001))
    list3.append(FApproximation(listPrice, x, 0.003))
    list4.append(FApproximation(listPrice, x, 0.01))
# création de fichier output de la même manière que le précédent
df = pandas.DataFrame(list(zip(list1, list2,list3, list4)), columns=['returns', 'CDF for h = 0.001', 'CDF for h = 0.003', 'CDF for h = 0.01'])
df = df.sort_values(by=['returns'])
print(df)
df.to_csv('TD1 outputPricesRandom.csv')

# Avec le même raisonnement, créer les Values at risk
for i in range (len(listQuantile)) :
    VaRlisth1.append(-1 * df[df['CDF for h = 0.001'] >= 1 - listQuantile[i]].iloc[0,0])
    VaRlisth2.append(-1 * df[df['CDF for h = 0.003'] >= 1 - listQuantile[i]].iloc[0,0])
    VaRlisth3.append(-1 * df[df['CDF for h = 0.01'] >= 1 - listQuantile[i]].iloc[0,0])

dfVar = pandas.DataFrame(list(zip(listQuantile, VaRlisth1, VaRlisth2, VaRlisth3)), columns=['Quantile', 'VaR for h = 0.001', 'VaR for h = 0.003', 'VaR for h = 0.01'])
dfVar.to_csv('TD1 outputVaRRandom.csv')