from cmath import sqrt, sinh
import pandas

# Le cours de l'action a dernier moment sk
def Sk(listSkMinusOne, dataset):
    listSk = []
    for i in range (len(listSkwave) - 1):
        listSk.append(listSkMinusOne[i+1])
    listSk.append(dataset.iloc[len(dataset) - 1, 4])
    return listSk

# Sk vague Sk-1 - h mais également somme des nkipi/nk où nki et pi sont les prix et volume instantannés
def Skwave(x, dataset):
    count = 0
    sumnip = 0
    sumnk = 0
    while count < len(dataset) and dataset.iloc[count, 0] < x / 24:
        if dataset.iloc[count, 0] > (x - 1)/24 and not pandas.isnull(dataset.iloc[count, 2]):
            sumnip = sumnip + dataset.iloc[count, 2] * dataset.iloc[count, 4]
            sumnk = sumnk + dataset.iloc[count, 2]
        count = count + 1
    return sumnip / sumnk

# Calcul de la volatilité
def volatility(dataset):
    mean = 0
    meansquares = 0
    for i in range (len(dataset)-1):
        mean = dataset.iloc[i+1,4]/dataset.iloc[i,4] - 1 + mean
        meansquares = (dataset.iloc[i+1,4]/dataset.iloc[i,4] - 1)**2 + meansquares
    mean = mean/(len(dataset) - 1)
    meansquares = meansquares/(len(dataset) - 1)
    # moyenne des carrés moins le carré de la moyenne le tout à la racine
    return sqrt(meansquares - mean**2).real

# sk vague = sk-1 - h, h = sk-1 - sk vague
def computeH(listOfSkWaves, dataset):
    listH =  []
    listSkminusOne = []
    toBeAdded = 0
    for i in range (24):
        if i == 0:
            listSkminusOne.append(dataset.iloc[0, 4])
        else :
            for j in range(len(dataset)):
                if dataset.iloc[j, 0] <= i/24 :
                    toBeAdded = dataset.iloc[j, 4]
            listSkminusOne.append(toBeAdded)
    # Pour tous les éléments de la liste
    for i in range (24):
        listH.append(listSkminusOne[i] - listOfSkWaves[i])
    return listSkminusOne, listH

def listOfnk(dataset):
    listnk = []
    for i in range(25):
        toBeAdded = 0
        if i != 0:
            for j in range(len(dataset)):
                if dataset.iloc[j, 0] > (i - 1) / 24 and dataset.iloc[j, 0] <= i / 24 and not pandas.isnull(dataset.iloc[j, 2]):
                    # volume * signe de l'opération
                    toBeAdded = dataset.iloc[j,2]*dataset.iloc[j, 3] + toBeAdded
            # - car on veut le volume à liquider donc ce que l'on vend
            listnk.append(toBeAdded * -1)
    return listnk

# liste des tau g car sk = sk-1 + terme erreur - tau g
def listRandomArithmeticDynamic(Sk, Skminusone):
    list = []
    for i in range (len(Sk)):
        list.append(Skminusone[i] - Sk[i])
    # liste des tau g
    return list

# Calcul de l'espérance avec la formule du cours
def expectation(xks, nks, dynamics, hs):
    expectation = 0
    for i in range (len(xks)):
        expectation = expectation + xks[i]*dynamics[i] + nks[i]*hs[i]
    return expectation

# Idem pour la variance
def variance(std, xks, tau):
    sumOfxktau = 0
    for i in range (len(xks)):
        sumOfxktau = sumOfxktau + tau * xks[i] * xks[i]
    return std*std*sumOfxktau

# somme des tau g xk = gamma * 1/2 ( X^2 - somme des nk carrés), isoler gamma
def computegamma(dynamics, x, nks, xks):
    sumdyn = 0
    sumnksquare = 0
    for i in range (len(dynamics)):
        sumdyn = sumdyn + dynamics[i] * xks[i]
        sumnksquare = sumnksquare + nks[i]*nks[i]
    return 2*sumdyn/(x*x - sumnksquare)

# h1 et h2 = xi * signe de nk + eta * nk / tau, pour trouver eta, il suffit de trouver 2 nk à signe opposés et addtionner leurs h
def computeEta(listH, nks, tau):
    indexpositive = 0
    indexnegative = 0
    for i in range (len(nks)):
        if nks[i] < 0:
            indexnegative = i
        else:
            indexpositive = i
    return (listH[indexpositive] + listH[indexnegative])*tau/(nks[indexpositive]+nks[indexnegative])

# Etant donné qu'on a eta, tau, h et nk on peut facilement trouver xi
def computeXi(listH, nks, tau, eta):
    indexnegative = 0
    for i in range (len(nks) - 1):
        if nks[i] < 0:
            indexnegative = i
    return listH[indexnegative] - eta * nks[indexnegative]/tau

# Espérance pour la variane nulle
def expForVzero(eta, xi, x, tau):
    return xi*x + eta*x*x/tau

# Minimisation par étude du signe de la dérivée, extremum pour dérivée nulle
def minimize(tau, nks, gamma, lbda, sigma):
    listxks = []
    for i in range (len(nks)):
        listxks.append(gamma*nks[i]/(lbda*sigma*sigma*tau*2))
    return listxks

# xk avec la formule du sinus hyperbolique
def computeXk(X, lbda, sigma, eta, tau):
    listxks = []
    for i in range(24):
        listxks.append((sinh(sqrt(lbda * sigma * sigma/ eta)*(1 - ((i+1) - 1/2*tau)))*
                        X/sinh(sqrt(lbda * sigma * sigma/ eta))).real)
    return listxks

# lecture sur fichier
data = pandas.read_csv("Dataset TD4.csv")
# print(data)

listSkwave = []
listnk = listOfnk(data)
sumOfNk = 0
# Somme des nk
for i in range (24):
    listSkwave.append(Skwave(i + 1, data))
    sumOfNk = sumOfNk + listnk[i]

# trouver les sk moins 1 et les h
listSkMinusOne, listH = computeH(listSkwave, data)
listxk = []
toBeAdded = sumOfNk
# liste des xk = X - somme des nk
for i in range (24):
    # cummulé
    toBeAdded = toBeAdded - listnk[i]
    listxk.append(toBeAdded)

# Création des variables une par une
standardDev = volatility(data)
Sks = Sk(listSkMinusOne, data)
listdynamic = listRandomArithmeticDynamic(Sks, listSkMinusOne)
expected = expectation(listxk, listnk, listdynamic, listH)
var = variance(standardDev, listxk, 1/24)
gamma = computegamma(listdynamic, sumOfNk, listnk, listxk)
eta = computeEta(listH, listnk, 1/24)
xi = computeXi(listH, listnk, 1/24, eta)
expectedForV0 = expForVzero(eta, xi, sumOfNk, 1/24)
lbda = (expectedForV0 - expected)/(0 - var)

# Afffichage de vérification
print(listxk)
print(listSkwave)
print(Sks)
print(listdynamic)
print(expected)
print(var)
print(gamma)
print(eta)
print(xi)
print(expectedForV0)
print(lbda)
print(expectedForV0)

# Fichier avec toutes les listes utilisées
df1 = pandas.DataFrame(list(zip(listnk, listxk, Sks, listSkMinusOne, listSkwave, listdynamic, listH)), columns=['nk', 'xk', 'Sk', 'S(k-1)', '~Sk', 'Txkg(nk/T)', 'h(nk/T)'])
df1.to_csv('TD4 elements.csv')

# liste des paramètres estimés
df2 = pandas.DataFrame(list(zip(['standardDeviation', 'expectation', 'variance', 'gamma', 'eta', 'xi', 'expectation for V = 0', 'lambda'], [standardDev, expected, var, gamma, eta, xi, expectedForV0, lbda])), columns=['Parameters', 'Values'])
df2.to_csv('TD4 parameters.csv')

# xk par sinus hyperbolique
xkSinhyp = computeXk(sumOfNk, lbda, standardDev.real, eta, 1/24)
print(xkSinhyp)
expected2 = expectation(xkSinhyp, listnk, listdynamic, listH)
var2 = variance(standardDev, xkSinhyp, 1/24)
lbda2 = (expected - expected2)/(var - var2)
print(lbda2)

# estimation des xk avec les nouveau lambdas
estimatedxks = minimize(1/24, listnk, gamma, lbda2, standardDev)
print(estimatedxks)

# Output avec les estimation des xk à chaque heure
df3 = pandas.DataFrame(list(zip(listnk, listxk, estimatedxks, xkSinhyp, Sks)), columns=['Original nk', 'To be liquidated (xk)', 'Efficient xk', 'Hyperbolic sinus', 'Price'])
df3.to_csv('TD4 Estimation xk.csv')

