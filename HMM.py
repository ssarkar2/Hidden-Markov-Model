import numpy as np
class HMM():
    #note in class the states ent from 1 to 3. here it goes from 0 to 2
    #Similarly for observations, indexing starts from 0
    def __init__(self, A, B, Pi, O):
        self.A = A; self.B = B; self.Pi = Pi; self.O = O
        self.numstates = len(self.A); self.numobs = len(self.B[0]);  self.lengthobs = len(O)

    def calculateAlpha(self):
        A = self.A; B = self.B; O = self.O; Pi = self.Pi; numstates = self.numstates; numobs = self.numobs; lengthobs = self.lengthobs
        #this function iteratively generates alphas
        #initialize alpha(1,i), i = 0...numstates-1
        alpha = []
        alphatemp = []
        #for timestep 0  (note in class this is t = 1
        string = 'Printing Alpha calculation\n'
        for j in range(0,numstates): #[0,1,2,...numstates-1]
            string = string + str(Pi[j]) + ' * ' + str(B[j][O[0]]) + '=' + str(Pi[j] * B[j][O[0]]) + '\n'
            alphatemp.append(Pi[j] * B[j][O[0]])
        alpha.append(alphatemp)
        #print alpha

        string = string + '\n'
        for t in range(1,lengthobs):  #timesteps 1 to lengthobs-1. In class notation this is timesteps 2 to lengthobs
            alphatemp = []
            for j in range(0,numstates):
                sum = 0
                string = string + '('
                for i in range(0,numstates):
                    sum = sum + alpha[t-1][i]*A[i][j]
                    string = string + ' +' + str(alpha[t-1][i]) + '*' + str(A[i][j])
                string = string + ') * ' + str(B[j][O[t]]) + ' = ' + str(sum*B[j][O[t]]) + '\n'
                alphatemp.append(sum*B[j][O[t]])
            string = string + '\nxxxxxxxxxxxx\n'
            alpha.append(alphatemp)
            #print alphatemp
        self.alpha = alpha
        print string
        return alpha

    def findProbObservation(self):
        self.probObs = sum(self.alpha[-1])
        return self.probObs

    def calculateBeta(self):
        A = self.A; B = self.B; O = self.O; Pi = self.Pi; numstates = self.numstates; numobs = self.numobs; lengthobs = self.lengthobs
        #this function iteratively generates betas
        #initialize beta(lengthobs-1,i), i = 1...numstates
        beta = []
        string = 'Printing Beta calculation\n'
        for i in range(0,lengthobs):
            betatemp = []
            for j in range(0, numstates): #[0,1,2,...numstates-1]
                if i == lengthobs-1:
                    betatemp.append(1)
                    string = string + '1, '
                else:
                    betatemp.append(0)
                    string = string + '0, '
            beta.append(betatemp)
            string = string + '\n'

        for t in list(reversed(range(0,lengthobs-1))):
            betatemp = []
            for i in range(0,numstates):
                sum = 0
                for j in range(0,numstates):
                    string = string + ' + ' + str(A[i][j]) + '*' + str(B[j][O[t+1]]) + '*' + str(beta[t+1][j])
                    sum = sum + A[i][j]*B[j][O[t+1]]*beta[t+1][j]
                string = string + ' = ' + str(sum) + '\n'
                beta[t][i] = sum
            string = string + '\n yyyyyyyyyyyyyyyyyyyyyyy \n'
        #print beta
        self.beta = beta
        print string
        return beta

    def calculateGamma(self):  #not tested
        gamma = []
        string = 'Printing Gamma calculation\n'
        for t in range(0,self.lengthobs):
            gammatemp = [];
            for i in range(self.numstates):
                string = string + str(self.alpha[t][i]) + '*' + str(self.beta[t][i]) + '/' + str(self.probObs) + ' = '  +str(self.alpha[t][i]*self.beta[t][i]/self.probObs) + '   '
                gammatemp.append(self.alpha[t][i]*self.beta[t][i]/self.probObs)
            gamma.append(gammatemp)
            string = string + '\n'
            #print sum(gammatemp)  #should sum to 1
        self.gamma = gamma
        print string
        return gamma

    def calculateQ(self):  #not tested
        self.q = [self.gamma[t].index(max(self.gamma[t])) for t in range(0,self.lengthobs)]
        return self.q

    def viterbiForward(self):
        numstates = self.numstates; lengthobs = self.lengthobs
        delta = []
        deltatemp = []
        psi = []
        psitemp = []
        stringdelta = 'Printing Delta/Psi calculation\n'
        for j in range(0,numstates):
            deltatemp.append(Pi[j] * B[j][O[0]])
            psitemp.append(0)
            stringdelta = stringdelta + str(Pi[j]) + '*' + str(B[j][O[0]]) + '=' + str(Pi[j] * B[j][O[0]]) + '  maxindex 0\n '
        stringdelta = stringdelta + 'zzzzzzzzzzzzzzzzzz \n'
        delta.append(deltatemp); psi.append(psitemp)


        for t in range(1,lengthobs):
            deltatemp = []
            psitemp = []
            for j in range(0,numstates):
                maxtemp = []
                stringdelta = stringdelta + 'max('
                for i in range(0,numstates):
                    stringdelta = stringdelta + str(delta[t-1][i]) + '*' + str(A[i][j]) + ', '
                    maxtemp.append(delta[t-1][i]*A[i][j])
                stringdelta = stringdelta + ') * ' + str(B[j][O[t]]) + ' = ' + str(max(maxtemp)*B[j][O[t]]) + ', maxindex = ' + str(maxtemp.index(max(maxtemp)) + 1) +  '\n'
                psitemp.append(maxtemp.index(max(maxtemp)))
                deltatemp.append(max(maxtemp)*B[j][O[t]])
            stringdelta = stringdelta + '\n zzzzzzzzzzzz \n'
            delta.append(deltatemp); psi.append(psitemp)
        self.psi = psi
        self.delta = delta
        self.pstar = max(delta[lengthobs-1])
        self.qstar = delta[lengthobs-1].index(max(delta[lengthobs-1]))
        print stringdelta
        return [[psi], [delta], [self.pstar], [self.qstar]]

    def viterbiBackward(self):
        lengthobs = self.lengthobs
        qstarall = [-1] * lengthobs
        qstarall[lengthobs-1] = self.qstar
        for t in list(reversed(range(0,lengthobs-1))):
            qstarall[t] = self.psi[t+1][qstarall[t+1]]

        self.qstarall = qstarall#[::-1]
        return self.qstarall

    def calculateZeta(self):
        numstates = self.numstates; lengthobs = self.lengthobs
        zeta = [-1] * numstates; zeta = [zeta] * numstates; zeta = [zeta] * lengthobs
        alpha = self.alpha; beta = self.beta; probObs = self.probObs
        A = self.A; B = self.B
        zeta = []
        string = 'Printing Zeta calculation\n'
        for t in range(lengthobs-1):
            tmpi = []
            for i in range(numstates):
                tmpj = []
                for j in range(numstates):
                    string = string + '+' + str(alpha[t][i]) + '*' + str(A[i][j]) + '*' + str(B[j][O[t+1]]) + '*' + str(beta[t+1][j]) + '/' + str(probObs) + ' = ' + str((alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j])/probObs) + ', '
                    tmpj = tmpj + [(alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j])/probObs]
                string = string + '\n'
                tmpi = tmpi + [tmpj]
            string = string + '\n aaaaaaaaaaaaaaaaaaaaa \n'
            zeta = zeta + [tmpi]
        self.zeta = zeta
        print string
        return zeta


    def estimateModel(self, iterations, initA, initB, initPi):
        currA = initA; currB = initB; currPi = initPi
        numstates = self.numstates; numobs = self.numobs; lengthobs = self.lengthobs
        for iter in range(iterations):
            hmmmodeltemp = HMM(currA, currB, currPi, O)
            alphaCurr = hmmmodeltemp.calculateAlpha(); probObsCurr = hmmmodeltemp.findProbObservation()
            betaCurr = hmmmodeltemp.calculateBeta()
            gammaCurr = hmmmodeltemp.calculateGamma()
            zetaCurr = hmmmodeltemp.calculateZeta()

            #calculate Pi for this iteration
            newPi = gammaCurr[0]

            #calculate A for this iteration
            newA = []
            stringA = 'Estimate A\n'
            for i in range(numstates):
                tempA = []
                for j in range(numstates):
                    stringA = stringA + '('+'+'.join([str(zetaCurr[t][i][j]) for t in range(lengthobs-1)]) + ')/(' + '+'.join([str(gammaCurr[t][i]) for t in range(lengthobs-1)]) +') = ' + str(sum([zetaCurr[t][i][j] for t in range(lengthobs-1)])/sum([gammaCurr[t][i] for t in range(lengthobs-1)])) + '\n'
                    tempA = tempA + [sum([zetaCurr[t][i][j] for t in range(lengthobs-1)])/sum([gammaCurr[t][i] for t in range(lengthobs-1)])]
                newA = newA + [tempA]
                stringA = stringA + 'xxxxxxxxxxxxxxx\n'
                print stringA

            #calculate B for this iteration
            newB = []
            print 'Estimate B\n'
            stringB = ''; stringB2 = ''
            for j in range(numstates):
                tempB = []
                for k in range(numobs):
                    sumGammaNum = 0.0; sumGammaDenom = 0.0
                    for t in range(lengthobs):
                        if O[t] == k:
                            stringB = stringB + '+' + str(gammaCurr[t][j])
                            sumGammaNum = sumGammaNum + gammaCurr[t][j]
                        stringB2 = stringB2 + '+' + str(gammaCurr[t][j])
                        sumGammaDenom = sumGammaDenom + gammaCurr[t][j]
                    print '(' + stringB + ')/(' + stringB2 + ') = ' + str(sumGammaNum/sumGammaDenom) + '\n'
                    stringB = ''; stringB2 = ''
                    tempB = tempB + [sumGammaNum/sumGammaDenom]
                newB = newB + [tempB]

            currA = newA; currB = newB; currPi = newPi

        return [currA, currB, currPi]




A = [[0.4, 0.3, 0.3], [0.2, 0.7, 0.1], [0.3, 0.2, 0.5]]
B = [[0.5, 0.2, 0.2, 0.1], [0.1, 0.6, 0.2, 0.1], [0.4, 0.1, 0.2, 0.3]]
#Pi = [0.3, 0.4, 0.3]
#Pi = [0.2, 0.53, 0.27]
Pi = [0.34, 0.54, 0.12]
#R=0 G=1 B=2 Y=3
O = [1,0,0,2,3,1]

hmmmodel = HMM(A, B, Pi, O)
alpha = hmmmodel.calculateAlpha() #tested
print alpha
probObs = hmmmodel.findProbObservation() #tested
print probObs
beta = hmmmodel.calculateBeta() #tested
print beta
gamma = hmmmodel.calculateGamma()  #not tested
print gamma
q = hmmmodel.calculateQ()  #not tested
print q
[psi, delta, pstar, qstar] = hmmmodel.viterbiForward() #tested
print psi
print delta
print pstar
print qstar
qstarall = hmmmodel.viterbiBackward() #tested
print qstarall
zeta = hmmmodel.calculateZeta()
print np.array(zeta)  #printing as numpy array for prettyprint
[Aest, Best, Piest] = hmmmodel.estimateModel(1, A, B, Pi)

print np.array(Aest)
print np.array(Best)
print np.array(Piest)



##
##print 'xxxxxxxxxxxxxxHWxxxxxxxxxxx'
###homework
##A = [[0.8, 0.2], [0.2, 0.8]]
##B = [[0.4,0.1,0.4,0.1], [0.1,0.4,0.1,0.4]]
##Pi = [0.5, 0.5]
###A = 0, C = 1, G = 2, T = 3
##O = [1,2,3,1,0,2]
##hmmhw = HMM(A, B, Pi, O)
##alpha = hmmhw.calculateAlpha()
##print alpha
##probObs = hmmhw.findProbObservation()
##print probObs
##beta = hmmhw.calculateBeta()
##print beta
##gamma = hmmhw.calculateGamma()
##print gamma
##q = hmmhw.calculateQ()
##print q
##[psi, delta, pstar, qstar] = hmmhw.viterbiForward()
##print psi
##print delta
##print pstar
##print qstar
##qstarall = hmmhw.viterbiBackward() #tested
##print qstarall


