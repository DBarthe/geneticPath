#!/usr/bin/python

import sys
from random import randint, random
from itertools import permutations
from math import sqrt

PLACES_NBR = 50

class Geo:
    def __init__(self):
        self._distances = []
        self._placesNbr = 0

    @property
    def distances(self):
        return self._distances

    @property
    def placesNbr(self):
        return self._placesNbr

    def __initDistances(self, placesNbr):
        self._placesNbr = placesNbr
        self._distances = [[ 0 for b in range(self._placesNbr)] for a in range(self._placesNbr)]

    def distance(self,a,b):
        return self._distances[a][b]

    def makeFakesDistances(self, placesNbr):
        self.__initDistances(placesNbr)
        for a in range(self._placesNbr):
            for b in range(a+1,self._placesNbr): 
                d = randint(10,1000)
                self._distances[a][b] = d
                self._distances[b][a] = d

    def makeFakesPositions(self, placesNbr):
        self.__initDistances(placesNbr)
        positions = [ (randint(0,1000),randint(0,1000)) for i in range(placesNbr)]
        for a in range(self._placesNbr):
            for b in range(a+1,self._placesNbr): 
                (xa,ya),(xb,yb) = positions[a], positions[b]
                d = int(sqrt((xa - xb)**2 + (ya - yb)**2))
                self._distances[a][b] = d
                self._distances[b][a] = d

    def printDistances(self):
        for a,sub in enumerate(self._distances):
            for b,dist in enumerate(sub):
                print (a,b,dist)

class Chromosome:
    def __init__(self, code=[]):
        self._code = code
        self._fitness = 0

    def __str__(self):
        return str(self._code) + (" => %d" % self._fitness)

    @property
    def code(self):
        return self._code

    @property
    def fitness(self):
        return self._fitness

    def updateFitness(self, geo):
        self._fitness = 0
        for i in range(len(self._code) - 1):
            self._fitness += geo.distance(self._code[i], self._code[i+1])
        return self

    def genRandomCode(self, placesNbr):
        self._code = []
        for i in range(placesNbr):
            p = randint(0,placesNbr-1)
            while p in self._code:
                p = randint(0,placesNbr-1)
            self._code.append(p)
        return self

    def mutation(self, n):
        for i in range(n):
            i1 = randint(0,len(self._code)-1)
            i2 = randint(0,len(self._code)-1)
            self._code[i1], self._code[i2] = self._code[i2], self._code[i1]

class Population:
    def __init__(self, geo, crossoverFreq=0.8, mutationFreq=0.2):
        print "[crossover=%.2f mutation=%.2f]" % (crossoverFreq,mutationFreq)
        self._geo = geo
        self._chromosomes = []
        self._chrmNbr = 0
        self._crossoverFreq = crossoverFreq
        self._mutationFreq = mutationFreq
        
    def __str__(self):
        s = "Population nbr = %d\n" % self._chrmNbr
        for ch in self._chromosomes:
            s += str(ch)+'\n'
        return s

    def __iter__(self):
        while True:
            self.doIter()
            yield self._chromosomes[0] 

    def __selectChrm(self, wheel):
        r = random()
        step = self._chrmNbr / 2
        i = step
        inf = (r < wheel[i])
        while  inf or r >= wheel[i+1]:
            if step > 1:
                step /= 2
            if inf:
                i -= step
            else:
                i += step
            inf = (r < wheel[i])
        return self._chromosomes[i]

    def __selectChrmCouple(self, wheel):
        return (self.__selectChrm(wheel), self.__selectChrm(wheel))

    def __crossover(self, chrm1, chrm2):
        if random() < self._crossoverFreq:
            parent1, parent2 = chrm1.code, chrm2.code
            point = randint(1, self._geo.placesNbr-1)
            chrm1.code = parent1[0:point] + [p for p in parent2 if p not in parent1[0:point]]
            chrm2.code = parent2[0:point] + [p for p in parent1 if p not in parent2[0:point]]
        return (chrm1, chrm2)

    def __mutation(self, chrm):
        if random() < self._mutationFreq:
            chrm.mutation(randint(1, self._geo.placesNbr / 2))
        return chrm

    def __makeWeightList(self):
        weightList = [ float(self._chromosomes[-1].fitness) / float(chrm.fitness) for chrm in self._chromosomes ]
        s = sum(weightList)
        return [ w/s for w in weightList]
    
    def __makeWheel(self):
        weightList = self.__makeWeightList()
        wheel = [0]
        for w in weightList:
            wheel.append(w + wheel[-1])
        wheel.append(2)
        return wheel


    def __selectPart(self, wheel, n, newChrmList):
        for i in range(n / 2):
            (chrm1,chrm2) = self.__selectChrmCouple(wheel)
            (chrm1,chrm2) = self.__crossover(chrm1, chrm2)
            chrm1 = self.__mutation(chrm1)
            chrm2 = self.__mutation(chrm2)
            newChrmList.append(chrm1)
            newChrmList.append(chrm2)
        return newChrmList

    def __selection(self):
        newChrmList = []
        wheel = self.__makeWheel()
        self.__selectPart(wheel, self._chrmNbr, newChrmList)
        newChrmList.append(self._chromosomes[0]);
        newChrmList.append(Chromosome().genRandomCode(self._geo.placesNbr).updateFitness(self._geo))
        return newChrmList

    def sum(self):
        sum = 0
        for c in self._chromosomes: sum += c.fitness
        return sum

    def sort(self):
        self._chromosomes.sort(cmp=lambda x,y: x.fitness - y.fitness)

    def best(self):
        return min(self._chromosomes, key=lambda x: x.fitness)

    def generate(self, chrmNbr):
        self._chrmNbr = chrmNbr
        for i in range(self._chrmNbr):
            chromosome = Chromosome()
            chromosome.genRandomCode(self._geo.placesNbr)
            chromosome.updateFitness(self._geo)
            self._chromosomes.append(chromosome)
        self.sort()

    def doIter(self):
        newChrmList = self.__selection()
        self._chromosomes = newChrmList
        self._chrmNbr = len(newChrmList)
        self.sort()


def evalPath(path, geo):
    score = 0
    for i in range(geo.placesNbr-1):
        score += geo.distance(path[i], path[i+1])
    return score

def searchPathIt(geo):
    minDist = 424242424224242
    for path in permutations(range(geo.placesNbr)):
        minDist = min(evalPath(path, geo), minDist)
        yield minDist

def main(argv):

    print "Usage : " + argv[0] + " [places-nbr population crossover-freq mutation-freq]"
    global PLACES_NBR
    geo = Geo()
    if len(argv) >= 5:
        geo.makeFakesPositions(int(argv[1]))
        pop = Population(geo, float(argv[3]), float(argv[4]))
        pop.generate(int(argv[2]))
    else:
        geo.makeFakesPositions(PLACES_NBR)
        pop = Population(geo)
        pop.generate(100)
    cnt = 0
    fittest = 424242424224242
    for c in pop:
        if c.fitness < fittest:
            fittest = c.fitness
            print "iteration %d: new fittest => %d" % (cnt, c.fitness)
        elif cnt % 100 == 0:
            print "iteration %d" % (cnt)
        cnt += 1
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
