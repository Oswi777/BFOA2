import copy
import math
import random
import numpy
from multiprocessing import Manager, Pool
from evaluadorBlosum import evaluadorBlosum
import concurrent.futures

class bacteria():
    def __init__(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))

    def resetListas(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))

    def cuadra(self, numSec, poblacion):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])[:numSec]
            maxLen = max(len(seq) for seq in bacterTmp)
            for t in range(numSec):
                gap_count = maxLen - len(bacterTmp[t])
                if gap_count > 0:
                    bacterTmp[t].extend(["-"] * gap_count)
            poblacion[i] = tuple(bacterTmp)

    def tumbo(self, numSec, poblacion, numGaps):
        for i in range(len(poblacion)):
            bacterTmp = list(poblacion[i])
            for j in range(numGaps):
                seqnum = random.randint(0, len(bacterTmp)-1)
                pos = random.randint(0, len(bacterTmp[seqnum]))
                bacterTmp[seqnum] = bacterTmp[seqnum][:pos] + ["-"] + bacterTmp[seqnum][pos:]
            poblacion[i] = tuple(bacterTmp)

    def creaGranListaPares(self, poblacion):
        for i in range(len(poblacion)):
            pares = []
            bacterTmp = list(poblacion[i])
            for j in range(len(bacterTmp[0])):
                column = self.getColumn(bacterTmp, j)
                pares += self.obtener_pares_unicos(column)
            self.granListaPares[i] = pares

    def evaluaFila(self, fila, num):
        evaluador = evaluadorBlosum()
        score = sum(evaluador.getScore(a, b) for a, b in fila)
        self.blosumScore[num] = score

    def evaluaBlosum(self):
        with Pool() as pool:
            args = [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]
            pool.starmap(self.evaluaFila, args)

    def getColumn(self, bacterTmp, colNum):
        return [bacterTmp[i][colNum] for i in range(len(bacterTmp))]

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)

    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, blosumScore, d, w = args
        diff = (blosumScore[indexBacteria] - otherBlosumScore) ** 2.0
        self.NFE[indexBacteria] += 1
        return d * numpy.exp(w * diff)

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        with Pool() as pool:
            args = [(indexBacteria, otherScore, self.blosumScore, d, w) for otherScore in self.blosumScore]
            results = pool.map(self.compute_diff, args)

        total = sum(results)
        if atracTrue:
            self.tablaAtract[indexBacteria] = total
        else:
            self.tablaRepel[indexBacteria] = total

    def creaTablaAtract(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, True)

    def creaTablaRepel(self, poblacion, d, w):
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria, d, w, False)

    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)

    def creaTablaInteraction(self):
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        for i in range(len(self.tablaInteraction)):
            self.tablaFitness[i] = self.blosumScore[i] + self.tablaInteraction[i]

    def getNFE(self):
        return sum(self.NFE)

    def obtieneBest(self, globalNFE):
        bestIdx = max(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        print("-------------------   Best: ", bestIdx, " Fitness: ", self.tablaFitness[bestIdx],
              "BlosumScore ", self.blosumScore[bestIdx], "Interaction: ", self.tablaInteraction[bestIdx],
              "NFE: ", globalNFE)
        return bestIdx, self.tablaFitness[bestIdx]

    # ✅ MEJORA 1: Reemplazo elitista con perturbación ligera
    def replaceWorst(self, poblacion, best):
        worst = min(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        mejor_bacteria = list(copy.deepcopy(poblacion[best]))
        numSec = len(mejor_bacteria)
        seq_idx = random.randint(0, numSec - 1)
        pos = random.randint(0, len(mejor_bacteria[seq_idx]))
        mejor_bacteria[seq_idx] = mejor_bacteria[seq_idx][:pos] + ["-"] + mejor_bacteria[seq_idx][pos:]
        poblacion[worst] = tuple(mejor_bacteria)

    # ✅ MEJORA 2: Mutación fina dirigida
    def mutacion_fina_dirigida(self, poblacion, top_n=1):
        evaluador = evaluadorBlosum()
        for i in range(top_n):
            bacter = list(poblacion[i])
            worst_cols = self.identificar_peores_columnas(bacter, evaluador)
            for col in worst_cols:
                for s in range(len(bacter)):
                    if len(bacter[s]) > col:
                        bacter[s].insert(col, "-")
            poblacion[i] = tuple(bacter)

    def identificar_peores_columnas(self, bacter, evaluador):
        col_scores = []
        for col in range(len(bacter[0])):
            columna = [bacter[fila][col] for fila in range(len(bacter))]
            pares = self.obtener_pares_unicos(columna)
            score = sum(evaluador.getScore(a, b) for (a, b) in pares)
            col_scores.append((col, score))
        col_scores.sort(key=lambda x: x[1])
        return [idx for idx, _ in col_scores[:2]]  # Top 2 peores columnas
