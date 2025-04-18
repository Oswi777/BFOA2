from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
from fastaReader import fastaReader
import numpy
import copy
import csv
from datetime import datetime

if __name__ == "__main__":
    numeroDeBacterias = 4
    iteraciones = 3
    tumbo = 2
    secuencias = fastaReader().seqs
    names = fastaReader().names

    for i in range(len(secuencias)):
        secuencias[i] = list(secuencias[i])

    globalNFE = 0

    dAttr = 0.1
    wAttr = 0.002
    hRep = dAttr
    wRep = 0.001

    manager = Manager()
    numSec = len(secuencias)
    print("numSec: ", numSec)

    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))

    def poblacionInicial():
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)

    operadorBacterial = bacteria(numeroDeBacterias)
    veryBest = [None, None, None]

    start_time = time.time()
    poblacionInicial()

    for it in range(iteraciones):
        print(f"Iteración {it+1} - Tumbo...")
        operadorBacterial.tumbo(numSec, poblacion, tumbo)

        print("Cuadrando...")
        operadorBacterial.cuadra(numSec, poblacion)

        print("Generando gran lista de pares...")
        operadorBacterial.creaGranListaPares(poblacion)

        print("Evaluando BLOSUM...")
        operadorBacterial.evaluaBlosum()

        print("Calculando interacción atractiva y repulsiva...")
        operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, hRep, wRep)

        operadorBacterial.creaTablaInteraction()
        operadorBacterial.creaTablaFitness()

        globalNFE += operadorBacterial.getNFE()

        bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)

        if (veryBest[0] is None) or (bestFitness > veryBest[1]):
            veryBest[0] = bestIdx
            veryBest[1] = bestFitness
            veryBest[2] = copy.deepcopy(poblacion[bestIdx])

        # ✅ Reemplazo elitista con perturbación ligera
        operadorBacterial.replaceWorst(poblacion, veryBest[0])

        # ✅ Mutación fina dirigida al mejor individuo
        operadorBacterial.mutacion_fina_dirigida(poblacion, top_n=1)

        # Limpiar para siguiente iteración
        operadorBacterial.resetListas(numeroDeBacterias)

    print("Very Best: ", veryBest)
    total_time = time.time() - start_time
    print(f"--- Tiempo total: {total_time:.2f} segundos ---")

    # Guardar resultados en CSV
    output_file = "resultados_bfoa.csv"
    try:
        with open(output_file, 'x', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Run", "Fitness", "BlosumScore", "Interacción", "Tiempo"])
    except FileExistsError:
        pass

    fitness = veryBest[1]
    blosum = operadorBacterial.blosumScore[veryBest[0]]
    interaccion = operadorBacterial.tablaInteraction[veryBest[0]]
    run_id = datetime.now().strftime("%Y%m%d%H%M%S")

    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([run_id, fitness, blosum, interaccion, total_time])
