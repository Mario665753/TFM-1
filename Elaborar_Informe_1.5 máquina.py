#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:00:54 2024

@author: alumno
"""
import re
import sys
import pandas as pd
import argparse
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

#Declarando el parser y cogiendo argumentos

parser = argparse.ArgumentParser()
parser.add_argument('--f', action="store", dest='file', type=str, nargs='*',
                    help="Introduce el nombre de los archivos tsv para los que desea generar el informe entre comillas. Ejemplo: \"Archivo 1.txt\", \"Archivo 2.tsv\", etc.")
parser.add_argument('--g', action="store", dest='gold_standard', type=str, nargs=1,
                    help="Introduzca el nombre del archivo gold_standard en formato tsv")
args = parser.parse_args()

#Cogiendo el input
d={}
print('Leyendo archivos...')
i = 0
for x in args.file:
    print(i, '/', len(args.file))
    a=pd.read_csv(x, sep='\t', dtype=str) #Leo el df
    d.update({x:a}) #Lo meto en el diccionario con su nombre como clave
    d[x]['Score'] = d[x]['Score'].astype(float) #Hago que los scores sean float
    i = i + 1

print('¿Desea formatear los archivos para poder trabajar con ellos?', '\n', '1) Sí', '\n', '2) No')
cond = int(input())
if cond == 1:
    print('Poniendo los archivos en formato para trabajar con ellos...')
    #Creo una nueva columna para analizar los genes del tirón con el operador 'in'
    for x in args.file:
        print('\t', 'Formateando archivo ', x, '...')
        gene1_2=[]
        gene2_1=[]
        for i in range(len(d[x])):
            if i % 30000 == 0:
                print(i, '/', len(d[x]))
            new_gene1_2 = d[x]['Gene1'][i] + '-' + d[x]['Gene2'][i]
            gene1_2.append(new_gene1_2)
            new_gene2_1 = d[x]['Gene2'][i] + '-' + d[x]['Gene1'][i]
            gene2_1.append(new_gene2_1)
        d[x]['Gene1-2'] = gene1_2
        d[x]['Gene2-1'] = gene2_1
        a = re.split('\.', x) #Estaría guay hacer un diccionario para esto, la verdad
        name = a[0] + '_processed.txt'
        print('\t', '¡Completado! Guardando archivo procesado a ', name)
        d[x].to_csv(name, index=False, sep='\t')
print('¿Quiere procesar el gold standard?', '\n', '1) Sí', '\n', '2) No')
cond = int(input())
if cond == 1:
    #Ahora formateo y quito redundancias del gold_standard
    gold_standard = pd.read_csv(args.gold_standard[0], sep='\t')
    print('\t', 'Formateando el gold_standard...')
    gene1_2=[]
    gene2_1=[]
    for i in range(len(gold_standard)):
        if i % 30000 == 0:
            print(i, '/', len(d[x]))
        new_gene1_2 = gold_standard['Gene1'][i] + '-' + gold_standard['Gene2'][i]
        gene1_2.append(new_gene1_2)
        new_gene2_1 = gold_standard['Gene2'][i] + '-' + gold_standard['Gene1'][i]
        gene2_1.append(new_gene2_1)
    gold_standard['Gene1-2'] = gene1_2
    gold_standard['Gene2-1'] = gene2_1
    b = re.split('\.', args.gold_standard[0])
    name = b[0] + '_processed.txt'
    gold_standard.to_csv(name, index=False, sep='\t')
    print('\t', '¡Completado! Guardando gold_standard procesado a ', name)
'''
    print('\t', 'Eliminando redundancias...')
    for i in range(len(gold_standard)):
        if i % 30000 == 0:
            print(i, '/', len(d[x]))
        e1 = sum(gold_standard['Gene1-2'][i+1:len(gold_standard)] == gold_standard['Gene1-2'][i])
        e2 = sum(gold_standard['Gene2-1'][i+1:len(gold_standard)] == gold_standard['Gene1-2'][i])
        e = e1 + e2
        if e>0:
            gold_standard = gold_standard.drop(i)
    gold_standard.index = range(len(gold_standard))
    b = re.split('\.', args.gold_standard[0])
    name = b[0] + '_processed.txt'
    gold_standard.to_csv(name, index=False, sep='\t')
    print('\t', '¡Completado! Guardando gold_standard procesado a ', name)
'''
print('EL PROGRAMA SOLO TRABAJA CON LOS ARCHIVOS DADOS. PARA TRABAJAR CON LOS ARCHIVOS PROCESADOS DEBE INTRODUCIRLOS')

while 1==1:
    #Seleccionando opción
    print('\n', 'Este script tiene como finalidad realizar las siguientes funciones:',
          '\n', '0) Salir',
          '\n', '1) Hacer filtrado de genes redundantes quedándose con los de mayor score, por archivo',
          '\n', '2) Cálculo de coeficientes de correlación de Spearman entre métodos y generación de gráficas',
          '\n', '3) Obtener valores de precisión y recall y calcular curva precisión/recall',
          '\n', '4) Obtener las interacciones totales no redundantes de todos los métodos',
          '\n', '5) Obtener las interacciones exclusivas no redundantes de cada método y comunes entre métodos',
          '\n', '6) Ver cuántos del n-top de los comunes son gold-standard',
          '\n', '7) Hacer filtrado estadístico de los datos más relevantes',
          '\n', 'Inserte la opción que quiera realizar:', '\n')
    #4 Y 5 SE PODRÍAN HACER A LA VEZ
    case = int(input())
    if case == 0:
        sys.exit("Vuelva cuando quiera")
###################################################  ELIMINACIÓN DE REDUNDANCIAS  #################################################
    elif case == 1:
        for x in args.file:
            pairs=dict()
            old_score=0
            i = 0
            with open(x, "r") as file:
                print('Procesando el archivo ', x)
                for line in file:
                    g1,g2,score=line.strip().split("\t")
                    first,second=sorted([g1,g2])
                    if(first in pairs.keys()):
                        if(second in pairs[first].keys()):
                            old_score=pairs[first][second]
                            if (score>old_score):
                                pairs[first][second]=score
                        else:
                            pairs[first][second]=score
                    else:
                        pairs[first]=dict()
                        pairs[first][second]=score
                    i = i+1
                    print('Lines processed: ', i)
            name = 'no_redundancies/' + x
            with open(name, 'w') as file:
                for f in pairs.keys():
                    for s in pairs[f].keys():
                        file.write(f + "\t" + s + "\t" + pairs[f][s] + '\n')
        cond = 0
        print('\n', '¿Desea procesar el gold standard?', '\n', '1) Sí', '\n', '2) No')
        cond = int(input())
        if cond == 1:
            with open(args.gold_standard[0], "r") as file:
                print('Procesando el archivo ', args.gold_standard[0])
                pairs=dict()
                i = 0
                for line in file:
                    g1,g2=line.strip().split("\t")
                    first,second=sorted([g1,g2])
                    if(first in pairs.keys()):
                        if(second in pairs[first].keys()) == False:
                            pairs[first][second]=0
                    else:
                        pairs[first] = dict()
                        pairs[first][second] = 0
                    i = i+1
                    print('Lines processed: ', i)
        name = 'no_redundancies/' + args.gold_standard[0]
        with open(name, 'w') as file:
            for f in pairs.keys():
                for s in pairs[f].keys():
                    file.write(f + "\t" + s + '\n')
    elif case == 2:
#########################################################  SPEARMAN  #######################################################################
        e=0
        while e==0:
            print("Se va a calcular el coeficiente de Spearman y las gráficas con los archivos que ha indicado al principio. Si quiere cambiar de archivos, salga del programa y vuelva a iniciarlo indicando dichos archivos")
            print("\n", "1) Continuar")
            print("\n", "2) Salir")
            e = int(input())
        if e==2: #Salir
            break
        '''
        if e==1:
            e = 0
            print('Puede elegir entre obtener el coeficiente de Spearman para todos los pares del método o solo para los pares del gold standard: ')
            print('\n', '1) Todos los pares')
            print('\n', '2) Solo gold_standard')
            e = int(input())
            i = 0
            #Creo un diccionario con las keys que sean el gold standard
            with open(args.gold_standard[0], 'r') as file:
                gs = {}
                for line in file:
                    g1,g2=line.strip().split('\t')
                    if (first in gs.keys()):
                        if (second in gs[first].keys()) == False:
                            gs[first][second]=0
                    else:
                        gs[first] = dict()
                        gs[first][second]=0
            for i in range(len(args.file)):
                for j in range(i+1, len(args.file)): #Esto será para hacer comparaciones combinatorias
                    with open(args.file[i], 'r') as file1:
                        with open(args.file[i+1], 'r') as file2:
                            for line1 in file1:
                                m1g1,m1g2,m1score=line.strip().split("\t")
'''                
        if e==1: 
            e = 0
            print('Puede elegir entre obtener el coeficiente de Spearman para todos los pares del método o solo para los pares del gold standard: ')
            print('\n', '1) Todos los pares')
            print('\n', '2) Solo gold_standard')
            e = int(input())
            a={}
            i=0
            z=pd.DataFrame(columns=['Method comparison', 'Spearman coefficient', 'p-value'])
            gold_standard = pd.read_csv(args.gold_standard[0], sep='\t')
            z=pd.DataFrame(columns=['Method comparison', 'Spearman coefficient', 'p-value'])
            for x in args.file:#Recorro los archivos
                a[i]=d[x] #Paso los índices a números, será más fácil
                i=i+1
            if e==1:
                for i in range(len(a)-1):
                    for j in range(i+1, len(a)): #De esta forma hace comparaciones combinatorias
                        GS_method_1 = pd.DataFrame() #Gold_standard que ha capturado el primer método CAMBIA ESTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
                        GS_method_2 = pd.DataFrame()# Gold_standard que ha capturado el segundo método CAMBIA ESTOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
                        name_1 = re.split('\.', args.file[i])[0] #Nombre del primer método
                        name_2 = re.split('\.', args.file[j])[0] #Nombre del segundo método
                        print('Calculando Spearman entre los métodos ', name_1, ' y ', name_2, '...')
                        print('\t', 'Obteniendo pares comunes a ambos métodos...')
                        #Ahora voy a crear un solo dataframe con los comunes
                        #Primero renombro las columnas 'Score' para que tengan el nombre del método de donde vienen
                        name_1_score = name_1 + '_Score'
                        name_2_score = name_2 + '_Score'
                        a[i].rename(columns={'Score':name_1_score}, inplace = True)
                        a[j].rename(columns={'Score':name_2_score}, inplace = True)
                        #Ahora los ordeno y les pongo el orden, incluyendo en el nombre de la columna donde está el orden el nombre del método del que proviene
                        name_1_rank = 'Rank_' + name_1
                        name_2_rank = 'Rank_' + name_2
                        a[i] = a[i].sort_values(by=name_1_score, ascending=False)
                        a[i][name_1_rank] = range(len(a[i]))
                        a[j] = a[j].sort_values(by=name_2_score, ascending=False)
                        a[j][name_2_rank] = range(len(a[j]))
                        a[i] = a[i].sort_values(by='Gene1-2')
                        a[i] = a[i].sort_values(by='Gene1-2') #Ordeno los valores por su nombre, de esa forma el orden relativo de los pares es el mismo
                        common = a[i][a[i]['Gene1-2'].isin(a[j]['Gene1-2'])][['Gene1-2', name_1_score, name_1_rank]] #Hago intersección de los que son comunes
                        d = a[j][a[j]['Gene1-2'].isin(a[i]['Gene1-2'])][['Gene1-2', name_2_score, name_2_rank]] #Cojo la información del score y el rank del segundo método
                        d = d.reset_index(drop = True)#Reseteo los índices
                        common = common.reset_index(drop = True)#Reseteo los índices
                        common[[name_2_score, name_2_rank]] = d[[name_2_score, name_2_rank]]#Ahora todo cabe perfectamente
                ##################################################IMPRIMO LOS PARES COMUNES A LOS DOS MÉTODOS###########################
                    name = 'Spearman/Common/common_' + name_1 + '-' + name_2 + '_Spearman.txt'
                    print('\t', '¡Completado! Guardando los pares del gold_standard comunes a ambos métodos en ', name)
                    #name = 'DeepRig-PIDC_GS_Spearman.txt'
                    common.to_csv(name, sep='\t', index=False) #Imprimo los resultados en el archivo
                ##################################################################################################################################
                    print('\t', 'Generando gráfica ', name_1, ' vs ', name_2) #Genero la gráfica, los 1000 mejores del primer método contra el segundo
                    common = common.sort_values(by=name_1_rank)
                    plt.scatter(common[name_1_rank][0:1000], common[name_2_rank][0:1000])
                    plt.xlabel(name_1)
                    plt.ylabel(name_2)
                    title = 'Pair ranking in ' + name_1 + ' vs ' + name_2
                    plt.title(title)
                    name = 'Spearman/graphs/Common/common_' + name_1 + '_vs_' + name_2 + '.pdf'
                    plt.savefig(name, format='pdf')
                    plt.clf()
                    print('\t', '¡Completado! Guardando a ', name)
                    print('\t', 'Generando gráfica ', name_2, ' vs ', name_1) #Genero la gráfica, los 1000 mejores del segundo método contra el primero
                    common = common.sort_values(by=name_2_rank)
                    plt.scatter(common[name_2_rank][0:1000], common[name_1_rank][0:1000])
                    plt.xlabel(name_2)
                    plt.ylabel(name_1)
                    title = 'Pair ranking of gold standard in ' + name_2 + ' vs ' + name_1
                    plt.title(title)
                    name = 'Spearman/graphs/Common/common_' + name_2 + '_vs_' + name_1 + '.pdf'
                    plt.savefig(name, format='pdf')
                    plt.clf()
                    print('\t', '¡Completado! Guardando a ', name)
                    #Genero informe de coeficiente y p-valor asociado
                    print('\t', 'Generando informe de coeficiente y p-valor asociado para ', name_1, name_2)
                    coef, p = spearmanr(common[name_1_score], common[name_2_score]) #Calculo coeficiente y p-valor
                    mc= name_1 + '-' + name_2 #Hago nombre conjunto de métodos
                    z1={'Method comparison':[mc], 'Spearman coefficient':[coef], 'p-value':[p]} #Meto la información en df
                    z1=pd.DataFrame.from_dict(z1) 
                    z=pd.concat([z,z1]) #Lo concateno con lo que tenía
                    print('¡Completado! Pasando al siguiente método...')
                    ###########################################
                    plt.scatter(common[name_2_rank], common[name_1_rank])
                    plt.xlabel(name_2)
                    plt.ylabel(name_1)
                    title = 'Pair ranking in ' + name_2 + ' vs pair ranking in ' + name_1
                    #plt.xlabel('DeepRig')
                    #plt.ylabel('PIDC')
                    #title = 'Pair ranking in DeepRirg vs Pair ranking in PIDC'
                    plt.title(title)
                    name = 'Spearman/graphs/Common/common_' + name_2 + ' _vs_ ' + name_1 + '.pdf'
                    plt.savefig(name, format='pdf')
                    plt.clf()
                    #########################################################Todo repetición para que la última figura salga
                print('Guardando informe de coeficientes en Spearman/Informe_coeficientes_common.txt')
                z.to_csv('Spearman/Informe_coeficientes_common.txt', sep='\t', index=False)
            else:###############################################  SPEARMAN CON GOLD STANDARD  ######################################
                for i in range(len(a)-1):
                    for j in range(i+1, len(a)): #De esta forma hace comparaciones combinatorias
                        GS_method_1 = pd.DataFrame() #Gold_standard que ha capturado el primer método
                        GS_method_2 = pd.DataFrame()# Gold_standard que ha capturado el segundo método
                        name_1 = re.split('\.', args.file[i])[0] #Nombre del primer método
                        name_2 = re.split('\.', args.file[j])[0] #Nombre del segundo método
                        print('Calculando Spearman entre los métodos ', name_1, ' y ', name_2, '...')
                        print('\t', 'Obteniendo pares del gold_standard de ambos métodos...')
                        for k in range(len(gold_standard)):
                            if k % 30000 == 0:
                                print(k, '/', len(gold_standard))
                                #print(gold_standard)
                            e1 = sum(gold_standard['Gene1-2'][k] == a[i]['Gene1-2'])#Miro si está en el método uno como G1-G2
                            e2 = sum(gold_standard['Gene1-2'][k] == a[i]['Gene2-1'])#Miro si está en el método dos como G2-G1
                            if e1:
                                GS_method_1 = pd.concat([GS_method_1, a[i][gold_standard['Gene1-2'][k] == a[i]['Gene1-2']][['Gene1-2', 'Score']]]) #Lo cojo
                            elif e2:
                                b = a[i][gold_standard['Gene1-2'][k] == a[i]['Gene2-1']][['Gene2-1', 'Score']]#Cojo el score y el nombre
                                b = b.rename(columns={'Gene2-1':'Gene1-2'})#Cambio el nombre de la columna
                                GS_method_1 = pd.concat([GS_method_1, b]) #Lo incluyo
                            e1 = sum(gold_standard['Gene1-2'][k] == a[j]['Gene1-2'])#Miro si está en el método uno como G1-G2
                            e2 = sum(gold_standard['Gene1-2'][k] == a[j]['Gene2-1'])#Miro si está en el método dos como G2-G1
                            if e1:
                                GS_method_2 = pd.concat([GS_method_2, a[j][gold_standard['Gene1-2'][k] == a[j]['Gene1-2']][['Gene1-2', 'Score']]]) #Lo cojo
                            elif e2:
                                c=a[j][gold_standard['Gene1-2'][k] == a[j]['Gene2-1']][['Gene2-1', 'Score']] #Cojo el Score y el nombre
                                c = c.rename(columns={'Gene2-1':'Gene1-2'})#Cambio el nombre de la columna
                                GS_method_2 = pd.concat([GS_method_2, c]) #Lo incluyo
                            GS_method_1.index = range(len(GS_method_1))
                            GS_method_2.index = range(len(GS_method_2))
                            gs_method_1 = 'Spearman/GS/' + 'GS_' + name_1 + '.txt' 
                            gs_method_2 = 'Spearman/GS/' + 'GS_' + name_2 + '.txt'
                        print('\t', '¡Completado! Guardando los pares del gold_standard encontrados de ambos métodos en ', gs_method_1, ' y en ', gs_method_2)
                ###############################IMPRIMO LOS GOLD_STANDARD EN CADA UNO DE LOS MÉTODOS#############################
                        GS_method_1.to_csv(gs_method_1, sep='\t', index=False) #Imprimo resultados
                        GS_method_2.to_csv(gs_method_2, sep='\t', index=False)
                ##########################################################################################################
                        print('\t', 'Obteniendo pares del gold_standard que son comunes a ambos métodos...')
                        #Ahora voy a crear un solo dataframe con los comunes
                        #Primero renombro las columnas 'Score' para que tengan el nombre del método de donde vienen
                        name_1_score = name_1 + '_Score'
                        name_2_score = name_2 + '_Score'
                        GS_method_1.rename(columns={'Score':name_1_score}, inplace = True)
                        GS_method_2.rename(columns={'Score':name_2_score}, inplace = True)
                        #Ahora los ordeno y les pongo el orden, incluyendo en el nombre de la columna donde está el orden el nombre del método del que proviene
                        name_1_rank = 'Rank_' + name_1
                        name_2_rank = 'Rank_' + name_2
                        GS_method_1 = GS_method_1.sort_values(by=name_1_score, ascending=False)
                        GS_method_1[name_1_rank] = range(len(GS_method_1))
                        GS_method_2 = GS_method_2.sort_values(by=name_2_score, ascending=False)
                        GS_method_2[name_2_rank] = range(len(GS_method_2))
                        GS_method_1 = GS_method_1.sort_values(by='Gene1-2')
                        GS_method_2 = GS_method_2.sort_values(by='Gene1-2') #Ordeno los valores por su nombre, de esa forma el orden relativo de los pares es el mismo
                        GS_common = GS_method_1[GS_method_1['Gene1-2'].isin(GS_method_2['Gene1-2'])][['Gene1-2', name_1_score, name_1_rank]] #Hago intersección de los que son comunes
                        d = GS_method_2[GS_method_2['Gene1-2'].isin(GS_method_1['Gene1-2'])][['Gene1-2', name_2_score, name_2_rank]] #Cojo la información del score y el rank del segundo método
                        d = d.reset_index(drop = True)#Reseteo los índices
                        GS_common = GS_common.reset_index(drop = True)#Reseteo los índices
                        GS_common[[name_2_score, name_2_rank]] = d[[name_2_score, name_2_rank]]#Ahora todo cabe perfectamente
                ##################################################IMPRIMO LOS GOLD_STANDARD COMUNES A LOS DOS MÉTODOS###########################
                        name = 'Spearman/GS/GS_' + name_1 + '-' + name_2 + '_Spearman.txt'
                        print('\t', '¡Completado! Guardando los pares del gold_standard comunes a ambos métodos en ', name)
                        #name = 'DeepRig-PIDC_GS_Spearman.txt'
                        GS_common.to_csv(name, sep='\t', index=False) #Imprimo los resultados en el archivo
                ##################################################################################################################################
                        print('\t', 'Generando gráfica ', name_1, ' vs ', name_2) #Genero la gráfica, los 1000 mejores del primer método contra el segundo
                        GS_common = GS_common.sort_values(by=name_1_rank)
                        plt.scatter(GS_common[name_1_rank][0:1000], GS_common[name_2_rank][0:1000])
                        plt.xlabel(name_1)
                        plt.ylabel(name_2)
                        title = 'Pair ranking of gold standard in ' + name_1 + ' vs ' + name_2
                        plt.title(title)
                        name = 'Spearman/graphs/GS/GS_' + name_1 + '_vs_' + name_2 + '.pdf'
                        plt.savefig(name, format='pdf')
                        plt.clf()
                        print('\t', '¡Completado! Guardando a ', name)
                        print('\t', 'Generando gráfica ', name_2, ' vs ', name_1) #Genero la gráfica, los 1000 mejores del segundo método contra el primero
                        GS_common = GS_common.sort_values(by=name_2_rank)
                        plt.scatter(GS_common[name_2_rank][0:1000], GS_common[name_1_rank][0:1000])
                        plt.xlabel(name_2)
                        plt.ylabel(name_1)
                        title = 'Pair ranking of gold standard in ' + name_2 + ' vs ' + name_1
                        plt.title(title)
                        name = 'Spearman/graphs/GS/GS_' + name_2 + '_vs_' + name_1 + '.pdf'
                        plt.savefig(name, format='pdf')
                        plt.clf()
                        print('\t', '¡Completado! Guardando a ', name)
                        #Genero informe de coeficiente y p-valor asociado
                        print('\t', 'Generando informe de coeficiente y p-valor asociado para ', name_1, name_2)
                        coef, p = spearmanr(GS_common[name_1_score], GS_common[name_2_score]) #Calculo coeficiente y p-valor
                        mc= name_1 + '-' + name_2 #Hago nombre conjunto de métodos
                        z1={'Method comparison':[mc], 'Spearman coefficient':[coef], 'p-value':[p]} #Meto la información en df
                        z1=pd.DataFrame.from_dict(z1) 
                        z=pd.concat([z,z1]) #Lo concateno con lo que tenía
                        print('¡Completado! Pasando al siguiente método...')
                        ###########################################
                        plt.scatter(GS_common[name_1_rank], GS_common[name_2_rank])
                        plt.xlabel(name_2)
                        plt.ylabel(name_1)
                        title = 'Pair ranking in ' + name_2 + ' vs pair ranking in ' + name_1
                        #plt.xlabel('DeepRig')
                        #plt.ylabel('PIDC')
                        #title = 'Pair ranking in DeepRirg vs Pair ranking in PIDC'
                        plt.title(title)
                        #name = 'Spearman/graphs/' + name_1[0] + ' _vs_ ' + name_2[0] + '.pdf'
                        plt.savefig(name, format='pdf')
                        plt.clf()
                        #########################################################Todo repetición para que la última figura salga
                print('Guardando informe de coeficientes en Spearman/Informe_coeficientes_GS.txt')
                z.to_csv('Spearman/Informe_coeficientes_GS.txt', sep='\t', index=False)
    elif case == 3:
 ################################################################  VALORES DE PRECISIÓN Y RECALL Y CURVAS ROC  ####################################################################
        print('Obteniendo curvas ROC...')
        z=pd.DataFrame(columns=['Threshold <=', 'true_positives', 'false_positives', 'true_negatives', 'false_negatives', 'Recall', 'Precision'])
        gold_standard = pd.read_csv(args.gold_standard[0], sep='\t')
        for x in args.file:
            print('\t', 'Procesando el archivo ', x)
            print('\n', 'Para obtener la curva ROC es necesario establecer de forma iterativa qué interacciones consideramos como significativas y cuáles no lo son. Mientras menos pares se requieran para calcular un nuevo punto de la curva, más tiempo tardará. Introduzca cada cuántos pares quiere calcular un punto para la curva ROC: ')
            jump = int(input())
            print('\t', 'Obteniendo puntos de la curva ROC...')
            z=pd.DataFrame(columns=['Threshold <=', 'true_positives', 'false_positives', 'true_negatives', 'false_negatives', 'Recall', 'Precision'])
            for i in range(0, len(d[x]), jump): #Este será el threshold que se vaya moviendo.
                if i % 30000 == 0:
                    print(i, '/', len(d[x]))
                neg = len(d[x]['Gene1-2']) - i - 1 #Esto es el número de instancias que hemos clasificado como positivas
                pos = i + 1#Esto es el número de instancias que hemos clasificado como negativas
                false_neg = len(set(d[x]['Gene1-2'][i+1:len(d[x])]) & set(gold_standard['Gene1-2'])) + len(set(d[x]['Gene2-1'][i+1:len(d[x])]) & set(gold_standard['Gene1-2'])) #Auellos pares que estén presentes en el gold standard ya sea como g1-2 o g2-1
                true_pos = len(set(d[x]['Gene1-2'][0:i+1]) & set(gold_standard['Gene1-2'])) + len(set(d[x]['Gene2-1'][0:i+1]) & set(gold_standard['Gene1-2'])) #Aquellas que coinciden con el gold_standard y han sido clasificadas como negativas
                false_pos = pos-true_pos #El número de instancias clasificadas como positivas menos las que son correctas
                true_neg = neg - false_neg
                Precision = true_pos/pos
                recall = true_pos/(true_pos + false_neg)
                z1={'Threshold <=':[i],'true_positives':[true_pos], 'false_positives':[false_pos], 'true_negatives':[true_neg], 'false_negatives':[false_neg],
                    'Recall':[recall], 'Precision':[Precision]}
                z1=pd.DataFrame.from_dict(z1)
                z=pd.concat([z,z1])
            print('¡Completado!')
            print('Generando gráfica')
            plt.plot(z['Recall'], z['Precision'])
            plt.xlabel('Recall')
            plt.ylabel('Precision')
            title = 'ROC curve'
            plt.title(title)
            a = re.split('\.', x)
            name = 'ROC/' + a[0] + '_ROC.pdf'
            print('Exportando gráfica a ', name)
            plt.savefig(name, format='pdf')
            plt.clf()
            print('¡Completado!')
            name = 'ROC/' + a[0] + '_ROC points.txt'
            print('Exportando puntos de la curva a ', name)
            z.to_csv(name, sep='\t', index=False)
            print('¡Completado!')
    elif case == 4:
############################################################################  INTERACCIONES TOTALES NO REDUNDANTES  ##############################################################
        print('Lo que se va a hacer es concatenar los 3 archivos que ha dado. Después debe salir del programa, volver a iniciarlo y pasarlo por el primer apartado para eliminar redundancias.')
        all_red = pd.DataFrame(columns=['Gene1-2', 'Gene2-1', 'Score' ])
        for x in args.file:
            print('\t', 'Generando columna \'Gene2-1\' para el archivo ', x, '...')
            d[x]['Gene2-1']=""
            for i in range(len(d[x])):
                if i % 30000 == 0:
                    print(i, '/', len(d[x]))
                a = re.split('-', d[x]['Gene1-2'][i])
                d[x]['Gene2-1'][i] = a[1] + '-' + a[0]
            print('\t', 'Concatenando dataframes...')
            all_red = pd.concat([all_red, d[x]])
        print('¡Completado! Exportando todas las interacciones a Todas_interacciones.tsv')
        all_red.to_csv('Todas_interacciones.tsv', sep='\t', index=False)
    elif case == 5:
##################################################################################  INTERACCIONES COMUNES Y EXCLUSIVAS  ############################################################
        ################################################  INTERACCIONES COMUNES  ##############################################
        i = 0
        a = {}
        for x in args.file:#Recorro los archivos
                a[i]=d[x] #Paso los índices a números, será más fácil
                i=i+1
        print('Obteniendo interacciones comunes...')
        common = set(a[0]['Gene1-2']).intersection(set(a[1]['Gene1-2'])) #Hago intersección entre los pares G1-G2 del primer y el segundo método
        common = common | set(a[0]['Gene1-2']).intersection(set(a[1]['Gene2-1'])) #A esto añado también los G2-G1
        for i in range(2,len(a)): #Con esto obtengo los pares comunes a todos los métodos
            print('\t', 'Sets procesados: ', i, '/', len(a))
            common = common.intersection(set(a[i]['Gene1-2'])) #Repito la intersección para todos los métodos tanto G1-G2 como G2-G1
            common = common | common.intersection(set(a[i]['Gene2-1']))
        common = pd.DataFrame({'Gene1-2':list(common)}) #Lo convierto a DataFrame
        common = common.sort_values(by='Gene1-2').reset_index(drop=True)# Y los ordeno por nombre, así el orden relativo de los comunes será el mismo, también en el resto de métodos si están ordenados
        print('Obteniendo Scores...')
        for i in range(len(a)):#Con esto ahora voy a coger los scores
            print('Sets procesados: ', i, '/', len(a))
            name = re.split('\.', args.file[i])[0] #Cojo el nombre y lo convierto en el nombre de la columna donde va a ir el score del método que toque
            name_score = name + '_Score'
            a[i] = a[i].sort_values(by='Gene1-2').reset_index(drop = True) #Ordeno los pares del método para que el orden relativo de los pares sea el mismo
            common[name_score] = range(len(common))#Esto es para que no me de error en la siguiente línea
            if sum(common['Gene1-2'].isin(a[i]['Gene1-2'])) > 0:
                print('He entradooooooooooooooooooooooooooooooooo')
                common[name_score].loc[common['Gene1-2'].isin(a[i]['Gene1-2'])] = list(a[i]['Score'][a[i]['Gene1-2'].isin(common['Gene1-2'])]) #Localizo los Scores de G1-G2 en el método y lo añado a los comunes
            #Hay errores en la línea de arriba en la máquina, míralos
            if sum(common['Gene1-2'].isin(a[i]['Gene2-1'])) > 0:
                print(len(common[name_score].loc[common['Gene1-2'].isin(a[i]['Gene2-1'])]))
                print(len(list(a[i]['Score'][a[i]['Gene2-1'].isin(common['Gene1-2'])])))
                print(len(a[i]['Gene2-1'])!=len(set(a[i]['Gene2-1'])))
                print(set(a[i]['Gene2-1'][a[i]['Gene2-1'].isin(common['Gene1-2'])]) - set(common['Gene1-2'].loc[common['Gene1-2'].isin(a[i]['Gene2-1'])]))
                common[name_score].loc[common['Gene1-2'].isin(a[i]['Gene2-1'])] = list(a[i]['Score'][a[i]['Gene2-1'].isin(common['Gene1-2'])])#Hago lo mismo con G2-G1
        #DA PROBLEMAS EN LA LÍNEA ANTERIOR TIENES QUE MIRARLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        print('Exportando interacciones comunes no redundantes a Interacciones comunes no redundantes.txt')
        common.to_csv('Interacciones comunes no redundantes.txt', sep='\t', index = False)
        '''
        ##################################################  INTERACCIONES EXCLUSIVAS  ################################################
        print('Obteniendo interacciones exclusivas...')
        for i in range(len(a)):
            print('\t', 'Sets procesados: ', i, '/', len(a))
            exclusive = set(a[i]['Gene1-2']) #Establezco el set de exclusivos que va a partir de todos los pares completos del método
            for j in range(len(a)): #Esto es para ir comparando el método que sea con los siguientes y si coincide consigo mismo que no haga la comparación
                if j != i:
                    exclusive = exclusive - set(a[j]['Gene1-2']) #Le quita al método los pares G1-G2 del otro
                    exclusive = exclusive - set(a[j]['Gene2-1']) #Hace lo mismo con los G2-G1
            exclusive = pd.DataFrame({'Gene1-2':list(exclusive)}) #Lo transformo en DataFrame
            exclusive = exclusive.sort_values(by='Gene1-2').reset_index(drop=True) #Los ordeno por nombre
            print('Obteniendo Scores...')
            name = re.split('.', args.file[i])[0] #Obtengo el nombre del método para convertirlo en la columna 'Score'
            name_score = name + '_Score'
            a[i] = a[i].sort_values(by='Gene1-2').reset_index(drop = True)
            exclusive[name_score] = range(len(exclusive)) #Esto es para que no me de error en la siguiente línea
            exclusive[name_score].loc[exclusive['Gene1-2'].isin(a[i]['Gene1-2'])] = list(a[i]['Score'][a[i]['Gene1-2'].isin(exclusive['Gene1-2'])]) #Localizo los Scores de los genes exclusivos y los añado al dataframe
            print('Exportando interacciones exclusivas a Interacciones exclusivas/' + name + '.txt')
            exclusive.to_csv('Interacciones exclusivas/' + name + '.txt', sep='\t', index = False)
        '''
    elif case == 6:
        comunes = pd.read_csv('Interacciones comunes no redundantes.txt', sep='\t')
        comunes = comunes.sort_values(by='results_aracne_2_processed_Score', ascending=False)
        comunes['Rank_DeepSEM']=range(len(comunes))
        comunes = comunes.sort_values(by='results_SCODE_processed_Score', ascending=False)
        comunes['Rank_DeepRig']=range(len(comunes))
        comunes = comunes.sort_values(by='results_PIDC_processed_Score', ascending=False)
        comunes['Rank_PORTIA']=range(len(comunes))
        comunes.to_csv('Interacciones comunes IA ranked.txt', index=False, sep='\t')
        with open('Interacciones comunes IA ranked.txt', "r") as file:
            with open('Interacciones comunes IA added rank.txt', 'w') as file2:
                for line in file:
                    g12, score_1, score_2, score_3, rank_DR, rank_SCODE, rank_PIDC = line.strip().split("\t")
                    if g12 == 'Gene1-2':
                        file2.write('Gene1' + '\t' + 'Gene2' + '\t' + 'Gene1-2' + '\t' + score_1+ '\t'+ score_2+ '\t'+ score_3+ '\t'+ rank_DR+ '\t'+ rank_SCODE+ '\t'+ rank_PIDC+ '\t'+ 'Added_rank'+ '\n')
                        continue
                    add_rank = rank_DR + rank_SCODE + rank_PIDC
                    #print(add_rank)
                    #print(g12.split('-', maxsplit=1))
                    g1, g2=g12.split('-', maxsplit=1)
                    file2.write(g1+ '\t'+ g2+ '\t'+ score_1+ '\t'+ g12 + '\t' + score_2+ '\t'+ score_3+ '\t'+ rank_DR+ '\t'+ rank_SCODE+ '\t'+ rank_PIDC+ '\t'+ add_rank + '\n')
        comunes_rank = pd.read_csv('Interacciones comunes IA added rank.txt', sep='\t')
        comunes_rank = comunes_rank.sort_values(by='Added_rank')
        comunes_rank.to_csv('Interacciones comunes IA added rank.txt', sep='\t', index=False)
        '''
        print('Introduzca n-top: ')
        a = int(input())
        with open(args.gold_standard[0], "r") as file1:
            with open('Interacciones comunes added rank.txt', 'r') as file2:
                with open('n-top_common.txt', 'w') as file3:
                    i = 0
                    for line2 in file2:
                        if i == a:
                            break
                        print(line2.strip().split('\t'))
                        g1,g2,g12,score_1,score_2,score_3,rank_DR,rank_SCODE,rank_PIDC,added_rank=line2.strip().split("\t")
                        g1,g2=sorted([g1,g2])
                        for line1 in file1:
                            G1,G2,s,d = line1.strip().split('\t')
                            G1,G2 = sorted([G1,G2])
                            if (g1==G1) and (g2==G2):
                                file3.write(g1+ '\t'+ g2+ '\t'+ score_1+ '\t'+ score_2+ '\t'+ score_3+ '\t'+ rank_DR+ '\t'+ rank_SCODE+ '\t'+ rank_PIDC+ '\t'+ added_rank + '\n')
                        i = i+1
                        print('n-top processed: ', i)
#Hay algo que no funciona bien porque no te está detectando los comunes con el gold standard
                        '''