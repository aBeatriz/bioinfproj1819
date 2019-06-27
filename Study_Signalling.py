# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:10:29 2019

@author: abeat
"""

from __future__ import print_function
import cobra
from file_utils import read_pickle, pickle_object

'''First: import_data
Second: run
Later: get_database, get_model, get_metabolites_model, get_path_ends, get_pathways, save_initial_data, save_data
'''

class Study_Signalling ():
    
    def __init__(self):
        self.final_pathway = {} #cascata de sinalizacao
        self.database_data = {}
        self.model = None
        self.path_end = [] #todos os targets que nao servem de tf
        self.nodes = []
    
    def import_data(self, data_from = 'NULL'): #import dos dados da base de dados trrust e reg
        if data_from == 'study_and_database':
            self.model = cobra.io.read_sbml_model('Modelo_case.xml')
            trrust = open('trrust_rawdata.human.tsv', 'r')
            lines_t = trrust.readlines()
            for line in lines_t:
                l = line.replace('\n', '').split('\t')
                tf = l[0]
                target = l[1]
                if tf not in self.database_data:
                    self.database_data[tf] = []
                if target not in self.database_data[tf]:
                    self.database_data[tf].append(target) 
            trrust.close()
            reg = open('new_kegg.human.reg.direction.txt', 'r')
            lines_r = reg.readlines()
            for i in range(1, len(lines_r)):
                l = lines_r[i].split(' ')
                tf = l[0]
                target = l[2]
                if tf not in self.database_data:
                    self.database_data[tf] = []
                if target not in self.database_data[tf]:
                    self.database_data[tf].append(target) 
            reg.close()
        elif data_from == 'NULL':
            self.database_data = read_pickle('database.pkl')
            self.model = read_pickle('model.pkl')
     
    def get_database(self):
        return self.database_data
   
    def get_model(self):
        return self.model
    
    def get_metabolites_model(self):
        met = []
        for k in range(len(self.model.metabolites)):
            metabolites = self.model.metabolites[k].name.replace('p', '').replace('n', '').split('-')
            for i in metabolites:
                if i not in met:
                    met.append(i)
        return met
    
    def __add_tf(self, tf, target): #seq1: base de dados; target_seq1: target correspondente na base de dados
        if tf not in self.final_pathway:
            self.final_pathway[tf] = []
            if tf not in self.nodes:
                self.nodes.append(tf)
        if target not in self.final_pathway[tf]: #evitar targets repetidos
            self.final_pathway[tf].append(target)
            if target not in self.nodes:
                self.nodes.append(target)
    
    def __build_pathway(self, metabolites):
        to_visit = []
        for metabolite in metabolites: #comecar com os que aparecem no modelo
            if metabolite in self.database_data:
                for target in self.database_data[metabolite]:
                    self.__add_tf(metabolite, target)
                    if target not in to_visit:
                        to_visit.append(target)
                        
        for tf in to_visit: #continuar a cascata
            if tf in self.database_data: #quer dizer que esta presente na base de dados
                targets = self.database_data[tf]
                for target in targets:
                    self.__add_tf(tf, target)
                    if target not in to_visit:
                        to_visit.append(target)
            elif tf not in self.database_data and tf not in self.path_end:
                self.path_end.append(tf)
            
    def save_initial_data(self):
        pickle_object(self.database_data, 'database.pkl')
        pickle_object(self.model, 'model.pkl')
        
    def save_data(self):
        pickle_object(self.final_pathway, 'final_pathway.pkl')
        pickle_object(self.path_end, 'possible_genes.pkl')         
    
    def run(self):
        self.__build_pathway(self.get_metabolites_model())
        self.save_data()
        
    def print_final_pathway(self):
        for tf in self.final_pathway.keys():
            print (tf, " -> ", self.final_pathway[tf])

    def get_path_ends(self):
        return sorted(self.path_end)
    
    def get_pathways(self):
        return self.final_pathway

    def get_targets(self, tf):
        return self.final_pathway[tf]
    
    def get_tfs(self, target):
        lista = []
        for k in self.final_pathway.keys():
            if target in self.final_pathway[k]:
                lista.append(k)
        return lista

    def active_proteins(self):
        lista = []
        with open('active_proteins.txt', 'r') as f:
            lines = f.readlines()
            for i in range(1, len(lines)):
                lista.append(lines[i].replace('"', '').replace('\n', '').split(',')[1])
        metabolites = []
        for k in range(len(lista)):
            met = lista[k].split('-')
            for i in met:
                if i not in metabolites:
                    metabolites.append(i)
        self.final_pathway = {}
        self.path_end = []
        self.__build_pathway(metabolites)
        

if __name__ == "__main__":
    Pathway = Study_Signalling()
    Pathway.import_data('study_and_database')   
    Pathway.save_initial_data()
    Pathway.run()
