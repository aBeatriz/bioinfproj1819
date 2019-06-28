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
Later: get_database, get_model, get_metabolites_model, get_path_ends, get_pathways, save_initial_data, save_results, print_pathway, get_targets, get_rfs, get_metabolites_path_ends
'''

class Study_Signalling ():
    
    def __init__(self):
        self.final_pathway = {}
        self.database_data = {}
        self.model = None
        self.path_end = [] #possible genes
        self.nodes = []
        self.metabolites_path_ends = {} #pahtways starting points and the genes associated
    
    def import_data(self, data_from = 'NULL'): #to import the necessary data - the data can be empty or 'study_and_database'
        if data_from == 'study_and_database':
            self.model = cobra.io.read_sbml_model('Modelo_case.xml') #study case model
            trrust = open('trrust_rawdata.human.tsv', 'r') #database TRRUST
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
            reg = open('new_kegg.human.reg.direction.txt', 'r') #database RegNetwork
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
        elif data_from == 'NULL': #faster -- database and model are already loaded in these files
            self.database_data = read_pickle('database.pkl')
            self.model = read_pickle('model.pkl')
     
    def get_database(self): #return the repository of regulatory factors and their targets
        return self.database_data
   
    def get_model(self): #return the model
        return self.model
    
    def get_metabolites_model(self): #returns the pool of metabolites from the model
        met = []
        for k in range(len(self.model.metabolites)):
            metabolites = self.model.metabolites[k].name.replace('p', '').replace('n', '').split('-')
            for i in metabolites:
                if i not in met:
                    met.append(i)
        return met
    
    def __add_rf(self, rf, target): #add the regulatory factor and its target to the pathway
        if rf not in self.final_pathway:
            self.final_pathway[rf] = []
            if rf not in self.nodes:
                self.nodes.append(rf)
        if target not in self.final_pathway[rf]:
            self.final_pathway[rf].append(target)
            if target not in self.nodes:
                self.nodes.append(target)
    
    def __build_pathway(self, metabolites): #build the pathway and the list of possible genes
        to_visit = []
        for metabolite in metabolites: #find the pathways' starting points
            if metabolite in self.database_data:
                for target in self.database_data[metabolite]:
                    self.__add_rf(metabolite, target)
                    if target not in to_visit:
                        to_visit.append(target)
                self.metabolites_path_ends[metabolite] = []
                        
        for rf in to_visit: #continue to build the pathways
            if rf in self.database_data:
                targets = self.database_data[rf]
                for target in targets:
                    self.__add_rf(rf, target)
                    if target not in to_visit:
                        to_visit.append(target)
            elif rf not in self.database_data and rf not in self.path_end:
                self.path_end.append(rf)
            
    def save_initial_data(self): #save the initial data already loaded
        pickle_object(self.database_data, 'database.pkl')
        pickle_object(self.model, 'model.pkl')
        
    def save_results(self): #save the results
        pickle_object(self.final_pathway, 'final_pathway.pkl')
        pickle_object(self.path_end, 'possible_genes.pkl')         
    
    def run(self): #build the pathways and save the results
        self.__build_pathway(self.get_metabolites_model())
        self.save_results()
        
    def print_pathway(self): #prints the final pathway
        for tf in self.final_pathway.keys():
            print (tf, " -> ", self.final_pathway[tf])

    def get_path_ends(self): #returns the list with the possible genes
        return sorted(self.path_end)
    
    def get_pathways(self): #returns the final pathway
        return self.final_pathway

    def get_targets(self, rf): #return the target of the regulatory factor in the final pathway
        return self.final_pathway[rf]
    
    def get_rfs(self, target):
        lista = []
        for k in self.final_pathway.keys():
            if target in self.final_pathway[k]:
                lista.append(k)
        return lista
        
    def __associate_metabolites_with_path_ends(self): #find the relationships between the pathways' starting points with the possible genes
        association = {}
        for k in self.metabolites_path_ends:
            association[k] = []
        for rf in self.final_pathway:
            for target in self.final_pathway[rf]:
                for k in association:
                    if rf == k:
                        if target in self.path_end and target not in self.metabolites_path_ends[k]:
                            self.metabolites_path_ends[k].append(target)
                        else:
                            association[k].append(target)
                    else:
                        if rf in association[k]:
                            if target in self.path_end and target not in self.metabolites_path_ends[k]:
                                self.metabolites_path_ends[k].append(target)
                            else:
                                association[k].append(target)
    
    def get_metabolites_path_ends(self): #returns the dictionary with the pathways' starting points and the respective possible genes
        self.__associate_metabolites_with_path_ends()
        return self.metabolites_path_ends

if __name__ == "__main__":
    Pathway = Study_Signalling()
    Pathway.import_data('study_and_database')   
    Pathway.save_initial_data()
    Pathway.run()
    Pathway.save_results()
    print('database: ', Pathway.get_database())
    print('model: ', Pathway.get_model())
    print('metabolites_model: ', Pathway.get_metabolites_model())
    print('path_ends: ', Pathway.get_path_ends())
    print('pathways: ', Pathway.get_pathways())
    print('targets(EGFR): ', Pathway.get_targets('EGFR'))
    print('rfs(STAT3): ', Pathway.get_rfs('STAT3'))
    print('metabolites_path_ends: ', Pathway.get_metabolites_path_ends())
    Pathway.print_pathway()