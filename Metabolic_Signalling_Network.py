# -*- coding: utf-8 -*-
"""
Created on Mon May 27 01:05:46 2019

@author: abeat
"""
'''First:  run_single_ko or run_multiple_ko
Later: get_matched_genes, get_results, get_results_csv, show_biomass, find_blocked_essential_reactions, blocked_essential_reactions_to_csv, get_different_blocked_essentiall_reactions
'''

from Study_Metabolism import Study_Metabolism
from Study_Signalling import Study_Signalling
import cobra
from file_utils import pickle_object, read_pickle
import pandas as pd
import numpy


class Metabolic_Signalling_Network():
    
    def __init__(self): 
        self.signalling_path = Study_Signalling()
        self.metabolism_data = Study_Metabolism()
        self.matched_genes = {}
        self.model = None
        self.res_fba = {}
        self.res_fva = {}
        self.blocked_reactions = {}
        self.essential_reactions = {}
            
    def __initialize_study_signalling(self): #initialize the study_signalling object, imports the necessary data and build the pathways
        self.signalling_path.import_data() #can be empty or 'study_and_database'
        self.signalling_path.run()
    
    def __initialize_study_metabolism(self): #initialize the study_metabolis object, imports the recon and do the conversion
        self.metabolism_data.import_recon() #pode ser vazio ou 'recon'
        self.metabolism_data.run()
        self.model = self.metabolism_data.get_recon()
    
    def initialize(self): #initialize both the objects from the other classes
        self.__initialize_study_signalling()
        self.__initialize_study_metabolism()
        
    def __match_genes(self, genes_list): #match the genes from a given list with the genes in Recon3D
        model_genes = self.metabolism_data.get_genes_with_ids()
        for gene in genes_list:
            if gene in model_genes.keys():
                self.matched_genes[gene] = model_genes[gene]
    
    def get_matched_genes(self): #return the matched genes
        return self.matched_genes

    def simulate_nc(self): #performs the FBA and FVA simulations
        self.res_fba['NC'] = self.model.optimize()
        self.res_fva['NC'] = cobra.flux_analysis.flux_variability_analysis(self.model)

    def simulate_single_ko(self): #performs the FBA and FVA simulations for single KO
        for gene in self.matched_genes:
            gene_id = self.matched_genes[gene]
            with self.model as temp_model:
               temp_model.genes.get_by_id(gene_id).knock_out()
               self.res_fba[gene] = temp_model.optimize()
    
    def simulate_multiple_ko(self): #performs the FBA and FVA simulations for multiple KO
        dic_genes = self.signalling_path.get_metabolites_path_ends()
        for k in dic_genes:
            list_genes = dic_genes[k]
            self.matched_genes = {}
            self.__match_genes(list_genes)
            name = k
            with self.model as temp_model:
                for gene in self.matched_genes:
                    gene_id = self.matched_genes[gene]
                    temp_model.genes.get_by_id(gene_id).knock_out()
                self.res_fva[name] = cobra.flux_analysis.flux_variability_analysis(temp_model)
                self.res_fba[name] = temp_model.optimize()
    
    def run_multiple_ko(self): #runs the simulations for multiple KO and NC
        self.initialize()
        self.simulate_nc()
        self.simulate_multiple_ko()
        pickle_object(self.res_fba, 'res_fba_multiple_ko.pkl')
        pickle_object(self.res_fva, 'res_fva_multiple_ko.pkl')
    
    def run_single_ko (self): #runs the simulations for single KO and NC
        self.initialize()
        self.__match_genes(self.signalling_path.get_path_ends())
        self.simulate_nc()
        self.simulate_single_ko()
        pickle_object(self.res_fba, 'res_fba.pkl')
    
    def __import_results_fva(self): #import the results from fva simulation
        self.res_fva = read_pickle('res_fva_multiple_ko.pkl')
            
    def get_results(self, res_from = 'fba'): #return the results from the simulations as data frames
        if res_from == 'fba':
            data =  pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()})
            return data
        elif res_from == 'fva':
            self.__import_results_fva()
            return self.res_fva
    
    def get_results_csv(self, results_from = 'fba'): #create a csv file with the results
        if results_from == 'fba':
            data = self.get_results()
            data.to_csv('Results_fba.csv')
        elif results_from == 'fva':
            self.__import_results_fva()
            lista = list(self.res_fva.keys())
            data = self.res_fva[lista[0]]
            data.columns = (str(lista[0])+'_minimum', str(lista[0])+'_maximum')
            for i in range(1, len(lista)):
                data[str(lista[i])+'_minimum'] = self.res_fva[lista[i]]['minimum']
                data[str(lista[i])+'_maximum'] = self.res_fva[lista[i]]['maximum']
            data.to_csv('Results_fva.csv')
    
    def show_biomass(self): #return the biomass production for each condition
        return pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()}).loc['biomass_reaction',:]
        
    def find_blocked_essential_reactions(self): #determine the blocked and essential reactions
        self.__import_results_fva()
        for k in self.res_fva:
            self.blocked_reactions[k] = []
            self.essential_reactions[k] = []
            data = self.res_fva[k]
            for i in range(len(data)):
                minimum = data.iloc[i]['minimum']
                maximum = data.iloc[i]['maximum']
                if minimum == 0 and maximum == 0:
                    self.blocked_reactions[k].append(data.iloc[i].name)
                elif minimum < 0 and maximum < 0:
                    self.essential_reactions[k].append(data.iloc[i].name)
                elif minimum > 0 and maximum > 0:
                    self.essential_reactions[k].append(data.iloc[i].name)    
    
    def blocked_essential_reactions_to_csv(self): #create a csv file with the blocked and essential reactions
        blocked, essential = self.find_blocked_essential_reactions()
        r_blocked = {}
        for k in blocked:
            r_blocked[k] = []
            for r in blocked[k]:
                r_function = self.model.reactions.get_by_id(r).subsystem.split("'")[1]
                r_blocked[k].append(r_function)
            filename = 'fva_multiple_blocked_' + str(k) + '.txt'
            with open(filename, 'w') as file:
                file.write('\n'.join(r_blocked[k]))
        r_essential = {}
        for k in essential:
            r_essential[k] = []
            for r in essential[k]:
                r_function = self.model.reactions.get_by_id(r).subsystem.split("'")[1]
                r_essential[k].append(r_function)
            filename = 'fva_multiple_essential_' + str(k) + '.txt'
            with open(filename, 'w') as file:
                file.write('\n'.join(r_essential[k]))
    
    
    def __get_differences(self): #return the differences between conditions
        b = list(self.blocked_reactions.keys())
        e = list(self.essential_reactions.keys())
        diff_b = {}
        for i in range(len(b)):
            for j in range(len(b)):
                if i != j:
                    name = str(b[i]) + '-' + str(b[j])
                    diff_b[name] = list(numpy.setdiff1d(self.blocked_reactions[b[i]], self.blocked_reactions[b[j]]))
        diff_e = {}
        for i in range(len(e)):
            for j in range(len(b)):
                if i != j:
                    name = str(e[i]) + '-' + str(e[j])
                    diff_e[name] = list(numpy.setdiff1d(self.essential_reactions[e[i]], self.essential_reactions[e[j]]))
        return diff_b, diff_e
    
    def get_different_blocked_essential_reactions(self): ##create a csv file with the results of the different blocked and essential reactions
        blocked, essential = self.get_differences()
        r_blocked = {}
        r_essential = {}
        for k in blocked:
            reactions = blocked[k]
            r_blocked[k] = [[], []] #[name], [function]
            for r in reactions:
                r_name = self.model.reactions.get_by_id(r).name
                r_function = self.model.reactions.get_by_id(r).subsystem.split("'")[1]
                if r_name not in r_blocked[k][0]:
                    r_blocked[k][0].append(r_name)
                    r_blocked[k][1].append(r_function)
            filename = 'fva_multiple_blocked_' + str(k) + '.txt'
            with open(filename, 'w') as file:
                file.write('\n'.join(r_blocked[k][1]))
        for k in essential:
            reactions = essential[k]
            r_essential[k] = [[], []] #[name], [function]
            for r in reactions:
                r_name = self.model.reactions.get_by_id(r).name
                r_function = self.model.reactions.get_by_id(r).subsystem.split("'")[1]
                if r_name not in r_essential[k][0]:
                    r_essential[k][0].append(r_name)
                    r_essential[k][1].append(r_function)
            filename = 'fva_multiple_essential_' + str(k) + '.txt'
            with open(filename, 'w') as file:
                file.write('\n'.join(r_essential[k][1]))
        return r_blocked, r_essential
        
    
if __name__ == "__main__":
    Study = Metabolic_Signalling_Network()
    Study.run_single_ko()
    Study.run_multiple_ko()