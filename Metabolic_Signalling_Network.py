# -*- coding: utf-8 -*-
"""
Created on Mon May 27 01:05:46 2019

@author: abeat
"""

from Study_Metabolism import Study_Metabolism
from Study_Signalling import Study_Signalling
import cobra
from file_utils import pickle_object, read_pickle
import pandas as pd


class Metabolic_Signalling_Network():
    
    def __init__(self): 
        self.signalling_path = Study_Signalling()
        self.metabolism_data = Study_Metabolism()
        self.matched_genes = {}
        self.model = None
        self.res_fba = {}
        self.res_fva = {}
        self.res_fva_final = {}
            
    def __initialize_study_signalling(self):
        self.signalling_path.import_data() #pode ser vazio ou 'study_and_database'
        self.signalling_path.run()
    
    def __initialize_study_metabolism(self):
        self.metabolism_data.import_recon() #pode ser vazio ou 'recon'
        self.metabolism_data.run()
        self.model = self.metabolism_data.get_recon()
    
    def initialize(self):
        self.__initialize_study_signalling()
        self.__initialize_study_metabolism()
        
    def match_genes(self):
        path_ends = self.signalling_path.get_path_ends()
        model_genes = self.metabolism_data.get_genes_with_ids()
        for gene in path_ends:
            if gene in model_genes.keys():
                self.matched_genes[gene] = model_genes[gene]
    
    def get_matched_genes(self):
        return self.matched_genes
        
    def simulate_ko(self):
        for gene in self.matched_genes:
            gene_id = self.matched_genes[gene]
            with self.model as temp_model:
               temp_model.genes.get_by_id(gene_id).knock_out()
               self.res_fba[gene] = temp_model.optimize()
#                   self.res_fva[gene] = cobra.flux_analysis.flux_variability_analysis(temp_model)
    
    def simulate_nc(self):
        self.res_fba['NC'] = self.model.optimize()
#        self.res_fva['NC'] = cobra.flux_analysis.flux_variability_analysis(self.model)
    
    def run (self):
        self.initialize()
        self.match_genes()
#        self.simulate_nc()
#        self.simulate_ko()
        self.save_simulation()
    
    def save_simulation(self):
#        pickle_object(self.res_fba, 'res_fba.pkl')
        pickle_object(self.matched_genes, 'matched.pkl')
#        pickle_object(self.res_fva, 'res_fva.pkl')
            
    def get_results(self, res_from = 'fba'):
        if res_from == 'fba':
            data =  pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()})
            data.to_csv('Results_fba.csv')
            return data
        elif res_from == 'fva':
            data =  pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()})
            data.to_csv('Results_fva.csv')
            return data
        elif res_from == 'matched':
            with open('matched.txt', 'w') as f: #Para o go
                lista = [s.split('.')[0] for s in list(self.matched_genes.values())]
                f.write('\n'.join(lista))
    
    def show_item(self, res_from = 'fba', item = 'biomass_reaction'):
        if res_from == 'fba':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()}).loc['biomass_reaction',:]
        elif res_from == 'fva':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()}).loc['biomass_reaction',:]
    
    def get_genes_associated_bc_literature(self):
        genes_literature = []
        with open('gene_breast_cancer.txt', 'r') as file:
            lines = file.readlines()
            for l in range(1, len(lines)):
                line = lines[l]
                genes_literature.append(line.split('\t')[5])
        res = []
        for gene in self.matched_genes:
            if gene in genes_literature:
                res.append(self.matched_genes[gene].split('.')[0])
        with open('genes_bc_macthed.txt', 'w') as file:
            file.write('\n'.join(res))
    
    def fva_results(self):
        res = {'blocked_r':None, 'essential_r':None, 'essential_g':None}
        for k in self.res_fva:
            fva = self.res_fva[k]
            res['blocked_r'] = fva.find_blocked
            pass
    
    def alter_reactions_bounds(self):
        pass
            
    
if __name__ == "__main__":
    Study = Metabolic_Signalling_Network()
    Study.run()
#    Study.get_results()
#    Study.get_results('matched')
#    print((Study.show_item()))
#    print(Study.fva_results())
    print(Study.get_matched_genes().keys())