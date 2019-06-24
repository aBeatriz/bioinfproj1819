# -*- coding: utf-8 -*-
"""
Created on Mon May 27 01:05:46 2019

@author: abeat
"""

from Estudo_recon import Study_Metabolism
import cobra
from file_utils import pickle_object, read_pickle
import pandas as pd


class Metabolic_Signalling_Network():
    
    def __init__(self): 
        self.metabolism_data = Study_Metabolism()
        self.macthed_genes = {}
        self.model = None
#        self.res_fba = {}
        self.res_fva = {}
    
    def __initialize_study_metabolism(self):
        self.metabolism_data.import_recon() #pode ser vazio ou 'recon'
        self.metabolism_data.run()
        self.model = self.metabolism_data.get_recon()
    
    def initialize(self):
        self.__initialize_study_metabolism()
        
    def match_restantes_genes(self):
        path_ends = read_pickle('restantes.pkl')
        model_genes = self.metabolism_data.get_genes_with_ids()
        for gene in path_ends:
            if gene in model_genes.keys():
                self.macthed_genes[gene] = model_genes[gene]
    
    def get_matched_genes(self):
        return self.macthed_genes
        
    def simulate_ko(self):
        for gene in self.macthed_genes:
            gene_id = self.macthed_genes[gene]
            with self.model as temp_model:
               temp_model.genes.get_by_id(gene_id).knock_out()
               self.res_fva[gene] = cobra.flux_analysis.flux_variability_analysis(temp_model)
#               self.res_fba[gene] = temp_model.optimize()
    
    def simulate_wt(self):
        self.res_fva['WT'] = cobra.flux_analysis.flux_variability_analysis(self.model)
#        self.res_fba['WT'] = self.model.optimize()
    
    def run (self):
        self.initialize()
        self.match_restantes_genes()
        self.simulate_wt()
        self.simulate_ko()
        self.save_simulation()
    
    def save_simulation(self):
#        pickle_object(self.res_fba, 'res_fba.pkl')
        pickle_object(self.res_fva, 'res_fva.pkl')
            
    def get_results(self, res_from = 'fva'):
        if res_from == 'fba':
            data =  pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()})
            data.to_csv('Results_fba_r.csv')
            return data
        elif res_from == 'fva':
            data =  pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()})
            data.to_csv('Results_fva.csv')
            return data
        elif res_from == 'matched':
            with open('matched.txt', 'w') as f: #Para o go
                lista = [s.split('.')[0] for s in list(self.macthed_genes.values())]
                f.write('\n'.join(lista))
    
    def show_item(self, res_from = 'fva', item = 'biomass_reaction'):
        if res_from == 'fba':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()}).loc['biomass_reaction',:]
        elif res_from == 'fva':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()}).loc['biomass_reaction',:]
    
if __name__ == "__main__":
    Study = Metabolic_Signalling_Network()
    Study.run()
    print((Study.get_results()))
