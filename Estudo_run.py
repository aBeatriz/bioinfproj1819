from Estudo_recon import Study_Metabolism
from Estudo_database import Study_Signalling
import cobra
from file_utils import pickle_object
import pandas as pd


class Metabolic_Signalling_Network():
    
    def __init__(self): 
        self.signalling_path = Study_Signalling()
        self.metabolism_data = Study_Metabolism()
        self.macthed_genes = {}
        self.model = None
        self.res_fba = {}
        self.res_fva = {}
            
    def __initialize_study_signalling(self):
        self.signalling_path.import_data('study_and_database') #pode ser vazio ou 'study_and_database'
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
        model_genes = self.metabolism_data.get_genes_with_index()
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
               self.res_fba[gene] = temp_model.optimize()
               self.res_fva[gene] = cobra.flux_analysis.flux_variability_analysis(temp_model)
    
    def simulate_wt(self):
        self.res_fba['WT'] = self.model.optimize()
        self.res_fva['WT'] = cobra.flux_analysis.flux_variability_analysis(self.model)
    
    def run (self):
        self.initialize()
        self.match_genes()
        self.simulate_wt()
        self.simulate_ko()
        self.save_simulation()
    
    def save_simulation(self):
        pickle_object(self.res_fba, 'res_fba.pkl')
        pickle_object(self.res_fva, 'res_fva.pkl')
            
    def get_results(self, res_from = 'fba'):
        if res_from == 'fba':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()})
        elif res_from == 'fva':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()})
    
    def show_item(self, res_from = 'fba', item = 'biomass_reaction'):
        if res_from == 'fba':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fba.items()}).loc['biomass_reaction',:]
        elif res_from == 'fva':
            return pd.DataFrame({k:v.fluxes for k,v in self.res_fva.items()}).loc['biomass_reaction',:]
    
if __name__ == "__main__":
    Study = Metabolic_Signalling_Network()
    Study.run()
    print(Study.get_results())
