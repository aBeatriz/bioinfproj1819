from Estudo_recon import Study_Metabolism
from Estudo_database import Study_Signalling
import cobra
from file_utils import pickle_object, read_pickle
import pandas as pd


class Metabolic_Signalling_Network():
    
    def __init__(self): 
        self.signalling_path = Study_Signalling()
        self.metabolism_data = Study_Metabolism()
        self.macthed_genes = {}
        self.model = None
        self.res_wt_ko = {}
            
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
        model_genes = self.metabolism_data.get_genes_with_index()
        for gene in path_ends:
            if gene in model_genes.keys():
                self.macthed_genes[gene] = model_genes[gene]
    
    def get_matched_genes(self):
        return self.macthed_genes
        
    def simulate_ko(self):
        for gene in self.macthed_genes.keys():
            gene_id = self.macthed_genes[gene]
            with self.model as temp_model:
               temp_model.genes.get_by_id(gene_id).knock_out()
               self.res_wt_ko[gene] = temp_model               
    
    def simulate_wt(self):
        self.res_wt_ko['wt'] = self.model
    
    def run (self):
        self.initialize()
        self.match_genes()
        self.simulate_wt()
        self.simulate_ko()
    
    def save_simulation_fba(self):
        pickle_object(self.res_wt_ko, 'res_wt_ko.pkl')
            
    def get_results_simulation_fba(self):
        return pd.DataFrame({k:v.optimize().fluxes for k,v in self.res_wt_ko.items()})
    
    def show_item(self, item = 'biomass_reaction'):
        return pd.DataFrame({k:v.optimize().fluxes for k,v in self.res_wt_ko.items()}).loc['biomass_reaction',:]
    
    
if __name__ == "__main__":
    Study = Metabolic_Signalling_Network()
    Study.run()
    print(Study.get_results_simulation_fba())
