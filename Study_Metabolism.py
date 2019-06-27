# -*- coding: utf-8 -*-
"""
Created on Mon May 27 01:05:21 2019

@author: abeat
"""

from __future__ import print_function
import cobra
from file_utils import read_pickle, pickle_object

'''First: import_data
Second: run
Later: get_recon, get_genes_with_ids, save_recon, save_data
'''

class Study_Metabolism():
    
    def __init__(self):
        self.converter = {}
        self.recon = None
        self.genes = {}
         
    def __import_converter(self):
        file = open('conversion_genes.txt', 'r') #doc from https://www.genenames.org/download/custom/
        lines = file.readlines()
        for line in lines:
            l = line.split('\t')
            symbol = l[0]
            ncbi_id = l[1]
            if ncbi_id not in self.converter:
                self.converter[ncbi_id] = symbol
        file.close()
    
    def import_recon(self, data_from = 'NULL'):
        if data_from == 'recon':
            self.recon = cobra.io.load_matlab_model('Recon3D_301.mat')
        else:
            self.recon = read_pickle('Recon3D.pkl')
    
    def __conversion_genes_ids(self):
        for i in range(len(self.recon.genes)):
            gene_id = self.recon.genes[i].id
            ncbi_id = gene_id.split('.')[0] #primeira parte corresponde ao numero no ncbi
            if ncbi_id in self.converter: #nao tem todos
                symbol = self.converter[ncbi_id]
                if symbol not in self.genes: #existem mais que uma entrada para cada gene
                    self.genes[symbol] = gene_id #atribui-se cada simbolo ao id respetivo que se encontra no recon
    
    def get_recon(self):
        return self.recon
    
    def get_genes_with_ids(self):
        return self.genes
    
    def save_recon(self):
        pickle_object(self.recon, 'Recon3D.pkl')
    
    def save_data(self):
        pickle_object(self.genes, 'genes_symbols_recon.pkl')
        
    def run(self):
        self.__import_converter()
        self.__conversion_genes_ids()
        self.save_data()
        
           
if __name__ == "__main__":
    Metabolic_data = Study_Metabolism()
    Metabolic_data.save_recon()
    Metabolic_data.run()    
