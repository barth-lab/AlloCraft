import re
import os
import sys
import pandas as pd
#import numpy as np
import urllib
from datetime import datetime
from collections import Counter
#import requests
#from bs4 import BeautifulSoup as bs
#import seaborn as sb
from time import sleep
#import matplotlib.pylab as plot

#import functools
from functools import wraps

#pd.options.display.max_rows    = 2000
#pd.options.display.max_columns = 40
#pd.options.display.max_colwidth= 80




def get_id_tuple(f, args, kwargs):
    """ 
    Some quick'n'dirty way to generate a unique key for an specific call.
    """
    l = [id(f)]
    for arg in args:
        l.append(arg)

    for k, v in kwargs.items():
        l.append(k)
        l.append(v)

    return l

def memoize(f):
    """ 
    Some basic memoizer
    """

    _memoized = {}

    @wraps(f)
    def memoized(*args, **kwargs):
        key = get_id_tuple(f, args, kwargs)

        if str(key) not in _memoized:
            _memoized[str(key)] = f(*args, **kwargs)
        return _memoized[str(key)]
    return memoized



#@functools.lru_cache()
#@memoize
class pdbPandas:
    
    def __init__(self, pdb_file, load_only=False):

        start = datetime.now()

            
        if not os.path.exists(pdb_file) and \
           not self.is_path_to_file(pdb_file):
            self.pd_pdb = self.read_data(self.download_pdb_data(pdb_file))
                
        else: self.pd_pdb = self.read_data(pdb_file)
            
        #self.pd_pdb = self.pd_pdb.replace(r'', np.nan, regex=False)
        
        
        #needed to make sure the 'res_seq' column can be converted to int
        #if there is no value for TER lines for this column
        if len(self.pd_pdb[self.pd_pdb['type'] == 'TER']):
            for ter in self.pd_pdb[self.pd_pdb['type'] == 'TER'].index:
                if self.pd_pdb.loc[ter, 'res_seq'] == '':
                    self.pd_pdb.loc[ter, 'res_seq'] = self.pd_pdb.loc[ter-1, 'res_seq']
                    
        self.pd_pdb['res_seq'] = self.pd_pdb['res_seq'].astype(int)
        self.pd_pdb = self.pd_pdb.sort_values(by='res_seq')


        self.distance_memo = dict()






        #print(end-start, 'loaded pdb')
        
        
        
    def is_path_to_file(self, string):
        if ('/' in string or '\\' in string) and \
           not re.search('\.[\w\d]{2,4}$', string) is None:
            return True
        else: return False
        

    def download_pdb_data(self, query):
        
        url = ('http://www.rcsb.org/pdb/download/downloadFile.do'
               '?fileFormat=pdb&compression=NO&structureId=%s' 
                % query.upper()
              )

        try:
            print('downloading %s...' % query.upper(), end='')
            pdb_data = urllib.request.urlopen(url)

            print('done.')

            return pdb_data.read().decode()

        except Exception as e: 
            print(e)
            return None
        
        
    
    def read_data(self, pdb_file):
        
        file_info = list()
        
        try:
            for line in pdb_file.split('\n') if not os.path.exists(pdb_file) else open(pdb_file):
                result = re.search('^(ATOM|TER|HETATM)', line)
                
                #fix for short TER lines from Rosetta
                if 'TER' in line and len(line) < 80:
                    line += ' ' * (80 - len(line))

                if not result is None:
                    line_info = [line[0:6], line[6:11], line[12:16], 
                                 line[16], line[17:20], line[21], 
                                 line[22:26], line[26], line[30:38], 
                                 line[38:47], line[46:54], line[54:60], 
                                 line[60:66], line[72:76], 
                                 line[76:78], line[78:80]]

                    line_info = [info.strip() for info in line_info]
                    file_info.append(line_info)
        
        #windows path is too long error
        except ValueError:
            
            for line in pdb_file.split('\n'):
                result = re.search('^(ATOM|TER|HETATM)', line)

                if not result is None:
                    line_info = [line[0:6], line[6:11], line[12:16], 
                                 line[16], line[17:20], line[21], 
                                 line[22:26], line[26], line[30:38], 
                                 line[38:47], line[46:54], line[54:60], 
                                 line[60:66], line[72:76], 
                                 line[76:78], line[78:80]]

                    line_info = [info.strip() for info in line_info]
                    file_info.append(line_info)
                    
            
                
                
        pd_pdb = pd.DataFrame(file_info, columns=
                              ['type', 'serial', 'atom_name', 'alt_loc', 
                               'res_name', 'chain_id', 'res_seq', 'i_code',
                               'x', 'y', 'z', 'occupancy', 'temp_factor', 
                               'seg_id', 'element', 'charge'], dtype=int)
        
        return pd_pdb.sort_values(['res_seq', 'serial', 'type'])
        
        
    
    def extract_first_chain(self, pd_pdb):
        
        term = None
        
        try: term = pd_pdb[(pd_pdb['type'] == 'TER')].index.tolist()[0]
        except: term = pd_pdb[(pd_pdb['type'] == 'TER')].index
        
        try: return pd_pdb[pd_pdb.index <= term]
        except: None
            
            
    def create_new_numbering(self):
        
        count = 1
        
        numbers = sorted(list(set(self.pd_pdb.res_seq.astype(int))))
        
        if 'new_nums' in self.pd_pdb.columns: 
            self.pd_pdb.drop('new_nums', axis=1, inplace=True)
        
        for number in numbers:
        
            self.pd_pdb.loc[self.pd_pdb['res_seq'] == number, 'new_nums'] = count
            count += 1
            
        self.pd_pdb['new_nums'] = self.pd_pdb['new_nums'].astype('int64')
            
            
            
                
    def save_pdb(self, pdb_file, renumber=False):
        
        outfile = open(pdb_file, mode='w+')
        out_data = None
        
        if not renumber: out_data = self.pd_pdb.drop(['new_nums'], axis=1)
        else: 
            out_data = self.pd_pdb.copy()
            out_data['res_seq'] = out_data['new_nums']
            out_data.drop(['new_nums'], axis=1, inplace=True)
        
        
        for line in out_data.values:
            
            outfile.write(('{:6s}{:>5n}  {:<3s}{:1s}'\
                           '{:3s} {:1s}{:>4} {:1s}  '\
                           '{:>8s}{:>8s}{:>8s}{:>6s}'\
                           '{:>6s}       {:4s}{:2s}{:2s}\n').format(*line))
            

    def get_residue_name(self, position, one_letter=False):


        aa_trans = {"ALA" : "A", "ARG" : "R", "ASN" : "N", 
                    "ASP" : "D", "CYS" : "C", "GLN" : "Q", 
                    "GLU" : "E", "GLY" : "G", "HIS" : "H", 
                    "ILE" : "I", "LEU" : "L", "LYS" : "K", 
                    "MET" : "M", "PHE" : "F", "PRO" : "P", 
                    "SER" : "S", "THR" : "T", "TRP" : "W", 
                    "TYR" : "Y", "VAL" : "V"}

        residue = self.pd_pdb[(self.pd_pdb['alt_loc'] != 'B') & 
                              (self.pd_pdb['atom_name'] == 'CA') &
                              (self.pd_pdb['res_seq'] == position)]['res_name'].iloc[0]

        if one_letter:

            return aa_trans[residue.upper()]

        else:

            return residue.upper()


            
    def get_sequence(self, one_letter=False):
        
        aa_trans = {"ALA" : "A", "ARG" : "R", "ASN" : "N", 
                    "ASP" : "D", "CYS" : "C", "GLN" : "Q", 
                    "GLU" : "E", "GLY" : "G", "HIS" : "H", 
                    "ILE" : "I", "LEU" : "L", "LYS" : "K", 
                    "MET" : "M", "PHE" : "F", "PRO" : "P", 
                    "SER" : "S", "THR" : "T", "TRP" : "W", 
                    "TYR" : "Y", "VAL" : "V"}
        
        ca_atoms = self.pd_pdb[(self.pd_pdb.type == 'ATOM') & 
                               (self.pd_pdb.atom_name == 'CA') &
                               ((self.pd_pdb.alt_loc == '') | 
                                (self.pd_pdb.alt_loc == 'A'))]
        
        
        if one_letter: 
            return ''.join(aa_trans[aa] for aa in ca_atoms.res_name)
        
        else:
            return list(ca_atoms.res_name)
            
    def print_aa_frequency(self):
        aa_stat = Counter(self.get_sequence())
        for name, aa in aa_stat.most_common():
            print('{:3} {:3} {:6.2f}%'.format(name, aa, aa/len(self.get_sequence())*100))
            
    def remove_hetatms(self):
        self.pd_pdb = self.pd_pdb[self.pd_pdb['type'] != 'HETATM']
        
        self.create_new_numbering()
            
                   

    def get_adj_residues(self, residue, distance):


        start = datetime.now()

        try:
            return self.distance_memo[residue]
        except:
            pass

        residue_atoms = self.pd_pdb[self.pd_pdb.res_seq == residue]

        adj_residues = list()

        for atom in residue_atoms.index:
            #print(atom[['x', 'y', 'z']])

            if self.pd_pdb.iloc[atom].type == 'TER':
                continue

            atom_x = float(self.pd_pdb.iloc[atom]['x'])
            atom_y = float(self.pd_pdb.iloc[atom]['y'])
            atom_z = float(self.pd_pdb.iloc[atom]['z'])

            atom_positions = self.pd_pdb[['res_seq', 'x', 'y', 'z']]
            atom_positions = atom_positions.replace('', pd.np.nan)
            atom_positions = atom_positions.dropna(how='any')

            atom_positions['x'] = atom_positions['x'].astype(float)
            atom_positions['y'] = atom_positions['y'].astype(float)
            atom_positions['z'] = atom_positions['z'].astype(float)


            atom_positions['distance'] = ((atom_positions['x'] - atom_x) ** 2 +
                                          (atom_positions['y'] - atom_y) ** 2 +
                                          (atom_positions['z'] - atom_z) ** 2) ** 0.5



            #print(set(atom_positions[atom_positions['distance'] <= distance]['res_seq']))

            adj_residues += list(atom_positions[atom_positions['distance'] <= distance]['res_seq'])

        end = datetime.now()


        self.distance_memo[residue] = sorted(set(adj_residues))

        #print(end-start, 'distances calculated')

        return sorted(set(adj_residues))

        #print(atom_positions)
    
#    def map_contacts(self):
#        
#        ca_atoms = self.pd_pdb[(self.pd_pdb['alt_loc'] != 'B') & 
#                               (self.pd_pdb['atom_name'] == 'CA')]
#        
#        distances = list()
#        
#        #print(ca_atoms)
#        
#        for atom in ca_atoms.index:
#            
#            atom_x = float(ca_atoms.loc[atom]['x'])
#            atom_y = float(ca_atoms.loc[atom]['y'])
#            atom_z = float(ca_atoms.loc[atom]['z'])
#            
#            
#            dist = list(((ca_atoms.x.astype(float) - atom_x) ** 2 + 
#                         (ca_atoms.y.astype(float) - atom_y) ** 2 +
#                         (ca_atoms.z.astype(float) - atom_z) ** 2) ** 0.5)
#
#            dist = [x if x <= 10 else 15 for x  in dist ]
#            
#            distances.append(dist)
#            
#        #print(distances)
#        
#        fig, ax = plot.subplots(figsize=(30,20))
#
#        dist_map = sb.heatmap(distances, ax=ax, annot=False, fmt='.0f',)
#        dist_map_fig = dist_map.get_figure()
#        dist_map_fig.savefig('c:/users/danny/documents/contact_hm_5FTN_chainA.png', format='png')
#
#        #dist_map.get_figure()


if __name__ == '__main__':

    start = datetime.now()
    #pdb_path = "C:/Users/Danny/Documents/Research/Scripts/mc4r_from_4iar.pdb"
    #pdb = pdbPandas(pdb_path)

    pdb_path = "C:\\Users\\Danny\\Documents\\Research\\Collaborations\\KayLE-2016\\Structures\\Modified\\5FTN_chainA.pdb"
    pdb = pdbPandas(pdb_path)

    #pdb = pdbPandas('5iu4')
    

    #pdb.pd_pdb[(pdb.pd_pdb['alt_loc'] != 'B') & 
    #           (pdb.pd_pdb['atom_name'] == 'CA')]


    #print(pdb.get_adj_residues(150, 5))
    #print(pdb.get_adj_residues(150, 7))

    adj_r = pdb.get_adj_residues(150, 10)
    #print(len(pdb.get_adj_residues(150, 10)))

    

    end = datetime.now()

    pdb.map_contacts()

    print(end-start)

