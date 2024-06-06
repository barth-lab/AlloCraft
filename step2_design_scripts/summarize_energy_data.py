
import os
import re
import sys
from random import random
from glob import glob
import pandas as pd
from collections import OrderedDict
import click

pd.set_option('max_colwidth', 100)


DB_PATH  = '/path/to/ism/sims_db.csv'
SIM_PATH = '/path/to/ism/'

EXCLUDE_LIST = ['input', 'input_file', 'input_files']

pd_data  = None




#Load or Create sims_db.csv
if not os.path.exists(DB_PATH):
    db_columns = ['receptor', 'structure', 'mutations', 
                  'energy', 'state', 'sims', 'directory', 
                  'path', 'min_energy_file']

    pd_data = pd.DataFrame(columns=db_columns)

    pd_data.to_csv(DB_PATH, sep='|', )

    print('the DB does not exist.')

else:
    pd_data = pd.read_csv(DB_PATH, sep='|', index_col=0)

@click.command()
@click.argument('path', default=SIM_PATH, type=click.Path(exists=True))
@click.option('--ignore/--heed', default=False, help='Ignore "STARTED" and "FINISHED" status files.')
@click.option('--threshold', '-t', default=200, help='Summarize energies with at least n number of pdbs')
@click.option('--save/--test', default=True, help='Save parsed data to DB or test run only.')
def get_parameters(path=SIM_PATH, save=False, ignore=False, threshold=200):

    walk_sim_path(the_path=path, save_db=save, ignore=ignore, threshold=threshold)


def get_pdb_energies(the_path):

    path    = the_path
    pdb_dir = OrderedDict()


    if not os.path.exists(path):
            sys.exit('The path "{}" does not exist.'.format(path))

    if path[-1] == '*':
            path = path[0:-1]

    if not path[-1] == '/':
            path = path + '/'

    #print(path)

    pdb_files = glob(path + '*.pdb')

    for file in (pdb for pdb in pdb_files if re.search('.*?_\d{4}\.pdb$', pdb)):

            with open(file) as infile:
                
                filename = os.path.split(file)[1]
                
                for line in infile:
                    
                    #the bk_tot line always starts with a space
                    if line[0] == ' ':
                        
                        result = re.search('bk_tot:\s(.*)', line)
                    else:
                        result = None
                    
                    if result is not None:

                            pdb_dir[filename] = float(result.group(1))

                            break

    sorted_pdb_dict = OrderedDict(sorted(pdb_dir.items(), key=lambda x: x[1]))

    #print(sorted_pdb_dict)

    pd_data = pd.DataFrame(pd.DataFrame({'file' : list(sorted_pdb_dict.keys()), 
                                         'energy' : list(sorted_pdb_dict.values())}))

    return pd_data


def walk_sim_path(the_path='/raid1/dkeri/sims/11-21-2018/', save_db=False, ignore=False, threshold=200):

    global pd_data
    
    #print(the_path)


    for path_info in os.walk(the_path):
        

        print(path_info[0])

        dir_name    = path_info[0].rsplit('/')[-1]
        parent_name = path_info[0].rsplit('/')[-2]


        if dir_name in EXCLUDE_LIST:
                continue

        if 'STARTED' in path_info[2]:

            if ignore:
                pass

            elif not 'FINISHED' in path_info[2]:

                continue

        
        pdb_files   = [file for file in path_info[2] 
                       if re.search('.*?_\d{4}\.pdb$', file) is not None]
       
        energy_file = path_info[0] + '/energies.csv'
 
        if ((len(pdb_files) >= 1 and not os.path.exists(energy_file))
         or (len(pdb_files) >= 1 and sum([1 for x in open(energy_file)])-1 < len(pdb_files))):


            if len(pdb_files) < threshold:

                continue

            print('files in ', path_info[0])

            e_data = get_pdb_energies(path_info[0])

            e_data.to_csv(path_info[0] + '/energies.csv', index=False)
	    

        else:
            #for folder in path_info[1]:
            #    walk_sim_path(folder, save_db=save_db, ignore=ignore, threshold=threshold)

            pass


    if save_db:
        pd_data.to_csv(DB_PATH, sep='|')


if __name__ == '__main__':

    '''if len(sys.argv) > 1:
                    if os.path.exists(sys.argv[1]):
                        walk_sim_path(sys.argv[1], save_db=True)
            
                    else:
                        sys.exit('{} does not exist.'.format(sys.argv[1]))
            
                else:
            
                    walk_sim_path(SIM_PATH, save_db=True)'''


    get_parameters()

#pd_data.to_csv(DB_PATH, sep='|')
