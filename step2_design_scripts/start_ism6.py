
import re
import os
import sys
import click
import shutil
import subprocess

from pdbPandas import pdbPandas
from datetime  import datetime
from pathlib   import Path

#structures, span files
INPUT_ROOT_PATH    = '/path/to/demo/step2/input/'
ROSETTA_PATH       = '/path/ro/Rosetta_Intel/bin/rosetta.intel'
ROSETTA_DB_PATH    = '/path/to/Rosetta_Intel/rosetta_database'


TODAY_DATE         = datetime.strftime(datetime.now(), "%Y-%m-%d")
OUTPUT_ROOT_PATH   = '/path/to/ism/'
#OUTPUT_ROOT_PATH   = '/home/dk/sims/'
#SIM_OUTPUT_PATH    = click.Path(OUTPUT_ROOT_PATH / TODAY_DATE)

CLUSTER_QUEUE_PATH =  '/path/to/ism/cluster_queue.txt'



print('Using ISM v6')

@click.command()
@click.argument('structs', nargs=-1, ) #type=click.Path(exists=True))
@click.argument('mutants', )#help='Mutations to introduce'
@click.option('--hours', default=3, help='length of the simulation run in hours.')#
@click.option('--nstruct', default=200, help='number of structures to generate.')
@click.option('--radius', default=5.0, help='radius (in angstroms) of NATAA shell.')
@click.option('--cytosolic/--membrane', default=False)
@click.option('--test/--no-test', default=False)
@click.option('--natro/--minimize', default=False, help='disregard mutations and radius and only score structure')
@click.option('--inplace/--paths', default=False, help='output files will be placed in current directory.')
@click.option('--rosetta_path', default=ROSETTA_PATH, type=click.Path(exists=True))
@click.option('--rosetta_db_path', default=ROSETTA_DB_PATH, type=click.Path(exists=True))
@click.option('--input_root', default=INPUT_ROOT_PATH, type=click.Path(exists=True))
@click.option('--output_root', default=OUTPUT_ROOT_PATH, type=click.Path(exists=True))
@click.option('--resfile', default=None, type=click.Path(exists=True), help='path to the resfile')
@click.option('--msdresfile', default=None, type=click.Path(exists=True), help='path to a MSD resfile, which will be converted to the old format')
@click.option('--spanfile', default=None, type=click.Path(exists=True), help='path to the spanfile')
@click.option('--queue/--launch', default=True, help='Add to queue or launch directly with sbatch')
def get_parameters(structs, 
                   mutants, 
                   nstruct, 
                   radius, 
                   hours, 
                   cytosolic, 
                   test, 
                   natro,
                   inplace,
                   rosetta_path,
                   rosetta_db_path,
                   input_root,
                   output_root,
                   resfile,
                   msdresfile,
                   spanfile,
                   queue,
                   ):


    set_up_sim(structs=structs, 
               mutants=mutants.strip(), 
               nstruct=nstruct, 
               radius=radius, 
               hours=hours, 
               is_cytosolic=cytosolic,
               test=test, 
               inplace=inplace, 
               rosetta_path=Path(rosetta_path),
               rosetta_db_path=Path(rosetta_db_path),
               input_root=Path(input_root),
               output_root=Path(output_root),
               resfile=resfile,
               msdresfile=msdresfile,
               spanfile=spanfile,
               queue=queue,
               natro=natro,
              )




def set_up_sim(structs, 
               mutants, 
               nstruct, 
               radius, 
               hours=3, 
               is_cytosolic=False,
               test=False, 
               inplace=False, 
               rosetta_path=ROSETTA_PATH,
               rosetta_db_path=ROSETTA_DB_PATH,
               input_root=INPUT_ROOT_PATH,
               output_root=OUTPUT_ROOT_PATH,
               resfile=None,
               msdresfile=None,
               spanfile=None,
               queue=True,
               natro=False,

              ):


    mutants_list    = parse_mutants(mutants)
    #destination     = None


    #recursively run this function if a mutation is to X
    for mut in mutants_list:
        if 'X' in mut:
            for new_aa in 'ACDEFGHIKLMNPQRSTVWY':

                new_mutants_list = mutants_list.copy()
                new_mutants_list[mutants_list.index(mut)] = mut[:-1] + new_aa
                #print(mutants_list, new_mutants_list)

                set_up_sim(structs, 
                           '-'.join(new_mutants_list),
                           nstruct=nstruct, 
                           radius=radius, 
                           hours=hours, 
                           is_cytosolic=is_cytosolic,
                           test=test, 
                           natro=natro,
                           inplace=inplace, 
                           rosetta_path=rosetta_path,
                           rosetta_db_path=rosetta_db_path,
                           input_root=input_root,
                           output_root=output_root,
                           spanfile=spanfile,
                           queue=True,
                          )

            #Return to end this function
            #Let the recursive functions do the rest
            return None

    print()
    print('Loading pdb\'s...')
    print()

    pd_struct_dict = load_pdb_list(pdb_list=structs, input_root=input_root)


    validate_mutations(mutants_list, pd_struct_dict)

    #to-do:
    #done - set up directories
    #copy pdb and span file (if necessary  to output dir)
    #create resfile in output dir
    #create paths file in output dir
    #create bash file
    #create a parameter file in the output directory

    
    for struct in structs:

        print('Setting up output directory for ', struct)

        setup_output_directory(struct=struct, 
                               mutants=mutants, 
                               radius=radius, 
                               nstruct=nstruct,
                               hours=hours,
                               pd_struct_dict=pd_struct_dict,
                               natro=natro,
                               inplace=inplace, 
                               rosetta_db_path=rosetta_db_path,
                               input_root=input_root, 
                               output_root=output_root,
                               resfile=resfile,
                               msdresfile=msdresfile,
                               spanfile=spanfile,
                               queue=queue,
                               is_cytosolic=is_cytosolic,
                              )



    print('done.')
    print()



def validate_mutations(mutants_list, pd_struct_dict):


    output_str = 'Residue mismatch: {}\n{} is {}'

    for mut in mutants_list:

        for pdb_name, pdb in pd_struct_dict.items():

            res_in_pdb = pdb.get_residue_name(int(re.search('(\d+)', mut).group(1)), True)


            if not mut[0] == res_in_pdb:
                sys.exit(output_str .format(mut, pdb_name, res_in_pdb))

                #or raise error instead?



def load_pdb_list(pdb_list, input_root=INPUT_ROOT_PATH):

    pdb_path_list  = [Path(x) for x in pdb_list]

    pd_struct_dict = dict()

    for pdb_file in pdb_path_list:



        if not pdb_file.exists():

            if pdb_file.is_absolute():
                sys.exit('{} does not exist.'.format(pdb_file))

            else:

                root_pdb_path = input_root/pdb_file.name

                print('WARNING: ./{} does not exist.'.format(pdb_file))
                print('Using {} instead.'.format(root_pdb_path))
                print('')

                pd_struct_dict[str(root_pdb_path)] = pdbPandas(str(root_pdb_path))

        else:

            pd_struct_dict[str(pdb_file)] = pdbPandas(str(pdb_file))


    return pd_struct_dict



def setup_output_directory(struct, 
                           mutants, 
                           radius, 
                           nstruct,
                           hours,
                           pd_struct_dict, 
                           natro=False,
                           inplace=False, 
                           rosetta_db_path=ROSETTA_DB_PATH,
                           input_root=INPUT_ROOT_PATH, 
                           output_root=OUTPUT_ROOT_PATH, 
                           resfile=None,
                           msdresfile=None,
                           spanfile=None,
                           queue=True,
                           is_cytosolic=False,

                          ):

    global TODAY_DATE

    cwd      = Path.cwd()
    pdb_file = Path(struct)


    if not is_cytosolic:

        if not spanfile is None:

            span_file = Path(spanfile)
        else:

            span_file = pdb_file.parent/('{}.span'.format(pdb_file.stem))


        print(span_file)

        print(span_file.exists())




    mutants_list = parse_mutants(mutants)


    struct_var_dir_name = '{}_{}_{}A'.format(pdb_file.stem, '-'.join(mutants_list), radius)


    if inplace:
        sim_output_path = cwd/pdb_file.stem/struct_var_dir_name/TODAY_DATE

    else:
        sim_output_path = output_root/pdb_file.stem/struct_var_dir_name/TODAY_DATE


    print(sim_output_path.name)

    while sim_output_path.exists():

            if not '_' in str(sim_output_path.name):
                
                sim_output_path = Path(str(sim_output_path) + '_1')

            else:

                n = int(str(sim_output_path)[-1])

                sim_output_path = Path(str(sim_output_path)[:-1] + str(n+1))


    sim_out_input_path = sim_output_path/'input'


    print('Creating sim directory: {}'.format(sim_output_path))
    print()


    sim_output_path.mkdir(parents=True)
    sim_out_input_path.mkdir()

    if not pdb_file.exists():

        if pdb_file.is_absolute() or not pdb_file.parent == Path('.'):
            sys.exit('{} does not exist.'.format(pdb_file))

        else:

            root_pdb_path  = input_root/pdb_file.name
            root_span_path = input_root/'{}.span'.format(pdb_file.stem)

            print('WARNING: ./{} does not exist.'.format(pdb_file))
            print('Using {} instead.'.format(root_pdb_path))
            print('')

            root_pdb_path.copy(sim_out_input_path/pdb_file.name)
            root_span_path.copy(sim_out_input_path/'{}.span'.format(pdb_file.stem))



    else:

        root_pdb_path = pdb_file

        pdb_file.copy(sim_out_input_path/pdb_file.name)

        if not is_cytosolic:

            span_file.copy(sim_out_input_path/'{}.span'.format(pdb_file.stem))





    if not resfile is None:

        resfile_path      = Path(resfile)
        dest_resfile_path = sim_out_input_path/'{}.resfile'.format(pdb_file.stem)

        resfile_path.copy(dest_resfile_path)

    elif not msdresfile is None:

        convert_msd2ism_res(resfile=msdresfile, 
                            struct_name=pdb_file.stem, 
                            length=len(pd_struct_dict[str(root_pdb_path)].get_sequence(1)), 
                            destination=sim_out_input_path
                           )

    else:

        create_resfile(mutants=mutants, 
                       nataa_list=get_nataa_list(mutants, radius, pd_struct_dict), 
                       radius=radius, 
                       length=len(pd_struct_dict[str(root_pdb_path)].get_sequence(1)), 
                       struct_name=pdb_file.stem, 
                       destination=sim_out_input_path,
                       natro=natro,
                      )

        dest_resfile_path = sim_out_input_path/'{}_{}_{}A.resfile'.format(pdb_file.stem, '-'.join(mutants_list), radius)




    create_paths_file(destination=sim_out_input_path, 
                      files_path=sim_out_input_path, 
                      db_path=rosetta_db_path, 
                      alt_files_path=sim_out_input_path,
                     )


    create_bash_file(structure_path=sim_out_input_path/pdb_file.name, 
                     destination=sim_output_path, 
                     nstruct=nstruct, 
                     mutants=mutants, 
                     radius=radius, 
                     hours=hours, 
                     is_cytosolic=is_cytosolic,
                     queue=queue
                    )



def get_nataa_list(mutants, radius, pd_struct_dict):
    '''
    test
    '''

    mutants_list    = parse_mutants(mutants)

    mut_pos_list    = [int(re.search('(\d+)', x).group(1)) for x in mutants_list]

    temp_nataa_list = [x.get_adj_residues(y, radius) 
                        for x in pd_struct_dict.values() 
                        for y in mut_pos_list]

    temp_cons_list = list()

    for adj_list in temp_nataa_list:
        temp_cons_list += adj_list

    nataa_list = set(temp_cons_list)
    print(nataa_list)
    return list(nataa_list)



def convert_msd2ism_res(resfile, struct_name, length, destination=OUTPUT_ROOT_PATH):

    res_path = Path(resfile)

    if not res_path.exists():
        sys.exit('{} does not exist.'.format(res_path))


    nataa_list = list()

    for line in res_path.read_text().split('\n'):

        result = re.search('(\d+)\s+\w\s(\w+)', line)

            
        if not result is None:
            
            nataa_list.append(int(result.group(1)))




    resfile    = ' start\n'

    natro_line = '   {:>4} {:>4} NATRO\n'
    nataa_line = '   {:>4} {:>4} NATAA\n'


    for n in range(1,length+1):
        
        if n in nataa_list:
            resfile += nataa_line.format(n, n)
            
        else:
            resfile += natro_line.format(n, n)


    resfile_name     = '{}.resfile'.format(struct_name)


    new_resfile_path = Path(destination/resfile_name)


    new_resfile_path.write_text(resfile)



def create_resfile(mutants, nataa_list, radius, length, struct_name, destination=OUTPUT_ROOT_PATH, natro=False):
    
    mutants_list    = parse_mutants(mutants)
    mutant_pos_list = [int(re.search('(\d+)', x).group(1)) for x in mutants_list]

    #print(mutants_list)

    resfile    = ' start\n'

    natro_line = '   {:>4} {:>4} NATRO\n'
    nataa_line = '   {:>4} {:>4} NATAA\n'
    pikaa_line = '   {:>4} {:>4} PIKAA  {}\n'

    #print(mutants_list)


    if not natro:
        for pos in range(1, length+1):
            if pos in mutant_pos_list:

                to_res   = re.search('(\d+)(\w+)', mutants_list[mutant_pos_list.index(pos)]).group(2)

                resfile += pikaa_line.format(pos, pos, to_res)

            elif pos in nataa_list:
                resfile += nataa_line.format(pos, pos)

            else:
                resfile += natro_line.format(pos, pos)


        resfile_name = '{}_{}_{}A.resfile'.format(struct_name, 
                                                  '-'.join(mutants_list), 
                                                  radius)

    else:

        for pos in range(1, length+1):

            resfile += natro_line.format(pos, pos)



        resfile_name = '{}_{}_{}A.resfile'.format(struct_name, 
                                                  '-'.join(mutants_list), 
                                                  radius)

    

    resfile_path = Path(destination/resfile_name)

    if not resfile_path.exists():
        resfile_path.write_text(resfile)

    else:
        sys.exit('resfile already exists in {}'.format(destination))



def create_paths_file(destination, 
                      files_path, 
                      alt_files_path,
                      db_path=ROSETTA_DB_PATH,
                      #overwrite=False,
                     ):

    paths_path            = destination / 'paths.txt'
    modded_db_path        = Path(db_path)
    modded_files_path     = Path(files_path)
    modded_alt_files_path = Path(alt_files_path)



    if not (destination).exists():
        sys.exit('{} does not exist.'.format(paths_path))

    #if not overwrite and (paths_path).exists():
    if (paths_path).exists():
        sys.exit('{} already exists.'.format(paths_path))
    

    if not (modded_db_path).exists():
        sys.exit('{} does not exist.'.format(modded_db_path))

    if not (modded_files_path).exists():
        sys.exit('{} does not exist.'.format(modded_files_path))

    if not (modded_alt_files_path).exists():
        sys.exit('{} does not exist.'.format(modded_alt_files_path))


    paths_file = ["Rosetta Input/Output Paths (order essential)",
                  "path is first '/', './',or  '../' to next whitespace, must end with '/'",
                  "INPUT PATHS:",
                  "pdb1                            {}/".format(modded_files_path),
                  "pdb2                            {}/".format(modded_alt_files_path),
                  "membrane spanning regions       {}/".format(modded_files_path),
                  "alternate data files            {}/".format(modded_db_path),
                  "fragments                       {}/".format(modded_files_path),
                  "structure dssp,ssa (dat,jones)  {}/".format(modded_files_path),
                  "sequence fasta,dat,jones        {}/".format(modded_files_path),
                  "constraints                     {}/".format(modded_files_path),
                  "starting structure              {}/".format(modded_files_path),
                  "data files                      {}/".format(modded_db_path),
                  "OUTPUT PATHS:",
                  "movie                           ./",
                  "pdb path                        ./",
                  "score                           ./",
                  "status                          ./",
                  "user                            ./",
                  "FRAGMENTS: (use '*****' in place of pdb name and chain)",
                  "2                                      number of valid fragment files",
                  "3                                      frag file 1 size",
                  "aa*****03_05.200_v1_3                                name",
                  "9                                      frag file 2 size",
                  "aa*****09_05.200_v1_3                                name",
                 ]


    paths_path.write_text('\n'.join(paths_file))



def create_bash_file(structure_path, 
                     destination, 
                     nstruct, 
                     mutants, 
                     radius, 
                     hours, 
                     is_cytosolic=False,
                     queue=True
                    ):

    global ROSETTA_PATH

    pikaa         = parse_mutants(mutants)

    pdb_filename  = structure_path.name
    pdb_name      = structure_path.stem
 
    res_name      = '{}_{}_{}A.resfile'.format(pdb_name,'-'.join(pikaa), radius)
    span_name     = '{}.span'.format(pdb_name)
 
    res_path      = None
    span_path     = None
    paths_path    = None
 
    bash_file     = None
    rosetta_flags = None


    #print(destination, res_name)


    #START CHECK RES FILE LOCATION
    if (destination / res_name).exists():
        res_path = destination / res_name

    elif (destination / 'input' / res_name).exists():
        res_path = destination / 'input' / res_name

    elif (destination / '{}.resfile'.format(pdb_name)).exists():
        res_path = destination / '{}.resfile'.format(pdb_name)

    elif (destination /'input/{}.resfile'.format(pdb_name)).exists():
        res_path = destination / 'input/{}.resfile'.format(pdb_name)

    else:
        print(destination / res_name)
        print(destination / 'input' / res_name)
        sys.exit('res file not found.')
    #END CHECK RES FILE LOCATION


    #START CHECK SPAN FILE LOCATION
    if not is_cytosolic:
        if (destination / span_name).exists():
            span_path = destination / span_name

        elif (destination / 'input' / span_name).exists():
            span_path = destination / 'input/' / span_name

        else:
            sys.exit('span file not found.')
    #END CHECK SPAN FILE LOCATION


    #START CHECK PATHS FILE LOCATION
    if (destination / 'paths.txt').exists():
        paths_path = destination / 'paths.txt'

    elif (destination / 'input/paths.txt').exists():
        paths_path = destination / 'input/paths.txt'

    else:
        sys.exit('paths file not found.')
    #END CHECK PATHS FILE LOCATION

    bash_file =   ['#!/bin/bash',
                   '',
                   '#SBATCH --chdir {}'.format(destination),
                   '#SBATCH --job-name ISM',
#		   '#SBATCH --partition=serial',
                   '#SBATCH -n 1                    # number of cores',
                   '#SBATCH -t 0-{0:02d}:00              # time (D-HH:MM)'.format(hours),
                   '#SBATCH -o slurm.%N.%j.out      # STDOUT',
                   '#SBATCH -e slurm.%N.%j.err      # STDERR',
                   '',
                   '',
                   'touch STARTED'
                   ]


    if not is_cytosolic:
        rosetta_flags = [str(ROSETTA_PATH),
                         '-s {}'.format(pdb_filename),
                         '-nstruct {}'.format(nstruct),
                         '-resfile {}'.format(res_path),
                         '-pose1', '-cst_mode', '-cst_design', '-fa_input', 
                         '-pose_memb', '-memb_solv', '-memb_env', '-Wmbenv 0.55',
                         '-memb_hb', '-thickness 12.5', '-steepness 10', 
                         '-spanfile {}'.format(span_path),
                         '-cst_min', '-chi_move', '-bb_move', 
                         '-paths {}'.format(paths_path), 
                         '-ex1', '-ex1aro', 
                         '-ex2', '-ex3', '-ex4', '-extrachi_cutoff 0', 
                         '-seed_offset $$', 
                         #'-out:file:silent {}_{}_{}A.silent'.format(pdb_name,'-'.join(pikaa), radius),

                         #remove cst_min, bb_move for non-repacked alanine scanning
                        ]
    else:
        rosetta_flags = [str(ROSETTA_PATH),
                         '-s {}'.format(pdb_filename),
                         '-nstruct {}'.format(nstruct),
                         '-resfile {}'.format(res_path),
                         '-pose1', '-cst_mode', '-cst_design', '-fa_input', 
                         '-cst_min', '-chi_move', '-bb_move', 
                         '-paths {}'.format(paths_path), 
                         '-ex1', '-ex1aro', 
                         '-ex2', '-ex3', '-ex4', '-extrachi_cutoff 0', 
                         '-seed_offset $$', 
                         #'-out:file:silent {}_{}_{}A.silent'.format(pdb_name,'-'.join(pikaa), radius),

                         #remove cst_min, bb_move for non-repacked alanine scanning
                        ]

    bash_file.append(' '.join(rosetta_flags))


    #so I can update sim_db while other sims are still running
    bash_file.append('\n')
    bash_file.append('touch FINISHED')


    outfile_name = 'sim_{}_{}.sh'.format(pdb_name, '-'.join(pikaa))
    outfile_path = Path(destination / outfile_name)

    outfile_path.write_text('\n'.join(bash_file))

    if queue:

        add_to_cluster_queue(outfile_path)

    else:

        launch_on_cluster(outfile_path)



def parse_mutants(mut_string):

    mut_list     = mut_string.upper().split('-')
    mut_out_list = list()


    for mut in mut_list:
        results  = re.search('(\w)(\d+)(\w+)', mut)

        if results is None:
            sys.exit('Mutation {} is not properly formatted.'.format(mut))

        #if mut[0].upper() not in 'ACDEFGHIKLMNPQRSTVWYX':
        #    sys.exit('Mutation {} is invalid.'.format(mut))

        
        init_res = mut[0]
        position = results.group(2)
        final_res= results.group(3)

        #print(final_res, len(final_res))


        if len(final_res) == 1 and final_res not in 'ACDEFGHIKLMNPQRSTVWYX':
            sys.exit('Mutation {} is invalid.'.format(final_res))


        elif len(final_res) > 1:
    
            for res in final_res:
                if res not in 'ACDEFGHIKLMNPQRSTVWY':
                    sys.exit('Mutation {} is invalid.'.format(mut))



        mut_out_list.append(mut.upper())

        #print(mut_out_list)

    #print(sorted(mut_out_list, key=lambda x: int(x[1:-1])) )

    return sorted(mut_out_list, key=lambda x: int(re.search('(\d+)', x).group(1)))


def add_to_cluster_queue(path, queue_path=CLUSTER_QUEUE_PATH):
    
    file_to_add_path = Path(path)

    queue_file_path  = Path(queue_path)


    if not file_to_add_path.exists():
        sys.exit('{} does not exist, can\'t add to cluster queue.')

    
    if queue_file_path.exists():
        
        queue_file_text = queue_file_path.read_text()

    else:
        queue_file_text = ''

    queue_file_path.write_text(queue_file_text + 
                               str(file_to_add_path) + '\n')


def launch_on_cluster(path):
    
    file_to_run_path = Path(path)

    if not file_to_run_path.exists():
    
        sys.exit('{} does not exist, run on cluster.')

    
    
    command = 'sbatch {}'.format(str(file_to_run_path))

    command_process  = subprocess.run(command, stdout=subprocess.PIPE, shell=True)

    print(command_process.stdout.decode().strip(), )




#"monkey-patching" Path()
def _copy(self, target):
    import shutil
    assert self.is_file()
    # str() only there for Python < (3, 6)
    shutil.copy(str(self), str(target))  

Path.copy = _copy



if __name__ == '__main__':

    get_parameters()
