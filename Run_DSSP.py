#!/usr/bin/env python
# coding: utf-8

# In[8]:


###############################################
##Dmitry Sutormin, Konstantin Gilep, 2020##
##N-to-C-terminus assymetry in protein domain secondary structure elements composition##

#Takes representative PDB structures linked to PFAM families and domains (prepared with Pfam_to_pdb.py).
#Runs DSSP of defined regions of the representative structures and collects output:
#elements of secondary structure, Phi, and Psi angles.
###############################################

#######
#Packages to be imported.
#######

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


#######
#Variables to be defined.
#######

#Path to representative structures.
Representative_structures_path='/home/d.sutormin/Project/PFAM_representative_structures/Representative_pdb_domains.txt'

#Path for output data.
DSSP_output_path='/home/d.sutormin/Project/DSSP/DSSP_Representative_pdb_domains.txt'

#Path for lost structures.
Lost_structures_path='/home/d.sutormin/Project/DSSP/DSSP_Representative_pdb_domains_lost_structures.txt'


###########
#This function takes pdb_id, chain_id, start and stop of residues according to the PDB numeration.
#It finds file in the server pdb directory and if no file found with such pdb_id returns None. 
#Function iterates throughout residue IDs INCLUDING!!! start and stop.
#It can not deal with residues with hetero-flag and insertion code. 
#In this case and in case of absence of some residue returns None
###########

def SecStr(pdb_id, chain_id, start, stop):
    
    #Change pdb_id to lower cases - as in local pdb db. 
    pdb_id = pdb_id.lower()
    
    #Read pdb structure if it exists.
    p = PDBParser()
    try:
        structure = p.get_structure(pdb_id, f'/home/m.pak/pdb/pdb{pdb_id}.pdb')
    except FileNotFoundError:
        print(f'File not found, proceed...  {pdb_id}')
        return None, None, None
    model = structure[0]
    
    #Run DSSP.
    try:
        dssp = DSSP(model, f'/home/m.pak/pdb/pdb{pdb_id}.pdb')
    except:
        print(f'DSSP unable to process the structure {pdb_id}, proceed...')
        return None, None, None
    
    #Keep annotation of secondaty structure elements, Phi and Psi angles for defined region of structure.
    sec_str = ''
    phi_lst = []
    psi_lst = []
    
    #INCLUDES STOP!!!!
    for num in range(start, stop+1):
        try: 
            res_key = (chain_id, (' ', num, ' ')) #Can not deal with hetero-flag and insertion code
            res = dssp[res_key]
        except:
            print(f'{res_key} not found in {pdb_id}, proceed...')
            continue
        
        sec_str += res[2]
        phi_lst.append(res[4])
        psi_lst.append(res[5])
        
    return sec_str, phi_lst, psi_lst


###########
#Read file containing representative PDB structures linked to PFAM database.
#File is prepared using Pfam_to_pdb.py script.
###########

def read_representative_structures_list(pathin):
    filein=open(pathin, 'r')
    
    #Keep structures info in dictionary.
    Repr_structures_dict={}
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['PFAM_id']:
            PFAM_ID=line[0]
            PDB_id=line[1]
            Chain_id=line[2]
            Start=int(line[3])
            Stop=int(line[4])
            Resolution=float(line[6])
            Repr_structures_dict[PFAM_ID]=[PDB_id, Chain_id, Start, Stop, Resolution]
    
    filein.close()
    return Repr_structures_dict


###########
#Iterate over representative structures. Run DSSp and collect information on 
#secondary structure elements, phi and psi angles.
###########

def iterate_structures(Repr_structures_dict):
    
    DSSP_data_dict={}
    Dict_of_lost_structures={}
    for pfam_id, pdb_info in Repr_structures_dict.items():
        pdb_id=pdb_info[0]
        chain_id=pdb_info[1]
        start=pdb_info[2]
        stop=pdb_info[3]
        resolution=pdb_info[4]
        print(f'Now working with {pdb_id}:{chain_id}:{start}-{stop}')
        sec_str, phi_lst, psi_lst=SecStr(pdb_id, chain_id, start, stop)
        if sec_str is not None:
            DSSP_data_dict[pfam_id]=[pdb_id, chain_id, start, stop, resolution, sec_str, phi_lst, psi_lst]
        else:
            Dict_of_lost_structures[pfam_id]=[pdb_id, chain_id, start, stop, resolution, sec_str, phi_lst, psi_lst]
        
    return DSSP_data_dict, Dict_of_lost_structures


###########
#Write DSSP data.
###########

def write_dssp_data(DSSP_data_dict, pathout):
    
    fileout=open(pathout, 'w')
    fileout.write('PFAM_id\tPDB_id\tChain\tStart\tEnd\tResolution\tSecondary_structures\tPhi_list\tPsi_list\n')
    for pfam_id, dssp_info in DSSP_data_dict.items():
        fileout.write(f'{pfam_id}\t{dssp_info[0]}\t{dssp_info[1]}\t{dssp_info[2]}\t{dssp_info[3]}\t{dssp_info[4]}\t{dssp_info[5]}\t{dssp_info[6]}\t{dssp_info[7]}\n')
    
    fileout.close()
    return


###########
#Wrapper function.
###########

def wrapper(PFAM_rep_structures_inpath, DSSP_rep_structures_outpath, Lost_structures_outpath):
    #Take representative structures.
    Repr_structures_dict=read_representative_structures_list(PFAM_rep_structures_inpath)
    #Iterate over representative structures, run DSSP, collect data: secondary structures, Phi, Psi.
    DSSP_data_dict, Dict_of_lost_structures=iterate_structures(Repr_structures_dict)
    #Write DSSP-generated data.
    write_dssp_data(DSSP_data_dict, DSSP_rep_structures_outpath)
    #Write lost structures.
    write_dssp_data(Dict_of_lost_structures, Lost_structures_outpath)
    
    return

wrapper(Representative_structures_path,DSSP_output_path, Lost_structures_path)

