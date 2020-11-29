###############################################
##Dmitry Sutormin, 2020##
##N-to-C-terminus assymetry in protein domain secondary structure elements composition##

#Pareses pfam file connecting pfam domains with pdb entities (pdb_pfamA_reg).
#Computes distribution of length of structure elements to pick up the most complete pdb structure, defined as representative structure.
#If there are several structures satisfying the definition, parses pdb file from pfam db to select structure with the highest resolution. The structure is now called representative structure.
#Returns a list of pdb id of representative structures for pfam domains with sequence coordinates.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import matplotlib.pyplot as plt

#######
#Variables to be defined.
#######

#Path to pdb_pfamA_reg file.
pdb_pfamA_reg="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Documents\Skoltech_PhD\Structural_bioinformatics\Project\Databases\pdb_pfamA_reg.txt"


#Path to pdb file.
pdb="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Documents\Skoltech_PhD\Structural_bioinformatics\Project\Databases\pdb.txt"


#Path to output file.
Representative_pdb_domains_info="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Documents\Skoltech_PhD\Structural_bioinformatics\Project\Representative_pdb_domains\Representative_pdb_domains.txt"


#######
#Filter representative structures by length.
#######

def representative_by_length(family_ar):
    #Sort by domain length.
    family_ar_sorted=sorted(family_ar, key=lambda x: x[4], reverse=True)
    #Take all the longest domains of the same length.
    representative_length=family_ar_sorted[0][4]
    representative_structures_length_ar=[]
    for structure in family_ar_sorted:
        if structure[4]==representative_length:
            representative_structures_length_ar.append(structure)    
    return representative_structures_length_ar


#######
#Filter representative structures by resolution.
#######

def representative_by_resolution(representative_structures_length_ar, pdb_structures_dict):
    #Sort by structure resolution.
    representative_structures_resolution_sorted=sorted(representative_structures_length_ar, key=lambda x: pdb_structures_dict[x[0]], reverse=False)
    #Take structure with the highest resolution.
    representative_structure=representative_structures_resolution_sorted[0]
    representative_structure_resolution=pdb_structures_dict[representative_structure[0]]
    #Add resolution to structure info.
    representative_structure_info=representative_structure+[representative_structure_resolution]
    return representative_structure_info


#######
#Reads pdb file with description of pfam-related pdb structures.
#######

def read_pdb(pathin):
    #Read pdb file line by line, make a dictionary of structures.
    filein=open(pathin, 'r')
    
    pdb_structures_dict={}
    for line in filein:
        line=line.rstrip().split('\t')    
        pdb_id=line[0]
        #If resolution is not defined (Solution NMR structure), set resolution to infinity.
        if line[4]=='NULL':
            pdb_resolution=100
        else:
            pdb_resolution=float(line[4])
        pdb_structures_dict[pdb_id]=pdb_resolution
    
    filein.close()
    return pdb_structures_dict


#######
#Reads pdb_pfamA_reg, constructs a dictionary of pfam domains.
#Some infographics of pfam database.
#######

def read_pdb_pfamA_reg(pathin):
    #Read pdb_pfamA_reg file line by line, make a dictionary of pfam families.
    filein=open(pathin, 'r')
    
    pfam_family_dict={}
    for line in filein:
        line=line.rstrip().split('\t')
        pdb_id=line[2]
        pfam_family_id=line[3]
        chain_name=line[5]
        pdb_domain_start=int(line[6])
        pdb_domain_end=int(line[8])
        #Create new family.
        if pfam_family_id not in pfam_family_dict:
            pfam_family_dict[pfam_family_id]=[[pdb_id, chain_name, pdb_domain_start, pdb_domain_end, pdb_domain_end-pdb_domain_start]]
        #Add pdb entry to existing family.
        else:
            pfam_family_dict[pfam_family_id].append([pdb_id, chain_name, pdb_domain_start, pdb_domain_end, pdb_domain_end-pdb_domain_start])
            
    print(f'Number of pfam families: {len(pfam_family_dict)}')
    
    filein.close()
    return pfam_family_dict


#######
#Write representative structures list.
#######

def write_representative_structures(Representative_structures, pathout):
    #Write file with information about representative pdb structures, selected by length and resolution.
    fileout=open(pathout, 'w')
    
    fileout.write('PFAM_id\tPDB_id\tChain\tStart\tEnd\tLength\tResolution\n')
    
    for family_name, representative_structure in Representative_structures.items():
        fileout.write(f'{family_name}\t{representative_structure[0]}\t{representative_structure[1]}\t{representative_structure[2]}\t{representative_structure[3]}\t{representative_structure[4]}\t{representative_structure[5]}\n')
    
    fileout.close()
    return


#######
#Find representative structures.
#######

def find_representative_structures(pfam_family_dict, pdb_structures_dict):
    Number_of_structures_ar=[]
    
    Number_of_repr_structures_len_ar=[]
    
    Representative_structures={}
    
    Representative_structures_resolution=[]
    
    Representative_structures_length=[]
    
    for pfam_family_name, family_ar in pfam_family_dict.items():
        Number_of_structures_ar.append(len(family_ar))
        #Select representative structures by length (the longest).
        representative_structures_length_ar=representative_by_length(family_ar)
        Number_of_repr_structures_len_ar.append(len(representative_structures_length_ar))
        #Select representative structures by resolution (the highest resolution) from the longest structures selected.
        representative_structure=representative_by_resolution(representative_structures_length_ar, pdb_structures_dict)
        Representative_structures[pfam_family_name]=representative_structure
        Representative_structures_resolution.append(representative_structure[5])
        Representative_structures_length.append(representative_structure[4])
        
    print(f'Total number of pfam-associated structures: {sum(Number_of_structures_ar)}')
    print(f'Total number of pfam-associated representative (by length) structures: {sum(Number_of_repr_structures_len_ar)}')
    print(f'Total number of pfam-associated representative (by resolution) structures: {len(Representative_structures)}')
    
    Representative_structures_length_high=[length for length in Representative_structures_length if length >= 125]
    Representative_structures_length_low=[length for length in Representative_structures_length if length < 125]
    
    fig, plot=plt.subplots(1,4,figsize=(11,2.5), dpi=100)
    plot[0].hist(Number_of_structures_ar, rwidth=0.85, edgecolor='black', linewidth=0.5)
    plot[0].set_xlabel('Number of PDB structures')
    plot[0].set_ylabel('Number of PFAM families')
    plot[0].set_yscale('log')
    plot[0].set_title('Full PFAM')
    plot[0].annotate(f'Number of\nstructures\n{sum(Number_of_structures_ar)}', (0.4,0.6), xycoords='axes fraction')
    plot[1].hist(Number_of_repr_structures_len_ar, rwidth=0.85, color='#ffb566', edgecolor='black', linewidth=0.5)
    plot[1].set_xlabel('Number of PDB structures')
    plot[1].set_ylabel('Number of PFAM families') 
    plot[1].set_yscale('log')
    plot[1].set_title('Representative by length')
    plot[1].annotate(f'Number of\nstructures\n{sum(Number_of_repr_structures_len_ar)}', (0.4,0.6), xycoords='axes fraction')
    plot[2].hist(Representative_structures_resolution, rwidth=0.85, color='#ab658c', edgecolor='black', linewidth=0.5)
    plot[2].set_xlabel(r'Resolution, $\AA$')
    plot[2].set_ylabel('Number of structures') 
    plot[2].set_yscale('log')
    plot[2].set_title('Representative structures')
    plot[2].annotate(f'Number of\nstructures\n{len(Representative_structures)}', (0.4,0.6), xycoords='axes fraction')
    plot[3].hist(Representative_structures_length, bins=500, rwidth=0.85, color='#94fff1', edgecolor='black', linewidth=0.1)
    plot[3].axvline(x=125, ls='--', linewidth=0.7)
    plot[3].set_xlabel('Domain length, aa')
    plot[3].set_ylabel('Number of structures') 
    plot[3].set_yscale('log')
    plot[3].set_title('Representative structures')
    plot[3].annotate(f'Number of\nstructures\n{len(Representative_structures)}', (0.4,0.6), xycoords='axes fraction') 
    plot[3].annotate(f'>=125 : {len(Representative_structures_length_high)}', (0.4,0.4), xycoords='axes fraction') 
    plot[3].annotate(f'<125 : {len(Representative_structures_length_low)}', (0.4,0.3), xycoords='axes fraction')
    plt.tight_layout()
    plt.show()    
        
    return Representative_structures


#######
#Wrapper function.
#######

def wrapper(pdb_pfamA_reg_in, pdb_in, representative_structures_out):
    #Read input files: pdb_pfamA_reg and pdb.
    pfam_family_dict=read_pdb_pfamA_reg(pdb_pfamA_reg_in)
    pdb_structures_dict=read_pdb(pdb_in)
    
    #Find representative structures.
    Representative_structures=find_representative_structures(pfam_family_dict, pdb_structures_dict)
    #Write representative structures info.
    write_representative_structures(Representative_structures, representative_structures_out)
    
    return

wrapper(pdb_pfamA_reg, pdb, Representative_pdb_domains_info)