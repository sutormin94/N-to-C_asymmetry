###############################################
##Dmitry Sutormin, 2020##
##N-to-C-terminus assymetry in protein domain secondary structure elements composition##

#Pareses DSSP output generated for representative structures with Run_DSSP.py
#
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
from sklearn.datasets.samples_generator import make_blobs

#######
#Variables to be defined.
#######

#Path to DSSP data file.
DSSP_data_inpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Documents\Skoltech_PhD\Structural_bioinformatics\Project\DSSP\DSSP\DSSP_Representative_pdb_domains.txt"


#Path to output file.
DSSP_data_outpath="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Documents\Skoltech_PhD\Structural_bioinformatics\Project\DSSP\\"


#######
#Read DSSP output.
#######

def read_dssp_data(DSSP_inpath):
    
    #Read DSSP data keep in dictionary.
    DSSP_data_dict={}
    filein=open(DSSP_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')  
        if line[0] not in ['PFAM_id']:
            PFAM_id=line[0]
            PDB_id=line[1]
            Chain_id=line[2]
            Start=int(line[3])
            End=int(line[4])
            Resolution=float(line[5])
            SSE_string=line[6]
            Phi_list=[float(x) for x in line[7].lstrip('[').rstrip(']').split(', ')]
            Psi_list=[float(x) for x in line[8].lstrip('[').rstrip(']').split(', ')]
            DSSP_data_dict[PFAM_id]=[PDB_id, Chain_id, Start, End, Resolution, SSE_string, Phi_list, Psi_list]
    
    filein.close()
    return DSSP_data_dict


#######
#Filter data, split data by domain length.
#######

def define_length_groups(DSSP_data_dict, min_len, thr_len):
    
    Domain_length_ar=[]
    Discarded_short_structures={}
    Short_structures={}
    Long_structures={}
    for pfam_id, dssp_data in DSSP_data_dict.items():
        domain_len=len(dssp_data[5])
        Domain_length_ar.append(domain_len)
        
        #Classify domains by length.
        if domain_len<min_len:
            Discarded_short_structures[pfam_id]=dssp_data
        elif min_len<=domain_len<thr_len:
            Short_structures[pfam_id]=dssp_data
        elif domain_len>=thr_len:
            Long_structures[pfam_id]=dssp_data    
    
    #Plot domain length distribution.
    fig, plot=plt.subplots(1,1,figsize=(3,3), dpi=100)
    plot.hist(Domain_length_ar, bins=50, rwidth=0.85, color='#94fff1', edgecolor='black', linewidth=0.1)
    plot.axvline(x=min_len, ls='--', linewidth=0.7, color='black')
    plot.axvline(x=thr_len, ls='--', linewidth=0.7)
    plot.set_xlabel('Domain length, aa')
    plot.set_ylabel('Number of structures') 
    plot.set_yscale('log')
    plot.set_title('Representative structures\nafter DSSP')
    plot.annotate(f'Number of\nstructures\n{len(Domain_length_ar)}', (0.5,0.7), xycoords='axes fraction') 
    plot.annotate(f'x<{min_len} : {len(Discarded_short_structures)}', (0.5,0.6), xycoords='axes fraction', size=7)
    plot.annotate(f'{min_len}<=x<{thr_len} : {len(Short_structures)}', (0.5,0.5), xycoords='axes fraction', size=7)
    plot.annotate(f'x>={thr_len} : {len(Long_structures)}', (0.5,0.4), xycoords='axes fraction', size=7) 
    
    plt.tight_layout()
    plt.show()
    
    return Short_structures, Long_structures


#######
#Take phi, psi angles for N- and C-termini.
#######

def phi_psi_N_to_C(structures_dict, window_width):
        
    phi_N=[]
    phi_C=[]
    psi_N=[]
    psi_C=[] 
    phi=[]
    psi=[]
    for pfam_id, dssp_data in structures_dict.items():
        phi_list=dssp_data[6]
        psi_list=dssp_data[7]
        phi_N+=phi_list[:window_width]
        phi_C+=phi_list[-window_width:]
        psi_N+=psi_list[:window_width]
        psi_C+=psi_list[-window_width:] 
        phi+=phi_list
        psi+=psi_list
        
    return phi_N, phi_C, psi_N, psi_C, phi, psi


#######
#Create secondary structure element frequency matrix.
#######

def ss_element_frequency_matrix(structures_dict, window_width):
    
    #Create positioned matrix of secondary structure elements.
    ss_matrix_N=[]
    ss_matrix_C=[]
    for i in range(window_width):
        column_N=[]
        column_C=[]
        for pfam_id, dssp_data in structures_dict.items():
            SSE_string=dssp_data[5]
            column_N.append(SSE_string[i])
            column_C.append(SSE_string[-window_width+i])
        ss_matrix_N.append(column_N)
        ss_matrix_C.append(column_C)
            
    #Create position frequency matrix of secondary structure elements.
    DSSP_alphabet=['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']
    ss_pfm_N={}
    ss_pfm_C={}
    ss_pfm_conf_N={}
    ss_pfm_conf_C={}
    for letter_code in DSSP_alphabet:
        #Keep frequences of ss elements.
        frequency_row_N=[]
        frequency_row_C=[]
        #Keep boundaries of confident interval.
        conf_upper_N=[]
        conf_upper_C=[]
        conf_lower_N=[]
        conf_lower_C=[]        
        for i in range(len(ss_matrix_N)):
            column_letter_freq_N=ss_matrix_N[i].count(letter_code)/float(len(ss_matrix_N[i]))
            confident_interval_N=st.binom.interval(0.95, len(ss_matrix_N[i]), column_letter_freq_N, loc=0)
            lower_N=confident_interval_N[0]/len(ss_matrix_N[i])
            upper_N=confident_interval_N[1]/len(ss_matrix_N[i])
            conf_lower_N.append(lower_N)
            conf_upper_N.append(upper_N)
            frequency_row_N.append(column_letter_freq_N)
            
            column_letter_freq_C=ss_matrix_C[i].count(letter_code)/float(len(ss_matrix_C[i]))
            confident_interval_C=st.binom.interval(0.95, len(ss_matrix_C[i]), column_letter_freq_C, loc=0)
            lower_C=confident_interval_C[0]/len(ss_matrix_C[i])
            upper_C=confident_interval_C[1]/len(ss_matrix_C[i])
            conf_lower_C.append(lower_C)
            conf_upper_C.append(upper_C)
            frequency_row_C.append(column_letter_freq_C)  
          
        ss_pfm_N[letter_code]=frequency_row_N
        ss_pfm_C[letter_code]=frequency_row_C
        
        ss_pfm_conf_N[letter_code+'_upper']=conf_upper_N
        ss_pfm_conf_N[letter_code+'_lower']=conf_lower_N
        ss_pfm_conf_C[letter_code+'_upper']=conf_upper_C
        ss_pfm_conf_C[letter_code+'_lower']=conf_lower_C        
       
    print(ss_pfm_N) 
    print(ss_pfm_C)   
    return ss_matrix_N, ss_matrix_C, ss_pfm_N, ss_pfm_C, ss_pfm_conf_N, ss_pfm_conf_C


#######
#Enrichment of secondary structure elements N- over C-terminus.
#######

def ss_ele_enrichment(ss_pfm_N, ss_pfm_C):
    
    enrichment_N_to_C_dict={}
    for letter_code, frequency_row_N in ss_pfm_N.items():
        frequency_N=np.array(frequency_row_N)
        frequency_C=np.array(ss_pfm_C[letter_code][::-1])
        enrichment_N_to_C=frequency_N/frequency_C
        enrichment_N_to_C_dict[letter_code]=enrichment_N_to_C
    
    return enrichment_N_to_C_dict


#######
#Terminal beta-strands relations: simultaneous occurence or independant?.
#######

def termini_dependance(ss_matrix_N, ss_matrix_C, ss_pfm_N, ss_pfm_C, local_window_width):
    
    DSSP_alphabet=['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']
    
    #Calculate observed frequences of elements co-occurence coordinate-wise.
    Observed_matrices_dict={}
    for letter_code_1 in DSSP_alphabet:
        for letter_code_2 in DSSP_alphabet:
            Ni_Cj_freq_2d_ar=[]
            for Ni in range(local_window_width):
                Cj_freq_1d_ar=[]
                for Cj in range(local_window_width):
                    N_column_i=ss_matrix_N[Ni]
                    C_column_j=ss_matrix_C[-Cj-1]
                    Ni_Cj_counts=0
                    for structure_k in range(len(N_column_i)):
                        if (N_column_i[structure_k]==letter_code_1) and (C_column_j[structure_k]==letter_code_2):
                            Ni_Cj_counts+=1
                    Ni_Cj_frequency=Ni_Cj_counts/float(len(N_column_i))
                    Cj_freq_1d_ar.append(Ni_Cj_frequency)
            
                Ni_Cj_freq_2d_ar.append(Cj_freq_1d_ar)
            
            Ni_Cj_freq_2d_ar_np=np.array(Ni_Cj_freq_2d_ar)
            Observed_matrices_dict[letter_code_1+letter_code_2]=Ni_Cj_freq_2d_ar_np
        
    #Calculate expected frequences of elements co-occurence coordinate-wise.
    Expected_matrices_dict={}
    for letter_code_1 in DSSP_alphabet:
        for letter_code_2 in DSSP_alphabet:
            Expected_Ni_Cj_freq_2d_ar=[]
            for Ni in range(local_window_width):
                Expected_Cj_freq_1d_ar=[]
                for Cj in range(local_window_width):  
                    Ni_frequency=ss_pfm_N[letter_code_1][Ni]
                    Cj_frequency=ss_pfm_C[letter_code_2][-Cj-1]
                    Expected_Ni_Cj_frequency=Ni_frequency*Cj_frequency
                    Expected_Cj_freq_1d_ar.append(Expected_Ni_Cj_frequency)
                Expected_Ni_Cj_freq_2d_ar.append(Expected_Cj_freq_1d_ar)
                
            Expected_Ni_Cj_freq_2d_ar_np=np.array(Expected_Ni_Cj_freq_2d_ar)
            Expected_matrices_dict[letter_code_1+letter_code_2]=Expected_Ni_Cj_freq_2d_ar_np

    #Calculate observed over expected ratio of frequences of elements co-occurence coordinate-wise.
    Obs_over_exp_matrices_dict={}
    for letter_code_1 in DSSP_alphabet:
        for letter_code_2 in DSSP_alphabet:
            Obs_over_exp_freq_matrix=np.divide(Observed_matrices_dict[letter_code_1+letter_code_2], Expected_matrices_dict[letter_code_1+letter_code_2])
            Obs_over_exp_matrices_dict[letter_code_1+letter_code_2]=Obs_over_exp_freq_matrix
            
            print(letter_code_1, letter_code_2, Obs_over_exp_freq_matrix)
    
    return Obs_over_exp_matrices_dict


#######
#Analyse N-to-C asymmetry in secondary structure elements frequency.
#######

def distribution_of_frequences(ss_pfm_N, ss_pfm_C):
    ss_pfm=ss_pfm_N+ss_pfm_C
    print(len(ss_pfm))
    
    fig, plot=plt.subplots(1,1,figsize=(4,4), dpi=100)
    weights=np.ones_like(ss_pfm)/(len(ss_pfm)) #Taken from https://stackoverflow.com/questions/42481698/probability-density-histogram-with-matplotlib-doesnt-make-sense     
    plot.hist(ss_pfm, bins=10, weights=weights)
    plot.set_xlabel('Frequency of ss element')
    plot.set_ylabel('Fraction of positions')
    plt.show()

    return


#######
#Analyse N-to-C asymmetry in secondary structure elements frequency.
#######

def N_to_C_asymmetry(ss_pfm_N_short, ss_pfm_C_short, ss_pfm_conf_N_short, ss_pfm_conf_C_short, ss_pfm_N_long, ss_pfm_C_long, ss_pfm_conf_N_long, ss_pfm_conf_C_long, window_width):
    
    ##Plot distribution of frequences, get confidential interval.
    ##Distribution are far from being normal, estimation is not good -> depricated.
    #distribution_of_frequences(ss_pfm_N_short['H'], ss_pfm_C_short['H'])
    #distribution_of_frequences(ss_pfm_N_long['H'], ss_pfm_C_long['H'])
    
    #Plot frequency of ss elements as a function of distance from N- and C-termini.
    X_N=range(window_width)
    X_C=range(60, 60+window_width)
    xticks_ar=list(range(0,111,10))
    xticklabels_ar=[0,10,20,30,40,50,-50,-40,-30,-20,-10,0]
    xticks_spec=[-4, 114]
    xticklabels_spec=['N-', '-C']
    fig, plot=plt.subplots(2,2,figsize=(12,6), dpi=100)
    plot[0,0].plot(X_N, ss_pfm_N_short['H'], color='#94fff1', linewidth=2, label=r'$\alpha$-helix short domains')
    plot[0,0].fill_between(X_N, ss_pfm_conf_N_short['H_lower'], ss_pfm_conf_N_short['H_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    plot[0,0].plot(X_C, ss_pfm_C_short['H'], color='#94fff1', linewidth=2)
    plot[0,0].fill_between(X_C, ss_pfm_conf_C_short['H_lower'], ss_pfm_conf_C_short['H_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    
    plot[0,0].plot(X_N, ss_pfm_N_long['H'], color='#ab658c', linewidth=2, label=r'$\alpha$-helix long domains')
    plot[0,0].fill_between(X_N, ss_pfm_conf_N_long['H_lower'], ss_pfm_conf_N_long['H_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    plot[0,0].plot(X_C, ss_pfm_C_long['H'], color='#ab658c', linewidth=2) 
    plot[0,0].fill_between(X_C, ss_pfm_conf_C_long['H_lower'], ss_pfm_conf_C_long['H_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    
    plot[0,0].plot(X_N, ss_pfm_N_short['E'], color='#59acff', linewidth=2, label=r'$\beta$-strand short domains')
    plot[0,0].fill_between(X_N, ss_pfm_conf_N_short['E_lower'], ss_pfm_conf_N_short['E_upper'], color='#59acff', alpha=0.3, linewidth=0)
    plot[0,0].plot(X_C, ss_pfm_C_short['E'], color='#59acff', linewidth=2)
    plot[0,0].fill_between(X_C, ss_pfm_conf_C_short['E_lower'], ss_pfm_conf_C_short['E_upper'], color='#59acff', alpha=0.3, linewidth=0)
    
    plot[0,0].plot(X_N, ss_pfm_N_long['E'], color='#ffec59', linewidth=2, label=r'$\beta$-strand long domains')
    plot[0,0].fill_between(X_N, ss_pfm_conf_N_long['E_lower'], ss_pfm_conf_N_long['E_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    plot[0,0].plot(X_C, ss_pfm_C_long['E'], color='#ffec59', linewidth=2)    
    plot[0,0].fill_between(X_C, ss_pfm_conf_C_long['E_lower'], ss_pfm_conf_C_long['E_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    
    plot[0,0].set_xticks(xticks_ar)
    plot[0,0].set_xticklabels(xticklabels_ar)
    plot[0,0].set_xticks(xticks_spec, minor=True)
    plot[0,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,0].tick_params(axis='x', which='minor', length=0)
    plot[0,0].set_xlabel('Distance, aa')
    plot[0,0].set_ylabel('Frequency') 
    plot[0,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')
    
    plot[0,1].plot(X_N, ss_pfm_N_short['B'], color='#94fff1', linewidth=2, label=r'$\beta$-bridge short domains')
    plot[0,1].fill_between(X_N, ss_pfm_conf_N_short['B_lower'], ss_pfm_conf_N_short['B_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    plot[0,1].plot(X_C, ss_pfm_C_short['B'], color='#94fff1', linewidth=2)
    plot[0,1].fill_between(X_C, ss_pfm_conf_C_short['B_lower'], ss_pfm_conf_C_short['B_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    
    plot[0,1].plot(X_N, ss_pfm_N_long['B'], color='#ab658c', linewidth=2, label=r'$\beta$-bridge long domains')
    plot[0,1].fill_between(X_N, ss_pfm_conf_N_long['B_lower'], ss_pfm_conf_N_long['B_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    plot[0,1].plot(X_C, ss_pfm_C_long['B'], color='#ab658c', linewidth=2) 
    plot[0,1].fill_between(X_C, ss_pfm_conf_C_long['B_lower'], ss_pfm_conf_C_long['B_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    
    plot[0,1].plot(X_N, ss_pfm_N_short['G'], color='#59acff', linewidth=2, label=r'3-10-helix short domains')
    plot[0,1].fill_between(X_N, ss_pfm_conf_N_short['G_lower'], ss_pfm_conf_N_short['G_upper'], color='#59acff', alpha=0.3, linewidth=0)
    plot[0,1].plot(X_C, ss_pfm_C_short['G'], color='#59acff', linewidth=2)
    plot[0,1].fill_between(X_C, ss_pfm_conf_C_short['G_lower'], ss_pfm_conf_C_short['G_upper'], color='#59acff', alpha=0.3, linewidth=0)
    
    plot[0,1].plot(X_N, ss_pfm_N_long['G'], color='#ffec59', linewidth=2, label=r'3-10-helix long domains')
    plot[0,1].fill_between(X_N, ss_pfm_conf_N_long['G_lower'], ss_pfm_conf_N_long['G_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    plot[0,1].plot(X_C, ss_pfm_C_long['G'], color='#ffec59', linewidth=2)    
    plot[0,1].fill_between(X_C, ss_pfm_conf_C_long['G_lower'], ss_pfm_conf_C_long['G_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    
    plot[0,1].set_xticks(xticks_ar)
    plot[0,1].set_xticklabels(xticklabels_ar)
    plot[0,1].set_xticks(xticks_spec, minor=True)
    plot[0,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,1].tick_params(axis='x', which='minor', length=0)
    plot[0,1].set_xlabel('Distance, aa')
    plot[0,1].set_ylabel('Frequency')  
    plot[0,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center')
    
    plot[1,0].plot(X_N, ss_pfm_N_short['I'], color='#94fff1', linewidth=2, label=r'$\pi$-helix short domains')
    plot[1,0].fill_between(X_N, ss_pfm_conf_N_short['I_lower'], ss_pfm_conf_N_short['I_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    plot[1,0].plot(X_C, ss_pfm_C_short['I'], color='#94fff1', linewidth=2)
    plot[1,0].fill_between(X_C, ss_pfm_conf_C_short['I_lower'], ss_pfm_conf_C_short['I_upper'], color='#94fff1', alpha=0.3, linewidth=0)

    plot[1,0].plot(X_N, ss_pfm_N_long['I'], color='#ab658c', linewidth=2, label=r'$\pi$-helix long domains')
    plot[1,0].fill_between(X_N, ss_pfm_conf_N_long['I_lower'], ss_pfm_conf_N_long['I_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    plot[1,0].plot(X_C, ss_pfm_C_long['I'], color='#ab658c', linewidth=2) 
    plot[1,0].fill_between(X_C, ss_pfm_conf_C_long['I_lower'], ss_pfm_conf_C_long['I_upper'], color='#ab658c', alpha=0.3, linewidth=0)

    plot[1,0].plot(X_N, ss_pfm_N_short['T'], color='#59acff', linewidth=2, label=r'turn short domains')
    plot[1,0].fill_between(X_N, ss_pfm_conf_N_short['T_lower'], ss_pfm_conf_N_short['T_upper'], color='#59acff', alpha=0.3, linewidth=0)
    plot[1,0].plot(X_C, ss_pfm_C_short['T'], color='#59acff', linewidth=2)
    plot[1,0].fill_between(X_C, ss_pfm_conf_C_short['T_lower'], ss_pfm_conf_C_short['T_upper'], color='#59acff', alpha=0.3, linewidth=0)

    plot[1,0].plot(X_N, ss_pfm_N_long['T'], color='#ffec59', linewidth=2, label=r'turn long domains')
    plot[1,0].fill_between(X_N, ss_pfm_conf_N_long['T_lower'], ss_pfm_conf_N_long['T_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    plot[1,0].plot(X_C, ss_pfm_C_long['T'], color='#ffec59', linewidth=2)  
    plot[1,0].fill_between(X_C, ss_pfm_conf_C_long['T_lower'], ss_pfm_conf_C_long['T_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    
    plot[1,0].set_xticks(xticks_ar)
    plot[1,0].set_xticklabels(xticklabels_ar)
    plot[1,0].set_xticks(xticks_spec, minor=True)
    plot[1,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,0].tick_params(axis='x', which='minor', length=0)
    plot[1,0].set_xlabel('Distance, aa')
    plot[1,0].set_ylabel('Frequency')  
    plot[1,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center')
    
    plot[1,1].plot(X_N, ss_pfm_N_short['S'], color='#94fff1', linewidth=2, label=r'bend short domains')
    plot[1,1].fill_between(X_N, ss_pfm_conf_N_short['S_lower'], ss_pfm_conf_N_short['S_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    plot[1,1].plot(X_C, ss_pfm_C_short['S'], color='#94fff1', linewidth=2)
    plot[1,1].fill_between(X_C, ss_pfm_conf_C_short['S_lower'], ss_pfm_conf_C_short['S_upper'], color='#94fff1', alpha=0.3, linewidth=0)
    
    plot[1,1].plot(X_N, ss_pfm_N_long['S'], color='#ab658c', linewidth=2, label=r'bend long domains')
    plot[1,1].fill_between(X_N, ss_pfm_conf_N_long['S_lower'], ss_pfm_conf_N_long['S_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    plot[1,1].plot(X_C, ss_pfm_C_long['S'], color='#ab658c', linewidth=2) 
    plot[1,1].fill_between(X_C, ss_pfm_conf_C_long['S_lower'], ss_pfm_conf_C_long['S_upper'], color='#ab658c', alpha=0.3, linewidth=0)
    
    plot[1,1].plot(X_N, ss_pfm_N_short['-'], color='#59acff', linewidth=2, label=r'unstructured short domains')
    plot[1,1].fill_between(X_N, ss_pfm_conf_N_short['-_lower'], ss_pfm_conf_N_short['-_upper'], color='#59acff', alpha=0.3, linewidth=0)
    plot[1,1].plot(X_C, ss_pfm_C_short['-'], color='#59acff', linewidth=2)
    plot[1,1].fill_between(X_C, ss_pfm_conf_C_short['-_lower'], ss_pfm_conf_C_short['-_upper'], color='#59acff', alpha=0.3, linewidth=0)
        
    plot[1,1].plot(X_N, ss_pfm_N_long['-'], color='#ffec59', linewidth=2, label=r'unstructured long domains')
    plot[1,1].fill_between(X_N, ss_pfm_conf_N_long['-_lower'], ss_pfm_conf_N_long['-_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    plot[1,1].plot(X_C, ss_pfm_C_long['-'], color='#ffec59', linewidth=2)  
    plot[1,1].fill_between(X_C, ss_pfm_conf_C_long['-_lower'], ss_pfm_conf_C_long['-_upper'], color='#ffec59', alpha=0.3, linewidth=0)
    
    plot[1,1].set_xticks(xticks_ar)
    plot[1,1].set_xticklabels(xticklabels_ar)
    plot[1,1].set_xticks(xticks_spec, minor=True)
    plot[1,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,1].tick_params(axis='x', which='minor', length=0)
    plot[1,1].set_xlabel('Distance, aa')
    plot[1,1].set_ylabel('Frequency')  
    plot[1,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center')    
    
    
    plt.tight_layout()
    plt.show()    
    
    return


#######
#Enrichment of ss elements at N-terminus over the elements at C-terminus.
#######

def N_to_C_enrichment(enrichment_N_to_C_dict_short, enrichment_N_to_C_dict_long, window_width):
    
    #Plot enrichment of ss elements as a function of distance from termini.
    X=range(window_width)
    xticks_ar=list(range(0,51,10))
    xticklabels_ar=[0,10,20,30,40,50]
    xticks_spec=[-4]
    xticklabels_spec=['N-\nC-', ]
    fig, plot=plt.subplots(2,2,figsize=(12,6), dpi=100)
    plot[0,0].plot(X, enrichment_N_to_C_dict_short['H'], color='#94fff1', linewidth=2, label=r'$\alpha$-helix short domains')
    plot[0,0].plot(X, enrichment_N_to_C_dict_long['H'], color='#ab658c', linewidth=2, label=r'$\alpha$-helix long domains')
    plot[0,0].plot(X, enrichment_N_to_C_dict_short['E'], color='#59acff', linewidth=2, label=r'$\beta$-strand short domains')
    plot[0,0].plot(X, enrichment_N_to_C_dict_long['E'], color='#ffec59', linewidth=2, label=r'$\beta$-strand long domains')   
    plot[0,0].set_xticks(xticks_ar)
    plot[0,0].set_xticklabels(xticklabels_ar)
    plot[0,0].set_xticks(xticks_spec, minor=True)
    plot[0,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,0].tick_params(axis='x', which='minor', length=0)
    plot[0,0].set_xlabel('Distance, aa')
    plot[0,0].set_ylabel('Enrichment p(N)/p(C)') 
    plot[0,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')
    
    plot[0,1].plot(X, enrichment_N_to_C_dict_short['B'], color='#94fff1', linewidth=2, label=r'$\beta$-bridge short domains')
    plot[0,1].plot(X, enrichment_N_to_C_dict_long['B'], color='#ab658c', linewidth=2, label=r'$\beta$-bridge long domains')
    plot[0,1].plot(X, enrichment_N_to_C_dict_short['G'], color='#59acff', linewidth=2, label=r'3-10-helix short domains')
    plot[0,1].plot(X, enrichment_N_to_C_dict_long['G'], color='#ffec59', linewidth=2, label=r'3-10-helix long domains')   
    plot[0,1].set_xticks(xticks_ar)
    plot[0,1].set_xticklabels(xticklabels_ar)
    plot[0,1].set_xticks(xticks_spec, minor=True)
    plot[0,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,1].tick_params(axis='x', which='minor', length=0)
    plot[0,1].set_xlabel('Distance, aa')
    plot[0,1].set_ylabel('Enrichment p(N)/p(C)')  
    plot[0,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')
    
    plot[1,0].plot(X, enrichment_N_to_C_dict_short['I'], color='#94fff1', linewidth=2, label=r'$\pi$-helix short domains')
    plot[1,0].plot(X, enrichment_N_to_C_dict_long['I'], color='#ab658c', linewidth=2, label=r'$\pi$-helix long domains')
    plot[1,0].plot(X, enrichment_N_to_C_dict_short['T'], color='#59acff', linewidth=2, label=r'turn short domains')
    plot[1,0].plot(X, enrichment_N_to_C_dict_long['T'], color='#ffec59', linewidth=2, label=r'turn long domains')   
    plot[1,0].set_xticks(xticks_ar)
    plot[1,0].set_xticklabels(xticklabels_ar)
    plot[1,0].set_xticks(xticks_spec, minor=True)
    plot[1,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,0].tick_params(axis='x', which='minor', length=0)
    plot[1,0].set_xlabel('Distance, aa')
    plot[1,0].set_ylabel('Enrichment p(N)/p(C)')  
    plot[1,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')
    
    plot[1,1].plot(X, enrichment_N_to_C_dict_short['S'], color='#94fff1', linewidth=2, label=r'bend short domains')
    plot[1,1].plot(X, enrichment_N_to_C_dict_long['S'], color='#ab658c', linewidth=2, label=r'bend long domains') 
    plot[1,1].plot(X, enrichment_N_to_C_dict_short['-'], color='#59acff', linewidth=2, label=r'unstructured short domains')
    plot[1,1].plot(X, enrichment_N_to_C_dict_long['-'], color='#ffec59', linewidth=2, label=r'unstructured long domains')  
    plot[1,1].set_xticks(xticks_ar)
    plot[1,1].set_xticklabels(xticklabels_ar)
    plot[1,1].set_xticks(xticks_spec, minor=True)
    plot[1,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,1].tick_params(axis='x', which='minor', length=0)
    plot[1,1].set_xlabel('Distance, aa')
    plot[1,1].set_ylabel('Enrichment p(N)/p(C)')  
    plot[1,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')    
    
    
    plt.tight_layout()
    plt.show()        
    
    return

#######
#Plot co-occurence of ss elements at termini.
#######

def plot_co_occurence(Obs_over_exp_matrices_short, Obs_over_exp_matrices_long, local_window_width):
    
    fig, plot=plt.subplots(2,2,figsize=(7,7), dpi=100)
    ticks_ar=list(range(-1, local_window_width, int(local_window_width/10)))
    ticks_ar[0]=0
    print(list(ticks_ar))
    ticklabels_ar=np.array(list(ticks_ar))+1
    plot00=plot[0,0].imshow(Obs_over_exp_matrices_short['EH'], cmap='gnuplot', vmin=0.3, vmax=2.1, interpolation='nearest')
    plot[0,0].set_title(r'$\beta-\alpha$ short domains')
    plo00_cbar=plot[0,0].figure.colorbar(plot00, ax=plot[0,0], shrink=0.7)
    plot[0,0].set_xticks(ticks_ar)
    plot[0,0].set_xticklabels(ticklabels_ar) 
    plot[0,0].set_yticks(ticks_ar)
    plot[0,0].set_yticklabels(ticklabels_ar)    
    plot[0,0].set_xlabel(r'Distance from C-terminus, aa ($\alpha$)')
    plot[0,0].set_ylabel(r'Distance from N-terminus, aa ($\beta$)')     
    plo00_cbar.ax.set_ylabel('p(Obs)/p(Exp)', rotation=-90, va="bottom")
    
    plot01=plot[0,1].imshow(Obs_over_exp_matrices_long['EH'], cmap='gnuplot', vmin=0.3, vmax=2.1, interpolation='nearest')
    plot[0,1].set_title(r'$\beta-\alpha$-helix long domains')
    plo01_cbar=plot[0,1].figure.colorbar(plot01, ax=plot[0,1], shrink=0.7)  
    plot[0,1].set_xticks(ticks_ar)
    plot[0,1].set_xticklabels(ticklabels_ar) 
    plot[0,1].set_yticks(ticks_ar)
    plot[0,1].set_yticklabels(ticklabels_ar)    
    plot[0,1].set_xlabel(r'Distance from C-terminus, aa ($\alpha$)')
    plot[0,1].set_ylabel(r'Distance from N-terminus, aa ($\beta$)')     
    plo01_cbar.ax.set_ylabel('p(Obs)/p(Exp)', rotation=-90, va="bottom")
    
    plot10=plot[1,0].imshow(Obs_over_exp_matrices_short['HE'], cmap='gnuplot', vmin=0.3, vmax=2.1, interpolation='nearest')
    plot[1,0].set_title(r'$\alpha-\beta$-strand short domains')
    plo10_cbar=plot[0,0].figure.colorbar(plot10, ax=plot[1,0], shrink=0.7)
    plot[1,0].set_xticks(ticks_ar)
    plot[1,0].set_xticklabels(ticklabels_ar) 
    plot[1,0].set_yticks(ticks_ar)
    plot[1,0].set_yticklabels(ticklabels_ar)    
    plot[1,0].set_xlabel(r'Distance from C-terminus, aa ($\beta$)')
    plot[1,0].set_ylabel(r'Distance from N-terminus, aa ($\alpha$)')     
    plo10_cbar.ax.set_ylabel('p(Obs)/p(Exp)', rotation=-90, va="bottom")
    
    plot11=plot[1,1].imshow(Obs_over_exp_matrices_long['HE'], cmap='gnuplot', vmin=0.3, vmax=2.1, interpolation='nearest')
    plot[1,1].set_title(r'$\alpha-\beta$-strand long domains')
    plo11_cbar=plot[1,1].figure.colorbar(plot11, ax=plot[1,1], shrink=0.7)  
    plot[1,1].set_xticks(ticks_ar)
    plot[1,1].set_xticklabels(ticklabels_ar) 
    plot[1,1].set_yticks(ticks_ar)
    plot[1,1].set_yticklabels(ticklabels_ar)    
    plot[1,1].set_xlabel(r'Distance from C-terminus, aa ($\beta$)')
    plot[1,1].set_ylabel(r'Distance from N-terminus, aa ($\alpha$)')     
    plo11_cbar.ax.set_ylabel('p(Obs)/p(Exp)', rotation=-90, va="bottom")    
    
    plt.tight_layout()
    plt.show()    
    return


#######
#Create Ramachandran plots.
#######

def plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi):
    
    fig, plot=plt.subplots(2,3,figsize=(11,7), dpi=100)
    ticks_ar=list(range(-180, 181, 60))
    ticklabels_ar=ticks_ar
    plot[0,0].scatter(sphi_N, spsi_N, s=0.01, c='black')
    plot[0,0].set_title(r'Short domains, N-terminus')
    plot[0,0].set_xticks(ticks_ar)
    plot[0,0].set_xticklabels(ticklabels_ar) 
    plot[0,0].set_yticks(ticks_ar)
    plot[0,0].set_yticklabels(ticklabels_ar)    
    plot[0,0].set_xlabel(r'$\phi$')
    plot[0,0].set_ylabel(r'$\psi$')   
    plot[0,0].set_xlim([-180, 180])
    plot[0,0].set_ylim([-180, 180])
    
    plot[0,1].scatter(sphi_C, spsi_C, s=0.01, c='black')
    plot[0,1].set_title(r'Short domains, C-terminus')
    plot[0,1].set_xticks(ticks_ar)
    plot[0,1].set_xticklabels(ticklabels_ar) 
    plot[0,1].set_yticks(ticks_ar)
    plot[0,1].set_yticklabels(ticklabels_ar)    
    plot[0,1].set_xlabel(r'$\phi$')
    plot[0,1].set_ylabel(r'$\psi$')
    plot[0,1].set_xlim([-180, 180])
    plot[0,1].set_ylim([-180, 180])    
    
    plot[0,2].scatter(sphi, spsi, s=0.01, c='black')
    plot[0,2].set_title(r'Short domains, all')
    plot[0,2].set_xticks(ticks_ar)
    plot[0,2].set_xticklabels(ticklabels_ar) 
    plot[0,2].set_yticks(ticks_ar)
    plot[0,2].set_yticklabels(ticklabels_ar)    
    plot[0,2].set_xlabel(r'$\phi$')
    plot[0,2].set_ylabel(r'$\psi$') 
    plot[0,2].set_xlim([-180, 180])
    plot[0,2].set_ylim([-180, 180])    
    
    plot[1,0].scatter(lphi_N, lpsi_N, s=0.01, c='black')
    plot[1,0].set_title(r'Long domains, N-terminus')
    plot[1,0].set_xticks(ticks_ar)
    plot[1,0].set_xticklabels(ticklabels_ar) 
    plot[1,0].set_yticks(ticks_ar)
    plot[1,0].set_yticklabels(ticklabels_ar)    
    plot[1,0].set_xlabel(r'$\phi$')
    plot[1,0].set_ylabel(r'$\psi$')   
    plot[1,0].set_xlim([-180, 180])
    plot[1,0].set_ylim([-180, 180])    
    
    plot[1,1].scatter(lphi_C, lpsi_C, s=0.01, c='black')
    plot[1,1].set_title(r'Long domains, C-terminus')
    plot[1,1].set_xticks(ticks_ar)
    plot[1,1].set_xticklabels(ticklabels_ar) 
    plot[1,1].set_yticks(ticks_ar)
    plot[1,1].set_yticklabels(ticklabels_ar)    
    plot[1,1].set_xlabel(r'$\phi$')
    plot[1,1].set_ylabel(r'$\psi$') 
    plot[1,1].set_xlim([-180, 180])
    plot[1,1].set_ylim([-180, 180])    
    
    plot[1,2].scatter(lphi, lpsi, s=0.01, c='black')
    plot[1,2].set_title(r'Long domains, all')
    plot[1,2].set_xticks(ticks_ar)
    plot[1,2].set_xticklabels(ticklabels_ar) 
    plot[1,2].set_yticks(ticks_ar)
    plot[1,2].set_yticklabels(ticklabels_ar)    
    plot[1,2].set_xlabel(r'$\phi$')
    plot[1,2].set_ylabel(r'$\psi$')   
    plot[1,2].set_xlim([-180, 180])
    plot[1,2].set_ylim([-180, 180])    
    
    plt.tight_layout()
    plt.show()    
    return


#######
#Create Ramachandran plots using KDE.
#######

def plot_Ramachandran_KDE(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi):
    
    #Define the borders. Taken from https://towardsdatascience.com/simple-example-of-2d-density-plots-in-python-83b83b934f67
    deltaX=(180+180)/36
    deltaY=(180+180)/36
    xmin=-180-deltaX
    xmax=180+deltaX
    ymin=-180-deltaY
    ymax=180+deltaY
    print(xmin, xmax, ymin, ymax)   #Create meshgrid
    xx, yy = np.mgrid[xmin:xmax:360j, ymin:ymax:360j] 
    
    positions=np.vstack([xx.ravel(), yy.ravel()])
    values=np.vstack([sphi_N, spsi_N])
    kernel=st.gaussian_kde(values)
    f=np.reshape(kernel(positions).T, xx.shape)  
    
    fig=plt.figure(figsize=(8,8))
    ax=fig.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    cfset=ax.contourf(xx, yy, f, cmap='coolwarm')
    ax.imshow(np.rot90(f), cmap='coolwarm', extent=[xmin, xmax, ymin, ymax])
    cset=ax.contour(xx, yy, f, colors='k')
    ax.clabel(cset, inline=1, fontsize=10)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\psi$')
    #plt.title('2D Gaussian Kernel density estimation')    
    
    plt.tight_layout()
    plt.show() 
    
    return
    

#######
#Wrapper function.
#######

def wrapper(DSSP_inpath):
    #Define domain length trhresholds.
    min_len=50
    thr_len=130
    
    #Define distance from termini to analyse.
    window_width=50
    local_window_width=50
    
    #Read DSSP data.
    DSSP_data_dict=read_dssp_data(DSSP_inpath)
    #Classify domains by length.
    Short_structures, Long_structures=define_length_groups(DSSP_data_dict, min_len, thr_len)
    
    #Get phi, psi angles for N- and C-termini.
    sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi=phi_psi_N_to_C(Short_structures, window_width)
    lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi=phi_psi_N_to_C(Long_structures, window_width)
    
    #Create Ramachandran plots.
    ##plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi)
    ##plot_Ramachandran_KDE(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi)
    
    #Compute position frequency matrices.
    ss_matrix_N_short, ss_matrix_C_short, ss_pfm_N_short, ss_pfm_C_short, ss_pfm_conf_N_short, ss_pfm_conf_C_short=ss_element_frequency_matrix(Short_structures, window_width)
    ss_matrix_N_long, ss_matrix_C_long, ss_pfm_N_long, ss_pfm_C_long, ss_pfm_conf_N_long, ss_pfm_conf_C_long=ss_element_frequency_matrix(Long_structures, window_width)
    
    #Plot frequency of ss elements as a function of a distance from termini.
    N_to_C_asymmetry(ss_pfm_N_short, ss_pfm_C_short, ss_pfm_conf_N_short, ss_pfm_conf_C_short, ss_pfm_N_long, ss_pfm_C_long, ss_pfm_conf_N_long, ss_pfm_conf_C_long, window_width)
    
    #Enrichment of ss elements N- over C-terminus.
    enrichment_N_to_C_dict_short=ss_ele_enrichment(ss_pfm_N_short, ss_pfm_C_short)
    enrichment_N_to_C_dict_long=ss_ele_enrichment(ss_pfm_N_long, ss_pfm_C_long)
    
    #Plot enrichment of frequency of ss elements at N-terminus over C-terminus.
    N_to_C_enrichment(enrichment_N_to_C_dict_short, enrichment_N_to_C_dict_long, window_width)
    
    #Analyse co-occurence of secondary structure elements at protein termini.
    ##Obs_over_exp_matrices_short=termini_dependance(ss_matrix_N_short, ss_matrix_C_short, ss_pfm_N_short, ss_pfm_C_short, local_window_width)
    ##Obs_over_exp_matrices_long=termini_dependance(ss_matrix_N_long, ss_matrix_C_long, ss_pfm_N_long, ss_pfm_C_long, local_window_width)
    
    #Plot co-occurence of secondary structure elements at protein termini.
    ##plot_co_occurence(Obs_over_exp_matrices_short, Obs_over_exp_matrices_long, local_window_width)
    
    return

wrapper(DSSP_data_inpath)