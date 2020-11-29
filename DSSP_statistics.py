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
            Phi_list=list(line[7])
            Psi_list=list(line[8])
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
    for letter_code in DSSP_alphabet:
        frequency_row_N=[]
        frequency_row_C=[]
        for i in range(len(ss_matrix_N)):
            column_letter_freq_N=ss_matrix_N[i].count(letter_code)/float(len(ss_matrix_N[i]))
            frequency_row_N.append(column_letter_freq_N)
            column_letter_freq_C=ss_matrix_C[i].count(letter_code)/float(len(ss_matrix_C[i]))
            frequency_row_C.append(column_letter_freq_C)            
        ss_pfm_N[letter_code]=frequency_row_N
        ss_pfm_C[letter_code]=frequency_row_C
       
    print(ss_pfm_N) 
    print(ss_pfm_C)   
    return ss_pfm_N, ss_pfm_C


#######
#Analyse N-to-C asymmetry in secondary structure elements frequency.
#######

def N_to_C_asymmetry(ss_pfm_N_short, ss_pfm_C_short, ss_pfm_N_long, ss_pfm_C_long, window_width):
    
    #Plot frequency of ss elements as a function of distance from N- and C-termini.
    X_N=range(window_width)
    X_C=range(60, 60+window_width)
    xticks_ar=list(range(0,111,10))
    xticklabels_ar=[0,10,20,30,40,50,-50,-40,-30,-20,-10,0]
    xticks_spec=[-4, 114]
    xticklabels_spec=['N-', '-C']
    fig, plot=plt.subplots(2,2,figsize=(12,6), dpi=100)
    plot[0,0].plot(X_N, ss_pfm_N_short['H'], color='#94fff1', linewidth=2, label=r'$\alpha$-helix short domains')
    plot[0,0].plot(X_C, ss_pfm_C_short['H'], color='#94fff1', linewidth=2)
    plot[0,0].plot(X_N, ss_pfm_N_long['H'], color='#ab658c', linewidth=2, label=r'$\alpha$-helix long domains')
    plot[0,0].plot(X_C, ss_pfm_C_long['H'], color='#ab658c', linewidth=2) 
    plot[0,0].plot(X_N, ss_pfm_N_short['E'], color='#59acff', linewidth=2, label=r'$\beta$-strand short domains')
    plot[0,0].plot(X_C, ss_pfm_C_short['E'], color='#59acff', linewidth=2)
    plot[0,0].plot(X_N, ss_pfm_N_long['E'], color='#ffec59', linewidth=2, label=r'$\beta$-strand long domains')
    plot[0,0].plot(X_C, ss_pfm_C_long['E'], color='#ffec59', linewidth=2)    
    plot[0,0].set_xticks(xticks_ar)
    plot[0,0].set_xticklabels(xticklabels_ar)
    plot[0,0].set_xticks(xticks_spec, minor=True)
    plot[0,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,0].tick_params(axis='x', which='minor', length=0)
    plot[0,0].set_xlabel('Distance, aa')
    plot[0,0].set_ylabel('Frequency') 
    plot[0,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='upper center')
    
    plot[0,1].plot(X_N, ss_pfm_N_short['B'], color='#94fff1', linewidth=2, label=r'$\beta$-bridge short domains')
    plot[0,1].plot(X_C, ss_pfm_C_short['B'], color='#94fff1', linewidth=2)
    plot[0,1].plot(X_N, ss_pfm_N_long['B'], color='#ab658c', linewidth=2, label=r'$\beta$-bridge long domains')
    plot[0,1].plot(X_C, ss_pfm_C_long['B'], color='#ab658c', linewidth=2) 
    plot[0,1].plot(X_N, ss_pfm_N_short['G'], color='#59acff', linewidth=2, label=r'3-10-helix short domains')
    plot[0,1].plot(X_C, ss_pfm_C_short['G'], color='#59acff', linewidth=2)
    plot[0,1].plot(X_N, ss_pfm_N_long['G'], color='#ffec59', linewidth=2, label=r'3-10-helix long domains')
    plot[0,1].plot(X_C, ss_pfm_C_long['G'], color='#ffec59', linewidth=2)    
    plot[0,1].set_xticks(xticks_ar)
    plot[0,1].set_xticklabels(xticklabels_ar)
    plot[0,1].set_xticks(xticks_spec, minor=True)
    plot[0,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[0,1].tick_params(axis='x', which='minor', length=0)
    plot[0,1].set_xlabel('Distance, aa')
    plot[0,1].set_ylabel('Frequency')  
    plot[0,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center center')
    
    plot[1,0].plot(X_N, ss_pfm_N_short['I'], color='#94fff1', linewidth=2, label=r'$\pi$-helix short domains')
    plot[1,0].plot(X_C, ss_pfm_C_short['I'], color='#94fff1', linewidth=2)
    plot[1,0].plot(X_N, ss_pfm_N_long['I'], color='#ab658c', linewidth=2, label=r'$\pi$-helix long domains')
    plot[1,0].plot(X_C, ss_pfm_C_long['I'], color='#ab658c', linewidth=2) 
    plot[1,0].plot(X_N, ss_pfm_N_short['T'], color='#59acff', linewidth=2, label=r'turn short domains')
    plot[1,0].plot(X_C, ss_pfm_C_short['T'], color='#59acff', linewidth=2)
    plot[1,0].plot(X_N, ss_pfm_N_long['T'], color='#ffec59', linewidth=2, label=r'turn long domains')
    plot[1,0].plot(X_C, ss_pfm_C_long['T'], color='#ffec59', linewidth=2)    
    plot[1,0].set_xticks(xticks_ar)
    plot[1,0].set_xticklabels(xticklabels_ar)
    plot[1,0].set_xticks(xticks_spec, minor=True)
    plot[1,0].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,0].tick_params(axis='x', which='minor', length=0)
    plot[1,0].set_xlabel('Distance, aa')
    plot[1,0].set_ylabel('Frequency')  
    plot[1,0].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center center')
    
    plot[1,1].plot(X_N, ss_pfm_N_short['S'], color='#94fff1', linewidth=2, label=r'bend short domains')
    plot[1,1].plot(X_C, ss_pfm_C_short['S'], color='#94fff1', linewidth=2)
    plot[1,1].plot(X_N, ss_pfm_N_long['S'], color='#ab658c', linewidth=2, label=r'bend long domains')
    plot[1,1].plot(X_C, ss_pfm_C_long['S'], color='#ab658c', linewidth=2) 
    plot[1,1].plot(X_N, ss_pfm_N_short['-'], color='#59acff', linewidth=2, label=r'unstructured short domains')
    plot[1,1].plot(X_C, ss_pfm_C_short['-'], color='#59acff', linewidth=2)
    plot[1,1].plot(X_N, ss_pfm_N_long['-'], color='#ffec59', linewidth=2, label=r'unstructured long domains')
    plot[1,1].plot(X_C, ss_pfm_C_long['-'], color='#ffec59', linewidth=2)    
    plot[1,1].set_xticks(xticks_ar)
    plot[1,1].set_xticklabels(xticklabels_ar)
    plot[1,1].set_xticks(xticks_spec, minor=True)
    plot[1,1].set_xticklabels(xticklabels_spec, minor=True)
    plot[1,1].tick_params(axis='x', which='minor', length=0)
    plot[1,1].set_xlabel('Distance, aa')
    plot[1,1].set_ylabel('Frequency')  
    plot[1,1].legend(fontsize=8.5, ncol=2, handlelength=0.7, frameon=False, columnspacing=0.7, loc='center center')    
    
    
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
    
    #Read DSSP data.
    DSSP_data_dict=read_dssp_data(DSSP_inpath)
    #Classify domains by length.
    Short_structures, Long_structures=define_length_groups(DSSP_data_dict, min_len, thr_len)
    #Compute position frequency matrices.
    ss_pfm_N_short, ss_pfm_C_short=ss_element_frequency_matrix(Short_structures, window_width)
    ss_pfm_N_long, ss_pfm_C_long=ss_element_frequency_matrix(Long_structures, window_width)
    
    #Plot frequency of ss elements as a function of a distance from termini.
    N_to_C_asymmetry(ss_pfm_N_short, ss_pfm_C_short, ss_pfm_N_long, ss_pfm_C_long, window_width)
    
    return

wrapper(DSSP_data_inpath)