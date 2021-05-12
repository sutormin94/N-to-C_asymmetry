###############################################
##Dmitry Sutormin, 2020##
##N-to-C-terminus assymetry in protein domain secondary structure elements composition##

#Pareses DSSP output generated for representative structures with Run_DSSP.py
#Makes Ramachandran plots for N- and C-termini regions of protein domains.
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

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
    
    
#####
#calculates density on Ramachandran plot (returns numpy array with desity values)
#####

def density_Ramachandran(lst_phi, lst_psi, plot_resol=1):
    #plot_resol is plot resolution, every value in array correspond to the amount of phi,psi dots
    #inside the plot_resol*plot_resol (default 1*1) square

    #amount of dots on Ramachandran plot with such resolution
    dots = int(360//plot_resol)

    ram_density_matrix = np.zeros((dots,dots), dtype='int32')
    for phi, psi in zip(lst_phi, lst_psi):
        if phi >= 180 or psi >= 180:
            continue
        #calculates coordinates of plot_resol*plot_resol square to which this dot belongs to,
        #shift all values by 180 to make all indices positive
        coord_phi = int((phi+180)//plot_resol)
        coord_psi = int((psi+180)//plot_resol)
        ram_density_matrix[coord_psi, coord_phi] += 1
    return ram_density_matrix


#####
#makes log(density) Ramachandran plot
#####
def density_plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi, plot_resol=1):
    #plot_resol is plot resolution, every value on plot correspond to the amount of phi,psi dots
    #inside the plot_resol*plot_resol (default 1*1) square

    dots = int(360//plot_resol)
    fig, plot=plt.subplots(2,3,figsize=(11,7), dpi=100)
    ticks_ar=list(range(0, dots+1, int(60/plot_resol)))
    ticklabels_ar=range(-180, 181, 60)

    #log of density array are taken to make plot

    plot00=plot[0,0].imshow(np.log(density_Ramachandran(sphi_N, spsi_N, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo01_cbar=plot[0,0].figure.colorbar(plot00, ax=plot[0,0], shrink=0.7)
    plot[0,0].set_title(r'Short domains, N-terminus')
    plot[0,0].set_xticks(ticks_ar)
    plot[0,0].set_xticklabels(ticklabels_ar)
    plot[0,0].set_yticks(ticks_ar)
    plot[0,0].set_yticklabels(ticklabels_ar)
    plot[0,0].set_xlabel(r'$\phi$')
    plot[0,0].set_ylabel(r'$\psi$')
    plot[0,0].set_xlim([-0.5, dots-0.5])
    plot[0,0].set_ylim([-0.5, dots-0.5])

    plot01=plot[0,1].imshow(np.log(density_Ramachandran(sphi_C, spsi_C, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo01_cbar=plot[0,1].figure.colorbar(plot01, ax=plot[0,1], shrink=0.7)
    plot[0,1].set_title(r'Short domains, C-terminus')
    plot[0,1].set_xticks(ticks_ar)
    plot[0,1].set_xticklabels(ticklabels_ar)
    plot[0,1].set_yticks(ticks_ar)
    plot[0,1].set_yticklabels(ticklabels_ar)
    plot[0,1].set_xlabel(r'$\phi$')
    plot[0,1].set_ylabel(r'$\psi$')
    plot[0,1].set_xlim([-0.5, dots-0.5])
    plot[0,1].set_ylim([-0.5, dots-0.5])

    plot02=plot[0,2].imshow(np.log(density_Ramachandran(sphi, spsi, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo02_cbar=plot[0,2].figure.colorbar(plot02, ax=plot[0,2], shrink=0.7)
    plot[0,2].set_title(r'Short domains')
    plot[0,2].set_xticks(ticks_ar)
    plot[0,2].set_xticklabels(ticklabels_ar)
    plot[0,2].set_yticks(ticks_ar)
    plot[0,2].set_yticklabels(ticklabels_ar)
    plot[0,2].set_xlabel(r'$\phi$')
    plot[0,2].set_ylabel(r'$\psi$')
    plot[0,2].set_xlim([-0.5, dots-0.5])
    plot[0,2].set_ylim([-0.5, dots-0.5])

    plot10=plot[1,0].imshow(np.log(density_Ramachandran(lphi_N, lpsi_N, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo10_cbar=plot[1,0].figure.colorbar(plot10, ax=plot[1,0], shrink=0.7)
    plot[1,0].set_title(r'Long domains, N-terminus')
    plot[1,0].set_xticks(ticks_ar)
    plot[1,0].set_xticklabels(ticklabels_ar)
    plot[1,0].set_yticks(ticks_ar)
    plot[1,0].set_yticklabels(ticklabels_ar)
    plot[1,0].set_xlabel(r'$\phi$')
    plot[1,0].set_ylabel(r'$\psi$')
    plot[1,0].set_xlim([-0.5, dots-0.5])
    plot[1,0].set_ylim([-0.5, dots-0.5])

    plot11=plot[1,1].imshow(np.log(density_Ramachandran(lphi_C, lpsi_C, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo11_cbar=plot[1,1].figure.colorbar(plot11, ax=plot[1,1], shrink=0.7)
    plot[1,1].set_title(r'Long domains, C-terminus')
    plot[1,1].set_xticks(ticks_ar)
    plot[1,1].set_xticklabels(ticklabels_ar)
    plot[1,1].set_yticks(ticks_ar)
    plot[1,1].set_yticklabels(ticklabels_ar)
    plot[1,1].set_xlabel(r'$\phi$')
    plot[1,1].set_ylabel(r'$\psi$')
    plot[1,1].set_xlim([-0.5, dots-0.5])
    plot[1,1].set_ylim([-0.5, dots-0.5])

    plot12=plot[1,2].imshow(np.log(density_Ramachandran(lphi, lpsi, plot_resol)+1), cmap='hot', interpolation='nearest')
    plo12_cbar=plot[1,2].figure.colorbar(plot12, ax=plot[1,2], shrink=0.7)
    plot[1,2].set_title(r'Long domains')
    plot[1,2].set_xticks(ticks_ar)
    plot[1,2].set_xticklabels(ticklabels_ar)
    plot[1,2].set_yticks(ticks_ar)
    plot[1,2].set_yticklabels(ticklabels_ar)
    plot[1,2].set_xlabel(r'$\phi$')
    plot[1,2].set_ylabel(r'$\psi$')
    plot[1,2].set_xlim([-0.5, dots-0.5])
    plot[1,2].set_ylim([-0.5, dots-0.5])

    plt.tight_layout()
    #plt.savefig(DSSP_data_outpath+'phi_psi_density.png', dpi=100, bbox_inches='tight')
    plt.show()
    return

#####
#makes plot with log(ratio of Ramachandran plots density)
#####
def ratio_plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi, plot_resol=1):
    dots = int(360//plot_resol)
    fig, plot=plt.subplots(2,3,figsize=(11,7), dpi=100)
    ticks_ar=list(range(0, dots+1, int(60/plot_resol)))
    ticklabels_ar=range(-180, 181, 60)

    ratio_matrix_short_NC = np.divide(density_Ramachandran(sphi_N, spsi_N, plot_resol)+1, density_Ramachandran(sphi_C, spsi_C, plot_resol)+1)
    plot00=plot[0,0].imshow(np.log(ratio_matrix_short_NC), cmap='hot', interpolation='nearest', vmin=0, vmax=2)
    plot00_cbar=plot[0,0].figure.colorbar(plot00, ax=plot[0,0], shrink=0.7)
    plot[0,0].set_title(r'Short domains, lg(N/C)')
    plot[0,0].set_xticks(ticks_ar)
    plot[0,0].set_xticklabels(ticklabels_ar)
    plot[0,0].set_yticks(ticks_ar)
    plot[0,0].set_yticklabels(ticklabels_ar)
    plot[0,0].set_xlabel(r'$\phi$')
    plot[0,0].set_ylabel(r'$\psi$')
    plot[0,0].set_xlim([-0.5, dots-0.5])
    plot[0,0].set_ylim([-0.5, dots-0.5])

    ratio_matrix_short_CN = np.divide( density_Ramachandran(sphi_C, spsi_C, plot_resol)+1, density_Ramachandran(sphi_N, spsi_N, plot_resol)+1)
    plot10=plot[1,0].imshow(np.log(ratio_matrix_short_CN), cmap='hot', interpolation='nearest', vmin=0, vmax=2)
    plot10_cbar=plot[1,0].figure.colorbar(plot10, ax=plot[1,0], shrink=0.7)
    plot[1,0].set_title(r'Short domains, lg(C/N)')
    plot[1,0].set_xticks(ticks_ar)
    plot[1,0].set_xticklabels(ticklabels_ar)
    plot[1,0].set_yticks(ticks_ar)
    plot[1,0].set_yticklabels(ticklabels_ar)
    plot[1,0].set_xlabel(r'$\phi$')
    plot[1,0].set_ylabel(r'$\psi$')
    plot[1,0].set_xlim([-0.5, dots-0.5])
    plot[1,0].set_ylim([-0.5, dots-0.5])

    ratio_matrix_long_NC = np.divide(density_Ramachandran(lphi_N, lpsi_N, plot_resol)+1, density_Ramachandran(lphi_C, lpsi_C, plot_resol)+1)
    plot01=plot[0,1].imshow(np.log(ratio_matrix_long_NC), cmap='hot', interpolation='nearest', vmin=0, vmax=2)
    plot01_cbar=plot[0,1].figure.colorbar(plot01, ax=plot[0,1], shrink=0.7)
    plot[0,1].set_title(r'Long domains, lg(N/C)')
    plot[0,1].set_xticks(ticks_ar)
    plot[0,1].set_xticklabels(ticklabels_ar)
    plot[0,1].set_yticks(ticks_ar)
    plot[0,1].set_yticklabels(ticklabels_ar)
    plot[0,1].set_xlabel(r'$\phi$')
    plot[0,1].set_ylabel(r'$\psi$')
    plot[0,1].set_xlim([-0.5, dots-0.5])
    plot[0,1].set_ylim([-0.5, dots-0.5])

    ratio_matrix_long_CN = np.divide( density_Ramachandran(lphi_C, lpsi_C, plot_resol)+1, density_Ramachandran(lphi_N, lpsi_N, plot_resol)+1)
    plot11=plot[1,1].imshow(np.log(ratio_matrix_long_CN), cmap='hot', interpolation='nearest', vmin=0, vmax=2)
    plot11_cbar=plot[1,1].figure.colorbar(plot11, ax=plot[1,1], shrink=0.7)
    plot[1,1].set_title(r'Long domains, lg(C/N)')
    plot[1,1].set_xticks(ticks_ar)
    plot[1,1].set_xticklabels(ticklabels_ar)
    plot[1,1].set_yticks(ticks_ar)
    plot[1,1].set_yticklabels(ticklabels_ar)
    plot[1,1].set_xlabel(r'$\phi$')
    plot[1,1].set_ylabel(r'$\psi$')
    plot[1,1].set_xlim([-0.5, dots-0.5])
    plot[1,1].set_ylim([-0.5, dots-0.5])

    plot02=plot[0,2].imshow(np.log(ratio_matrix_long_NC)-np.log(ratio_matrix_short_NC), cmap='hot', interpolation='nearest', vmin=0)
    plot02_cbar=plot[0,2].figure.colorbar(plot02, ax=plot[0,2], shrink=0.7)
    plot[0,2].set_title(r'lg(long(N/C)/short(N/C))')
    plot[0,2].set_xticks(ticks_ar)
    plot[0,2].set_xticklabels(ticklabels_ar)
    plot[0,2].set_yticks(ticks_ar)
    plot[0,2].set_yticklabels(ticklabels_ar)
    plot[0,2].set_xlabel(r'$\phi$')
    plot[0,2].set_ylabel(r'$\psi$')
    plot[0,2].set_xlim([-0.5, dots-0.5])
    plot[0,2].set_ylim([-0.5, dots-0.5])

    plot12=plot[1,2].imshow(np.log(ratio_matrix_long_CN)-np.log(ratio_matrix_short_CN), cmap='hot', interpolation='nearest', vmin=0)
    plot12_cbar=plot[1,2].figure.colorbar(plot12, ax=plot[1,2], shrink=0.7)
    plot[1,2].set_title(r'lg(long(C/N)/short(C/N))')
    plot[1,2].set_xticks(ticks_ar)
    plot[1,2].set_xticklabels(ticklabels_ar)
    plot[1,2].set_yticks(ticks_ar)
    plot[1,2].set_yticklabels(ticklabels_ar)
    plot[1,2].set_xlabel(r'$\phi$')
    plot[1,2].set_ylabel(r'$\psi$')
    plot[1,2].set_xlim([-0.5, dots-0.5])
    plot[1,2].set_ylim([-0.5, dots-0.5])

    plt.tight_layout()
    #plt.savefig(DSSP_data_outpath+'phi_psi_ratio.png', dpi=100, bbox_inches='tight')
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
    window_width=25
    local_window_width=25
    
    #Define resolution of Ramachandran plot binning.
    plot_resol=20
    
    #Read DSSP data.
    DSSP_data_dict=read_dssp_data(DSSP_inpath)
    
    #Classify domains by length.
    Short_structures, Long_structures=define_length_groups(DSSP_data_dict, min_len, thr_len)
    
    #Get phi, psi angles for N- and C-termini.
    sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi=phi_psi_N_to_C(Short_structures, window_width)
    lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi=phi_psi_N_to_C(Long_structures, window_width)
    
    #Create Ramachandran plots.
    #plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi)
    #plot_Ramachandran_KDE(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi)
    
    #Create Ramachandran density plots.
    density_plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi, plot_resol)
    
    #Get ratio of Ramachandran densityy plots.
    ratio_plot_Ramachandran(sphi_N, sphi_C, spsi_N, spsi_C, sphi, spsi, lphi_N, lphi_C, lpsi_N, lpsi_C, lphi, lpsi, plot_resol)    
    
    return

wrapper(DSSP_data_inpath)