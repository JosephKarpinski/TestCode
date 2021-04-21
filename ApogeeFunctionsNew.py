

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ipywidgets as widgets
from IPython.display import display, clear_output
from mpl_toolkits import mplot3d
from termcolor import colored

import seaborn as sns
sns.set()

import astropy.io.fits as pyfits
from scipy import *
from pylab import *
import matplotlib as mpl
import os

# apogee_fields          - 36
# apogee_fields_RC       - 1087
# apogee_target2         - 2498
# apogee_target2RC       - 3429
# plotHistogramGC4       - 3754
# plotHistogramStreams   - 4134
# plotOverview2          - 5492
# plotCompare            - 5802
# plotMetric             - 6355
# plotOverviewGlon       - 6714
# plotMetricGlon         - 6961


### ***********************************************************************************************

def apogee_fields(fs1, rTarget, fs2, r2Target, rfile):
    

    #hdulist = pyfits.open('allStar-l31c.2.fits')
    hdulist = pyfits.open(rfile)

    starbad = 2**23 #bit flag for bad stars 
  
    gd = np.bitwise_and(fs1["ASPCAPFLAG"], starbad) == 0
    s1 = fs1[gd] 
    
    gd = np.bitwise_and(fs2["ASPCAPFLAG"], starbad) == 0
    s2 = fs2[gd] 
    
    
    print("Target: " + rTarget)
    print("\nTarget: " + r2Target)
    
    s1Label = rTarget.strip()
    s2Label = r2Target.strip()

    idx_BA  = where( (s1['ASPCAP_CLASS'] == 'BA') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx2_BA = where( (s2['ASPCAP_CLASS'] == 'BA') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    idx_Fd_a = where( (s1['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_b = where( (s1['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_c = where( (s1['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_d = where( (s1['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Fd_a = where( (s2['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_b = where( (s2['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_c = where( (s2['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_d = where( (s2['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_GKd_a = where( (s1['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_b = where( (s1['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_c = where( (s1['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_d = where( (s1['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKd_a = where( (s2['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_b = where( (s2['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_c = where( (s2['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_d = where( (s2['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]

    
    idx_GKg_a = where( (s1['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_b = where( (s1['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_c = where( (s1['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_d = where( (s1['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKg_a = where( (s2['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_b = where( (s2['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_c = where( (s2['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_d = where( (s2['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Md_a = where( (s1['ASPCAP_CLASS'] == 'Md_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_b = where( (s1['ASPCAP_CLASS'] == 'Md_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_c = where( (s1['ASPCAP_CLASS'] == 'Md_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_d = where( (s1['ASPCAP_CLASS'] == 'Md_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Md_a = where( (s2['ASPCAP_CLASS'] == 'Md_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_b = where( (s2['ASPCAP_CLASS'] == 'Md_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_c = where( (s2['ASPCAP_CLASS'] == 'Md_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_d = where( (s2['ASPCAP_CLASS'] == 'Md_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Mg_a = where( (s1['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_b = where( (s1['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_c = where( (s1['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_d = where( (s1['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Mg_a = where( (s2['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Mg_b = where( (s2['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_b.shape[0] == 1:
            #idx_Mg_b = idx_Mg_b.iloc[1:]
    idx2_Mg_c = where( (s2['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_c.shape[0] == 1:
            #idx_Mg_c = np.delete(idx_Mg_c, 0)
    idx2_Mg_d = where( (s2['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    

    print('\n' + s1Label + ': Stars '   + str(s1.shape[0]))
    
    print('\nBA' + ': Stars '   + str(s1[idx_BA].shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(s1[idx_Fd_a].shape[0]))
    print('Fd_b' + ': Stars '     + str(s1[idx_Fd_b].shape[0]))
    print('Fd_c' + ': Stars '     + str(s1[idx_Fd_c].shape[0]))
    print('Fd_d' + ': Stars '     + str(s1[idx_Fd_d].shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(s1[idx_GKd_a].shape[0]))
    print('GKd_b' + ': Stars '     + str(s1[idx_GKd_b].shape[0]))
    print('GKd_c' + ': Stars '     + str(s1[idx_GKd_c].shape[0]))
    print('GKd_d' + ': Stars '     + str(s1[idx_GKd_d].shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(s1[idx_GKg_a].shape[0]))
    print('GKg_b' + ': Stars '     + str(s1[idx_GKg_b].shape[0]))
    print('GKg_c' + ': Stars '     + str(s1[idx_GKg_c].shape[0]))
    print('GKg_d' + ': Stars '     + str(s1[idx_GKg_d].shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(s1[idx_Md_a].shape[0]))
    print('Md_b' + ': Stars '     + str(s1[idx_Md_b].shape[0]))
    print('Md_c' + ': Stars '     + str(s1[idx_Md_c].shape[0]))
    print('Md_d' + ': Stars '     + str(s1[idx_Md_d].shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(s1[idx_Mg_a].shape[0]))
    print('Mg_b' + ': Stars '     + str(s1[idx_Mg_b].shape[0]))
    print('Mg_c' + ': Stars '     + str(s1[idx_Mg_c].shape[0]))
    print('Mg_d' + ': Stars '     + str(s1[idx_Mg_d].shape[0]))
    
    
    print('\n' + s2Label + ': Stars '   + str(s2.shape[0]))
    
    print('\nBA' + ': Stars '   + str(s2[idx2_BA].shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(s2[idx2_Fd_a].shape[0]))
    print('Fd_b' + ': Stars '     + str(s2[idx2_Fd_b].shape[0]))
    print('Fd_c' + ': Stars '     + str(s2[idx2_Fd_c].shape[0]))
    print('Fd_d' + ': Stars '     + str(s2[idx2_Fd_d].shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(s2[idx2_GKd_a].shape[0]))
    print('GKd_b' + ': Stars '     + str(s2[idx2_GKd_b].shape[0]))
    print('GKd_c' + ': Stars '     + str(s2[idx2_GKd_c].shape[0]))
    print('GKd_d' + ': Stars '     + str(s2[idx2_GKd_d].shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(s2[idx2_GKg_a].shape[0]))
    print('GKg_b' + ': Stars '     + str(s2[idx2_GKg_b].shape[0]))
    print('GKg_c' + ': Stars '     + str(s2[idx2_GKg_c].shape[0]))
    print('GKg_d' + ': Stars '     + str(s2[idx2_GKg_d].shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(s2[idx2_Md_a].shape[0]))
    print('Md_b' + ': Stars '     + str(s2[idx2_Md_b].shape[0]))
    print('Md_c' + ': Stars '     + str(s2[idx2_Md_c].shape[0]))
    print('Md_d' + ': Stars '     + str(s2[idx2_Md_d].shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(s2[idx2_Mg_a].shape[0]))
    print('Mg_b' + ': Stars '     + str(s2[idx2_Mg_b].shape[0]))
    print('Mg_c' + ': Stars '     + str(s2[idx2_Mg_c].shape[0]))
    print('Mg_d' + ': Stars '     + str(s2[idx2_Mg_d].shape[0]))
    
    
    
    print('\nFields: ' + str(unique(s1['FIELD'])))
    
    print('\nFields: ' + str(unique(s2['FIELD'])))
    
    #################################################################################
    ################################################################################
    print("Plot 001")
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']
    rbins = linspace(-3,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('FeH \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('FeH', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 002")
    
    x1 = s1['LOGG']
    x2 = s2['LOGG']
    rbins = linspace(0,6,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('LogG \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('LogG', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()

    ################################################################################
    print("Plot 003")
    
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(2500,8500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 004")
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']
    rbins = linspace(-500,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Vhelio_Avg \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Vhelio_Avg', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 005")
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(-1,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 006")
    
    x1 = s1['H']
    x2 = s2['H']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('H \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('H', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 007")
    
    x1 = s1['J']
    x2 = s2['J']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('J \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('J', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 008")
    
    x1 = s1['K']
    x2 = s2['K']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('K', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 009")
    
    x1 = s1['J']
    x2 = s2['J']
    
    y1 = s1['K']
    y2 = s2['K']
    rbins = linspace(-0.5,1.5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1 - y1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2 - y2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('J-K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('J-K', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 010")
    
    x1 = s1['SNR']
    x2 = s2['SNR']
    rbins = linspace(0,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('SNR \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('SNR', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    
 
    
    
    ################################################################################
    ################################################################################
    print("Plot 011")
    
    x1 = s1['TEFF']
    y1 = s1['LOGG']
    x2 = s2['TEFF']
    y2 = s2['LOGG']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, y1,s=10, label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, y2,s=10, label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  Log G \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Log G', fontsize=15)
    ax.legend(fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(6, 0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    ###############################################################################
    print("Plot 012")
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()

 ###############################################################################
    print("Plot 013")
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Giants \nTeff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, 0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 014")
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_Fd_a], j1[idx_Fd_a] - k1[idx_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x1[idx_Fd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_Fd_b], j1[idx_Fd_b] - k1[idx_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x1[idx_Fd_b].shape[0]), color='b')
    ax.scatter(x1[idx_Fd_c], j1[idx_Fd_c] - k1[idx_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x1[idx_Fd_c].shape[0]), color='r')
    ax.scatter(x1[idx_Fd_d], j1[idx_Fd_d] - k1[idx_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x1[idx_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 015")
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_a], j2[idx2_Fd_a] - k2[idx2_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x2[idx2_Fd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_Fd_b], j2[idx2_Fd_b] - k2[idx2_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x2[idx2_Fd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_Fd_c], j2[idx2_Fd_c] - k2[idx2_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x2[idx2_Fd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_d], j2[idx2_Fd_d] - k2[idx2_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x2[idx2_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 016")
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKd_a], j1[idx_GKd_a] - k1[idx_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x1[idx_GKd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKd_b], j1[idx_GKd_b] - k1[idx_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x1[idx_GKd_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKd_c], j1[idx_GKd_c] - k1[idx_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x1[idx_GKd_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKd_d], j1[idx_GKd_d] - k1[idx_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x1[idx_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 017")
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_a], j2[idx2_GKd_a] - k2[idx2_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x2[idx2_GKd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKd_b], j2[idx2_GKd_b] - k2[idx2_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x2[idx2_GKd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKd_c], j2[idx2_GKd_c] - k2[idx2_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x2[idx2_GKd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_d], j2[idx2_GKd_d] - k2[idx2_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x2[idx2_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 018")
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKg_a], j1[idx_GKg_a] - k1[idx_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x1[idx_GKg_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKg_b], j1[idx_GKg_b] - k1[idx_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x1[idx_GKg_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKg_c], j1[idx_GKg_c] - k1[idx_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x1[idx_GKg_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKg_d], j1[idx_GKg_d] - k1[idx_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x1[idx_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 019")
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_a], j2[idx2_GKg_a] - k2[idx2_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x2[idx2_GKg_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKg_b], j2[idx2_GKg_b] - k2[idx2_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x2[idx2_GKg_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKg_c], j2[idx2_GKg_c] - k2[idx2_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x2[idx2_GKg_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_d], j2[idx2_GKg_d] - k2[idx2_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x2[idx2_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 2500)
    plt.ylim(1.5, -0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    
    ################################################################################
        
    
    
    ################################################################################
        
    
    ################################################################################
    print("Plot 020")
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(2500,8500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x1[idx_GKg_a].shape[0]))
    ax.hist(x1[idx_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x1[idx_GKg_b].shape[0]))
    ax.hist(x1[idx_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x1[idx_GKg_c].shape[0]))
    ax.hist(x1[idx_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x1[idx_GKg_d].shape[0]))
    ax.set_title(str(s1Label) +  '  \n ' + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 021")
    
    rbins = linspace(2500,8500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x2[idx2_GKg_a].shape[0]))
    ax.hist(x2[idx2_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x2[idx2_GKg_b].shape[0]))
    ax.hist(x2[idx2_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x2[idx2_GKg_c].shape[0]))
    ax.hist(x2[idx2_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x2[idx2_GKg_d].shape[0]))
    ax.set_title(str(s2Label) + ' \n '  + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 022")
    
    x1 = s1['TEFF']
    rbins = linspace(2500,8500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x1[idx_GKd_a].shape[0]))
    ax.hist(x1[idx_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x1[idx_GKd_b].shape[0]))
    ax.hist(x1[idx_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x1[idx_GKd_c].shape[0]))
    ax.hist(x1[idx_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x1[idx_GKd_d].shape[0]))
    ax.set_title(str(s1Label) + ' \n '  + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    ################################################################################
    print("Plot 023")
    
    x2 = s2['TEFF']
    rbins = linspace(2500,8500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x2[idx2_GKd_a].shape[0]))
    ax.hist(x2[idx2_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x2[idx2_GKd_b].shape[0]))
    ax.hist(x2[idx2_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x2[idx2_GKd_c].shape[0]))
    ax.hist(x2[idx2_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x2[idx2_GKd_d].shape[0]))
    ax.set_title(str(s2Label) + ' \n '  + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()
    
    
    
    ################################################################################
    print("Plot 024")
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s1['GLON'], s1['GLAT'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s2['GLON'], s2['GLAT'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    ################################################################################
    print("Plot 025")
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s1['RA'], s1['DEC'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    ################################################################################
    print("Plot 026")
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s2['RA'], s2['DEC'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    ################################################################################
    print("Plot 027")
    
    x1 = s1[(s1['GAIA_PMRA'].astype(float) > -9999)]
    x2 = s2[(s2['GAIA_PMRA'].astype(float) > -9999)]
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(x1['GAIA_PMRA'], x1['GAIA_PMDEC'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    ################################################################################
    print("Plot 028")
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(x2['GAIA_PMRA'], x2['GAIA_PMDEC'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    ################################################################################
    print("Plot 029")
    
    x11 = x1[(x1['LOGG'].astype(float) > -9999)]
    x22 = x2[(x2['LOGG'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x11['GAIA_PMRA'], x11['GAIA_PMDEC'], label=s1Label, marker='o', linewidth=0,
                c=x11.LOGG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('LOGG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 030")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x22['GAIA_PMRA'], x22['GAIA_PMDEC'], label=s2Label, marker='o', linewidth=0,
                c=x22.LOGG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('LOGG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 031")
    
    x111 = x1[(x1['FE_H'].astype(float) > -9999)]
    x222 = x2[(x2['FE_H'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x111['GAIA_PMRA'], x111['GAIA_PMDEC'], label=s1Label, marker='o', linewidth=0,
                c=x111.FE_H,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('FeH',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 032")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x222['GAIA_PMRA'], x222['GAIA_PMDEC'], label=s2Label, marker='o', linewidth=0,
                c=x222.FE_H,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('GAIA_PMRA', fontsize=15)
    plt.ylabel('GAIA_PMDEC', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('FeH',size=15)
    plt.show()
    
    
    
    ################################################################################
    print("Plot 033")
    
    #x1 = s1['TEFF']
    #x2 = s2['TEFF']
    
    x1 = s1[(s1['TEFF'].astype(float) > -9999)]
    
    x2 = s2[(s2['TEFF'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.TEFF,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('TEFF',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 034")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.TEFF,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('TEFF',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 035")
    
    #x1 = s1['LOGG']
    #x2 = s2['LOGG']
    
    x1 = s1[(s1['LOGG'].astype(float) > -9999)]
    
    x2 = s2[(s2['LOGG'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.LOGG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('LOGG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 036")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.LOGG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('LOGG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 037")
    
    
    x1 = s1[(s1['FE_H'].astype(float) > -9999)]
    
    x2 = s2[(s2['FE_H'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.FE_H,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('FeH',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 038")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.FE_H,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('FeH',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 039")
    
    
    x1 = s1[(s1['VHELIO_AVG'].astype(float) > -9999)]
    
    x2 = s2[(s2['VHELIO_AVG'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.VHELIO_AVG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('VHELIO_AVG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 040")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.VHELIO_AVG,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('VHELIO_AVG',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 041")
    
    
    x1 = s1[(s1['GAIA_RADIAL_VELOCITY'].astype(float) > -9999)]
    
    x2 = s2[(s2['GAIA_RADIAL_VELOCITY'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.GAIA_RADIAL_VELOCITY,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('GAIA_RADIAL_VELOCITY',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 042")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.GAIA_RADIAL_VELOCITY,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('GAIA_RADIAL_VELOCITY',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 043")
    
     
    x1 = s1[(s1['AL_FE'].astype(float) > -9999)]
    
    x2 = s2[(s2['AL_FE'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.AL_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('AL_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 044")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.AL_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('AL_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 045")
    
     
    x1 = s1[(s1['MG_FE'].astype(float) > -9999)]
    
    x2 = s2[(s2['MG_FE'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.MG_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('MG_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 046")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.MG_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('MG_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 047")
    
    x1 = s1[(s1['MN_FE'].astype(float) > -9999)]
    
    x2 = s2[(s2['MN_FE'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.MN_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('MN_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 048")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.MN_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('MN_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 049")
    
    x1 = s1[(s1['O_FE'].astype(float) > -9999)]
    
    x2 = s2[(s2['O_FE'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.O_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('O_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 050")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.O_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('O_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 051")
    
    x1 = s1[(s1['C_FE'].astype(float) > -9999)]
    
    x2 = s2[(s2['C_FE'].astype(float) > -9999)]
    
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x1['RA'], x1['DEC'], label=s1Label, marker='o', linewidth=0,
                c=x1.C_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('C_FE',size=15)
    plt.show()
    
    
    ################################################################################
    print("Plot 052")
    
    fig = plt.figure(figsize=(15, 12))
    
    plt.subplot(111)
    image = plt.scatter(x2['RA'], x2['DEC'], label=s2Label, marker='o', linewidth=0,
                c=x2.C_FE,edgecolor='k',cmap="jet_r")
    ax = plt.gca()
    plt.gca().invert_yaxis()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    bar = fig.colorbar(image,orientation="vertical",pad=0.01)
    bar.set_label('C_FE',size=15)
    plt.show()
    
    ################################################################################
    print("Plot 053")
    
    x1 = s1['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 054")
    
    x2 = s2['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 055")
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(0,0.5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    plt.show()
    
    ################################################################################
    print("Plot 056")

    
    x1 = s1['M_H']
    x2 = s2['M_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('M_H [dex]', fontsize=13)
    plt.title("M_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 057")
    
    x1 = s1['LOGG']
    x2 = s2['LOGG']

    rLow = 0
    rHi  = 6

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('LOGG [log (cgs)]', fontsize=13)
    plt.title("LOGG \nlog g from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 058")
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('FE_H [dex]', fontsize=13)
    plt.title("FE_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 059")
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']

    rLow = 2500
    rHi  = 8500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('TEFF [K]', fontsize=13)
    plt.title("TEFF \nTeff from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 060")
    
    x1 = s1['SNR']
    x2 = s2['SNR']

    rLow = 0
    rHi  = 500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('SNR', fontsize=13)
    plt.title("SNR \nmedian S/N per pixel in combined frame (at apStar sampling)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 061")
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VHELIO_AVG [km/s]', fontsize=13)
    plt.title("VHELIO_AVG \naverage radial velocity, weighted by S/N, using RVs ", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 062")
    
    x1 = s1['VSCATTER']
    x2 = s2['VSCATTER']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VSCATTER [km/s]', fontsize=13)
    plt.title("VSCATTER \nscatter of individual visit RVs around average", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 063")
    
    x1 = s1['VMACRO']
    x2 = s2['VMACRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMACRO (cgs)', fontsize=13)
    plt.title("VMACRO \nmacroturbulent velocity (f(log Teff, [M/H]) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 064")
    
    x1 = s1['VMICRO']
    x2 = s2['VMICRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMICRO (cgs)', fontsize=13)
    plt.title("VMICRO \nmicroturbulent velocity (fit for dwarfs, f(log g) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 065")
    
    x1 = s1['J']
    x2 = s2['J']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('J', fontsize=13)
    plt.title("J \n2MASS J mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 066")
    
    x1 = s1['H']
    x2 = s2['H']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('H', fontsize=13)
    plt.title("H \n2MASS H mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 067")
    
    
    x1 = s1['K']
    x2 = s2['K']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('K', fontsize=13)
    plt.title("K \n2MASS K mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 068")
    
    j1 = s1['J']
    j2 = s2['J']
    
    k1 = s1['K']
    k2 = s2['K']

    x1 = j1 - k1
    x2 = j2 - k2
    
    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('J-K', fontsize=13)
    plt.title("J-K \n2MASS J mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 069")
    
    x1 = s1['AL_FE']
    x2 = s2['AL_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('AL_FE', fontsize=13)
    plt.title("Aluminum AL_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 070")
    
    x1 = s1['C_FE']
    x2 = s2['C_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('C_FE', fontsize=13)
    plt.title("Carbon C_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 071")
    
    x1 = s1['CA_FE']
    x2 = s2['CA_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('CA_FE', fontsize=13)
    plt.title("Calcium CA_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 072")
    
    x1 = s1['CI_FE']
    x2 = s2['CI_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('CI_FE', fontsize=13)
    plt.title("CI_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 073")
    
    x1 = s1['CO_FE']
    x2 = s2['CO_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('CO_FE', fontsize=13)
    plt.title("Cobalt CO_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 074")
    
    x1 = s1['CR_FE']
    x2 = s2['CR_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('CR_FE', fontsize=13)
    plt.title("Chromium CR_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 075")
    
    x1 = s1['K_FE']
    x2 = s2['K_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('K_FE', fontsize=13)
    plt.title("Potassium K_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 076")
    
    x1 = s1['MG_FE']
    x2 = s2['MG_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('MG_FE', fontsize=13)
    plt.title("Magnesium MG_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 077")
    
    x1 = s1['MN_FE']
    x2 = s2['MN_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('MN_FE', fontsize=13)
    plt.title("Manganese MN_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 078")
    
    x1 = s1['N_FE']
    x2 = s2['N_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('N_FE', fontsize=13)
    plt.title("Nitrogen N_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 079")
    
    x1 = s1['NA_FE']
    x2 = s2['NA_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('NA_FE', fontsize=13)
    plt.title("Sodium NA_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    print("Plot 080")
    
    
    x1 = s1['NI_FE']
    x2 = s2['NI_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('NI_FE', fontsize=13)
    plt.title("Nickel NI_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 081")
    
    x1 = s1['O_FE']
    x2 = s2['O_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('O_FE', fontsize=13)
    plt.title("Oxygen O_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 082")
    
    x1 = s1['P_FE']
    x2 = s2['P_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('P_FE', fontsize=13)
    plt.title("Phosphorus P_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 083")
    
    x1 = s1['S_FE']
    x2 = s2['S_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('S_FE', fontsize=13)
    plt.title("Sulfur S_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 084")
    
    x1 = s1['SI_FE']
    x2 = s2['SI_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('SI_FE', fontsize=13)
    plt.title("Silicon SI_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 085")
    
    x1 = s1['TI_FE']
    x2 = s2['TI_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('TI_FE', fontsize=13)
    plt.title("Titanium TI_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 086")
    
    x1 = s1['TIII_FE']
    x2 = s2['TIII_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('TIII_FE', fontsize=13)
    plt.title("TIII_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    #print("Plot 087")
    
    x1 = s1['V_FE']
    x2 = s2['V_FE']

    rLow = -1
    rHi  = 1

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('V_FE', fontsize=13)
    plt.title("Vanadium V_FE", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    


### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************
    
def apogee_fields_RC(s1, rTarget, s2, r2Target, rfile):

    #hdulist = pyfits.open('allStar-l31c.2.fits')
    hdulist = pyfits.open(rfile)
    
    print("Target: " + rTarget)
    print("\nTarget: " + r2Target)
    
    s1Label = rTarget.strip()
    s2Label = r2Target.strip()

    idx_Fd_a = where( (s1['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_b = where( (s1['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_c = where( (s1['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_d = where( (s1['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Fd_a = where( (s2['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_b = where( (s2['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_c = where( (s2['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_d = where( (s2['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_GKd_a = where( (s1['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_b = where( (s1['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_c = where( (s1['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_d = where( (s1['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKd_a = where( (s2['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_b = where( (s2['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_c = where( (s2['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_d = where( (s2['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]

    
    idx_GKg_a = where( (s1['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_b = where( (s1['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_c = where( (s1['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_d = where( (s1['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKg_a = where( (s2['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_b = where( (s2['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_c = where( (s2['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_d = where( (s2['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Md_a = where( (s1['ASPCAP_CLASS'] == 'Md_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_b = where( (s1['ASPCAP_CLASS'] == 'Md_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_c = where( (s1['ASPCAP_CLASS'] == 'Md_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_d = where( (s1['ASPCAP_CLASS'] == 'Md_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Md_a = where( (s2['ASPCAP_CLASS'] == 'Md_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_b = where( (s2['ASPCAP_CLASS'] == 'Md_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_c = where( (s2['ASPCAP_CLASS'] == 'Md_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_d = where( (s2['ASPCAP_CLASS'] == 'Md_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Mg_a = where( (s1['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_b = where( (s1['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_c = where( (s1['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_d = where( (s1['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Mg_a = where( (s2['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Mg_b = where( (s2['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_b.shape[0] == 1:
            #idx_Mg_b = idx_Mg_b.iloc[1:]
    idx2_Mg_c = where( (s2['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_c.shape[0] == 1:
            #idx_Mg_c = np.delete(idx_Mg_c, 0)
    idx2_Mg_d = where( (s2['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    

    print('\n' + s1Label + ': Stars '   + str(s1.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx_Mg_d.shape[0]))
    
    
    print('\n' + s2Label + ': Stars '   + str(s2.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx2_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx2_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx2_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx2_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx2_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx2_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx2_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx2_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx2_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx2_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx2_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx2_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx2_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx2_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx2_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx2_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx2_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx2_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx2_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx2_Mg_d.shape[0]))
    
    
    
    print('\nFields: ' + str(unique(s1['FIELD'])))
    
    print('\nFields: ' + str(unique(s2['FIELD'])))
    
    #################################################################################
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']
    rbins = linspace(-3,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('FeH \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('FeH', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['LOGG']
    x2 = s2['LOGG']
    rbins = linspace(0,4,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('LogG \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('LogG', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()

    ################################################################################
    
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,8500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']
    rbins = linspace(-500,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Vhelio_Avg \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Vhelio_Avg', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(-1,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['H']
    x2 = s2['H']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('H \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('H', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['J']
    x2 = s2['J']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('J \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('J', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['K']
    x2 = s2['K']
    rbins = linspace(5,20,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('K', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['SNR']
    x2 = s2['SNR']
    rbins = linspace(0,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('SNR \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('SNR', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    
 
    
    
    ################################################################################
    ################################################################################
    
    x1 = s1['TEFF']
    y1 = s1['LOGG']
    x2 = s2['TEFF']
    y2 = s2['LOGG']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, y1,s=10, label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, y2,s=10, label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  Log G \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Log G', fontsize=15)
    ax.legend(fontsize=13)
    plt.xlim(6000, 3500)
    plt.ylim(0, 5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    show()
    
    ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()

 ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Giants \nTeff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_Fd_a], j1[idx_Fd_a] - k1[idx_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x1[idx_Fd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_Fd_b], j1[idx_Fd_b] - k1[idx_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x1[idx_Fd_b].shape[0]), color='b')
    ax.scatter(x1[idx_Fd_c], j1[idx_Fd_c] - k1[idx_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x1[idx_Fd_c].shape[0]), color='r')
    ax.scatter(x1[idx_Fd_d], j1[idx_Fd_d] - k1[idx_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x1[idx_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_a], j2[idx2_Fd_a] - k2[idx2_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x2[idx2_Fd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_Fd_b], j2[idx2_Fd_b] - k2[idx2_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x2[idx2_Fd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_Fd_c], j2[idx2_Fd_c] - k2[idx2_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x2[idx2_Fd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_d], j2[idx2_Fd_d] - k2[idx2_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x2[idx2_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKd_a], j1[idx_GKd_a] - k1[idx_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x1[idx_GKd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKd_b], j1[idx_GKd_b] - k1[idx_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x1[idx_GKd_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKd_c], j1[idx_GKd_c] - k1[idx_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x1[idx_GKd_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKd_d], j1[idx_GKd_d] - k1[idx_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x1[idx_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_a], j2[idx2_GKd_a] - k2[idx2_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x2[idx2_GKd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKd_b], j2[idx2_GKd_b] - k2[idx2_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x2[idx2_GKd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKd_c], j2[idx2_GKd_c] - k2[idx2_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x2[idx2_GKd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_d], j2[idx2_GKd_d] - k2[idx2_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x2[idx2_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKg_a], j1[idx_GKg_a] - k1[idx_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x1[idx_GKg_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKg_b], j1[idx_GKg_b] - k1[idx_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x1[idx_GKg_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKg_c], j1[idx_GKg_c] - k1[idx_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x1[idx_GKg_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKg_d], j1[idx_GKg_d] - k1[idx_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x1[idx_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_a], j2[idx2_GKg_a] - k2[idx2_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x2[idx2_GKg_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKg_b], j2[idx2_GKg_b] - k2[idx2_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x2[idx2_GKg_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKg_c], j2[idx2_GKg_c] - k2[idx2_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x2[idx2_GKg_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_d], j2[idx2_GKg_d] - k2[idx2_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x2[idx2_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + '\n Teff  vs  J - K ', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    ################################################################################
        
    
    
    ################################################################################
        
    
    ################################################################################
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x1[idx_GKg_a].shape[0]))
    ax.hist(x1[idx_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x1[idx_GKg_b].shape[0]))
    ax.hist(x1[idx_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x1[idx_GKg_c].shape[0]))
    ax.hist(x1[idx_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x1[idx_GKg_d].shape[0]))
    ax.set_title(str(s1Label) +  '  \n ' + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x2[idx2_GKg_a].shape[0]))
    ax.hist(x2[idx2_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x2[idx2_GKg_b].shape[0]))
    ax.hist(x2[idx2_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x2[idx2_GKg_c].shape[0]))
    ax.hist(x2[idx2_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x2[idx2_GKg_d].shape[0]))
    ax.set_title(str(s2Label) + ' \n '  + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x1[idx_GKd_a].shape[0]))
    ax.hist(x1[idx_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x1[idx_GKd_b].shape[0]))
    ax.hist(x1[idx_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x1[idx_GKd_c].shape[0]))
    ax.hist(x1[idx_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x1[idx_GKd_d].shape[0]))
    ax.set_title(str(s1Label) + ' \n '  + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x2[idx2_GKd_a].shape[0]))
    ax.hist(x2[idx2_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x2[idx2_GKd_b].shape[0]))
    ax.hist(x2[idx2_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x2[idx2_GKd_c].shape[0]))
    ax.hist(x2[idx2_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x2[idx2_GKd_d].shape[0]))
    ax.set_title(str(s2Label) + ' \n '  + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    
    ################################################################################
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s1['GLON'], s1['GLAT'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s2['GLON'], s2['GLAT'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    ################################################################################
    
    x1 = s1['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    x2 = s2['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(0,0.5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################

    x1 = s1['LOGG']
    x2 = s2['LOGG']

    rLow = 0
    rHi  = 4

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('LOGG [log (cgs)]', fontsize=13)
    plt.title("LOGG \nlog g from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['M_H']
    x2 = s2['M_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('M_H [dex]', fontsize=13)
    plt.title("M_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('FE_H [dex]', fontsize=13)
    plt.title("FE_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']

    rLow = 2500
    rHi  = 8500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('TEFF [K]', fontsize=13)
    plt.title("TEFF \nTeff from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['SNR']
    x2 = s2['SNR']

    rLow = 0
    rHi  = 500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('SNR', fontsize=13)
    plt.title("SNR \nmedian S/N per pixel in combined frame (at apStar sampling)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VHELIO_AVG [km/s]', fontsize=13)
    plt.title("VHELIO_AVG \naverage radial velocity, weighted by S/N, using RVs ", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VSCATTER']
    x2 = s2['VSCATTER']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VSCATTER [km/s]', fontsize=13)
    plt.title("VSCATTER \nscatter of individual visit RVs around average", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VMACRO']
    x2 = s2['VMACRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMACRO (cgs)', fontsize=13)
    plt.title("VMACRO \nmacroturbulent velocity (f(log Teff, [M/H]) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VMICRO']
    x2 = s2['VMICRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMICRO (cgs)', fontsize=13)
    plt.title("VMICRO \nmicroturbulent velocity (fit for dwarfs, f(log g) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['J']
    x2 = s2['J']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('J', fontsize=13)
    plt.title("J \n2MASS J mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['J0']
    x2 = s2['J0']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('J0', fontsize=13)
    plt.title("J0 \nExtinction-correct 2MASS J mag", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['H']
    x2 = s2['H']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('H', fontsize=13)
    plt.title("H \n2MASS H mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['H0']
    x2 = s2['H0']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('H0', fontsize=13)
    plt.title("H0 \nExtinction-correct 2MASS H mag", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['K']
    x2 = s2['K']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('K', fontsize=13)
    plt.title("K \n2MASS K mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['K0']
    x2 = s2['K0']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('K0', fontsize=13)
    plt.title("K0 \nExtinction-correct 2MASS K mag", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['RC_DIST']
    x2 = s2['RC_DIST']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('RC_DIST [kpc]', fontsize=13)
    plt.title("RC_DIST \nDistance", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['RC_DM']
    x2 = s2['RC_DM']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('RC_DM', fontsize=13)
    plt.title("RC_DM \nDistance modulus", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['RC_GALR']
    x2 = s2['RC_GALR']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('RC_GALR [kpc]', fontsize=13)
    plt.title("RC_GALR \nGalactocentric distance assuming R0 = 8 kpc, Z0 = 25 pc", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['RC_GALPHI']
    x2 = s2['RC_GALPHI']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('RC_GALPHI [radians]', fontsize=13)
    plt.title("RC_GALPHI \nGalactocentric azimuth assuming R0 = 8 kpc, Z0 = 25 pc", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['RC_GALZ']
    x2 = s2['RC_GALZ']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('RC_GALZ [kpc]', fontsize=13)
    plt.title("RC_GALZ \nDistance from the mid-plane assuming R0 = 8 kpc, Z0 = 25 pc", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['PMRA']
    x2 = s2['PMRA']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('PMRA [mas/yr]', fontsize=13)
    plt.title("PMRA \nProper motion from the UCAC-4 catalog", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['PMDEC']
    x2 = s2['PMDEC']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('PMDEC [mas/yr]', fontsize=13)
    plt.title("PMDEC \nProper motion from the UCAC-4 catalog", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVR']
    x2 = s2['GALVR']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVR [km/s]', fontsize=13)
    plt.title("GALVR \nGalactocentric radial velocity based on the UCAC-4 proper motion and assuming R0 = 8 kpc, \nZ0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame)",           fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVT']
    x2 = s2['GALVT']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVT [km/s]', fontsize=13)
    plt.title("GALVT \nGalactocentric tangential velocity based on the UCAC-4 proper motion and assuming R0 = 8 kpc, \nZ0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVZ']
    x2 = s2['GALVZ']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVZ [km/s]', fontsize=13)
    plt.title("GALVZ \nGalactocentric tangential velocity based on the UCAC-4 proper motion and assuming R0 = 8 kpc, \nZ0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['PMRA_HSOY']
    x2 = s2['PMRA_HSOY']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('PMRA_HSOY [mas/yr]', fontsize=13)
    plt.title("PMRA_HSOY \nProper motion from the HSOY catalog)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['PMDEC_HSOY']
    x2 = s2['PMDEC_HSOY']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('PMDEC_HSOY [mas/yr]', fontsize=13)
    plt.title("PMDEC_HSOY \nProper motion from the HSOY catalog)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVR_HSOY']
    x2 = s2['GALVR_HSOY']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVR_HSOY [km/s]', fontsize=13)
    plt.title("GALVR_HSOY \nGalactocentric radial velocity based on the HSOY proper motion and assuming R0 = 8 kpc,\n Z0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame))", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVT_HSOY']
    x2 = s2['GALVT_HSOY']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVT_HSOY [km/s]', fontsize=13)
    plt.title("GALVT_HSOY \nGalactocentric tangential velocity based on the HSOY proper motion and assuming R0 = 8 kpc,\n Z0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame))", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['GALVZ_HSOY']
    x2 = s2['GALVZ_HSOY']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('GALVZ_HSOY [km/s]', fontsize=13)
    plt.title("GALVZ_HSOY \nGalactocentric tangential velocity based on the HSOY proper motion and assuming R0 = 8 kpc,\n Z0 = 25 pc, and Vsun=[-11.1,30.24*8.,7.2] (left-handed coordinate frame))", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['METALS']
    x2 = s2['METALS']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('METALS [dex]', fontsize=13)
    plt.title("METALS \nBackward-compatible overall metallicity (from PARAM))", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['ALPHAFE']
    x2 = s2['ALPHAFE']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('ALPHAFE [dex]', fontsize=13)
    plt.title("ALPHAFE \nBackward-compatible overall alpha enhancement (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################

    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************
    
    
def apogee_target2(rTarget, rFlag, rBit, r2Target, r2Flag, r2Bit):

    #hdulist = pyfits.open('allStar-l31c.2.fits')
    hdulist = pyfits.open('allStar-r12-l33.fits')
    
    print("Target: " + rTarget)
    print("Flag:   "   + rFlag)
    print("Bit: "    + str(rBit))
    
    print("\nTarget: " + r2Target)
    print("Flag:   "   + r2Flag)
    print("Bit: "    + str(r2Bit))
    
    cut_Apogee_Target = (hdulist[1].data[rFlag] & 2**rBit)
    idx_Apogee_Target = where(cut_Apogee_Target)[0]
    rApogee_Target = hdulist[1].data[idx_Apogee_Target]
    s1 = rApogee_Target
    s1Label = rTarget
    
    cut_Apogee_Target2 = (hdulist[1].data[r2Flag] & 2**r2Bit)
    idx_Apogee_Target2 = where(cut_Apogee_Target2)[0]
    r2Apogee_Target = hdulist[1].data[idx_Apogee_Target2]
    s2 = r2Apogee_Target
    s2Label = r2Target
    
    idx_BA  = where( (s1['ASPCAP_CLASS'] == 'BA') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx2_BA = where( (s2['ASPCAP_CLASS'] == 'BA') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    idx_Fd_a = where( (s1['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_b = where( (s1['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_c = where( (s1['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_d = where( (s1['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Fd_a = where( (s2['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_b = where( (s2['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_c = where( (s2['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_d = where( (s2['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_GKd_a = where( (s1['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_b = where( (s1['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_c = where( (s1['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_d = where( (s1['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKd_a = where( (s2['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_b = where( (s2['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_c = where( (s2['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_d = where( (s2['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]

    
    idx_GKg_a = where( (s1['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_b = where( (s1['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_c = where( (s1['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_d = where( (s1['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKg_a = where( (s2['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_b = where( (s2['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_c = where( (s2['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_d = where( (s2['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Md_a = where( (s1['ASPCAP_CLASS'] == 'Md_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_b = where( (s1['ASPCAP_CLASS'] == 'Md_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_c = where( (s1['ASPCAP_CLASS'] == 'Md_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_d = where( (s1['ASPCAP_CLASS'] == 'Md_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Md_a = where( (s2['ASPCAP_CLASS'] == 'Md_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_b = where( (s2['ASPCAP_CLASS'] == 'Md_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_c = where( (s2['ASPCAP_CLASS'] == 'Md_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_d = where( (s2['ASPCAP_CLASS'] == 'Md_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Mg_a = where( (s1['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_b = where( (s1['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_c = where( (s1['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_d = where( (s1['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Mg_a = where( (s2['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Mg_b = where( (s2['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_b.shape[0] == 1:
            #idx_Mg_b = idx_Mg_b.iloc[1:]
    idx2_Mg_c = where( (s2['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_c.shape[0] == 1:
            #idx_Mg_c = np.delete(idx_Mg_c, 0)
    idx2_Mg_d = where( (s2['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    

    print('\n' + s1Label + ': Stars '   + str(s1.shape[0]))
    
    print('\nBA' + ': Stars '   + str(idx_BA.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx_Mg_d.shape[0]))
    
    
    print('\n' + s2Label + ': Stars '   + str(s2.shape[0]))
    
    print('\nBA' + ': Stars '   + str(idx2_BA.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx2_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx2_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx2_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx2_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx2_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx2_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx2_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx2_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx2_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx2_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx2_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx2_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx2_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx2_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx2_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx2_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx2_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx2_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx2_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx2_Mg_d.shape[0]))
    
    
    
    print('\nFields: ' + str(unique(s1['FIELD'])))
    
    print('\nFields: ' + str(unique(s2['FIELD'])))
    
    #################################################################################
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']
    rbins = linspace(-3,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('FeH \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('FeH', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['LOGG']
    x2 = s2['LOGG']
    rbins = linspace(0,6,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('LogG \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('LogG', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()

    ################################################################################
    
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,8500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']
    rbins = linspace(-500,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Vhelio_Avg \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Vhelio_Avg', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    y1 = s1['LOGG']
    x2 = s2['TEFF']
    y2 = s2['LOGG']
    rbins = linspace(0,6,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, y1,s=10, label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, y2,s=10, label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  Log G \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Log G', fontsize=15)
    ax.legend(fontsize=13)
    plt.xlim(8500, 3500)
    plt.ylim(0, 6)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    show()
    
    ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()

 ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Giants Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_Fd_a], j1[idx_Fd_a] - k1[idx_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x1[idx_Fd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_Fd_b], j1[idx_Fd_b] - k1[idx_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x1[idx_Fd_b].shape[0]), color='b')
    ax.scatter(x1[idx_Fd_c], j1[idx_Fd_c] - k1[idx_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x1[idx_Fd_c].shape[0]), color='r')
    ax.scatter(x1[idx_Fd_d], j1[idx_Fd_d] - k1[idx_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x1[idx_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_a], j2[idx2_Fd_a] - k2[idx2_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x2[idx2_Fd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_Fd_b], j2[idx2_Fd_b] - k2[idx2_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x2[idx2_Fd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_Fd_c], j2[idx2_Fd_c] - k2[idx2_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x2[idx2_Fd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_d], j2[idx2_Fd_d] - k2[idx2_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x2[idx2_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(8500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKd_a], j1[idx_GKd_a] - k1[idx_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x1[idx_GKd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKd_b], j1[idx_GKd_b] - k1[idx_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x1[idx_GKd_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKd_c], j1[idx_GKd_c] - k1[idx_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x1[idx_GKd_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKd_d], j1[idx_GKd_d] - k1[idx_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x1[idx_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_a], j2[idx2_GKd_a] - k2[idx2_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x2[idx2_GKd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKd_b], j2[idx2_GKd_b] - k2[idx2_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x2[idx2_GKd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKd_c], j2[idx2_GKd_c] - k2[idx2_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x2[idx2_GKd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_d], j2[idx2_GKd_d] - k2[idx2_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x2[idx2_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKg_a], j1[idx_GKg_a] - k1[idx_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x1[idx_GKg_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKg_b], j1[idx_GKg_b] - k1[idx_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x1[idx_GKg_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKg_c], j1[idx_GKg_c] - k1[idx_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x1[idx_GKg_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKg_d], j1[idx_GKg_d] - k1[idx_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x1[idx_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_a], j2[idx2_GKg_a] - k2[idx2_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x2[idx2_GKg_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKg_b], j2[idx2_GKg_b] - k2[idx2_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x2[idx2_GKg_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKg_c], j2[idx2_GKg_c] - k2[idx2_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x2[idx2_GKg_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_d], j2[idx2_GKg_d] - k2[idx2_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x2[idx2_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    ################################################################################
        
    
    
    ################################################################################
        
    
    ################################################################################
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x1[idx_GKg_a].shape[0]))
    ax.hist(x1[idx_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x1[idx_GKg_b].shape[0]))
    ax.hist(x1[idx_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x1[idx_GKg_c].shape[0]))
    ax.hist(x1[idx_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x1[idx_GKg_d].shape[0]))
    ax.set_title(str(s1Label) +  ' Sloan Digital Sky Survey Apogee DR15 \n ' + str(s1Label) + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x2[idx2_GKg_a].shape[0]))
    ax.hist(x2[idx2_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x2[idx2_GKg_b].shape[0]))
    ax.hist(x2[idx2_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x2[idx2_GKg_c].shape[0]))
    ax.hist(x2[idx2_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x2[idx2_GKg_d].shape[0]))
    ax.set_title(str(s2Label) + ' Sloan Digital Sky Survey Apogee DR15 \n ' + str(s2Label) + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x1[idx_GKd_a].shape[0]))
    ax.hist(x1[idx_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x1[idx_GKd_b].shape[0]))
    ax.hist(x1[idx_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x1[idx_GKd_c].shape[0]))
    ax.hist(x1[idx_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x1[idx_GKd_d].shape[0]))
    ax.set_title('Sloan Digital Sky Survey Apogee DR15 \n ' + str(s1Label) + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x2[idx2_GKd_a].shape[0]))
    ax.hist(x2[idx2_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x2[idx2_GKd_b].shape[0]))
    ax.hist(x2[idx2_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x2[idx2_GKd_c].shape[0]))
    ax.hist(x2[idx2_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x2[idx2_GKd_d].shape[0]))
    ax.set_title('Sloan Digital Sky Survey Apogee DR15 \n ' + str(s2Label) + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    
    ################################################################################
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s1['GLON'], s1['GLAT'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s2['GLON'], s2['GLAT'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    ################################################################################
    
    x1 = s1['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    x2 = s2['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(0,0.5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    ################################################################################

    x1 = s1['LOGG']
    x2 = s2['LOGG']

    rLow = 0
    rHi  = 6

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('LOGG [log (cgs)]', fontsize=13)
    plt.title("LOGG \nlog g from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['M_H']
    x2 = s2['M_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('M_H [dex]', fontsize=13)
    plt.title("M_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']

    rLow = -6
    rHi  = 2

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('FE_H [dex]', fontsize=13)
    plt.title("FE_H \n[Fe/H] from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']

    rLow = 2500
    rHi  = 8500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('TEFF [K]', fontsize=13)
    plt.title("TEFF \nTeff from ASPCAP analysis of combined spectrum (from PARAM)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['SNR']
    x2 = s2['SNR']

    rLow = 0
    rHi  = 500

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('SNR', fontsize=13)
    plt.title("SNR \nmedian S/N per pixel in combined frame (at apStar sampling)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VHELIO_AVG [km/s]', fontsize=13)
    plt.title("VHELIO_AVG \naverage radial velocity, weighted by S/N, using RVs ", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VSCATTER']
    x2 = s2['VSCATTER']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VSCATTER [km/s]', fontsize=13)
    plt.title("VSCATTER \nscatter of individual visit RVs around average", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VMACRO']
    x2 = s2['VMACRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMACRO (cgs)', fontsize=13)
    plt.title("VMACRO \nmacroturbulent velocity (f(log Teff, [M/H]) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['VMICRO']
    x2 = s2['VMICRO']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('VMICRO (cgs)', fontsize=13)
    plt.title("VMICRO \nmicroturbulent velocity (fit for dwarfs, f(log g) for giants)", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    x1 = s1['J']
    x2 = s2['J']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('J', fontsize=13)
    plt.title("J \n2MASS J mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    
    
    x1 = s1['H']
    x2 = s2['H']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('H', fontsize=13)
    plt.title("H \n2MASS H mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    ################################################################################
    
    x1 = s1['K']
    x2 = s2['K']

    rLow = -2000
    rHi  = 2000

    plt.rcParams["figure.figsize"] = (12, 6)
    sns.distplot(x1[(x1 >= rLow) & (x1 <= rHi)], color='r', hist=True,  
                    label=s1Label.strip() + " Stars " + str(x1.shape[0]))
    sns.distplot(x2[(x2 >= rLow) & (x2 <= rHi)], color='b', hist=True, 
                    label=s2Label.strip() + " Stars " + str(x2.shape[0]))
    plt.xlabel('K', fontsize=13)
    plt.title("K \n2MASS K mag [bad=99]", fontsize=13)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.legend(fontsize=13)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13);
    plt.show()
    
    ################################################################################
    
    
    
    return (rTarget, rFlag, rBit, r2Target, r2Flag, r2Bit)









### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************

def apogee_target2RC(rTarget, rFlag, rBit, r2Target, r2Flag, r2Bit):

    hdulist = pyfits.open('apogee-rc-DR14.fits')
    #hdulist = pyfits.open('allStar-l31c.2.fits')
    
    print("Target: " + rTarget)
    print("Flag:   "   + rFlag)
    print("Bit: "    + str(rBit))
    
    print("\nTarget: " + r2Target)
    print("Flag:   "   + r2Flag)
    print("Bit: "    + str(r2Bit))
    
    cut_Apogee_Target = (hdulist[1].data[rFlag] & 2**rBit)
    idx_Apogee_Target = where(cut_Apogee_Target)[0]
    rApogee_Target = hdulist[1].data[idx_Apogee_Target]
    s1 = rApogee_Target
    s1Label = rTarget
    
    cut_Apogee_Target2 = (hdulist[1].data[r2Flag] & 2**r2Bit)
    idx_Apogee_Target2 = where(cut_Apogee_Target2)[0]
    r2Apogee_Target = hdulist[1].data[idx_Apogee_Target2]
    s2 = r2Apogee_Target
    s2Label = r2Target
    
    
    idx_Fd_a = where( (s1['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_b = where( (s1['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_c = where( (s1['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Fd_d = where( (s1['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Fd_a = where( (s2['ASPCAP_CLASS'] == 'Fd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_b = where( (s2['ASPCAP_CLASS'] == 'Fd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_c = where( (s2['ASPCAP_CLASS'] == 'Fd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Fd_d = where( (s2['ASPCAP_CLASS'] == 'Fd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_GKd_a = where( (s1['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_b = where( (s1['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_c = where( (s1['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKd_d = where( (s1['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKd_a = where( (s2['ASPCAP_CLASS'] == 'GKd_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_b = where( (s2['ASPCAP_CLASS'] == 'GKd_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_c = where( (s2['ASPCAP_CLASS'] == 'GKd_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKd_d = where( (s2['ASPCAP_CLASS'] == 'GKd_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]

    
    idx_GKg_a = where( (s1['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_b = where( (s1['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_c = where( (s1['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_GKg_d = where( (s1['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_GKg_a = where( (s2['ASPCAP_CLASS'] == 'GKg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_b = where( (s2['ASPCAP_CLASS'] == 'GKg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_c = where( (s2['ASPCAP_CLASS'] == 'GKg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_GKg_d = where( (s2['ASPCAP_CLASS'] == 'GKg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Md_a = where( (s1['ASPCAP_CLASS'] == 'Md_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_b = where( (s1['ASPCAP_CLASS'] == 'Md_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_c = where( (s1['ASPCAP_CLASS'] == 'Md_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Md_d = where( (s1['ASPCAP_CLASS'] == 'Md_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Md_a = where( (s2['ASPCAP_CLASS'] == 'Md_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_b = where( (s2['ASPCAP_CLASS'] == 'Md_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_c = where( (s2['ASPCAP_CLASS'] == 'Md_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Md_d = where( (s2['ASPCAP_CLASS'] == 'Md_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    
    
    idx_Mg_a = where( (s1['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_b = where( (s1['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_c = where( (s1['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    idx_Mg_d = where( (s1['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s1['APOGEE_TARGET2'] & 2**9))[0]
    
    idx2_Mg_a = where( (s2['ASPCAP_CLASS'] == 'Mg_a') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    idx2_Mg_b = where( (s2['ASPCAP_CLASS'] == 'Mg_b') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_b.shape[0] == 1:
            #idx_Mg_b = idx_Mg_b.iloc[1:]
    idx2_Mg_c = where( (s2['ASPCAP_CLASS'] == 'Mg_c') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    #if idx_Mg_c.shape[0] == 1:
            #idx_Mg_c = np.delete(idx_Mg_c, 0)
    idx2_Mg_d = where( (s2['ASPCAP_CLASS'] == 'Mg_d') & logical_not (s2['APOGEE_TARGET2'] & 2**9))[0]
    

    print('\n' + s1Label + ': Stars '   + str(s1.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx_Mg_d.shape[0]))
    
    
    print('\n' + s2Label + ': Stars '   + str(s2.shape[0]))
    
    print('\nFd_a' + ': Stars '   + str(idx2_Fd_a.shape[0]))
    print('Fd_b' + ': Stars '     + str(idx2_Fd_b.shape[0]))
    print('Fd_c' + ': Stars '     + str(idx2_Fd_c.shape[0]))
    print('Fd_d' + ': Stars '     + str(idx2_Fd_d.shape[0]))
    
    #print('\nFd' + ': Stars '   + str(len(idx_Fd)))
    
    print('\nGKd_a' + ': Stars '   + str(idx2_GKd_a.shape[0]))
    print('GKd_b' + ': Stars '     + str(idx2_GKd_b.shape[0]))
    print('GKd_c' + ': Stars '     + str(idx2_GKd_c.shape[0]))
    print('GKd_d' + ': Stars '     + str(idx2_GKd_d.shape[0]))
    
    #print('\nGKd' + ': Stars '   + str(len(idx_GKd)))
    
    print('\nGKg_a' + ': Stars '   + str(idx2_GKg_a.shape[0]))
    print('GKg_b' + ': Stars '     + str(idx2_GKg_b.shape[0]))
    print('GKg_c' + ': Stars '     + str(idx2_GKg_c.shape[0]))
    print('GKg_d' + ': Stars '     + str(idx2_GKg_d.shape[0]))
    
    #print('\nGKg' + ': Stars '   + str(len(idx_GKg)))
    
    print('\nMd_a' + ': Stars '   + str(idx2_Md_a.shape[0]))
    print('Md_b' + ': Stars '     + str(idx2_Md_b.shape[0]))
    print('Md_c' + ': Stars '     + str(idx2_Md_c.shape[0]))
    print('Md_d' + ': Stars '     + str(idx2_Md_d.shape[0]))
    
    #print('\nMd' + ': Stars '   + str(len(idx_Md)))
    
    print('\nMg_a' + ': Stars '   + str(idx2_Mg_a.shape[0]))
    print('Mg_b' + ': Stars '     + str(idx2_Mg_b.shape[0]))
    print('Mg_c' + ': Stars '     + str(idx2_Mg_c.shape[0]))
    print('Mg_d' + ': Stars '     + str(idx2_Mg_d.shape[0]))
    
    
    
    print('\nFields: ' + str(unique(s1['FIELD'])))
    
    print('\nFields: ' + str(unique(s2['FIELD'])))
    
    #################################################################################
    
    x1 = s1['FE_H']
    x2 = s2['FE_H']
    rbins = linspace(-3,1,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('FeH \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('FeH', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['LOGG']
    x2 = s2['LOGG']
    rbins = linspace(0,4,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('LogG \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('LogG', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()

    ################################################################################
    
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,8500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['VHELIO_AVG']
    x2 = s2['VHELIO_AVG']
    rbins = linspace(-500,500,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Vhelio_Avg \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Vhelio_Avg', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.legend(fontsize=12)
    #L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    y1 = s1['LOGG']
    x2 = s2['TEFF']
    y2 = s2['LOGG']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, y1,s=10, label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, y2,s=10, label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  Log G \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Log G', fontsize=15)
    ax.legend(fontsize=13)
    plt.xlim(6000, 3500)
    plt.ylim(0, 5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    show()
    
    ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()

 ###############################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.scatter(x1, j1 - k1,s=10, 
               label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2, j2 - k2,s=10, 
               label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Giants Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.5)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_Fd_a], j1[idx_Fd_a] - k1[idx_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x1[idx_Fd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_Fd_b], j1[idx_Fd_b] - k1[idx_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x1[idx_Fd_b].shape[0]), color='b')
    ax.scatter(x1[idx_Fd_c], j1[idx_Fd_c] - k1[idx_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x1[idx_Fd_c].shape[0]), color='r')
    ax.scatter(x1[idx_Fd_d], j1[idx_Fd_d] - k1[idx_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x1[idx_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_a], j2[idx2_Fd_a] - k2[idx2_Fd_a],s=10, 
               label='Fd_a  ' +  ' Stars '  + str(x2[idx2_Fd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_Fd_b], j2[idx2_Fd_b] - k2[idx2_Fd_b],s=10, 
               label='Fd_b  ' +  ' Stars '  + str(x2[idx2_Fd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_Fd_c], j2[idx2_Fd_c] - k2[idx2_Fd_c],s=10, 
               label='Fd_c  ' +  ' Stars '  + str(x2[idx2_Fd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_Fd_d], j2[idx2_Fd_d] - k2[idx2_Fd_d],s=10, 
               label='Fd_d  ' +  ' Stars '  + str(x2[idx2_Fd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKd_a], j1[idx_GKd_a] - k1[idx_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x1[idx_GKd_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKd_b], j1[idx_GKd_b] - k1[idx_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x1[idx_GKd_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKd_c], j1[idx_GKd_c] - k1[idx_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x1[idx_GKd_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKd_d], j1[idx_GKd_d] - k1[idx_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x1[idx_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_a], j2[idx2_GKd_a] - k2[idx2_GKd_a],s=10, 
               label='GKd_a  ' +  ' Stars '  + str(x2[idx2_GKd_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKd_b], j2[idx2_GKd_b] - k2[idx2_GKd_b],s=10, 
               label='GKd_b  ' +  ' Stars '  + str(x2[idx2_GKd_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKd_c], j2[idx2_GKd_c] - k2[idx2_GKd_c],s=10, 
               label='GKd_c  ' +  ' Stars '  + str(x2[idx2_GKd_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKd_d], j2[idx2_GKd_d] - k2[idx2_GKd_d],s=10, 
               label='GKd_d  ' +  ' Stars '  + str(x2[idx2_GKd_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    j1 = s1['J']
    k1 = s1['K']
    x2 = s2['TEFF']
    j2 = s2['J']
    k2 = s2['K']
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x1[idx_GKg_a], j1[idx_GKg_a] - k1[idx_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x1[idx_GKg_a].shape[0]), color='orange')
    ax.scatter(x1[idx_GKg_b], j1[idx_GKg_b] - k1[idx_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x1[idx_GKg_b].shape[0]), color='b')
    ax.scatter(x1[idx_GKg_c], j1[idx_GKg_c] - k1[idx_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x1[idx_GKg_c].shape[0]), color='r')
    ax.scatter(x1[idx_GKg_d], j1[idx_GKg_d] - k1[idx_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x1[idx_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s1Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(0,5,)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    #ax.scatter(x1, j1 - k1,s=10, 
               #label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_a], j2[idx2_GKg_a] - k2[idx2_GKg_a],s=10, 
               label='GKg_a  ' +  ' Stars '  + str(x2[idx2_GKg_a].shape[0]), color='orange')
    ax.scatter(x2[idx2_GKg_b], j2[idx2_GKg_b] - k2[idx2_GKg_b],s=10, 
               label='GKg_b  ' +  ' Stars '  + str(x2[idx2_GKg_b].shape[0]), color='b')
    ax.scatter(x2[idx2_GKg_c], j2[idx2_GKg_c] - k2[idx2_GKg_c],s=10, 
               label='GKg_c  ' +  ' Stars '  + str(x2[idx2_GKg_c].shape[0]), color='r')
    ax.scatter(x2[idx2_GKg_d], j2[idx2_GKg_d] - k2[idx2_GKg_d],s=10, 
               label='GKg_d  ' +  ' Stars '  + str(x2[idx2_GKg_d].shape[0]), color='g')
    
    ax.set_title(str(s2Label) + ' Teff  vs  J - K \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('J-K', fontsize=15)
    ax.legend(frameon=0, fontsize=13)
    plt.xlim(6500, 3500)
    plt.ylim(1.5, 0.0)
    plt.xticks(size = 13)
    plt.yticks(size = 13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    ################################################################################
        
    
    
    ################################################################################
        
    
    ################################################################################
    
    x1 = s1['TEFF']
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x1[idx_GKg_a].shape[0]))
    ax.hist(x1[idx_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x1[idx_GKg_b].shape[0]))
    ax.hist(x1[idx_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x1[idx_GKg_c].shape[0]))
    ax.hist(x1[idx_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x1[idx_GKg_d].shape[0]))
    ax.set_title(str(s1Label) +  ' Sloan Digital Sky Survey Apogee DR15 \n ' + str(s1Label) + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKg_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKg_a  ' + str(x2[idx2_GKg_a].shape[0]))
    ax.hist(x2[idx2_GKg_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKg_b  ' + str(x2[idx2_GKg_b].shape[0]))
    ax.hist(x2[idx2_GKg_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKg_c  ' + str(x2[idx2_GKg_c].shape[0]))
    ax.hist(x2[idx2_GKg_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKg_d  ' + str(x2[idx2_GKg_d].shape[0]))
    ax.set_title(str(s2Label) + ' Sloan Digital Sky Survey Apogee DR15 \n ' + str(s2Label) + ' GK Giants', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    ################################################################################
    
    x1 = s1['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1[idx_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x1[idx_GKd_a].shape[0]))
    ax.hist(x1[idx_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x1[idx_GKd_b].shape[0]))
    ax.hist(x1[idx_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x1[idx_GKd_c].shape[0]))
    ax.hist(x1[idx_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x1[idx_GKd_d].shape[0]))
    ax.set_title('Sloan Digital Sky Survey Apogee DR15 \n ' + str(s1Label) + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    x2 = s2['TEFF']
    rbins = linspace(3500,6500,25)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2[idx2_GKd_a], bins=rbins,histtype='step',linewidth=2, color='red',
            label='GKd_a  ' + str(x2[idx2_GKd_a].shape[0]))
    ax.hist(x2[idx2_GKd_b], bins=rbins,histtype='step',linewidth=2, color='blue',
            label='GKd_b  ' + str(x2[idx2_GKd_b].shape[0]))
    ax.hist(x2[idx2_GKd_c], bins=rbins,histtype='step',linewidth=2, color='yellow',
            label='GKd_c  ' + str(x2[idx2_GKd_c].shape[0]))
    ax.hist(x2[idx2_GKd_d], bins=rbins,histtype='step',linewidth=2,  color='green',
            label='GKd_d  ' + str(x2[idx2_GKd_d].shape[0]))
    ax.set_title('Sloan Digital Sky Survey Apogee DR15 \n ' + str(s2Label) + ' GK Dwarfs', fontsize=15)
    ax.set_xlabel('Teff', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    show()
    
    
    
    ################################################################################
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s1['GLON'], s1['GLAT'], color='r', label=s1Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s1Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (12, 10)
    plt.subplot(111)
    plt.scatter(s2['GLON'], s2['GLAT'], color='r', label=s2Label, marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()
    plt.title(str(s2Label), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    ################################################################################
    
    x1 = s1['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    x2 = s2['SURVEY']
    rbins = linspace(0,5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='r')
    ax.set_title('Survey \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Survey', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    ################################################################################
    
    x1 = s1['ALPHA_M']
    x2 = s2['ALPHA_M']
    rbins = linspace(0,0.5,50)
    rf=figure(figsize=(12,6))
    ax=rf.add_subplot(111)
    ax.hist(x1, bins=rbins,histtype='step',linewidth=2, 
            label=s1Label +  ' Stars '  + str(x1.shape[0]), color='r')
    ax.hist(x2, bins=rbins,histtype='step',linewidth=2, 
            label=s2Label +  ' Stars '  + str(x2.shape[0]), color='b')
    ax.set_title('Alpha_M \nSloan Digital Sky Survey Apogee DR15', fontsize=15)
    ax.set_xlabel('Alpha_M', fontsize=15)
    ax.set_ylabel('Count', fontsize=15)
    ax.legend(fontsize=12, loc="upper left")
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=14)
    show()
    
    
    return (rTarget, rFlag, rBit, r2Target, r2Flag, r2Bit)





### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************


def plotHistogramGC4(rSelected, rLow, rHi, rBins, rTeffLow, rTeffHi, rHelioLow, rHelioHi, fehLow, fehHi):

    rTarget =  rSelected[(rSelected.logg   >= rLow) & (rSelected.logg < rHi)] 
    
    plt.subplot(2, 3, 1)
    plt.hist([rTarget.c_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True,  color='r');
    plt.xlabel('[C/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1 + '\n Logg ' + str(rLow) + ' - ' + str(rHi), fontsize=14)

    plt.subplot(2, 3, 2)
    plt.hist([rTarget.n_fe],  
        bins=rBins, range=(-.5, 0.5), stacked=True,  color='b');
    plt.xlabel('[N/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1 + '\n Logg ' + str(rLow) + ' - ' + str(rHi), fontsize=14)

    plt.subplot(2, 3, 3)
    plt.hist([rTarget.na_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True,  color='g');
    plt.xlabel('[Na/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    plt.title(selected1 + '\n Logg ' + str(rLow) + ' - ' + str(rHi), fontsize=14)
    
    plt.show()

    plt.subplot(2, 3, 1)
    plt.hist([rTarget.mg_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True,  color='y');
    plt.xlabel('[Mg/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()

    plt.subplot(2, 3, 2)
    plt.hist([rTarget.al_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True, color='orange');
    plt.xlabel('[Al/Fe]', fontsize=15)
    plt.grid()

    plt.subplot(2, 3, 3)
    plt.hist([rTarget.fe_h], 
        bins=rBins, range=(-2.75, .6), stacked=True, color='m');
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()

    plt.show()
    
    
    plt.subplot(2, 3, 1)
    plt.hist([rTarget.o_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True,  color='purple');
    plt.xlabel('[O/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()

    plt.subplot(2, 3, 2)
    plt.hist([rTarget.ni_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True, color='k');
    plt.xlabel('[Ni/Fe]', fontsize=15)
    plt.grid()
    
    plt.subplot(2, 3, 3)
    plt.hist([rTarget.mn_fe], 
        bins=rBins, range=(-.5, 0.5), stacked=True, color='r');
    plt.xlabel('[Mn/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    
    
    plt.show()
    
    
    plt.subplot(2, 3, 1)
    plt.hist([rTarget.logg], 
        bins = rBins, range = (rLow, rHi), stacked = True,  color = 'k');
    plt.xlabel('Logg', fontsize=14)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    
    plt.subplot(2, 3, 2)
    plt.hist([rTarget.teff], 
        bins = rBins, range = (rTeffLow, rTeffHi), stacked = True, color = 'r');
    plt.xlabel('Teff', fontsize=14)
    plt.grid()
    
    plt.subplot(2, 3, 3)
    plt.hist([rTarget.vhelio_avg], 
        bins = rBins, range = (rHelioLow, rHelioHi),stacked = True,  color = 'b');
    plt.xlabel('vhelio_avg', fontsize=14)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    
    plt.show()
    
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.c_fe, rTarget.n_fe,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.6, 0.5)
    ax.set_ylim(-.6, 1.25)
    plt.xlabel('[C/Fe]', fontsize=15)
    plt.ylabel('[N/Fe]', fontsize=15)
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.n_fe, rTarget.na_fe,  color='b', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-.6, 1.25)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[N/Fe]', fontsize=15)
    plt.ylabel('[Na/Fe]', fontsize=15)
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.scatter(rTarget.mg_fe, rTarget.al_fe,  color='g', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.25, 0.5)
    ax.set_ylim(-.5, 1.25)
    plt.xlabel('[Mg/Fe]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    plt.grid()
    
    plt.subplot(2, 2, 4)
    plt.scatter(rTarget.n_fe, rTarget.al_fe,  color='y', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.4, 1.2)
    ax.set_ylim(-.5, 1.25)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[N/Fe]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    plt.grid()

    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.na_fe, rTarget.mg_fe,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-.6, 1.25)
    plt.xlabel('[Na/Fe]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.o_fe, rTarget.mg_fe,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.5, 1.0)
    ax.set_ylim(-.6, 1.25)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[O/Fe]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.si_fe, rTarget.mg_fe,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.5, .6)
    ax.set_ylim(-.6, 1.0)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[Si/Fe]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.mn_fe, rTarget.mg_fe,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-0.5, 1.0)
    ax.set_ylim(-.6, 1.25)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Mn/Fe]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.m_h, rTarget.logg,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(-3, 1)
    ax.set_ylim(0, rHi)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[M/H]', fontsize=15)
    plt.ylabel('Log G', fontsize=15)
    plt.grid()
    
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.alpha_m, rTarget.logg,  color='r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, rHi)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[a/M]', fontsize=15)
    plt.ylabel('Log G', fontsize=15)
    plt.grid()
    
    
    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.teff, rTarget.logg,  c=rTarget.teff, cmap='jet_r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(0, rHi)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('Log G', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.teff, rTarget.fe_h,  c=rTarget.teff, cmap='jet_r', marker='.')  #, linewidths=0)
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-2.5, 1.0)
    ax.yaxis.set_label_position("right")
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('[Fe/H]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    
    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.ra, rTarget.dec,  color='b', marker='.', linewidths=0)
    ax = plt.gca()
    #ax.set_xlim(rTeffLow, rTeffHi)
    #ax.set_ylim(-2.5, 1.0)
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    #plt.gca().invert_xaxis()
    plt.grid()
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.glon, rTarget.glat,  color='b', marker='.', linewidths=0)
    ax = plt.gca()
    #ax.set_xlim(rTeffLow, rTeffHi)
    #ax.set_ylim(-2.5, 1.0)
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    ax.yaxis.set_label_position("right")
    #plt.gca().invert_xaxis()
    plt.grid()
    
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 3, 1)
    plt.scatter(rTarget.fe_h, rTarget.c_fe + rTarget.n_fe,  color='r', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .5)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[C + N/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    
    #plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 3, 2)
    plt.scatter(rTarget.fe_h, rTarget.o_fe,  color='b', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .5)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[O/Fe]', fontsize=15)
    #plt.title('[O/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 3, 3)
    plt.scatter(rTarget.fe_h, rTarget.mg_fe,  color='r', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.4, .5)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    #plt.title('[O/Fe]', fontsize=15)
    plt.grid()
    
    plt.show()
    
    
    plt.subplot(1, 3, 1)
    plt.scatter(rTarget.fe_h, rTarget.al_fe,  color='r', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-1, .5)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    #plt.title('[O/Fe]', fontsize=15)
    plt.grid()
    

    plt.subplot(1, 3, 2)
    plt.scatter(rTarget.fe_h, rTarget.mn_fe,  color='b', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .6)
    #ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Mn/Fe]', fontsize=15)
    #plt.title('[O/Fe]', fontsize=15)
    plt.grid()
    
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 3, 3)
    plt.scatter(rTarget.fe_h, rTarget.ni_fe,  color='r', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.4, .3)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Ni/Fe]', fontsize=15)
    #plt.title('[O/Fe]', fontsize=15)
    plt.grid()
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()
    
    
    plt.subplot(1, 2, 1)
    plt.hist([rTarget.ra], 
        bins = rBins, range = (0, 360), stacked = True, 
             label = str(r1_select_variable.value), color = 'r');
    plt.xlabel('Ra', fontsize=14)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    
    plt.subplot(1, 2, 2)
    plt.hist([rTarget.dec], 
        bins = rBins, range = (-360, 360), stacked = True, 
             label = str(r1_select_variable.value), color = 'b');
    plt.xlabel('Dec', fontsize=14)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    
    plt.show()
    
    
    print(str(selected1) + '\n shape ' + str(rSelected.shape))


    
    
    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************


def plotHistogramStreams(rSelected, selected1, rLow, rHi, rBins, rTeffLow, rTeffHi, rHelioLow, rHelioHi, fehLow, fehHi, rLoc, ra1, ra2):

    cData = rSelected
    
    cFd_a =  cData[(cData.aspcap_class == 'Fd_a')]
    if cFd_a.shape[0] == 1:
            cFd_a = cFd_a.iloc[1:]
    cFd_b =  cData[(cData.aspcap_class == 'Fd_b')]
    if cFd_b.shape[0] == 1:
            cFd_b = cFd_b.iloc[1:]
    cFd_c =  cData[(cData.aspcap_class == 'Fd_c')]
    if cFd_c.shape[0] == 1:
            cFd_c = cFd_c.iloc[1:]
    cFd_d =  cData[(cData.aspcap_class == 'Fd_d')]
    if cFd_d.shape[0] == 1:
            cFd_d = cFd_d.iloc[1:]
            
            
    cGKd_a =  cData[(cData.aspcap_class == 'GKd_a')]
    if cGKd_a.shape[0] == 1:
            cGKd_a = cGKd_a.iloc[1:]
    cGKd_b =  cData[(cData.aspcap_class == 'GKd_b')]
    if cGKd_b.shape[0] == 1:
            cGKd_b = cGKd_b.iloc[1:]
    cGKd_c =  cData[(cData.aspcap_class == 'GKd_c')]
    if cGKd_c.shape[0] == 1:
            cGKd_c = cGKd_c.iloc[1:]
    cGKd_d =  cData[(cData.aspcap_class == 'GKd_d')]
    if cGKd_d.shape[0] == 1:
            cGKd_d = cGKd_d.iloc[1:]

            
    cGKg_a =  cData[(cData.aspcap_class == 'GKg_a')]
    if cGKg_a.shape[0] == 1:
            cGKg_a = cGKg_a.iloc[1:]
    cGKg_b =  cData[(cData.aspcap_class == 'GKg_b')]
    if cGKg_b.shape[0] == 1:
            cGKg_b = cGKg_b.iloc[1:]
    cGKg_c =  cData[(cData.aspcap_class == 'GKg_c')]
    if cGKg_c.shape[0] == 1:
            cGKg_c = cGKg_c.iloc[1:]
    cGKg_d =  cData[(cData.aspcap_class == 'GKg_d')]
    if cGKg_d.shape[0] == 1:
            cGKg_d = cGKg_d.iloc[1:]

            
    cMd_a =  cData[(cData.aspcap_class == 'Md_a')]
    if cMd_a.shape[0] == 1:
            cMd_a = cMd_a.iloc[1:]
    cMd_b =  cData[(cData.aspcap_class == 'Md_b')]
    if cMd_b.shape[0] == 1:
            cMd_b = cMd_b.iloc[1:]
    cMd_c =  cData[(cData.aspcap_class == 'Md_c')]
    if cMd_c.shape[0] == 1:
            cMd_c = cMd_c.iloc[1:]
    cMd_d =  cData[(cData.aspcap_class == 'Md_d')]
    if cMd_d.shape[0] == 1:
            cMd_d = cMd_d.iloc[1:]

            
    cMg_a =  cData[(cData.aspcap_class == 'Mg_a')]
    if cMg_a.shape[0] == 1:
            cMg_a = cMg_a.iloc[1:]
    cMg_b =  cData[(cData.aspcap_class == 'Mg_b')]
    if cMg_b.shape[0] == 1:
            cMg_b = cMg_b.iloc[1:]
    cMg_c =  cData[(cData.aspcap_class == 'Mg_c')]
    if cMg_c.shape[0] == 1:
            cMg_c = cMg_c.iloc[1:]
    cMg_d =  cData[(cData.aspcap_class == 'Mg_d')]
    if cMg_d.shape[0] == 1:
            cMg_d = cMg_d.iloc[1:]

    print(str(selected1) + '\n' + '\n' + 'Stars: ' + str(rSelected.shape[0]) + '\n')

    frames_cFd = [cFd_a, cFd_b, cFd_c, cFd_d]
    cFd = pd.concat(frames_cFd)
    if cFd.shape[0] == 1:
            cFd = cFd.iloc[1:]
    print("\nFd: " + str(cFd.shape[0]))

    frames_cGKd = [cGKd_a, cGKd_b, cGKd_c, cGKd_d]
    cGKd = pd.concat(frames_cGKd)
    if cGKd.shape[0] == 1:
            cGKd = cGKd.iloc[1:]
    print("GKd: " + str(cGKd.shape[0]))

    frames_cGKg = [cGKg_a, cGKg_b, cGKg_c, cGKg_d]
    cGKg = pd.concat(frames_cGKg)
    if cGKg.shape[0] == 1:
            cGKg = cGKg.iloc[1:]
    print("GKg: " + str(cGKg.shape[0]))

    frames_cMd = [cMd_a, cMd_b, cMd_c, cMd_d]
    cMd = pd.concat(frames_cMd)
    if cMd.shape[0] == 1:
            cMd = cMd.iloc[1:]
    print("Md: " + str(cMd.shape[0]))

    frames_cMg = [cMg_a, cMg_b, cMg_c, cMg_d]
    cMg = pd.concat(frames_cMg)
    if cMg.shape[0] == 1:
            cMg = cMg.iloc[1:]
    print("Mg: " + str(cMg.shape[0]))

    
    print("\nFd_a:  " + str(cFd_a.shape[0]))
    print("Fd_b:  " + str(cFd_b.shape[0]))
    print("Fd_c:  " + str(cFd_c.shape[0]))
    print("Fd_d:  " + str(cFd_d.shape[0]))
    print("\nGKd_a: " + str(cGKd_a.shape[0]))
    print("GKd_b: " + str(cGKd_b.shape[0]))
    print("GKd_c: " + str(cGKd_c.shape[0]))
    print("GKd_d: " + str(cGKd_d.shape[0]))
    print("\nGKg_a: " + str(cGKg_a.shape[0]))
    print("GKg_b: " + str(cGKg_b.shape[0]))
    print("GKg_c: " + str(cGKg_c.shape[0]))
    print("GKg_d: " + str(cGKg_d.shape[0]))
    print("\nMd_a:  " + str(cMd_a.shape[0]))
    print("Md_b:  " + str(cMd_b.shape[0]))
    print("Md_c:  " + str(cMd_c.shape[0]))
    print("Md_d:  " + str(cMd_d.shape[0]))
    print("\nMg_a:  " + str(cMg_a.shape[0]))
    print("Mg_b:  " + str(cMg_b.shape[0]))
    print("Mg_c:  " + str(cMg_c.shape[0]))
    print("Mg_d:  " + str(cMg_d.shape[0]))



    plt.subplot(1, 2, 1)
    plt.hist([cGKg.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='r', label='Class GKg')
    plt.hist([cGKd.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='b', label='Class GKd')
    plt.hist([cFd.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='orange', label='Class Fd')
    plt.hist([cMd.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='y', label='Class Md')
    plt.hist([cMg.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='g', label='Class Mg')
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    #L=plt.legend(loc=rLoc, bbox_to_anchor=(0.5, 1.00), fontsize=13)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()

    plt.subplot(1, 2, 1)
    plt.hist([cGKg_a.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='y', label='Class GKg_a')
    plt.hist([cGKg_b.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='b', label='Class GKg_b')
    plt.hist([cGKg_c.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='r', label='Class GKg_c')
    plt.hist([cGKg_d.teff], bins=rBins, range=(3500, 6500), stacked=True,  color='orange', label='Class GKg_d')
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc=rLoc, bbox_to_anchor=(ra1, ra2), fontsize=13)
    #L=plt.legend(loc='center left', bbox_to_anchor=(ra1, ra2), fontsize=13)


    plt.subplot(1, 2, 2)
    plt.hist([cGKg_a.logg], bins=rBins, range=(0, 5), stacked=True,  color='y', label='Class GKg_a')
    plt.hist([cGKg_b.logg], bins=rBins, range=(0, 5), stacked=True,  color='b', label='Class GKg_b')
    plt.hist([cGKg_c.logg], bins=rBins, range=(0, 5), stacked=True,  color='r', label='Class GKg_c')
    plt.hist([cGKg_d.logg], bins=rBins, range=(0, 5), stacked=True,  color='orange', label='Class GKg_d')
    plt.xlabel('Logg', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc=rLoc, bbox_to_anchor=(ra1, ra2), fontsize=13)
    #L=plt.legend(loc=rLoc, bbox_to_anchor=(0.5, 1.00), fontsize=13)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)


    plt.show()

    
    rTarget =  rSelected[(rSelected.logg   >= rLow) & (rSelected.logg < rHi)] 
    
    rTarget00 =  rSelected[(rSelected.logg   >= 0) & (rSelected.logg < 0.5)]
    if rTarget00.shape[0] == 1:
        rTarget00 = rTarget00.iloc[1:]
    
    rTarget05 =  rSelected[(rSelected.logg   >= 0.5) & (rSelected.logg < 1.0)]
    if rTarget05.shape[0] == 1:
        rTarget05 = rTarget05.iloc[1:]
    
    rTarget10 =  rSelected[(rSelected.logg   >= 1.0) & (rSelected.logg < 1.5)]
    if rTarget10.shape[0] == 1:
        rTarget10 = rTarget10.iloc[1:]
    
    rTarget15 =  rSelected[(rSelected.logg   >= 1.5) & (rSelected.logg < 2.0)]
    if rTarget15.shape[0] == 1:
        rTarget15 = rTarget15.iloc[1:]
    
    rTarget20 =  rSelected[(rSelected.logg   >= 2.0) & (rSelected.logg < 2.5)]
    if rTarget20.shape[0] == 1:
        rTarget20 = rTarget20.iloc[1:]
    
    rTarget25 =  rSelected[(rSelected.logg   >= 2.5) & (rSelected.logg < 3.0)]
    if rTarget25.shape[0] == 1:
        rTarget25 = rTarget25.iloc[1:]
    
    rTarget30 =  rSelected[(rSelected.logg   >= 3.0) & (rSelected.logg < 3.5)]
    if rTarget30.shape[0] == 1:
        rTarget30 = rTarget30.iloc[1:]
    
    rTarget35 =  rSelected[(rSelected.logg   >= 3.5) & (rSelected.logg < 4.0)]
    if rTarget35.shape[0] == 1:
        rTarget35 = rTarget35.iloc[1:]
    
    rTarget99 =  rSelected[(rSelected.logg   == -9999)]
    if rTarget99.shape[0] == 1:
        rTarget99 = rTarget99.iloc[1:]
    
    print('Logg [0-0.5]: '   + str(rTarget00.shape[0]))
    print('Logg [0.5-1]: '   + str(rTarget05.shape[0]))
    print('Logg [1-1.5]: '   + str(rTarget10.shape[0]))
    print('Logg [1.5-2]: '   + str(rTarget15.shape[0]))
    print('Logg [2-2.5]: '   + str(rTarget20.shape[0]))
    print('Logg [2.5-3]: '   + str(rTarget25.shape[0]))
    print('Logg [3-3.5]: '   + str(rTarget30.shape[0]))
    print('Logg [3.5-4]: '   + str(rTarget35.shape[0]))
    print('Logg [-9999]: '   + str(rTarget99.shape[0]))
    
    
    
    rTeff3000 =  rSelected[(rSelected.teff > -9999) & (rSelected.teff < 3500)]
    rTeff3500 =  rSelected[(rSelected.teff >= 3500) & (rSelected.teff < 4000)]
    rTeff4000 =  rSelected[(rSelected.teff >= 4000) & (rSelected.teff < 4500)]
    rTeff4500 =  rSelected[(rSelected.teff >= 4500) & (rSelected.teff < 5000)]
    rTeff5000 =  rSelected[(rSelected.teff >= 5000) & (rSelected.teff < 5500)]
    rTeff5500 =  rSelected[(rSelected.teff >= 5500) & (rSelected.teff < 6000)]
    rTeff6000 =  rSelected[(rSelected.teff >= 6000) & (rSelected.teff < 6500)]
    rTeff6500 =  rSelected[(rSelected.teff >= 6500) ]
    
    rTeff9999 =  rSelected[(rSelected.teff   == -9999)]
                
    print('\nTeff[<3500]: '       + str(rTeff3000.shape[0]))
    print('Teff[3500-4000]: '   + str(rTeff3500.shape[0]))
    print('Teff[4000-4500]: '   + str(rTeff4000.shape[0]))
    print('Teff[4500-5000]: '   + str(rTeff4500.shape[0]))
    print('Teff[5000-5500]: '   + str(rTeff5000.shape[0]))
    print('Teff[5500-6000]: '   + str(rTeff5500.shape[0]))
    print('Teff[6000-6500]: '   + str(rTeff6000.shape[0]))
    print('Teff[>6500]: '       + str(rTeff6500.shape[0]))
    print('Teff[-9999]: '       + str(rTeff9999.shape[0]))
    
    
    rFeH25 =  rSelected[(rSelected.fe_h >= -2.6) & (rSelected.fe_h < -2.0)]
    if rFeH25.shape[0] == 1:
        rFeH25 = rFeH25.iloc[1:]
    rFeH20 =  rSelected[(rSelected.fe_h >= -2.0) & (rSelected.fe_h < -1.5)]
    if rFeH20.shape[0] == 1:
        rFeH20 = rFeH20.iloc[1:]
    rFeH15 =  rSelected[(rSelected.fe_h >= -1.5) & (rSelected.fe_h < -1.0)]
    if rFeH15.shape[0] == 1:
        rFeH15 = rFeH15.iloc[1:]
    rFeH10 =  rSelected[(rSelected.fe_h >= -1.0) & (rSelected.fe_h < -0.5)]
    if rFeH10.shape[0] == 1:
        rFeH10 = rFeH10.iloc[1:]
    rFeH05 =  rSelected[(rSelected.fe_h >= -0.5) & (rSelected.fe_h < -0.0)]
    if rFeH05.shape[0] == 1:
        rFeH05 = rFeH05.iloc[1:]
    rFeH00 =  rSelected[(rSelected.fe_h >= -0.0) & (rSelected.fe_h < 0.5)]
    if rFeH00.shape[0] == 1:
        rFeH00 = rFeH00.iloc[1:]
    
    rFeH99 =  rSelected[(rSelected.fe_h   == -9999)]
                
    print('\nFeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))
    
    
    rVhelio_Neg_500 =  rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg < -400)]
    if rVhelio_Neg_500.shape[0] == 1:
        rVhelio_Neg_500 = rVhelio_Neg_500.iloc[1:]
    rVhelio_Neg_400 =  rSelected[(rSelected.vhelio_avg >= -400) & (rSelected.vhelio_avg < -300)]
    if rVhelio_Neg_400.shape[0] == 1:
        rVhelio_Neg_400 = rVhelio_Neg_400.iloc[1:]
    rVhelio_Neg_300 =  rSelected[(rSelected.vhelio_avg >= -300) & (rSelected.vhelio_avg < -200)]
    if rVhelio_Neg_300.shape[0] == 1:
        rVhelio_Neg_300 = rVhelio_Neg_300.iloc[1:]
    rVhelio_Neg_200 =  rSelected[(rSelected.vhelio_avg >= -200) & (rSelected.vhelio_avg < -100)]
    if rVhelio_Neg_200.shape[0] == 1:
        rVhelio_Neg_200 = rVhelio_Neg_200.iloc[1:]
    rVhelio_Neg_100 =  rSelected[(rSelected.vhelio_avg >= -100) & (rSelected.vhelio_avg < -0)]
    if rVhelio_Neg_100.shape[0] == 1:
        rVhelio_Neg_100 = rVhelio_Neg_100.iloc[1:]
    rVhelio_Pos_000 =  rSelected[(rSelected.vhelio_avg >= 0) & (rSelected.vhelio_avg < 100)]
    if rVhelio_Pos_000.shape[0] == 1:
        rVhelio_Pos_000 = rVhelio_Pos_000.iloc[1:]
    rVhelio_Pos_100 =  rSelected[(rSelected.vhelio_avg >= 100) & (rSelected.vhelio_avg < 200)]
    if rVhelio_Pos_100.shape[0] == 1:
        rVhelio_Pos_100 = rVhelio_Pos_100.iloc[1:]
    rVhelio_Pos_200 =  rSelected[(rSelected.vhelio_avg >= 200) & (rSelected.vhelio_avg < 300)]
    if rVhelio_Pos_200.shape[0] == 1:
        rVhelio_Pos_200 = rVhelio_Pos_200.iloc[1:]
    rVhelio_Pos_300 =  rSelected[(rSelected.vhelio_avg >= 300) & (rSelected.vhelio_avg < 400)]
    if rVhelio_Pos_300.shape[0] == 1:
        rVhelio_Pos_300 = rVhelio_Pos_300.iloc[1:]
    rVhelio_Pos_400 =  rSelected[(rSelected.vhelio_avg >= 400)]
    if rVhelio_Pos_400.shape[0] == 1:
        rVhelio_Pos_400 = rVhelio_Pos_400.iloc[1:]
    
    rVhelio9999 =  rSelected[(rSelected.vhelio_avg   == -9999)]
    if rVhelio9999.shape[0] == 1:
        rVhelio9999 = rVhelio9999.iloc[1:]
    
    print('\nvhelio [< -400]:         ' + str(rVhelio_Neg_500.shape[0]))
    print('vhelio [>= -400< -300]:  ' + str(rVhelio_Neg_400.shape[0]))
    print('vhelio [>= -300 < -200]: ' + str(rVhelio_Neg_300.shape[0]))
    print('vhelio [>= -200 < -100]: ' + str(rVhelio_Neg_200.shape[0]))
    print('vhelio [>= -100 < -0]:   ' + str(rVhelio_Neg_100.shape[0]))
    print('vhelio [>= 0 < 100]:     ' + str(rVhelio_Pos_000.shape[0]))
    print('vhelio [>= 100 < 200]:   ' + str(rVhelio_Pos_100.shape[0]))
    print('vhelio [>= 200 < 300]:   ' + str(rVhelio_Pos_200.shape[0]))
    print('vhelio [>= 300 < 400]:   ' + str(rVhelio_Pos_300.shape[0]))
    print('vhelio [>= 400]:         ' + str(rVhelio_Pos_400.shape[0]))
    print('vhelio [-9999]:          ' + str(rVhelio9999.shape[0]))
    
    ##########################################################################################################  
    plt.subplot(1, 2, 1)
    plt.hist([rTarget99.c_fe + rTarget99.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='k')
    plt.hist([rTarget00.c_fe + rTarget00.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='y')
    plt.hist([rTarget05.c_fe + rTarget05.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='m')
    plt.hist([rTarget10.c_fe + rTarget10.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='g')
    plt.hist([rTarget15.c_fe + rTarget15.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='purple')
    plt.hist([rTarget20.c_fe + rTarget20.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='b')
    plt.hist([rTarget25.c_fe + rTarget25.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='gold')
    plt.hist([rTarget30.c_fe + rTarget30.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='orange')
    plt.hist([rTarget35.c_fe + rTarget35.n_fe], bins=rBins, range=(-.6, .5), stacked=True, color='r')
    #plt.hist([rTarget99.c_fe + rTarget99.n_fe], bins=rBins, range=(-.6, .5), stacked=True,normed = False, color='k')
    plt.xlabel('[(C+N)/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1 + '\n Logg ' + str(rLow) + ' - ' + str(rHi), fontsize=14)

    plt.subplot(1, 2, 2)
    plt.hist([rTarget99.o_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.hist([rTarget00.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='y', label='Logg [0-0.5]')
    plt.hist([rTarget05.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='m', label='Logg [0.5-1]')
    plt.hist([rTarget10.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='g', label='Logg [1-1.5]')
    plt.hist([rTarget15.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='purple', label='Logg [1.5-2]')
    plt.hist([rTarget20.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='b', label='Logg [2-2.5]')
    plt.hist([rTarget25.o_fe], bins=rBins, range=(-.6, .5), stacked=True, color='gold', label='Logg [2.5-3]')
    plt.hist([rTarget30.o_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='orange', label='Logg [3-3.5]')
    plt.hist([rTarget35.o_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='r', label='Logg [3.5-4]')
    #plt.hist([rTarget99.o_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.xlabel('[O/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    plt.title(selected1 + '\n Logg ' + str(rLow) + ' - ' + str(rHi), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    plt.subplot(1, 2, 1)
    plt.hist([rTarget99.mg_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.hist([rTarget00.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='y', label='Logg [0-0.5]')
    plt.hist([rTarget05.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='m', label='Logg [0.5-1]')
    plt.hist([rTarget10.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='g', label='Logg [1-1.5]')
    plt.hist([rTarget15.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='purple', label='Logg [1.5-2]')
    plt.hist([rTarget20.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='b', label='Logg [2-2.5]')
    plt.hist([rTarget25.mg_fe], bins=rBins, range=(-.6, .5), stacked=True, color='gold', label='Logg [2.5-3]')
    plt.hist([rTarget30.mg_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='orange', label='Logg [3-3.5]')
    plt.hist([rTarget35.mg_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='r', label='Logg [3.5-4]')
    #plt.hist([rTarget99.mg_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.xlabel('[Mg/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.hist([rTarget99.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True,  color='k', label='Logg [-9999]')
    plt.hist([rTarget00.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='y', label='Logg [0-0.5]')
    plt.hist([rTarget05.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='m', label='Logg [0.5-1]')
    plt.hist([rTarget10.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='g', label='Logg [1-1.5]')
    plt.hist([rTarget15.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='purple', label='Logg [1.5-2]')
    plt.hist([rTarget20.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='b', label='Logg [2-2.5]')
    plt.hist([rTarget25.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True, color='gold', label='Logg [2.5-3]')
    plt.hist([rTarget30.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True,  color='orange', label='Logg [3-3.5]')
    plt.hist([rTarget35.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True,  color='r', label='Logg [3.5-4]')
    #plt.hist([rTarget99.al_fe], bins=rBins, range=(-.6, 1.5), stacked=True,  color='k', label='Logg [-9999]')
    plt.xlabel('[Al/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.subplot(1, 2, 1)
    plt.hist([rTarget99.mn_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.hist([rTarget00.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='y', label='Logg [0-0.5]')
    plt.hist([rTarget05.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='m', label='Logg [0.5-1]')
    plt.hist([rTarget10.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='g', label='Logg [1-1.5]')
    plt.hist([rTarget15.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='purple', label='Logg [1.5-2]')
    plt.hist([rTarget20.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='blue', label='Logg [2-2.5]')
    plt.hist([rTarget25.mn_fe], bins=rBins, range=(-.6, .5), stacked=True, color='gold', label='Logg [2.5-3]')
    plt.hist([rTarget30.mn_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='orange', label='Logg [3-3.5]')
    plt.hist([rTarget35.mn_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='r', label='Logg [3.5-4]')
    #plt.hist([rTarget99.mn_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.xlabel('[Mn/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.hist([rTarget99.ni_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.hist([rTarget00.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='y', label='Logg [0-0.5]')
    plt.hist([rTarget05.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='m', label='Logg [0.5-1]')
    plt.hist([rTarget10.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='g', label='Logg [1-1.5]')
    plt.hist([rTarget15.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='purple', label='Logg [1.5-2]')
    plt.hist([rTarget20.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='blue', label='Logg [2-2.5]')
    plt.hist([rTarget25.ni_fe], bins=rBins, range=(-.6, .5), stacked=True, color='gold', label='Logg [2.5-3]')
    plt.hist([rTarget30.ni_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='orange', label='Logg [3-3.5]')
    plt.hist([rTarget35.ni_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='r', label='Logg [3.5-4]')
    #plt.hist([rTarget99.ni_fe], bins=rBins, range=(-.6, .5), stacked=True,  color='k', label='Logg [-9999]')
    plt.xlabel('[Ni/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.subplot(1, 2, 1)
    plt.hist([rTarget.logg], 
        bins=rBins, range=(0, 4), stacked=True, color='r');
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    plt.show()
    
    ##########################################################################################################  
    
    plt.subplot(1, 2, 1)
    plt.hist([rVhelio9999.c_fe + rVhelio9999.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist([rVhelio_Neg_100.c_fe + rVhelio_Neg_100.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='b', label='vhelio [+- 0 - 100]')
    plt.hist([rVhelio_Pos_000.c_fe + rVhelio_Pos_000.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist([rVhelio_Neg_300.c_fe + rVhelio_Neg_300.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio [+- 300 - 200]')
    plt.hist([rVhelio_Pos_200.c_fe + rVhelio_Pos_200.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist([rVhelio_Pos_100.c_fe + rVhelio_Pos_100.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='r', label='vhelio [+- 100 - 200]')
    plt.hist([rVhelio_Neg_200.c_fe + rVhelio_Neg_200.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist([rVhelio_Pos_300.c_fe + rVhelio_Pos_300.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio [+- 300 - 400]')
    plt.hist([rVhelio_Neg_400.c_fe + rVhelio_Neg_400.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist([rVhelio_Pos_400.c_fe + rVhelio_Pos_400.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio [+- 400]')
    plt.hist([rVhelio_Neg_500.c_fe + rVhelio_Neg_500.n_fe], bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[(C+N)/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)

    plt.subplot(1, 2, 2)
    plt.hist(rVhelio9999.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.o_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[O/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    plt.subplot(1, 2, 1)
    plt.hist(rVhelio9999.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
              color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.mg_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[Mg/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    #plt.title(selected1, fontsize=14)

    plt.subplot(1, 2, 2)
    plt.hist(rVhelio9999.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.al_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[Al/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    #plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.subplot(1, 2, 1)
    plt.hist(rVhelio9999.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
              color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
              color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
              color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.mn_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[Mn/Fe]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    #plt.title(selected1, fontsize=14)

    plt.subplot(1, 2, 2)
    plt.hist(rVhelio9999.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.ni_fe, bins=rBins, range=(-.6, .5), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[Ni/Fe]', fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    #plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    ##########################################################################################################  
 
    
    plt.subplot(1, 2, 1)
    plt.hist(rVhelio9999.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='k', label='vhelio [-9999]')
    plt.hist(rVhelio_Neg_100.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='b', label='vhelio +- [0-100]')
    plt.hist(rVhelio_Pos_000.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='b') # , label='vhelio [>= 0 < 100]')
    plt.hist(rVhelio_Neg_300.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='y', label='vhelio +- [200-300]')
    plt.hist(rVhelio_Pos_200.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='y') # , label='vhelio [>= 200 < 300]')
    plt.hist(rVhelio_Pos_100.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='r', label='vhelio +- [100-200]')
    plt.hist(rVhelio_Neg_200.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='r') # , label='vhelio [>= -200 < -100]')
    plt.hist(rVhelio_Pos_300.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='orange', label='vhelio +- [300-400]')
    plt.hist(rVhelio_Neg_400.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='orange') #, label='vhelio [>= -400< -300]')
    plt.hist(rVhelio_Pos_400.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='g', label='vhelio +- [400]')
    plt.hist(rVhelio_Neg_500.logg, bins=rBins, range=(0, 4), stacked=True, 
             color='g') # , label='vhelio [< -400]')
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    #ax.yaxis.set_label_position("right")
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.logg, rTarget.c_fe + rTarget.n_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[C + N/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.logg, rTarget.o_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.6, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[O/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.logg, rTarget.mg_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.4, .6)
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.logg, rTarget.al_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-1.5, 1.5)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    plt.grid()
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.logg, rTarget.mn_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[Mn/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.logg, rTarget.ni_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[Ni/Fe]', fontsize=15)
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.logg, rTarget.c_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[C/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.logg, rTarget.n_fe,  c=rTarget.teff, cmap='jet', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[N/Fe]', fontsize=15)
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH25.logg, rFeH25.c_fe, color='r', marker='o') 
    plt.scatter(rFeH20.logg, rFeH20.c_fe, color='orange', marker='o') 
    plt.scatter(rFeH15.logg, rFeH15.c_fe, color='y', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[C/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rFeH25.logg, rFeH25.n_fe, color='r', marker='o') 
    plt.scatter(rFeH20.logg, rFeH20.n_fe, color='orange', marker='o') 
    plt.scatter(rFeH15.logg, rFeH15.n_fe, color='y', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rLow, rHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Log G]', fontsize=15)
    plt.ylabel('[N/Fe]', fontsize=15)
    plt.grid()
    plt.show()
    
    
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 2, 1)
    plt.hist([rTarget.teff], 
        bins=rBins, range=(3500, 5500), stacked=True, color='r');
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.teff, rTarget.c_fe + rTarget.n_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[C + N/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.teff, rTarget.o_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-.6, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[O/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.teff, rTarget.mg_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-.4, .6)
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.teff, rTarget.al_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-1.5, 1.5)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.teff, rTarget.mn_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[Mn/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.teff, rTarget.ni_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Teff]', fontsize=15)
    plt.ylabel('[Ni/Fe]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.teff, rTarget.logg,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(0, rHi)
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('Log G', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.teff, rTarget.fe_h,  c=rTarget.teff, cmap='jet_r', marker='o')  
    ax = plt.gca()
    ax.set_xlim(rTeffLow, rTeffHi)
    ax.set_ylim(-2.5, 1.0)
    ax.yaxis.set_label_position("right")
    plt.xlabel('Teff', fontsize=15)
    plt.ylabel('[Fe/H]', fontsize=15)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 4)
    plt.subplot(1, 2, 1)
    plt.hist([rTarget.fe_h], 
        bins=rBins, range=(-2.5, .5), stacked=True, color='r');
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('Count', fontsize=15)
    ax = plt.gca()
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.fe_h, rTarget.c_fe + rTarget.n_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[C + N/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.fe_h, rTarget.o_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[O/Fe]', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.fe_h, rTarget.mg_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.4, .6)
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Mg/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.fe_h, rTarget.al_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-1.5, 1.5)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Al/Fe]', fontsize=15)
    plt.grid()
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.fe_h, rTarget.mn_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Mn/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.fe_h, rTarget.ni_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[Ni/Fe]', fontsize=15)
    plt.grid()
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.fe_h, rTarget.c_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.6, .6)
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[C/Fe]', fontsize=15)
    plt.grid()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.fe_h, rTarget.n_fe,  c=rTarget.teff, cmap='jet_r', marker='o') 
    ax = plt.gca()
    ax.set_xlim(fehLow, fehHi)
    ax.set_ylim(-.4, .6)
    ax.yaxis.set_label_position("right")
    plt.xlabel('[Fe/H]', fontsize=15)
    plt.ylabel('[N/Fe]', fontsize=15)
    plt.grid()
    plt.show()
    
    
    ####################
    
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget.o_fe, rTarget.na_fe,  color='r', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_ylim(-1, 1)
    ax.set_xlim(-0.5, .5)
    plt.xlabel('[O/Fe]', fontsize=15)
    plt.ylabel('[Na/Fe]', fontsize=15)
    plt.grid()
    plt.title(str(selected1), fontsize=13)
    
    plt.subplot(2, 2, 2)
    plt.scatter(rTarget.c_fe, rTarget.n_fe,  color='b', marker='.', linewidths=0)
    ax = plt.gca()
    ax.set_ylim(-0.5, 1)
    ax.set_xlim(-0.5, .5)
    plt.xlabel('[C/Fe]', fontsize=15)
    plt.ylabel('[N/Fe]', fontsize=15)
    plt.grid()
    plt.title(str(selected1), fontsize=13)
    
    plt.show()



    
    

    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************
    


def plotOverview2(r, r1):

    cData = r

    cFd_a =  cData[(cData.aspcap_class == 'Fd_a')]
    cFd_b =  cData[(cData.aspcap_class == 'Fd_b')]
    cFd_c =  cData[(cData.aspcap_class == 'Fd_c')]
    cFd_d =  cData[(cData.aspcap_class == 'Fd_d')]
    cGKd_a =  cData[(cData.aspcap_class == 'GKd_a')]
    cGKd_b =  cData[(cData.aspcap_class == 'GKd_b')]
    cGKd_c =  cData[(cData.aspcap_class == 'GKd_c')]
    cGKd_d =  cData[(cData.aspcap_class == 'GKd_d')]
    cGKg_a =  cData[(cData.aspcap_class == 'GKg_a')]
    cGKg_b =  cData[(cData.aspcap_class == 'GKg_b')]
    cGKg_c =  cData[(cData.aspcap_class == 'GKg_c')]
    cGKg_d =  cData[(cData.aspcap_class == 'GKg_d')]
    cMd_a =  cData[(cData.aspcap_class == 'Md_a')]
    cMd_b =  cData[(cData.aspcap_class == 'Md_b')]
    cMd_c =  cData[(cData.aspcap_class == 'Md_c')]
    cMd_d =  cData[(cData.aspcap_class == 'Md_d')]
    cMg_a =  cData[(cData.aspcap_class == 'Mg_a')]
    cMg_b =  cData[(cData.aspcap_class == 'Mg_b')]
    cMg_c =  cData[(cData.aspcap_class == 'Mg_c')]
    cMg_d =  cData[(cData.aspcap_class == 'Mg_d')]

    frames_cFd = [cFd_a, cFd_b, cFd_c, cFd_d]
    cFd = pd.concat(frames_cFd)

    frames_cGKd = [cGKd_a, cGKd_b, cGKd_c, cGKd_d]
    cGKd = pd.concat(frames_cGKd)

    frames_cGKg = [cGKg_a, cGKg_b, cGKg_c, cGKg_d]
    cGKg = pd.concat(frames_cGKg)

    frames_cMd = [cMd_a, cMd_b, cMd_c, cMd_d]
    cMd = pd.concat(frames_cMd)

    frames_cMg = [cMg_a, cMg_b, cMg_c, cMg_d]
    cMg = pd.concat(frames_cMg)

 
    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(cFd.ra, cFd.dec, color='b', label='Class Fd', marker='*', s=40)
    plt.scatter(cGKd.ra, cGKd.dec, color='y', label='Class GKd', marker='*', s=40)
    plt.scatter(cGKg.ra, cGKg.dec, color='r', label='Class GKg', marker='o', s=40)
    plt.scatter(cMd.ra, cMd.dec, color='g', label='Class Md', marker='*', s=40)
    plt.scatter(cMg.ra, cMg.dec, color='orange', label='Class Mg', marker='o', s=40)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(r1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()


    print(str(r1) + '\n' + '\n' + 'Stars: ' + str(r.shape[0]) + '\n')
    print("Fd: " + str(cFd.shape[0]))
    print("GKd: " + str(cGKd.shape[0]))
    print("GKg: " + str(cGKg.shape[0]))
    print("Md: " + str(cMd.shape[0]))
    print("Mg: " + str(cMg.shape[0]))



    rTarget00 =  r[(r.logg   >= 0) & (r.logg < 0.5)]
    if rTarget00.shape[0] == 1:
            rTarget00 = rTarget00.iloc[1:]
    
    rTarget05 =  r[(r.logg   >= 0.5) & (r.logg < 1.0)]
    if rTarget05.shape[0] == 1:
            rTarget05 = rTarget05.iloc[1:]
    
    rTarget10 =  r[(r.logg   >= 1.0) & (r.logg < 1.5)]
    if rTarget10.shape[0] == 1:
            rTarget10 = rTarget10.iloc[1:]

    rTarget15 =  r[(r.logg   >= 1.5) & (r.logg < 2.0)]
    if rTarget15.shape[0] == 1:
            rTarget15 = rTarget15.iloc[1:]

    rTarget20 =  r[(r.logg   >= 2.0) & (r.logg < 2.5)]
    if rTarget20.shape[0] == 1:
            rTarget20 = rTarget20.iloc[1:]

    rTarget25 =  r[(r.logg   >= 2.5) & (r.logg < 3.0)]
    if rTarget25.shape[0] == 1:
            rTarget25 = rTarget25.iloc[1:]

    rTarget30 =  r[(r.logg   >= 3.0) & (r.logg < 3.5)]
    if rTarget30.shape[0] == 1:
            rTarget30 = rTarget30.iloc[1:]

    rTarget35 =  r[(r.logg   >= 3.5) & (r.logg < 4.0)]
    if rTarget35.shape[0] == 1:
            rTarget35 = rTarget35.iloc[1:]

    rTarget99 =  r[(r.logg   == -9999)]
    if rTarget99.shape[0] == 1:
            rTarget99 = rTarget99.iloc[1:]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget99.ra, rTarget99.dec, color='k', label='Logg [-9999]', marker='x', s=20)
    plt.scatter(rTarget00.ra, rTarget00.dec, color='b', label='Logg [0-0.5]', marker='o', s=80)
    plt.scatter(rTarget05.ra, rTarget05.dec, color='m', label='Logg [0.5-1]', marker='o', s=80)
    plt.scatter(rTarget10.ra, rTarget10.dec, color='y', label='Logg [1-1.5]', marker='*', s=40)
    plt.scatter(rTarget15.ra, rTarget15.dec, color='purple', label='Logg [1.5-2]', marker='*', s=40)
    plt.scatter(rTarget20.ra, rTarget20.dec, color='blue', label='Logg [2-2.5]', marker='*', s=40)
    plt.scatter(rTarget25.ra, rTarget25.dec, color='orange', label='Logg [3-3.5]', marker='*', s=40)
    plt.scatter(rTarget35.ra, rTarget35.dec, color='r', label='Logg [3.5-4]', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(r1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()

    print(str(r1) + '\n' + '\n' + 'Stars: ' + str(r.shape[0]) + '\n')

    print('Logg [0-0.5]: '   + str(rTarget00.shape[0]))
    print('Logg [0.5-1]: '   + str(rTarget05.shape[0]))
    print('Logg [1-1.5]: '   + str(rTarget10.shape[0]))
    print('Logg [1.5-2]: '   + str(rTarget15.shape[0]))
    print('Logg [2-2.5]: '   + str(rTarget20.shape[0]))
    print('Logg [2.5-3]: '   + str(rTarget25.shape[0]))
    print('Logg [3-3.5]: '   + str(rTarget30.shape[0]))
    print('Logg [3.5-4]: '   + str(rTarget35.shape[0]))
    print('Logg [-9999]: '   + str(rTarget99.shape[0]))


    sFeh26 =  r[(r.fe_h >= -2.6) & (r.fe_h < -2.5)]
    sFeh25 =  r[(r.fe_h >= -2.5) & (r.fe_h < -2.4)]
    sFeh24 =  r[(r.fe_h >= -2.4) & (r.fe_h < -2.3)]
    sFeh23 =  r[(r.fe_h >= -2.3) & (r.fe_h < -2.2)]
    sFeh22 =  r[(r.fe_h >= -2.2) & (r.fe_h < -2.1)]
    sFeh21 =  r[(r.fe_h >= -2.1) & (r.fe_h < -2.0)]
    sFeh20 =  r[(r.fe_h >= -2.0) & (r.fe_h < -1.9)]
    sFeh19 =  r[(r.fe_h >= -1.9) & (r.fe_h < -1.8)]
    sFeh18 =  r[(r.fe_h >= -1.8) & (r.fe_h < -1.7)]
    sFeh17 =  r[(r.fe_h >= -1.7) & (r.fe_h < -1.6)]
    sFeh16 =  r[(r.fe_h >= -1.6) & (r.fe_h < -1.5)]
    sFeh15 =  r[(r.fe_h >= -1.5) & (r.fe_h < -1.4)]
    sFeh14 =  r[(r.fe_h >= -1.4) & (r.fe_h < -1.3)]
    sFeh13 =  r[(r.fe_h >= -1.3) & (r.fe_h < -1.2)]
    sFeh12 =  r[(r.fe_h >= -1.2) & (r.fe_h < -1.1)]
    sFeh11 =  r[(r.fe_h >= -1.1) & (r.fe_h < -1.0)]

    rFeH25 =  r[(r.fe_h >= -2.6) & (r.fe_h < -2.0)]
    rFeH20 =  r[(r.fe_h >= -2.0) & (r.fe_h < -1.5)]
    rFeH15 =  r[(r.fe_h >= -1.5) & (r.fe_h < -1.0)]
    rFeH10 =  r[(r.fe_h >= -1.0) & (r.fe_h < -0.5)]
    rFeH05 =  r[(r.fe_h >= -0.5) & (r.fe_h < -0.0)]
    rFeH00 =  r[(r.fe_h >= -0.0) & (r.fe_h < 0.5)]
    rFeH99 =  r[(r.fe_h   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH99.ra, rFeH99.dec, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.ra, rFeH25.dec, color='r', label='FeH [-2.6-2.0]', marker='o', s=80)
    plt.scatter(rFeH20.ra, rFeH20.dec, color='orange', label='FeH [-2.0-1.5]', marker='o', s=80)
    plt.scatter(rFeH15.ra, rFeH15.dec, color='y', label='FeH [-1.5-1.0]', marker='o', s=80)
    plt.scatter(rFeH10.ra, rFeH10.dec, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    plt.scatter(rFeH05.ra, rFeH05.dec, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    plt.scatter(rFeH00.ra, rFeH00.dec, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15), plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(r1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    print('\nFeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))

    print('\nFeH [-2.6-2.5]: ' + str(sFeh26.shape[0]))
    print('FeH [-2.5-2.4]: '   + str(sFeh25.shape[0]))
    print('FeH [-2.4-2.3]: '   + str(sFeh24.shape[0]))
    print('FeH [-2.3-2.2]: '   + str(sFeh23.shape[0]))
    print('FeH [-2.2-2.1]: '   + str(sFeh22.shape[0]))
    print('FeH [-2.1-2.0]: '   + str(sFeh21.shape[0]))
    print('FeH [-2.0-1.9]: '   + str(sFeh20.shape[0]))
    print('FeH [-1.9-1.8]: '   + str(sFeh19.shape[0]))
    print('FeH [-1.8-1.7]: '   + str(sFeh18.shape[0]))
    print('FeH [-1.7-1.6]: '   + str(sFeh17.shape[0]))
    print('FeH [-1.6-1.5]: '   + str(sFeh16.shape[0]))
    print('FeH [-1.5-1.4]: '   + str(sFeh15.shape[0]))
    print('FeH [-1.4-1.3]: '   + str(sFeh14.shape[0]))
    print('FeH [-1.3-1.2]: '   + str(sFeh13.shape[0]))
    print('FeH [-1.2-1.1]: '   + str(sFeh12.shape[0]))
    print('FeH [-1.1-1.0]: '   + str(sFeh11.shape[0]))


    rTeff3000 =  r[(r.teff > -9999) & (r.teff < 3500)]
    rTeff3500 =  r[(r.teff >= 3500) & (r.teff < 4000)]
    rTeff4000 =  r[(r.teff >= 4000) & (r.teff < 4500)]
    rTeff4500 =  r[(r.teff >= 4500) & (r.teff < 5000)]
    rTeff5000 =  r[(r.teff >= 5000) & (r.teff < 5500)]
    rTeff5500 =  r[(r.teff >= 5500) & (r.teff < 6000)]
    rTeff6000 =  r[(r.teff >= 6000) & (r.teff < 6500)]
    rTeff6500 =  r[(r.teff >= 6500) ]

    rTeff9999 =  r[(r.teff   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rTeff9999.ra, rTeff9999.dec, color='k', label='Teff[-9999]', marker='x', s=20)
    plt.scatter(rTeff3500.ra, rTeff3500.dec, color='r', label='Teff[3500-4000]', marker='o', s=80)
    plt.scatter(rTeff4000.ra, rTeff4000.dec, color='orange', label='Teff[4000-4500]', marker='o', s=80)
    plt.scatter(rTeff4500.ra, rTeff4500.dec, color='y', label='Teff[4500-5000]', marker='*', s=40)
    plt.scatter(rTeff5000.ra, rTeff5000.dec, color='purple', label='Teff[5000-5500]', marker='*', s=40)
    plt.scatter(rTeff5500.ra, rTeff5500.dec, color='m', label='Teff[5500-6000]', marker='*', s=40)
    plt.scatter(rTeff6000.ra, rTeff6000.dec, color='b', label='Teff[6000-6500]', marker='*', s=40)
    plt.scatter(rTeff6500.ra, rTeff6500.dec, color='b', label='Teff[>6500]', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15)
    plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(r1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()

    print('\nTeff[<3500]: '     + str(rTeff3000.shape[0]))
    print('Teff[3500-4000]: '   + str(rTeff3500.shape[0]))
    print('Teff[4000-4500]: '   + str(rTeff4000.shape[0]))
    print('Teff[4500-5000]: '   + str(rTeff4500.shape[0]))
    print('Teff[5000-5500]: '   + str(rTeff5000.shape[0]))
    print('Teff[5500-6000]: '   + str(rTeff5500.shape[0]))
    print('Teff[6000-6500]: '   + str(rTeff6000.shape[0]))
    print('Teff[>6500]: '       + str(rTeff6500.shape[0]))
    print('Teff[-9999]: '       + str(rTeff9999.shape[0]))


    #rVhelio_Draco =  rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg <= -200)]

    #rVhelio_Pal1 =  rSelected[(rSelected.vhelio_avg >= -90) & (rSelected.vhelio_avg <= -60)]
    #rVhelio_Pal2 =  rVhelio_Pal1[(rVhelio_Pal1.fe_h > -1.0) & (rVhelio_Pal1.fe_h < -0.2)]

    rVhelio_Neg_500 =  r[(r.vhelio_avg > -9999) & (r.vhelio_avg < -400)]
    rVhelio_Neg_400 =  r[(r.vhelio_avg >= -400) & (r.vhelio_avg < -300)]
    rVhelio_Neg_300 =  r[(r.vhelio_avg >= -300) & (r.vhelio_avg < -200)]
    rVhelio_Neg_200 =  r[(r.vhelio_avg >= -200) & (r.vhelio_avg < -100)]
    rVhelio_Neg_100 =  r[(r.vhelio_avg >= -100) & (r.vhelio_avg < -0)]
    rVhelio_Pos_000 =  r[(r.vhelio_avg >= 0) & (r.vhelio_avg < 100)]
    rVhelio_Pos_100 =  r[(r.vhelio_avg >= 100) & (r.vhelio_avg < 200)]
    rVhelio_Pos_200 =  r[(r.vhelio_avg >= 200) & (r.vhelio_avg < 300)]
    rVhelio_Pos_300 =  r[(r.vhelio_avg >= 300) & (r.vhelio_avg < 400)]
    rVhelio_Pos_400 =  r[(r.vhelio_avg >= 400)]
    rVhelio9999 =  r[(r.vhelio_avg   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rVhelio9999.ra, rVhelio9999.dec, color='k', label='vhelio [-9999]', marker='x', s=20)
    plt.scatter(rVhelio_Neg_500.ra, rVhelio_Neg_500.dec, color='r', label='vhelio [< -400]', marker='o', s=80)
    plt.scatter(rVhelio_Neg_400.ra, rVhelio_Neg_400.dec, color='orange', label='vhelio [>= -400< -300]', 
                marker='o', s=80)
    plt.scatter(rVhelio_Neg_300.ra, rVhelio_Neg_300.dec, color='y', label='vhelio [>= -300 < -200]', marker='o', s=80)
    plt.scatter(rVhelio_Neg_200.ra, rVhelio_Neg_200.dec, color='r', label='vhelio [>= -200 < -100]',  marker='*', s=60)
    plt.scatter(rVhelio_Neg_100.ra, rVhelio_Neg_100.dec, color='b', label='vhelio [>= -100 < -0]', marker='*', s=40)
    plt.scatter(rVhelio_Pos_000.ra, rVhelio_Pos_000.dec, color='b', label='vhelio [>= 0 < 100]', marker='*', s=40)
    plt.scatter(rVhelio_Pos_100.ra, rVhelio_Pos_100.dec, color='r', label='vhelio [>= 100 < 200]', marker='*', s=60)
    plt.scatter(rVhelio_Pos_200.ra, rVhelio_Pos_200.dec, color='y', label='vhelio [>= 200 < 300]', marker='o', s=80)
    plt.scatter(rVhelio_Pos_300.ra, rVhelio_Pos_300.dec, color='orange', label='vhelio [>= 300 < 400]', 
                    marker='o', s=80)
    plt.scatter(rVhelio_Pos_400.ra, rVhelio_Pos_400.dec, color='r', label='vhelio [>= 400]', marker='o', s=80)

    #plt.scatter(rVhelio_Pal1.ra, rVhelio_Pal1.dec, color='g', label='Pal1 vhelio [>= -90 <= -60]', 
                #marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15), plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(r1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()


    print('\nvhelio [< -400]:         ' + str(rVhelio_Neg_500.shape[0]))
    print('vhelio [>= -400< -300]:  ' + str(rVhelio_Neg_400.shape[0]))
    print('vhelio [>= -300 < -200]: ' + str(rVhelio_Neg_300.shape[0]))
    print('vhelio [>= -200 < -100]: ' + str(rVhelio_Neg_200.shape[0]))
    print('vhelio [>= -100 < -0]:   ' + str(rVhelio_Neg_100.shape[0]))
    print('vhelio [>= 0 < 100]:     ' + str(rVhelio_Pos_000.shape[0]))
    print('vhelio [>= 100 < 200]:   ' + str(rVhelio_Pos_100.shape[0]))
    print('vhelio [>= 200 < 300]:   ' + str(rVhelio_Pos_200.shape[0]))
    print('vhelio [>= 300 < 400]:   ' + str(rVhelio_Pos_300.shape[0]))
    print('vhelio [>= 400]:         ' + str(rVhelio_Pos_400.shape[0]))
    print('vhelio [-9999]:          ' + str(rVhelio9999.shape[0]))



    
    
    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************


def plotCompare(rX, rXlabel, rY, rYlabel, rS2, rFilter, rLogg, rFeH, rTeff, rVhelio, rBins, rLogMin, rLogMax):
    s2 = rS2
    
    if rFilter == "Vhelio":
        x = rX[(rX.vhelio_avg > rVhelio)]
    
    if rFilter == "FeH":
        x = rX[(rX.fe_h > rFeH)]
        
    if rFilter == "Teff":
        x = rX[(rX.teff > rTeff)]
        
    if rFilter == "Logg":
        x = rX[(rX.logg > rLogg)]

    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    
    if rFilter == "Vhelio":
        x = rX[(rX.vhelio_avg > -9999) & (rX.vhelio_avg <= rVhelio)]
    if rFilter == "FeH":
        x = rX[(rX.fe_h > -9999) & (rX.fe_h <= rFeH)]
    if rFilter == "Teff":
        x = rX[(rX.teff > -9999)  & (rX.teff < rTeff)]
    if rFilter == "Logg":
        x = rX[(rX.logg > -9999)  & (rX.logg < rLogg)]
            
    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    if rElementRatio_Neg_06.shape[0] == 1:
        rElementRatio_Neg_06 = rElementRatio_Neg_06.iloc[1:]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    if rElementRatio_Neg_05.shape[0] == 1:
        rElementRatio_Neg_05 = rElementRatio_Neg_05.iloc[1:]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    if rElementRatio_Neg_04.shape[0] == 1:
        rElementRatio_Neg_04 = rElementRatio_Neg_04.iloc[1:]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    if rElementRatio_Neg_03.shape[0] == 1:
        rElementRatio_Neg_03 = rElementRatio_Neg_03.iloc[1:]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    if rElementRatio_Neg_02.shape[0] == 1:
        rElementRatio_Neg_02 = rElementRatio_Neg_02.iloc[1:]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    if rElementRatio_Neg_01.shape[0] == 1:
        rElementRatio_Neg_01 = rElementRatio_Neg_01.iloc[1:]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    if rElementRatio_Pos_00.shape[0] == 1:
        rElementRatio_Pos_00 = rElementRatio_Pos_00.iloc[1:]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    if rElementRatio_Pos_01.shape[0] == 1:
        rElementRatio_Pos_01 = rElementRatio_Pos_01.iloc[1:]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    if rElementRatio_Pos_02.shape[0] == 1:
        rElementRatio_Pos_02 = rElementRatio_Pos_02.iloc[1:]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    if rElementRatio_Pos_03.shape[0] == 1:
        rElementRatio_Pos_03 = rElementRatio_Pos_03.iloc[1:]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    if rElementRatio_Pos_04.shape[0] == 1:
        rElementRatio_Pos_04 = rElementRatio_Pos_04.iloc[1:]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    if rElementRatio_Pos_05.shape[0] == 1:
        rElementRatio_Pos_05 = rElementRatio_Pos_05.iloc[1:]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    

    xElementRatio_Neg_04 = eval('rElementRatio_Neg_04.'+ str(s2))
    xElementRatio_Neg_03 = eval('rElementRatio_Neg_03.'+ str(s2))
    xElementRatio_Neg_02 = eval('rElementRatio_Neg_02.'+ str(s2))
    xElementRatio_Neg_01 = eval('rElementRatio_Neg_01.'+ str(s2))
    xElementRatio_Pos_00 = eval('rElementRatio_Pos_00.'+ str(s2))
    xElementRatio_Pos_01 = eval('rElementRatio_Pos_01.'+ str(s2))
    xElementRatio_Pos_02 = eval('rElementRatio_Pos_02.'+ str(s2))
    xElementRatio_Pos_03 = eval('rElementRatio_Pos_03.'+ str(s2))
    xElementRatio_Pos_04 = eval('rElementRatio_Pos_04.'+ str(s2))
    xElementRatio_Pos_05 = eval('rElementRatio_Pos_05.'+ str(s2))
    
    
    #####################################################################################################
    
    xElementRatio = eval('rX.'+ str(s2))
    yElementRatio = eval('rY.'+ str(s2))
    
    plt.rcParams["figure.figsize"] = (25, 10)
    plt.subplot(2, 3, 1)
    plt.hist([xElementRatio], range=(-.5, 1.0), stacked=True, color='r', bins=rBins,
         label=str(rXlabel))
    plt.xlabel(s2, fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.title(str(rXlabel), fontsize=15)
    plt.grid()
    #L=plt.legend(loc='upper center', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (25, 10)
    plt.subplot(2, 3, 2)
    plt.hist([yElementRatio], range=(-.5, 1.0), stacked=True, color='r', bins=rBins,
         label=str(rYlabel))
    plt.xlabel(s2, fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    plt.title(str(rYlabel), fontsize=15)
    #L=plt.legend(loc='upper center', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    #xElementRatio2 = xElementRatio > -9999
    #yElementRatio2 = yElementRatio > -9999
    
    #sns.kdeplot(xElementRatio2, shade=True, label=str(rXlabel), color='red');
    #sns.kdeplot(yElementRatio2, shade=True, label=str(rYlabel), color='blue');
    #plt.show()
    
    #sns.distplot(xElementRatio, color='red');
    #sns.distplot(yElementRatio, color='blue');
    #plt.show()
    
    
    
    #####################################################################################################
  
    
    rFeH25 =  rX[(rX.fe_h >= -2.6) & (rX.fe_h < -2.0)]
    rFeH20 =  rX[(rX.fe_h >= -2.0) & (rX.fe_h < -1.5)]
    rFeH15 =  rX[(rX.fe_h >= -1.5) & (rX.fe_h < -1.0)]
    rFeH10 =  rX[(rX.fe_h >= -1.0) & (rX.fe_h < -0.5)]
    rFeH05 =  rX[(rX.fe_h >= -0.5) & (rX.fe_h < -0.0)]
    rFeH00 =  rX[(rX.fe_h >= -0.0) & (rX.fe_h < 0.5)]
    rFeH99 =  rX[(rX.fe_h   == -9999)]
    
    sFeH25 =  rY[(rY.fe_h >= -2.6) & (rY.fe_h < -2.0)]
    sFeH20 =  rY[(rY.fe_h >= -2.0) & (rY.fe_h < -1.5)]
    sFeH15 =  rY[(rY.fe_h >= -1.5) & (rY.fe_h < -1.0)]
    sFeH10 =  rY[(rY.fe_h >= -1.0) & (rY.fe_h < -0.5)]
    sFeH05 =  rY[(rY.fe_h >= -0.5) & (rY.fe_h < -0.0)]
    sFeH00 =  rY[(rY.fe_h >= -0.0) & (rY.fe_h < 0.5)]
    sFeH99 =  rY[(rY.fe_h   == -9999)]

    xFeH25 = eval('rFeH25.'+ str(s2))
    xFeH20 = eval('rFeH20.'+ str(s2))
    xFeH15 = eval('rFeH15.'+ str(s2))
    xFeH10 = eval('rFeH10.'+ str(s2))
    xFeH05 = eval('rFeH05.'+ str(s2))
    xFeH00 = eval('rFeH00.'+ str(s2))
    xFeH99 = eval('rFeH99.'+ str(s2))
    
    yFeH25 = eval('sFeH25.'+ str(s2))
    yFeH20 = eval('sFeH20.'+ str(s2))
    yFeH15 = eval('sFeH15.'+ str(s2))
    yFeH10 = eval('sFeH10.'+ str(s2))
    yFeH05 = eval('sFeH05.'+ str(s2))
    yFeH00 = eval('sFeH00.'+ str(s2))
    yFeH99 = eval('sFeH99.'+ str(s2))
    

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, xFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 5500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(sFeH25.teff, yFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(sFeH20.teff, yFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(sFeH15.teff, yFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, xFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 5500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()
    
    
    #####################################################################################################
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.logg, xFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(sFeH99.logg, yxFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(sFeH25.logg, yFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(sFeH20.logg, yFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(sFeH15.logg, yFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(sFeH10.logg, yFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(sFeH05.logg, yFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(sFeH00.logg, yFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    
    plt.show()
    
    
    #####################################################################################################
    
    
    

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, rFeH05.logg, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(sFeH25.teff, sFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(sFeH20.teff, sFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(sFeH15.teff, sFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, rFeH05.logg, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()
    
    ####################################################################################################

    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    plt.scatter(sFeH99.teff, sFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(sFeH25.teff, sFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(sFeH20.teff, sFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(sFeH15.teff, sFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(sFeH10.teff, sFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(sFeH05.teff, sFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(sFeH00.teff, sFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    
    plt.show()
    
    
    ####################################################################################################


    print('\n' + str(rXlabel))
    print('FeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))
    
    print('\n' + str(rYlabel))
    print('FeH [-2.6-2.0]: '   + str(sFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(sFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(sFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(sFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(sFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(sFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(sFeH99.shape[0]))
    
    
    
    rTarget00 =  rX[(rX.logg   >= 0) & (rX.logg < 0.5)]
    rTarget05 =  rX[(rX.logg   >= 0.5) & (rX.logg < 1.0)]
    rTarget10 =  rX[(rX.logg   >= 1.0) & (rX.logg < 1.5)]
    rTarget15 =  rX[(rX.logg   >= 1.5) & (rX.logg < 2.0)]
    rTarget20 =  rX[(rX.logg   >= 2.0) & (rX.logg < 2.5)]
    rTarget25 =  rX[(rX.logg   >= 2.5) & (rX.logg < 3.0)]
    rTarget30 =  rX[(rX.logg   >= 3.0) & (rX.logg < 3.5)]
    rTarget35 =  rX[(rX.logg   >= 3.5) & (rX.logg < 4.0)]
    rTarget99 =  rX[(rX.logg   == -9999)]

    print('\n' + str(rXlabel))
    print('Logg [0-0.5]: '   + str(rTarget00.shape[0]))
    print('Logg [0.5-1]: '   + str(rTarget05.shape[0]))
    print('Logg [1-1.5]: '   + str(rTarget10.shape[0]))
    print('Logg [1.5-2]: '   + str(rTarget15.shape[0]))
    print('Logg [2-2.5]: '   + str(rTarget20.shape[0]))
    print('Logg [2.5-3]: '   + str(rTarget25.shape[0]))
    print('Logg [3-3.5]: '   + str(rTarget30.shape[0]))
    print('Logg [3.5-4]: '   + str(rTarget35.shape[0]))
    print('Logg [-9999]: '   + str(rTarget99.shape[0]))
    
    
    
    
    sTarget00 =  rY[(rY.logg   >= 0) & (rY.logg < 0.5)]
    sTarget05 =  rY[(rY.logg   >= 0.5) & (rY.logg < 1.0)]
    sTarget10 =  rY[(rY.logg   >= 1.0) & (rY.logg < 1.5)]
    sTarget15 =  rY[(rY.logg   >= 1.5) & (rY.logg < 2.0)]
    sTarget20 =  rY[(rY.logg   >= 2.0) & (rX.logg < 2.5)]
    sTarget25 =  rY[(rY.logg   >= 2.5) & (rX.logg < 3.0)]
    sTarget30 =  rY[(rY.logg   >= 3.0) & (rY.logg < 3.5)]
    sTarget35 =  rY[(rY.logg   >= 3.5) & (rY.logg < 4.0)]
    sTarget99 =  rY[(rY.logg   == -9999)]

    print('\n' + str(rYlabel))
    print('Logg [0-0.5]: '   + str(sTarget00.shape[0]))
    print('Logg [0.5-1]: '   + str(sTarget05.shape[0]))
    print('Logg [1-1.5]: '   + str(sTarget10.shape[0]))
    print('Logg [1.5-2]: '   + str(sTarget15.shape[0]))
    print('Logg [2-2.5]: '   + str(sTarget20.shape[0]))
    print('Logg [2.5-3]: '   + str(sTarget25.shape[0]))
    print('Logg [3-3.5]: '   + str(sTarget30.shape[0]))
    print('Logg [3.5-4]: '   + str(sTarget35.shape[0]))
    print('Logg [-9999]: '   + str(sTarget99.shape[0]))
    
    
    
    
    ####################################################################################################

    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 5500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(sFeH99.teff, yFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(sFeH25.teff, yFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(sFeH20.teff, yFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(sFeH15.teff, yFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(sFeH10.teff, yFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(sFeH05.teff, yFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(sFeH00.teff, yFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 5500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()
    
    
    ####################################################################################################


    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.logg, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(sFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(sFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(sFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(sFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(sFeH10.logg, yFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(sFeH05.logg, yFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(sFeH00.logg, yFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    #plt.xlim(0, 4)
    plt.xlim(rLogMax, rLogMin)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()

    ####################################################################################################


    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rXlabel), fontsize=14)
    #L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 2)
    #plt.scatter(sFeH99.teff, sFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(sFeH25.teff, sFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(sFeH20.teff, sFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(sFeH15.teff, sFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(sFeH10.teff, sFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(sFeH05.teff, sFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(sFeH00.teff, sFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.ylim(rLogMax, rLogMin)
    #plt.gca().set_ylim(bottom=0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(str(rYlabel), fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    
    plt.show()


    
    
    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************


def plotMetric(rSelected, selected1, s2, rFilter, rLogg, rFeH, rTeff, rVhelio, rBins):
    #s2 = r2_select_variable.value
    
    if rFilter == "Vhelio":
        x = rSelected[(rSelected.vhelio_avg > rVhelio)]
    
    if rFilter == "FeH":
        x = rSelected[(rSelected.fe_h > rFeH)]
        
    if rFilter == "Teff":
        x = rSelected[(rSelected.teff > rTeff)]
        
    if rFilter == "Logg":
        x = rSelected[(rSelected.logg > rLogg)]

    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    #plt.scatter(rElementRatio_Neg_9999.ra, rElementRatio_Neg_9999.dec, color='k', label='_nolegend_', 
            #marker='x', s=20)
    plt.scatter(rElementRatio_Neg_06.ra, rElementRatio_Neg_06.dec, color='orange', label='_nolegend_',      
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_05.ra, rElementRatio_Neg_05.dec, color='orange', label='_nolegend_',  
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_04.ra, rElementRatio_Neg_04.dec, color='gold', label='_nolegend_',     
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_03.ra, rElementRatio_Neg_03.dec, color='gold', label='_nolegend_',  
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_02.ra, rElementRatio_Neg_02.dec, color='b', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_01.ra, rElementRatio_Neg_01.dec, color='b', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_00.ra, rElementRatio_Pos_00.dec, color='r', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_01.ra, rElementRatio_Pos_01.dec, color='r', label='_nolegend_',     
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_02.ra, rElementRatio_Pos_02.dec, color='g', label='_nolegend_',   
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_03.ra, rElementRatio_Pos_03.dec, color='g', label='_nolegend_',    
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_04.ra, rElementRatio_Pos_04.dec, color='m', label='_nolegend_',    
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_05.ra, rElementRatio_Pos_05.dec, color='m', label='_nolegend_',
            marker='*', s=40)
    
    if rFilter == "Vhelio":
        x = rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg <= rVhelio)]
    if rFilter == "FeH":
        x = rSelected[(rSelected.fe_h > -9999) & (rSelected.fe_h <= rFeH)]
    if rFilter == "Teff":
        x = rSelected[(rSelected.teff > -9999)  & (rSelected.teff < rTeff)]
    if rFilter == "Logg":
        x = rSelected[(rSelected.logg > -9999)  & (rSelected.logg < rLogg)]
            
    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    if rElementRatio_Neg_06.shape[0] == 1:
        rElementRatio_Neg_06 = rElementRatio_Neg_06.iloc[1:]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    if rElementRatio_Neg_05.shape[0] == 1:
        rElementRatio_Neg_05 = rElementRatio_Neg_05.iloc[1:]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    if rElementRatio_Neg_04.shape[0] == 1:
        rElementRatio_Neg_04 = rElementRatio_Neg_04.iloc[1:]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    if rElementRatio_Neg_03.shape[0] == 1:
        rElementRatio_Neg_03 = rElementRatio_Neg_03.iloc[1:]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    if rElementRatio_Neg_02.shape[0] == 1:
        rElementRatio_Neg_02 = rElementRatio_Neg_02.iloc[1:]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    if rElementRatio_Neg_01.shape[0] == 1:
        rElementRatio_Neg_01 = rElementRatio_Neg_01.iloc[1:]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    if rElementRatio_Pos_00.shape[0] == 1:
        rElementRatio_Pos_00 = rElementRatio_Pos_00.iloc[1:]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    if rElementRatio_Pos_01.shape[0] == 1:
        rElementRatio_Pos_01 = rElementRatio_Pos_01.iloc[1:]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    if rElementRatio_Pos_02.shape[0] == 1:
        rElementRatio_Pos_02 = rElementRatio_Pos_02.iloc[1:]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    if rElementRatio_Pos_03.shape[0] == 1:
        rElementRatio_Pos_03 = rElementRatio_Pos_03.iloc[1:]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    if rElementRatio_Pos_04.shape[0] == 1:
        rElementRatio_Pos_04 = rElementRatio_Pos_04.iloc[1:]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    if rElementRatio_Pos_05.shape[0] == 1:
        rElementRatio_Pos_05 = rElementRatio_Pos_05.iloc[1:]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    #plt.scatter(rElementRatio_Neg_9999.ra, rElementRatio_Neg_9999.dec, color='k', label='[Ni/Fe = -9999]', 
            #marker='x', s=20)
    plt.scatter(rElementRatio_Neg_06.ra, rElementRatio_Neg_06.dec, color='orange', label='_nolegend_', 
                marker='o', s=80)
    plt.scatter(rElementRatio_Neg_05.ra, rElementRatio_Neg_05.dec, color='orange', label = s2 + ' [>= -0.5 < -0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_04.ra, rElementRatio_Neg_04.dec, color='gold', label= s2 + '[>= -0.4 < -0.2]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_03.ra, rElementRatio_Neg_03.dec, color='gold', label='_nolegend_', 
                marker='o', s=80)
    plt.scatter(rElementRatio_Neg_02.ra, rElementRatio_Neg_02.dec, color='b', label= s2 + '[>= -0.2 < -0.0]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_01.ra, rElementRatio_Neg_01.dec, color='b', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_00.ra, rElementRatio_Pos_00.dec, color='r', label= s2 + '[>= 0.0 < 0.2]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_01.ra, rElementRatio_Pos_01.dec, color='r', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_02.ra, rElementRatio_Pos_02.dec, color='g', label= s2 + '[>= 0.2 < 0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_03.ra, rElementRatio_Pos_03.dec, color='g', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_04.ra, rElementRatio_Pos_04.dec, color='m', label= s2 + '[>= 0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_05.ra, rElementRatio_Pos_05.dec, color='m', label='_nolegend_', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Ra', fontsize=15), plt.ylabel('Dec', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()

    print('\n' + s2 + ' [< -0.5]: ' + str(rElementRatio_Neg_06.shape[0]))
    print(s2 + ' [>= -0.5) & < -0.4]: '   + str(rElementRatio_Neg_05.shape[0]))
    print(colored(s2 +  ' [>= -0.4) & < -0.3]: '   + str(rElementRatio_Neg_04.shape[0]), 'yellow'))
    print(colored(s2 +  ' [>= -0.3) & < -0.2]: '   + str(rElementRatio_Neg_03.shape[0]), 'yellow'))
    print(colored(s2 +  ' [>= -0.2) & < -0.1]: '   + str(rElementRatio_Neg_02.shape[0]), 'blue'))
    print(colored(s2 +  ' [>= -0.1) & < -0.0]: '   + str(rElementRatio_Neg_01.shape[0]), 'blue'))
    print(colored(s2 +  ' [>=  0.0) & <  0.1]: '   + str(rElementRatio_Pos_00.shape[0]), 'red'))
    print(colored(s2 +  ' [>=  0.1) & <  0.2]: '   + str(rElementRatio_Pos_01.shape[0]), 'red'))
    print(colored(s2 +  ' [>=  0.2) & <  0.3]: '   + str(rElementRatio_Pos_02.shape[0]), 'green'))
    print(colored(s2 +  ' [>=  0.3) & <  0.4]: '   + str(rElementRatio_Pos_03.shape[0]), 'green'))
    print(colored(s2 +  ' [>=  0.4) & <  0.5]: '   + str(rElementRatio_Pos_04.shape[0]), 'magenta'))
    print(colored(s2 +  ' [>=  0.5]: '   + str(rElementRatio_Pos_05.shape[0]), 'magenta'))
    print(s2 +  ' [-9999]: '   + str(rElementRatio_Neg_9999.shape[0]))

    xElementRatio_Neg_04 = eval('rElementRatio_Neg_04.'+ str(s2))
    xElementRatio_Neg_03 = eval('rElementRatio_Neg_03.'+ str(s2))
    xElementRatio_Neg_02 = eval('rElementRatio_Neg_02.'+ str(s2))
    xElementRatio_Neg_01 = eval('rElementRatio_Neg_01.'+ str(s2))
    xElementRatio_Pos_00 = eval('rElementRatio_Pos_00.'+ str(s2))
    xElementRatio_Pos_01 = eval('rElementRatio_Pos_01.'+ str(s2))
    xElementRatio_Pos_02 = eval('rElementRatio_Pos_02.'+ str(s2))
    xElementRatio_Pos_03 = eval('rElementRatio_Pos_03.'+ str(s2))
    xElementRatio_Pos_04 = eval('rElementRatio_Pos_04.'+ str(s2))
    xElementRatio_Pos_05 = eval('rElementRatio_Pos_05.'+ str(s2))
    
    plt.rcParams["figure.figsize"] = (25, 10)
    plt.subplot(2, 3, 1)
    plt.hist([xElementRatio_Neg_04], range=(-.5, 1.0), stacked=True, color='y', bins=rBins,
         label=s2 +  ' [>= -0.4) & < -0.2]')
    plt.hist([xElementRatio_Neg_03], range=(-.5, 1.0), stacked=True, color='y', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Neg_02], range=(-.5, 1.0), stacked=True, color='b', bins=rBins,
         label=s2 +  ' [>= -0.2) & < -0.0]')
    plt.hist([xElementRatio_Neg_01], range=(-.5, 1.0), stacked=True, color='b', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_00], range=(-.5, 1.0), stacked=True, color='r', bins=rBins,
         label=s2 +  ' [>=  0.0) & <  0.2]')
    plt.hist([xElementRatio_Pos_01], range=(-.5, 1.0), stacked=True, color='r', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_02], range=(-.5, 1.0), stacked=True, color='g', bins=rBins,
         label=s2 +  ' [>=  0.2) & <  0.4]')
    plt.hist([xElementRatio_Pos_03], range=(-.5, 1.0), stacked=True, color='g', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_04], range=(-.5, 1.0), stacked=True, color='m', bins=rBins, label=s2 + ' [>=  0.4]')
    plt.hist([xElementRatio_Pos_05], range=(-.5, 1.0), stacked=True, color='m', bins=rBins, label='_nolegend_')
    plt.xlabel(s2, fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    rFeH25 =  rSelected[(rSelected.fe_h >= -2.6) & (rSelected.fe_h < -2.0)]
    rFeH20 =  rSelected[(rSelected.fe_h >= -2.0) & (rSelected.fe_h < -1.5)]
    rFeH15 =  rSelected[(rSelected.fe_h >= -1.5) & (rSelected.fe_h < -1.0)]
    rFeH10 =  rSelected[(rSelected.fe_h >= -1.0) & (rSelected.fe_h < -0.5)]
    rFeH05 =  rSelected[(rSelected.fe_h >= -0.5) & (rSelected.fe_h < -0.0)]
    rFeH00 =  rSelected[(rSelected.fe_h >= -0.0) & (rSelected.fe_h < 0.5)]
    rFeH99 =  rSelected[(rSelected.fe_h   == -9999)]

    xFeH25 = eval('rFeH25.'+ str(s2))
    xFeH20 = eval('rFeH20.'+ str(s2))
    xFeH15 = eval('rFeH15.'+ str(s2))
    xFeH10 = eval('rFeH10.'+ str(s2))
    xFeH05 = eval('rFeH05.'+ str(s2))
    xFeH00 = eval('rFeH00.'+ str(s2))
    xFeH99 = eval('rFeH99.'+ str(s2))

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, xFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 6500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.logg, xFeH05, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, rFeH05.logg, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    print('\nFeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 6500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.logg, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    



    
    
    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************


def plotOverviewGlon(rSelected, selected1):

    rTarget00 =  rSelected[(rSelected.logg   >= 0) & (rSelected.logg < 0.5)]
    if rTarget00.shape[0] == 1:
            rTarget00 = rTarget00.iloc[1:]
    
    rTarget05 =  rSelected[(rSelected.logg   >= 0.5) & (rSelected.logg < 1.0)]
    if rTarget05.shape[0] == 1:
            rTarget05 = rTarget05.iloc[1:]
    
    rTarget10 =  rSelected[(rSelected.logg   >= 1.0) & (rSelected.logg < 1.5)]
    if rTarget10.shape[0] == 1:
            rTarget10 = rTarget10.iloc[1:]

    rTarget15 =  rSelected[(rSelected.logg   >= 1.5) & (rSelected.logg < 2.0)]
    if rTarget15.shape[0] == 1:
            rTarget15 = rTarget15.iloc[1:]

    rTarget20 =  rSelected[(rSelected.logg   >= 2.0) & (rSelected.logg < 2.5)]
    if rTarget20.shape[0] == 1:
            rTarget20 = rTarget20.iloc[1:]

    rTarget25 =  rSelected[(rSelected.logg   >= 2.5) & (rSelected.logg < 3.0)]
    if rTarget25.shape[0] == 1:
            rTarget25 = rTarget25.iloc[1:]

    rTarget30 =  rSelected[(rSelected.logg   >= 3.0) & (rSelected.logg < 3.5)]
    if rTarget30.shape[0] == 1:
            rTarget30 = rTarget30.iloc[1:]

    rTarget35 =  rSelected[(rSelected.logg   >= 3.5) & (rSelected.logg < 4.0)]
    if rTarget35.shape[0] == 1:
            rTarget35 = rTarget35.iloc[1:]

    rTarget99 =  rSelected[(rSelected.logg   == -9999)]
    if rTarget99.shape[0] == 1:
            rTarget99 = rTarget99.iloc[1:]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rTarget99.glon, rTarget99.glat, color='k', label='Logg [-9999]', marker='x', s=20)
    plt.scatter(rTarget00.glon, rTarget00.glat, color='b', label='Logg [0-0.5]', marker='o', s=80)
    plt.scatter(rTarget05.glon, rTarget05.glat, color='m', label='Logg [0.5-1]', marker='o', s=80)
    plt.scatter(rTarget10.glon, rTarget10.glat, color='y', label='Logg [1-1.5]', marker='*', s=40)
    plt.scatter(rTarget15.glon, rTarget15.glat, color='purple', label='Logg [1.5-2]', marker='*', s=40)
    plt.scatter(rTarget20.glon, rTarget20.glat, color='blue', label='Logg [2-2.5]', marker='*', s=40)
    plt.scatter(rTarget25.glon, rTarget25.glat, color='orange', label='Logg [3-3.5]', marker='*', s=40)
    plt.scatter(rTarget35.glon, rTarget35.glat, color='r', label='Logg [3.5-4]', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()

    print(str(selected1) + '\n' + '\n' + 'Stars: ' + str(rSelected.shape[0]) + '\n')

    print('Logg [0-0.5]: '   + str(rTarget00.shape[0]))
    print('Logg [0.5-1]: '   + str(rTarget05.shape[0]))
    print('Logg [1-1.5]: '   + str(rTarget10.shape[0]))
    print('Logg [1.5-2]: '   + str(rTarget15.shape[0]))
    print('Logg [2-2.5]: '   + str(rTarget20.shape[0]))
    print('Logg [2.5-3]: '   + str(rTarget25.shape[0]))
    print('Logg [3-3.5]: '   + str(rTarget30.shape[0]))
    print('Logg [3.5-4]: '   + str(rTarget35.shape[0]))
    print('Logg [-9999]: '   + str(rTarget99.shape[0]))


    sFeh26 =  rSelected[(rSelected.fe_h >= -2.6) & (rSelected.fe_h < -2.5)]
    sFeh25 =  rSelected[(rSelected.fe_h >= -2.5) & (rSelected.fe_h < -2.4)]
    sFeh24 =  rSelected[(rSelected.fe_h >= -2.4) & (rSelected.fe_h < -2.3)]
    sFeh23 =  rSelected[(rSelected.fe_h >= -2.3) & (rSelected.fe_h < -2.2)]
    sFeh22 =  rSelected[(rSelected.fe_h >= -2.2) & (rSelected.fe_h < -2.1)]
    sFeh21 =  rSelected[(rSelected.fe_h >= -2.1) & (rSelected.fe_h < -2.0)]
    sFeh20 =  rSelected[(rSelected.fe_h >= -2.0) & (rSelected.fe_h < -1.9)]
    sFeh19 =  rSelected[(rSelected.fe_h >= -1.9) & (rSelected.fe_h < -1.8)]
    sFeh18 =  rSelected[(rSelected.fe_h >= -1.8) & (rSelected.fe_h < -1.7)]
    sFeh17 =  rSelected[(rSelected.fe_h >= -1.7) & (rSelected.fe_h < -1.6)]
    sFeh16 =  rSelected[(rSelected.fe_h >= -1.6) & (rSelected.fe_h < -1.5)]
    sFeh15 =  rSelected[(rSelected.fe_h >= -1.5) & (rSelected.fe_h < -1.4)]
    sFeh14 =  rSelected[(rSelected.fe_h >= -1.4) & (rSelected.fe_h < -1.3)]
    sFeh13 =  rSelected[(rSelected.fe_h >= -1.3) & (rSelected.fe_h < -1.2)]
    sFeh12 =  rSelected[(rSelected.fe_h >= -1.2) & (rSelected.fe_h < -1.1)]
    sFeh11 =  rSelected[(rSelected.fe_h >= -1.1) & (rSelected.fe_h < -1.0)]

    rFeH25 =  rSelected[(rSelected.fe_h >= -2.6) & (rSelected.fe_h < -2.0)]
    rFeH20 =  rSelected[(rSelected.fe_h >= -2.0) & (rSelected.fe_h < -1.5)]
    rFeH15 =  rSelected[(rSelected.fe_h >= -1.5) & (rSelected.fe_h < -1.0)]
    rFeH10 =  rSelected[(rSelected.fe_h >= -1.0) & (rSelected.fe_h < -0.5)]
    rFeH05 =  rSelected[(rSelected.fe_h >= -0.5) & (rSelected.fe_h < -0.0)]
    rFeH00 =  rSelected[(rSelected.fe_h >= -0.0) & (rSelected.fe_h < 0.5)]
    rFeH99 =  rSelected[(rSelected.fe_h   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH99.glon, rFeH99.glat, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.glon, rFeH25.glat, color='r', label='FeH [-2.6-2.0]', marker='o', s=80)
    plt.scatter(rFeH20.glon, rFeH20.glat, color='orange', label='FeH [-2.0-1.5]', marker='o', s=80)
    plt.scatter(rFeH15.glon, rFeH15.glat, color='y', label='FeH [-1.5-1.0]', marker='o', s=80)
    plt.scatter(rFeH10.glon, rFeH10.glat, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    plt.scatter(rFeH05.glon, rFeH05.glat, color='m', label='FeH [-0.5-0.0]', marker='*', s=40)
    plt.scatter(rFeH00.glon, rFeH00.glat, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15), plt.ylabel('Glat', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    print('\nFeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))

    print('\nFeH [-2.6-2.5]: ' + str(sFeh26.shape[0]))
    print('FeH [-2.5-2.4]: '   + str(sFeh25.shape[0]))
    print('FeH [-2.4-2.3]: '   + str(sFeh24.shape[0]))
    print('FeH [-2.3-2.2]: '   + str(sFeh23.shape[0]))
    print('FeH [-2.2-2.1]: '   + str(sFeh22.shape[0]))
    print('FeH [-2.1-2.0]: '   + str(sFeh21.shape[0]))
    print('FeH [-2.0-1.9]: '   + str(sFeh20.shape[0]))
    print('FeH [-1.9-1.8]: '   + str(sFeh19.shape[0]))
    print('FeH [-1.8-1.7]: '   + str(sFeh18.shape[0]))
    print('FeH [-1.7-1.6]: '   + str(sFeh17.shape[0]))
    print('FeH [-1.6-1.5]: '   + str(sFeh16.shape[0]))
    print('FeH [-1.5-1.4]: '   + str(sFeh15.shape[0]))
    print('FeH [-1.4-1.3]: '   + str(sFeh14.shape[0]))
    print('FeH [-1.3-1.2]: '   + str(sFeh13.shape[0]))
    print('FeH [-1.2-1.1]: '   + str(sFeh12.shape[0]))
    print('FeH [-1.1-1.0]: '   + str(sFeh11.shape[0]))


    rTeff3000 =  rSelected[(rSelected.teff > -9999) & (rSelected.teff < 3500)]
    rTeff3500 =  rSelected[(rSelected.teff >= 3500) & (rSelected.teff < 4000)]
    rTeff4000 =  rSelected[(rSelected.teff >= 4000) & (rSelected.teff < 4500)]
    rTeff4500 =  rSelected[(rSelected.teff >= 4500) & (rSelected.teff < 5000)]
    rTeff5000 =  rSelected[(rSelected.teff >= 5000) & (rSelected.teff < 5500)]
    rTeff5500 =  rSelected[(rSelected.teff >= 5500) & (rSelected.teff < 6000)]
    rTeff6000 =  rSelected[(rSelected.teff >= 6000) & (rSelected.teff < 6500)]
    rTeff6500 =  rSelected[(rSelected.teff >= 6500) ]

    rTeff9999 =  rSelected[(rSelected.teff   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rTeff9999.glon, rTeff9999.glat, color='k', label='Teff[-9999]', marker='x', s=20)
    plt.scatter(rTeff3500.glon, rTeff3500.glat, color='r', label='Teff[3500-4000]', marker='o', s=80)
    plt.scatter(rTeff4000.glon, rTeff4000.glat, color='orange', label='Teff[4000-4500]', marker='o', s=80)
    plt.scatter(rTeff4500.glon, rTeff4500.glat, color='y', label='Teff[4500-5000]', marker='*', s=40)
    plt.scatter(rTeff5000.glon, rTeff5000.glat, color='purple', label='Teff[5000-5500]', marker='*', s=40)
    plt.scatter(rTeff5500.glon, rTeff5500.glat, color='m', label='Teff[5500-6000]', marker='*', s=40)
    plt.scatter(rTeff6000.glon, rTeff6000.glat, color='b', label='Teff[6000-6500]', marker='*', s=40)
    plt.scatter(rTeff6500.glon, rTeff6500.glat, color='b', label='Teff[>6500]', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15)
    plt.ylabel('Glat', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)

    plt.show()

    print('\nTeff[<3500]: '     + str(rTeff3000.shape[0]))
    print('Teff[3500-4000]: '   + str(rTeff3500.shape[0]))
    print('Teff[4000-4500]: '   + str(rTeff4000.shape[0]))
    print('Teff[4500-5000]: '   + str(rTeff4500.shape[0]))
    print('Teff[5000-5500]: '   + str(rTeff5000.shape[0]))
    print('Teff[5500-6000]: '   + str(rTeff5500.shape[0]))
    print('Teff[6000-6500]: '   + str(rTeff6000.shape[0]))
    print('Teff[>6500]: '       + str(rTeff6500.shape[0]))
    print('Teff[-9999]: '       + str(rTeff9999.shape[0]))


    #rVhelio_Draco =  rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg <= -200)]

    #rVhelio_Pal1 =  rSelected[(rSelected.vhelio_avg >= -90) & (rSelected.vhelio_avg <= -60)]
    #rVhelio_Pal2 =  rVhelio_Pal1[(rVhelio_Pal1.fe_h > -1.0) & (rVhelio_Pal1.fe_h < -0.2)]

    rVhelio_Neg_500 =  rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg < -400)]
    rVhelio_Neg_400 =  rSelected[(rSelected.vhelio_avg >= -400) & (rSelected.vhelio_avg < -300)]
    rVhelio_Neg_300 =  rSelected[(rSelected.vhelio_avg >= -300) & (rSelected.vhelio_avg < -200)]
    rVhelio_Neg_200 =  rSelected[(rSelected.vhelio_avg >= -200) & (rSelected.vhelio_avg < -100)]
    rVhelio_Neg_100 =  rSelected[(rSelected.vhelio_avg >= -100) & (rSelected.vhelio_avg < -0)]
    rVhelio_Pos_000 =  rSelected[(rSelected.vhelio_avg >= 0) & (rSelected.vhelio_avg < 100)]
    rVhelio_Pos_100 =  rSelected[(rSelected.vhelio_avg >= 100) & (rSelected.vhelio_avg < 200)]
    rVhelio_Pos_200 =  rSelected[(rSelected.vhelio_avg >= 200) & (rSelected.vhelio_avg < 300)]
    rVhelio_Pos_300 =  rSelected[(rSelected.vhelio_avg >= 300) & (rSelected.vhelio_avg < 400)]
    rVhelio_Pos_400 =  rSelected[(rSelected.vhelio_avg >= 400)]
    rVhelio9999 =  rSelected[(rSelected.vhelio_avg   == -9999)]

    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    plt.scatter(rVhelio9999.glon, rVhelio9999.glat, color='k', label='vhelio [-9999]', marker='x', s=20)
    plt.scatter(rVhelio_Neg_500.glon, rVhelio_Neg_500.glat, color='r', label='vhelio [< -400]', marker='o', s=80)
    plt.scatter(rVhelio_Neg_400.glon, rVhelio_Neg_400.glat, color='orange', label='vhelio [>= -400< -300]', 
                marker='o', s=80)
    plt.scatter(rVhelio_Neg_300.glon, rVhelio_Neg_300.glat, color='y', label='vhelio [>= -300 < -200]', marker='o', s=80)
    plt.scatter(rVhelio_Neg_200.glon, rVhelio_Neg_200.glat, color='r', label='vhelio [>= -200 < -100]',  marker='*', s=60)
    plt.scatter(rVhelio_Neg_100.glon, rVhelio_Neg_100.glat, color='b', label='vhelio [>= -100 < -0]', marker='*', s=40)
    plt.scatter(rVhelio_Pos_000.glon, rVhelio_Pos_000.glat, color='b', label='vhelio [>= 0 < 100]', marker='*', s=40)
    plt.scatter(rVhelio_Pos_100.glon, rVhelio_Pos_100.glat, color='r', label='vhelio [>= 100 < 200]', marker='*', s=60)
    plt.scatter(rVhelio_Pos_200.glon, rVhelio_Pos_200.glat, color='y', label='vhelio [>= 200 < 300]', marker='o', s=80)
    plt.scatter(rVhelio_Pos_300.glon, rVhelio_Pos_300.glat, color='orange', label='vhelio [>= 300 < 400]', 
                    marker='o', s=80)
    plt.scatter(rVhelio_Pos_400.glon, rVhelio_Pos_400.glat, color='r', label='vhelio [>= 400]', marker='o', s=80)

    #plt.scatter(rVhelio_Pal1.ra, rVhelio_Pal1.dec, color='g', label='Pal1 vhelio [>= -90 <= -60]', 
                #marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15), plt.ylabel('Glat', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()


    print('\nvhelio [< -400]:         ' + str(rVhelio_Neg_500.shape[0]))
    print('vhelio [>= -400< -300]:  ' + str(rVhelio_Neg_400.shape[0]))
    print('vhelio [>= -300 < -200]: ' + str(rVhelio_Neg_300.shape[0]))
    print('vhelio [>= -200 < -100]: ' + str(rVhelio_Neg_200.shape[0]))
    print('vhelio [>= -100 < -0]:   ' + str(rVhelio_Neg_100.shape[0]))
    print('vhelio [>= 0 < 100]:     ' + str(rVhelio_Pos_000.shape[0]))
    print('vhelio [>= 100 < 200]:   ' + str(rVhelio_Pos_100.shape[0]))
    print('vhelio [>= 200 < 300]:   ' + str(rVhelio_Pos_200.shape[0]))
    print('vhelio [>= 300 < 400]:   ' + str(rVhelio_Pos_300.shape[0]))
    print('vhelio [>= 400]:         ' + str(rVhelio_Pos_400.shape[0]))
    print('vhelio [-9999]:          ' + str(rVhelio9999.shape[0]))



    
    
    
    
    
    
### ***********************************************************************************************   
### ***********************************************************************************************
### ***********************************************************************************************



def plotMetricGlon(rSelected, selected1, s2, rFilter, rLogg, rFeH, rTeff, rVhelio, rBins):
    #s2 = r2_select_variable.value
    
    if rFilter == "Vhelio":
        x = rSelected[(rSelected.vhelio_avg > rVhelio)]
    
    if rFilter == "FeH":
        x = rSelected[(rSelected.fe_h > rFeH)]
        
    if rFilter == "Teff":
        x = rSelected[(rSelected.teff > rTeff)]
        
    if rFilter == "Logg":
        x = rSelected[(rSelected.logg > rLogg)]

    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    plt.rcParams["figure.figsize"] = (25, 25)
    plt.subplot(2, 2, 1)
    #plt.scatter(rElementRatio_Neg_9999.glon, rElementRatio_Neg_9999.dec, color='k', label='_nolegend_', 
            #marker='x', s=20)
    plt.scatter(rElementRatio_Neg_06.glon, rElementRatio_Neg_06.glat, color='orange', label='_nolegend_',      
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_05.glon, rElementRatio_Neg_05.glat, color='orange', label='_nolegend_',  
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_04.glon, rElementRatio_Neg_04.glat, color='gold', label='_nolegend_',     
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_03.glon, rElementRatio_Neg_03.glat, color='gold', label='_nolegend_',  
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_02.glon, rElementRatio_Neg_02.glat, color='b', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Neg_01.glon, rElementRatio_Neg_01.glat, color='b', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_00.glon, rElementRatio_Pos_00.glat, color='r', label='_nolegend_',
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_01.glon, rElementRatio_Pos_01.glat, color='r', label='_nolegend_',     
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_02.glon, rElementRatio_Pos_02.glat, color='g', label='_nolegend_',   
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_03.glon, rElementRatio_Pos_03.glat, color='g', label='_nolegend_',    
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_04.glon, rElementRatio_Pos_04.glat, color='m', label='_nolegend_',    
            marker='*', s=40)
    plt.scatter(rElementRatio_Pos_05.glon, rElementRatio_Pos_05.glat, color='m', label='_nolegend_',
            marker='*', s=40)
    
    if rFilter == "Vhelio":
        x = rSelected[(rSelected.vhelio_avg > -9999) & (rSelected.vhelio_avg <= rVhelio)]
    if rFilter == "FeH":
        x = rSelected[(rSelected.fe_h > -9999) & (rSelected.fe_h <= rFeH)]
    if rFilter == "Teff":
        x = rSelected[(rSelected.teff > -9999)  & (rSelected.teff < rTeff)]
    if rFilter == "Logg":
        x = rSelected[(rSelected.logg > -9999)  & (rSelected.logg < rLogg)]
            
    x2 = eval('x.'+ str(s2))

    rElementRatio_Neg_06 =  x[(x2 > -9999) & (x2 < -0.5)]
    if rElementRatio_Neg_06.shape[0] == 1:
        rElementRatio_Neg_06 = rElementRatio_Neg_06.iloc[1:]
    rElementRatio_Neg_05 =  x[(x2 >= -0.5) & (x2 < -0.4)]
    if rElementRatio_Neg_05.shape[0] == 1:
        rElementRatio_Neg_05 = rElementRatio_Neg_05.iloc[1:]
    rElementRatio_Neg_04 =  x[(x2 >= -0.4) & (x2 < -0.3)]
    if rElementRatio_Neg_04.shape[0] == 1:
        rElementRatio_Neg_04 = rElementRatio_Neg_04.iloc[1:]
    rElementRatio_Neg_03 =  x[(x2 >= -0.3) & (x2 < -0.2)]
    if rElementRatio_Neg_03.shape[0] == 1:
        rElementRatio_Neg_03 = rElementRatio_Neg_03.iloc[1:]
    rElementRatio_Neg_02 =  x[(x2 >= -0.2) & (x2 < -0.1)]
    if rElementRatio_Neg_02.shape[0] == 1:
        rElementRatio_Neg_02 = rElementRatio_Neg_02.iloc[1:]
    rElementRatio_Neg_01 =  x[(x2 >= -0.1) & (x2 < -0.0)]
    if rElementRatio_Neg_01.shape[0] == 1:
        rElementRatio_Neg_01 = rElementRatio_Neg_01.iloc[1:]
    rElementRatio_Pos_00 =  x[(x2 >=  0.0) & (x2 <  0.1)]
    if rElementRatio_Pos_00.shape[0] == 1:
        rElementRatio_Pos_00 = rElementRatio_Pos_00.iloc[1:]
    rElementRatio_Pos_01 =  x[(x2 >=  0.1) & (x2 <  0.2)]
    if rElementRatio_Pos_01.shape[0] == 1:
        rElementRatio_Pos_01 = rElementRatio_Pos_01.iloc[1:]
    rElementRatio_Pos_02 =  x[(x2 >=  0.2) & (x2 <  0.3)]
    if rElementRatio_Pos_02.shape[0] == 1:
        rElementRatio_Pos_02 = rElementRatio_Pos_02.iloc[1:]
    rElementRatio_Pos_03 =  x[(x2 >=  0.3) & (x2 <  0.4)]
    if rElementRatio_Pos_03.shape[0] == 1:
        rElementRatio_Pos_03 = rElementRatio_Pos_03.iloc[1:]
    rElementRatio_Pos_04 =  x[(x2 >=  0.4) & (x2 <  0.5)]
    if rElementRatio_Pos_04.shape[0] == 1:
        rElementRatio_Pos_04 = rElementRatio_Pos_04.iloc[1:]
    rElementRatio_Pos_05 =  x[(x2 >=  0.5)]
    if rElementRatio_Pos_05.shape[0] == 1:
        rElementRatio_Pos_05 = rElementRatio_Pos_05.iloc[1:]
    rElementRatio_Neg_9999 =  x[(x2   == -9999)]
    
    #plt.scatter(rElementRatio_Neg_9999.ra, rElementRatio_Neg_9999.dec, color='k', label='[Ni/Fe = -9999]', 
            #marker='x', s=20)
    plt.scatter(rElementRatio_Neg_06.glon, rElementRatio_Neg_06.glat, color='orange', label='_nolegend_', 
                marker='o', s=80)
    plt.scatter(rElementRatio_Neg_05.glon, rElementRatio_Neg_05.glat, color='orange', label = s2 + ' [>= -0.5 < -0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_04.glon, rElementRatio_Neg_04.glat, color='gold', label= s2 + '[>= -0.4 < -0.2]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_03.glon, rElementRatio_Neg_03.glat, color='gold', label='_nolegend_', 
                marker='o', s=80)
    plt.scatter(rElementRatio_Neg_02.glon, rElementRatio_Neg_02.glat, color='b', label= s2 + '[>= -0.2 < -0.0]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Neg_01.glon, rElementRatio_Neg_01.glat, color='b', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_00.glon, rElementRatio_Pos_00.glat, color='r', label= s2 + '[>= 0.0 < 0.2]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_01.glon, rElementRatio_Pos_01.glat, color='r', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_02.glon, rElementRatio_Pos_02.glat, color='g', label= s2 + '[>= 0.2 < 0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_03.glon, rElementRatio_Pos_03.glat, color='g', label='_nolegend_', marker='o', s=80)
    plt.scatter(rElementRatio_Pos_04.glon, rElementRatio_Pos_04.glat, color='m', label= s2 + '[>= 0.4]', 
            marker='o', s=80)
    plt.scatter(rElementRatio_Pos_05.glon, rElementRatio_Pos_05.glat, color='m', label='_nolegend_', marker='o', s=80)
    ax = plt.gca()
    plt.xlabel('Glon', fontsize=15), plt.ylabel('Glat', fontsize=15)
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=15)
    plt.show()

    print('\n' + s2 + ' [< -0.5]: ' + str(rElementRatio_Neg_06.shape[0]))
    print(s2 + ' [>= -0.5) & < -0.4]: '   + str(rElementRatio_Neg_05.shape[0]))
    print(colored(s2 +  ' [>= -0.4) & < -0.3]: '   + str(rElementRatio_Neg_04.shape[0]), 'yellow'))
    print(colored(s2 +  ' [>= -0.3) & < -0.2]: '   + str(rElementRatio_Neg_03.shape[0]), 'yellow'))
    print(colored(s2 +  ' [>= -0.2) & < -0.1]: '   + str(rElementRatio_Neg_02.shape[0]), 'blue'))
    print(colored(s2 +  ' [>= -0.1) & < -0.0]: '   + str(rElementRatio_Neg_01.shape[0]), 'blue'))
    print(colored(s2 +  ' [>=  0.0) & <  0.1]: '   + str(rElementRatio_Pos_00.shape[0]), 'red'))
    print(colored(s2 +  ' [>=  0.1) & <  0.2]: '   + str(rElementRatio_Pos_01.shape[0]), 'red'))
    print(colored(s2 +  ' [>=  0.2) & <  0.3]: '   + str(rElementRatio_Pos_02.shape[0]), 'green'))
    print(colored(s2 +  ' [>=  0.3) & <  0.4]: '   + str(rElementRatio_Pos_03.shape[0]), 'green'))
    print(colored(s2 +  ' [>=  0.4) & <  0.5]: '   + str(rElementRatio_Pos_04.shape[0]), 'magenta'))
    print(colored(s2 +  ' [>=  0.5]: '   + str(rElementRatio_Pos_05.shape[0]), 'magenta'))
    print(s2 +  ' [-9999]: '   + str(rElementRatio_Neg_9999.shape[0]))

    xElementRatio_Neg_04 = eval('rElementRatio_Neg_04.'+ str(s2))
    xElementRatio_Neg_03 = eval('rElementRatio_Neg_03.'+ str(s2))
    xElementRatio_Neg_02 = eval('rElementRatio_Neg_02.'+ str(s2))
    xElementRatio_Neg_01 = eval('rElementRatio_Neg_01.'+ str(s2))
    xElementRatio_Pos_00 = eval('rElementRatio_Pos_00.'+ str(s2))
    xElementRatio_Pos_01 = eval('rElementRatio_Pos_01.'+ str(s2))
    xElementRatio_Pos_02 = eval('rElementRatio_Pos_02.'+ str(s2))
    xElementRatio_Pos_03 = eval('rElementRatio_Pos_03.'+ str(s2))
    xElementRatio_Pos_04 = eval('rElementRatio_Pos_04.'+ str(s2))
    xElementRatio_Pos_05 = eval('rElementRatio_Pos_05.'+ str(s2))
    
    plt.rcParams["figure.figsize"] = (25, 10)
    plt.subplot(2, 3, 1)
    plt.hist([xElementRatio_Neg_04], range=(-.5, 1.0), stacked=True, color='y', bins=rBins,
         label=s2 +  ' [>= -0.4) & < -0.2]')
    plt.hist([xElementRatio_Neg_03], range=(-.5, 1.0), stacked=True, color='y', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Neg_02], range=(-.5, 1.0), stacked=True, color='b', bins=rBins,
         label=s2 +  ' [>= -0.2) & < -0.0]')
    plt.hist([xElementRatio_Neg_01], range=(-.5, 1.0), stacked=True, color='b', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_00], range=(-.5, 1.0), stacked=True, color='r', bins=rBins,
         label=s2 +  ' [>=  0.0) & <  0.2]')
    plt.hist([xElementRatio_Pos_01], range=(-.5, 1.0), stacked=True, color='r', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_02], range=(-.5, 1.0), stacked=True, color='g', bins=rBins,
         label=s2 +  ' [>=  0.2) & <  0.4]')
    plt.hist([xElementRatio_Pos_03], range=(-.5, 1.0), stacked=True, color='g', bins=rBins, label='_nolegend_')
    plt.hist([xElementRatio_Pos_04], range=(-.5, 1.0), stacked=True, color='m', bins=rBins, label=s2 + ' [>=  0.4]')
    plt.hist([xElementRatio_Pos_05], range=(-.5, 1.0), stacked=True, color='m', bins=rBins, label='_nolegend_')
    plt.xlabel(s2, fontsize=15)
    ax = plt.gca()
    ax.yaxis.set_label_position("right")
    plt.grid()
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    rFeH25 =  rSelected[(rSelected.fe_h >= -2.6) & (rSelected.fe_h < -2.0)]
    rFeH20 =  rSelected[(rSelected.fe_h >= -2.0) & (rSelected.fe_h < -1.5)]
    rFeH15 =  rSelected[(rSelected.fe_h >= -1.5) & (rSelected.fe_h < -1.0)]
    rFeH10 =  rSelected[(rSelected.fe_h >= -1.0) & (rSelected.fe_h < -0.5)]
    rFeH05 =  rSelected[(rSelected.fe_h >= -0.5) & (rSelected.fe_h < -0.0)]
    rFeH00 =  rSelected[(rSelected.fe_h >= -0.0) & (rSelected.fe_h < 0.5)]
    rFeH99 =  rSelected[(rSelected.fe_h   == -9999)]
    
    xFeH25 = eval('rFeH25.'+ str(s2))
    xFeH20 = eval('rFeH20.'+ str(s2))
    xFeH15 = eval('rFeH15.'+ str(s2))
    xFeH10 = eval('rFeH10.'+ str(s2))
    xFeH05 = eval('rFeH05.'+ str(s2))
    xFeH00 = eval('rFeH00.'+ str(s2))
    xFeH99 = eval('rFeH99.'+ str(s2))

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 6500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.logg, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    #plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=40)
    #plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=40)
    #plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=40)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    print('\nFeH [-2.6-2.0]: '   + str(rFeH25.shape[0]))
    print('FeH [-2.0-1.5]: '   + str(rFeH20.shape[0]))
    print('FeH [-1.5-1.0]: '   + str(rFeH15.shape[0]))
    print('FeH [-1.0-0.5]: '   + str(rFeH10.shape[0]))
    print('FeH [-0.5-0.0]: '   + str(rFeH05.shape[0]))
    print('FeH [0.0-0.5]:  '   + str(rFeH00.shape[0]))
    print('FeH [-9999]:    '   + str(rFeH99.shape[0]))
    
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(3500, 6500)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    
    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.logg, xFeH99, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.logg, xFeH25, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.logg, xFeH20, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.logg, xFeH15, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.logg, xFeH10, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.logg, xFeH05, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.logg, xFeH00, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('LogG', fontsize=15), plt.ylabel(str(s2), fontsize=15)
    plt.xlim(0, 4)
    plt.ylim(-1, 1)
    #plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()

    plt.rcParams["figure.figsize"] = (15, 8)
    plt.subplot(2, 2, 1)
    #plt.scatter(rFeH99.teff, rFeH99.logg, color='k', label='FeH [-9999]', marker='x', s=20)
    #plt.scatter(rFeH25.teff, rFeH25.logg, color='r', label='FeH [-2.6-2.0]', marker='o', s=20)
    #plt.scatter(rFeH20.teff, rFeH20.logg, color='orange', label='FeH [-2.0-1.5]', marker='o', s=20)
    #plt.scatter(rFeH15.teff, rFeH15.logg, color='y', label='FeH [-1.5-1.0]', marker='o', s=20)
    plt.scatter(rFeH10.teff, rFeH10.logg, color='purple', label='FeH [-1.0-0.5]', marker='*', s=20)
    plt.scatter(rFeH05.teff, rFeH05.logg, color='g', label='FeH [-0.5-0.0]', marker='*', s=20)
    plt.scatter(rFeH00.teff, rFeH00.logg, color='b', label='FeH [0.0-0.5]', marker='*', s=20)
    ax = plt.gca()
    plt.xlabel('Teff', fontsize=15), plt.ylabel('Log (G)', fontsize=15)
    plt.xlim(4000, 5500)
    plt.ylim(3.5, 0.5)
    plt.gca().invert_xaxis()
    plt.grid()
    plt.title(selected1, fontsize=14)
    L=plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.show()
    



    
