import numpy as np
import math
import random as rnd
from scipy.linalg import null_space
from multiprocessing import Pool
from itertools import combinations
import itertools
from numpy import nan


class Simulation(object):

  # np.random.seed(0)
  def __init__(self, nD2DPairs, nCellUser=5, PowerCell_dB =10, PowerD2D_dB = 5, nEavesdroppers = 2):
    self.nSimBlocks = 10**3         # Number of simulated blocks
    self.nCellUsers = nCellUser     # Number of Cellular Users
    self.nD2DPairs = nD2DPairs      # Number of D2D Pairs
    self.nEavesdroppers = nEavesdroppers # Number of Eavesdroppers
    """
    Transmission Powers (Fixed Power Transmission)
    """     
    self.PowerCell_dB = PowerCell_dB
    self.PowerCell = 10**(self.PowerCell_dB/10)
    self.PowerD2D_dB = PowerD2D_dB
    self.PowerD2D = 10**(self.PowerD2D_dB/10)
    """
    Channel coefficients 
    """
    # Secret Transmission Links 
    self.H_c = (np.random.randn(self.nCellUsers,self.nSimBlocks)+1j*np.random.randn(self.nCellUsers,self.nSimBlocks))/math.sqrt(2)
    self.Gamma_c = abs(self.H_c)**2
    # Eavesdropping Links
    self.G_c = (np.random.randn(self.nEavesdroppers,self.nSimBlocks)+1j*np.random.randn(self.nEavesdroppers,self.nSimBlocks))/math.sqrt(2)
    self.Gamma_Ec = abs(self.G_c)**2

    # Secret Transmission Links D2D
    self.H_d = (np.random.randn(self.nD2DPairs,self.nSimBlocks)+1j*np.random.randn(self.nD2DPairs,self.nSimBlocks))/math.sqrt(2);
    self.Gamma_d = abs(self.H_d)**2;
    # Eavesdropping Links
    self.G_d = (np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nEavesdroppers)+1j*np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nEavesdroppers))/math.sqrt(2);
    self.Gamma_Ed = abs(self.G_d)**2;

    
    # Interference Links
    self.L_c  = (np.random.randn(self.nD2DPairs,self.nSimBlocks)+1j*np.random.randn(self.nD2DPairs,self.nSimBlocks))/math.sqrt(2/0.5);
    self.L_d  = (np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nCellUsers)+1j*np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nCellUsers))/math.sqrt(2/0.5)
    self.L_dd = (np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nD2DPairs)+1j*np.random.randn(self.nSimBlocks,self.nD2DPairs,self.nD2DPairs))/math.sqrt(2/0.5) 
    self.L_dd = ((1-np.eye(self.nD2DPairs)) * self.L_dd)  # diagonal elements=0(selfinterference=0)  vs = 1; ve = 1; vi = 0.5
    
    self.Gamma_Ic  = abs(self.L_c)**2;
    self.Gamma_Id  = abs(self.L_d)**2;
    self.Gamma_Idd = abs(self.L_dd)**2;



  def RcSumNoD2D(self, beta=2):
    """
    The Primary System's Secrecy Condition
      -arg: beta
      -returns: Rth
    """
    p = Pool(2)
    self.RC_sum_noD2D = 0;
    gamma_hc_max = p.map(max, (self.Gamma_c[:,l] for l in range(self.nSimBlocks)))
    gamma_gc_max = p.map(max, (self.Gamma_Ec[:,l] for l in range(self.nSimBlocks)))
    self.RC_sum_noD2D = np.sum((np.log2(1+np.array(gamma_hc_max)*self.PowerCell)-np.log2(1+np.array(gamma_gc_max)*self.PowerCell)))/self.nSimBlocks
    self.RC_sum_noD2D = max(0,self.RC_sum_noD2D)
    RC_sum_th =  beta*self.RC_sum_noD2D
    p.close()
    p.join()
    return RC_sum_th, self.RC_sum_noD2D

  def RcSum(self, Set_SD, Set_AN):
    """
    The Primary Secrecy Rate
    """

    self.RC_sum = 0;
    if not (Set_SD.size == 0 & Set_AN.size ==0):
      for l in range(self.nSimBlocks):
        SINR_hc_max = np.amax(self.Gamma_c[:,l][:,np.newaxis].T/(1+ np.sum(self.Gamma_Id[l,Set_SD,:], axis=0)[np.newaxis,:]*self.PowerD2D))
  
        index_hc = np.argmax(self.Gamma_c[:,l][:,np.newaxis].T/(1+ np.sum(self.Gamma_Id[l,Set_SD,:], axis=0)[np.newaxis,:]*self.PowerD2D))

        Z = null_space(np.vstack((np.squeeze(self.L_dd[l,Set_AN,Set_SD]).T, np.squeeze(self.L_d[l,Set_AN,index_hc]).T)))
        
        SINR_gc_max = np.amax(self.Gamma_Ec[:,l][:,np.newaxis].T/\
                              (1+ np.sum(self.Gamma_Ed[l,Set_SD,:], axis=0)[np.newaxis,:]*self.PowerD2D +\
                              np.diag(np.real(np.squeeze(self.G_d[l,Set_AN,:]).conj().T@(Z@Z.conj().T)@np.squeeze(self.G_d[l,Set_AN,:]))).T*self.PowerD2D))
        
        self.RC_sum = self.RC_sum + (np.log2(1.0+SINR_hc_max*self.PowerCell)-np.log2(1.0+SINR_gc_max*self.PowerCell))/self.nSimBlocks
    return max(0,self.RC_sum)      

  def RDk(self, Set_SD, Set_AN):
    """
    The D2D Secrecy Rate
      -arg: Set_SD, Set_AN
      -returns: RD_k
    """    
    self.RD = np.zeros((self.nD2DPairs,1))
    self.RD[Set_AN] = np.nan

    if not (Set_SD.size == 0 & Set_AN.size ==0):
      for l in range(self.nSimBlocks):

        Set_SD_R = Set_SD
        
        for k in list(Set_SD):

          Set_SD_R = list(set(Set_SD_R)-set([k]))

          index_hc = np.argmax(self.Gamma_c[:,l][:,np.newaxis].T/(1+ np.sum(self.Gamma_Id[l,Set_SD,:], axis=0)[np.newaxis,:]*self.PowerD2D))

          Z = null_space(np.vstack((np.squeeze(self.L_dd[l,Set_AN,Set_SD]).T, np.squeeze(self.L_d[l,Set_AN,index_hc]).T)))

          SINR_gd_max =  np.amax(self.Gamma_Ed[l,k,:]/\
                              (1+ np.sum(self.Gamma_Ed[l,Set_SD_R,:], axis=0)[np.newaxis,:]*self.PowerD2D +\
                              np.diag(np.real(np.squeeze(self.G_d[l,Set_AN,:]).conj().T@(Z@Z.conj().T)@np.squeeze(self.G_d[l,Set_AN,:]))).T*self.PowerD2D))
          
          SINR_hd = self.Gamma_d[k,l]/(1+self.Gamma_Ic[k,l]*self.PowerCell+np.sum(self.Gamma_Idd[l,Set_SD,k])*self.PowerD2D)

          # import pdb ; pdb.set_trace()

          # Gamma_d(k,l)/(1+Gamma_Ic(k,l)*Pc+sum(Gamma_Idd(Set_SD,k,l))*Pd)
          self.RD[k] = (self.RD[k]+(np.log2(1+SINR_hd*self.PowerD2D)-np.log2(1+SINR_gd_max*self.PowerD2D))/self.nSimBlocks)

    self.RD = (np.array(self.RD)).flatten()
    self.RD = [np.amax([0,rd]) for rd in self.RD]
    return np.array(self.RD)

  def getCombination(self, set_A, Exp=3):
    """
    Find all combination of Ks and Kj
    """
    K_A = len(set_A)
    if Exp==1:
        self.C = np.ones((1,self.nD2DPairs))*np.nan
        nC = 0;
        for nKs in range (1, (math.floor((self.nD2DPairs-1)/2)+1)):
            set_S = list(combinations(range(0,self.nD2DPairs),nKs))
            for ni in range(len(set_S)):
                self.C[nC,set_S[ni]]= set_S[ni]
                nC = nC + 1
                self.C = np.vstack((self.C, np.nan*np.ones((1,self.nD2DPairs))))
        return self.C  
    if Exp == 2:
        if K_A <= 2:
            return np.ones((1,self.nD2DPairs))*1j
        elif K_A == self.nD2DPairs:
            self.C = np.ones((1,self.nD2DPairs))*np.nan
            nC = 0
            for nKs in range (1, (math.floor((self.nD2DPairs-1)/2)+1)):
                set_S = list(combinations(range(0,self.nD2DPairs),nKs))
                for ni in range(len(set_S)):
                    self.C[nC,set_S[ni]]= set_S[ni]
                    nC = nC + 1
                    self.C = np.vstack((self.C, np.nan*np.ones((1,self.nD2DPairs))))
            self.C = self.C[:-1]
            return self.C 
        else:
            self.C = np.ones((1,self.nD2DPairs))*1j
            nC = 0
            for nKj in range(2, K_A+1):
                combJDevice = list(combinations(set_A,nKj))
                for nKs in range(1, min(nKj,self.nD2DPairs-nKj+1)):
                    for nj in range(len(combJDevice)):
                        set_S = list(set(range(0,self.nD2DPairs))-set(combJDevice[nj]))
                        combSDevice = list(combinations(set_S,nKs))
                        for ni in range (len(combSDevice)):
                            self.C[nC,combJDevice[nj]] = np.nan
                            self.C[nC,combSDevice[ni]]= combSDevice[ni]
                            nC = nC + 1
                            self.C = np.vstack((self.C, 1j*np.ones((1,self.nD2DPairs))))
            self.C = self.C[:-1]
            return self.C 

    if Exp == 3:
      if K_A <= 2:
        return np.ones((1,self.nD2DPairs))*1j
      else:
        self.C = []
        for seq in itertools.product([0,1,2], repeat=K_A):
            c_s = list(seq)
            if c_s.count(2) >= 2 and c_s.count(2) <= K_A-1 and c_s.count(1) >= 1 and c_s.count(1) <= c_s.count(2)-1:
                np_cs = np.array(c_s, dtype=complex) 
                np_cs[np_cs == 2] = np.nan
                np_cs[np_cs == 0] = 1j
                ks_idx = np.where(~np.isnan(np_cs) & np.isreal(np_cs))
                np_cs[ks_idx] = ks_idx
                self.C.append(np_cs)
        return np.array([element for (i,element) in enumerate(self.C)]) 
        
  def getEmin(self, nSubSlots):
    self.Emin = self.PowerD2D * (1/nSubSlots)
    return self.Emin

  def calAvailableEnergy(self,deviceEnergy,comb):
    '''
    Calculate Available Energy 
      Ek(t) = max(0, Ek(t-1) - Emin)
    '''
    txDevice = np.concatenate(np.where(np.isreal(comb) | np.isnan(comb)),axis=0)
    deviceAvailableEnergy = np.array(deviceEnergy)
    for idx_tx in txDevice:
      deviceAvailableEnergy[idx_tx] = max(0,deviceAvailableEnergy[idx_tx]- self.Emin)
    return deviceAvailableEnergy

  def getActiveDevices(self, deviceEnergy):
    '''
    Get Active Devices the satisfies Ek > Emin
    '''
    setA = np.argwhere(np.array(deviceEnergy) > self.Emin).flatten()
    return setA
