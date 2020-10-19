from sim import *
from addNode import *
from multiprocessing import Pool
import pandas as pd
from datetime import datetime
import operator
import csv
import argparse
import time

def Args():
    parser = argparse.ArgumentParser(description="Main handler for simulation")
    parser.add_argument("-start_d", "--start_device", type=int, help='starting D2D device range',required=True)
    parser.add_argument("-end_d", "--end_device", type=int, help='ending D2D device range (PS: end device excluded)')
    parser.add_argument("-ncellusr", "--n_cell_users", type=int, default=5, help='number of cellular user (default: 5)')
    parser.add_argument('-pcell',"--pwr_cell_db", type=int, default=10, help='Power of cellular user (default: 10)')
    parser.add_argument('-pd2d',"--pwr_d2d_db", type=int, default=5, help='Power of D2D user (default: 5)')
    parser.add_argument("-neve", "--n_eavesdropper", type=int, default=2, help=' number of eavesdroppers (default: 2)')
    parser.add_argument('-ncpu', "--ncpu", type=int, default=2, help='number of CPU (default: 2)', required=True)
    parser.add_argument('-nsubslot', "--nsubslot", type=int, default=1, help='number of sub-slots(default: 2)',required=True)
    parser.add_argument('-seed', "--seed", type=int, default=10, help='range of numpy random seed(default: 10)',required=True)
    parser.add_argument('-rdth', "--rdthreshold", type=int, help='Percentage of RD Max threshold')
    parser.add_argument('-exp', "--experiment", type=int, default=3, help='Experiment-1,2 or 3 (default:3)',required=True)

    args = parser.parse_args() 
    return args

def main(Combi):
    RD_sum = 0
    RC_sum = 0
    Set_SD = np.concatenate(np.where(np.isreal(Combi) & ~np.isnan(Combi)), axis=0)
    Set_AN = np.concatenate(np.where(np.isnan(Combi)), axis=0)[:,np.newaxis]
    RC_sum = Sim_class.RcSum(Set_SD, Set_AN)
    if RC_sum >= R_th:
        RD = Sim_class.RDk(Set_SD, Set_AN)
        RD_sum = np.nansum(RD)
    return RD_sum, RC_sum

def preprocessing(Combi_arr):
    if (np.sum(np.isnan(Combi_arr))+np.sum((~np.isnan(Combi_arr) & np.isreal(Combi_arr))) > 5):
        return Combi_arr

if __name__ == '__main__':
    '''
        Inilization
    '''
    args = Args()
    nCPU = args.ncpu
    nSubSlots = args.nsubslot
    rseeds = range(args.seed)

    if args.end_device is not None:
            rD2D = range(args.start_device, args.end_device)
    else:
            rD2D = [args.start_device]
    KS = np.zeros((len(rseeds),len(rD2D),nSubSlots))
    KJ = np.zeros((len(rseeds),len(rD2D),nSubSlots))
    R_throughput = np.zeros((len(rseeds),len(rD2D)))
    RC_sm = np.zeros((len(rseeds),len(rD2D)))
    time_taken = np.zeros((len(rseeds),len(rD2D),nSubSlots))

    for seed in rseeds:
        np.random.seed(seed)
        '''
        Energy harvest 
        '''
        Ek_h = 20.*np.random.rand(15)

        for nD2D_idx, nD2D in enumerate(rD2D):
            '''
            Inilize class 
            '''
            Sim_class = Simulation(nD2DPairs=nD2D, nCellUser=args.n_cell_users, PowerCell_dB =args.pwr_cell_db, PowerD2D_dB = args.pwr_d2d_db, nEavesdroppers = args.n_eavesdropper)
            Emin = Sim_class.getEmin(nSubSlots)
            R_th, Rcsum_noD2D = Sim_class.RcSumNoD2D()
            Ek = Ek_h[:nD2D]
            node_dict = []
            [node_dict.append({}) for slots in range(nSubSlots+1)]
            '''
            Inilize Root node
            '''
            node_dict[0]['rootNode'] = Add_node(name='root', RD_sum =0, RC_sum=Rcsum_noD2D, Rcsm_noD2D= Rcsum_noD2D, R_th=R_th, RD_acc = 0, Combi=range(nD2D), dEnergy = Ek)
            '''
            Build tree for each sub-slot
            '''
            for idx_dict, ndict in enumerate(node_dict[:-1]):
                start_time = time.time()
                print('Device: ', nD2D, ' | Sub-slot: ', (idx_dict+1),' | Seed: ', seed,' | Time: ', datetime.now().time().strftime("%I:%M:%S %p") )
                '''
                Combinations for each parent node
                '''
                for key in ndict:
                    set_A = Sim_class.getActiveDevices(node_dict[idx_dict][key].dEnergy)
                    '''
                    get all combination
                    '''
                    Combi = Sim_class.getCombination(set_A, args.experiment)
                    p = Pool(nCPU)
                    Combi = p.map(preprocessing, Combi)
                    p.close()
                    p.join()
                    Combi = np.array(list(filter(None.__ne__, Combi)))
                    '''
                    Divide the combinations into chunks for fast processing
                    '''
                    srt_lst = [np.count_nonzero(np.isnan(row)) for row in Combi]
                    srt_idx = np.argsort(-1*np.array(srt_lst))[:len(srt_lst)]
                    Combi = Combi[srt_idx]

                    combi_chunk = np.array_split(Combi,2)
                    for idx_combi_arr, combi_arr in enumerate(combi_chunk[:-1]):
                        p = Pool(nCPU)
                        '''
                        Check if Rcsum > R_th
                        '''
                        Rtrn_list = p.map(main, combi_arr)
                        p.close()
                        p.join()
                        Rtrn_np = np.array(Rtrn_list)
                        if not Rtrn_np.any():
                            RD_np = 0
                            RC_np = Rcsum_noD2D
                        else:
                            RD_np = Rtrn_np[:,0]
                            RC_np = Rtrn_np[:,1]

                        '''
                        Choose between optimal and sub-optimal solution
                        '''
                        if args.rdthreshold is not None:
                            max_rd=np.amax(RD_np)
                            min_rd_th = max_rd - (max_rd * (args.rdthreshold/100))
                            RD_Cidx  = np.concatenate(np.where(RD_np > min_rd_th),axis=0)
                            RD_del_Cidx = np.concatenate(np.where(RD_np == 0),axis=0)
                        else:
                            RD_Cidx = np.concatenate(np.where(RD_np > 0),axis=0)
                            RD_del_Cidx = np.concatenate(np.where(RD_np == 0),axis=0)
                        '''
                        Update Ek and save all the values of the combinations
                        '''
                        for r_i in RD_Cidx:
                            Ek = Sim_class.calAvailableEnergy(node_dict[idx_dict][key].dEnergy, combi_arr[r_i])
                            name_combi = str(node_dict[idx_dict][key].name)+'C'+str(r_i)
                            node_dict[idx_dict+1][name_combi] = Add_node(name=name_combi, Combi = combi_arr[r_i],Rcsm_noD2D= Rcsum_noD2D, R_th=R_th, RC_sum = RC_np[r_i], RD_sum = RD_np[r_i], RD_acc=(node_dict[idx_dict][key].RD_sum + RD_np[r_i]), dEnergy = Ek, parent= node_dict[idx_dict][key])

                        del_arr =[]
                        [del_arr.append(combi_arr[d_rth_i]) for d_rth_i  in RD_del_Cidx]
                        for _del_row in range(len(del_arr)):
                            del_row_member = np.concatenate(np.where(~np.isnan(del_arr[_del_row]) & np.isreal(del_arr[_del_row])),axis=0)
                            for ch_arr in combi_chunk[idx_combi_arr+1:]:
                                memb = [np.isin(del_row_member, ch_arr[l]) for l in range(ch_arr.shape[0])]
                                del_idx = [i for i, mem in enumerate(memb) if np.all(mem)]
                                combi_chunk[idx_combi_arr+1] = np.delete(ch_arr, del_idx, axis=0)

                end_time = time.time()
                print('--------------------------------------------------------------')        
                print('Time taken: ', end_time - start_time)
                print('--------------------------------------------------------------')  
                time_taken[seed,nD2D_idx,idx_dict] =  end_time - start_time

            '''
            Find max RD throughput and return KS and KJ.
            '''
            depth_tree = np.concatenate(np.where(node_dict[1:]),axis=0).size
            if depth_tree:
                leaf_list = list(node_dict[0]['rootNode'].leaves)
                leaf_dict =  {}
                for lf_idx, lf_val in enumerate(leaf_list):
                    leaf_dict.update({leaf_list[lf_idx].name : leaf_list[lf_idx].RD_acc})
                max_R = max(leaf_dict.items(), key=operator.itemgetter(1))[0]
                req_combi = list(node_dict[depth_tree][max_R].path)
                Final_comb = np.array([c.Combi for c in req_combi[1:]])
                for dpt_i in range(depth_tree):
                    KS[seed,nD2D_idx,dpt_i] = np.sum((~np.isnan(Final_comb[dpt_i]) & np.isreal(Final_comb[dpt_i])))
                    KJ[seed,nD2D_idx,dpt_i] = np.sum(np.isnan(Final_comb[dpt_i]))
                R_throughput[seed,nD2D_idx] = node_dict[depth_tree][max_R].RD_acc
                RC_sm[seed,nD2D_idx] = node_dict[depth_tree][max_R].RC_sum 

            '''
            Save to csv file
            '''
            # name_csv = 'csv_K_'+str(nD2D)+'_rdth_'+str(args.rdthreshold)+'_num_of_ss_'+str(args.nsubslot)+'_seed_'+str(seed)+'.csv'
            # with open(name_csv, 'w') as file:
            #     fieldnames = ['K','KS', 'KJ', 'R_throughput', 'RC_sm','Rcsm_noD2D', 'R_th', 'Time taken']
            #     writer = csv.DictWriter(file, fieldnames=fieldnames)
            #     writer.writeheader()
            #     writer.writerow({'K':rD2D[nD2D_idx], 'KS': KS[seed,nD2D_idx,:], 'KJ' : KJ[seed,nD2D_idx,:], 'R_throughput': R_throughput[seed,nD2D_idx],'RC_sm': RC_sm[seed,nD2D_idx], 'Rcsm_noD2D': Rcsum_noD2D, 'R_th':R_th,'Time taken': time_taken[seed,nD2D_idx,:]})
            # name_sim_pk = 'pkl_K_'+str(nD2D)+'_rdth_'+str(args.rdthreshold)+'_num_of_ss_'+str(args.nsubslot)+'_seed_'+str(seed)+'.pickle'
            # with open(name_sim_pk, 'wb') as handle:
            #     pkl.dump(node_dict, handle)
        
        name_sim = 'Acsv_K_'+str(args.start_device)+'_'+str(args.end_device)+'_rdth_'+str(args.rdthreshold)+'_num_of_ss_'+str(args.nsubslot)+'_seed_'+str(seed)+'.csv'
        with open(name_sim, 'w') as file:
            fieldnames = ['K','KS', 'KJ', 'R_throughput', 'RC_sm','Rcsm_noD2D', 'R_th','Time taken']
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in range(len(rD2D)):
                writer.writerow({'K':rD2D[row], 'KS': KS[seed,row,:], 'KJ' : KJ[seed,row,:], 'R_throughput': R_throughput[seed,row],'RC_sm': RC_sm[seed,row],'Rcsm_noD2D': Rcsum_noD2D, 'R_th':R_th, 'Time taken': time_taken[seed,row,:]})

    time_taken_E_seed = np.mean(time_taken, axis=0)
    KS_E_seed = np.round(np.mean(KS, axis=0))
    KJ_E_seed = np.round(np.mean(KJ, axis=0))
    R_thput_E_seed = np.mean(R_throughput, axis=0)
    RC_sm_E_seed = np.mean(RC_sm, axis=0)
    name_sim = 'Ecsv_K_'+str(args.start_device)+'_'+str(args.end_device)+'_rdth_'+str(args.rdthreshold)+'_num_of_ss_'+str(args.nsubslot)+'_seed_avg_.csv'
    with open(name_sim, 'w', newline='') as file:
            fieldnames = ['K','KS_E', 'KJ_E', 'R_throughput_E', 'RC_sm_E', 'Time taken']
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in range(len(rD2D)):
                writer.writerow({'K':rD2D[row], 'KS_E': KS_E_seed[row], 'KJ_E' : KJ_E_seed[row], 'R_throughput_E': R_thput_E_seed[row],'RC_sm_E': RC_sm_E_seed[row], 'Time taken': time_taken_E_seed[row]})
