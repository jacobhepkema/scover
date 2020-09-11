#!/usr/bin/env python

from __future__ import print_function

# Sys
import os
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import random
import argparse
import pickle

# Motif analysis
import logomaker 

# HDF5, math, plot, other
import h5py
import numpy as np
from sklearn.model_selection import KFold
from scipy.stats import spearmanr
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.Data import TensorDataset
from torch.utils.data.sampler import SubsetRandomSampler
import pandas as pd
from tqdm import trange 
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Custom functions
import scanem_utils as su
import scanem_model as sm

# argparse
parser = argparse.ArgumentParser(description='scanem pytorch version v0.2', 
formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('data', metavar='data', type=str,
                    help='path to tab-separated data file containing named columns "sequence" (containing DNA sequences) \
                        and "ind" (containing comma-separated values to be predicted)')
parser.add_argument('celldata', metavar='celldata', type=str,
                    help='path to tab-separated table containing information on cell types in dataset')
parser.add_argument('-name', metavar='EXPERIMENT_NAME', type=str, 
                    help='name of experiment that will be used as file prefix', 
                    default='exp_name')
parser.add_argument('-opt', metavar='OPTIMIZER', type=str, 
                    help='name of optimizer to use. Options: SGD, Adam, Adagrad', 
                    default='SGD')
parser.add_argument('-c', metavar='NUM_CALIBRATIONS', type=int,
                    help='amount of random parameter samplings to try', 
                    default=30)
parser.add_argument('-m', metavar='MOTIF_LEN', type=int,
                    help='Length of motif detectors', 
                    default=12)     
parser.add_argument('-d', metavar='NUMBER_OF_MOTIFS', type=int,
                    help='Number of motif detectors', 
                    default=300)               
parser.add_argument('-v', metavar='VAL_FACTOR', type=int, 
                    help='K in K-fold cross-validation', 
                    default=10)
parser.add_argument('-epochs', metavar='EPOCHS', type=int, 
                    help='number of epochs to train', 
                    default=30)
parser.add_argument('-model_output_dir', metavar='/PATH/TO/DIR', type=str, 
                    help='directory where model files will be stored. Will append experiment name to path', 
                    default='.')
parser.add_argument('-num_errtest', metavar='NUM_ERRTEST', type=int, 
                    help='every epochs / num_errtest learning steps validation error will be assessed', 
                    default=1)
parser.add_argument('-batch_size', metavar='BATCH_SIZE', type=int, 
                    help='batch size for SGD', 
                    default=128)
parser.add_argument('-num_candidates', metavar='NUM_CANDIDATES', type=int,
                    help='amount of runs to try with optimal initial parameters', 
                    default=10)                   
parser.add_argument('-epsilon_min', metavar='EPSILON_MIN', type=float,
                    help='the learning rate will be sampled from a loguniform distribution with this as the minimum value', 
                    default=5e-4)
parser.add_argument('-epsilon_max', metavar='EPSILON_MAX', type=float, 
                    help='the learning rate will be sampled from a loguniform distribution with this as the maximum value', 
                    default=5e-2)
parser.add_argument('-sigma_motifs_min', metavar='SIGMA_MOTIFS_MIN', type=float, 
                    help='the initial motif weights will be sampled from a loguniform distribution with this as the minimum value', 
                    default=1e-7)
parser.add_argument('-sigma_motifs_max', metavar='SIGMA_MOTIFS_MAX', type=float, 
                    help='the initial motif weights will be sampled from a loguniform distribution with this as the maximum value', 
                    default=1e-3)
parser.add_argument('-sigma_net_min', metavar='SIGMA_NET_MIN', type=float, 
                    help='the initial neural network weights will be sampled from a loguniform distribution with this as the minimum value', 
                    default=1e-5)
parser.add_argument('-sigma_net_max', metavar='SIGMA_NET_MAX', type=float, 
                    help='the initial neural network weights will be sampled from a loguniform distribution with this as the maximum value', 
                    default=1e-2)
parser.add_argument('-seed', metavar='SEED', type=int,
                    help='the seed for random number generators from torch and numpy',
                    default=42)


args = parser.parse_args()

name = args.name
data = args.data
celldata = args.celldata
model_output_dir = args.model_output_dir

val_factor = args.v
num_calibrations = args.c
num_candidates = args.num_candidates
epochs = args.epochs
num_errtest = args.num_errtest
batch_size = args.batch_size
optimizer = args.opt
criterion = "SGD"
m = args.m
d = args.d
epsilon_min = args.epsilon_min
epsilon_max = args.epsilon_max
sigma_motifs_min = args.sigma_motifs_min
sigma_motifs_max = args.sigma_motifs_max
sigma_net_min = args.sigma_net_min
sigma_net_max = args.sigma_net_max



# initialize random number generators
seed = args.seed
torch.cuda.manual_seed_all(seed)
torch.manual_seed(seed)
np.random.seed(seed)



print(file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
print("                     scanem pytorch version v0.2                     ", file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
print(file=sys.stderr)

print("Name \t\t\t=", name, file=sys.stderr)
print("Data path \t\t=", data, file=sys.stderr)
print("Celltype data path \t=", celldata, file=sys.stderr)
print("Model output dir \t=", model_output_dir, file=sys.stderr)
print("Validation factor \t=", str(val_factor), file=sys.stderr)
print("Calibrations \t\t=", str(num_calibrations), file=sys.stderr)
print("Candidates \t\t=", str(num_candidates), file=sys.stderr)
print("Epochs \t\t\t=", str(epochs), file=sys.stderr)
print("Num_errtest \t\t=", str(num_errtest), file=sys.stderr)
print("Batch size \t\t=", str(batch_size), file=sys.stderr)
print("Optimizer \t\t=", optimizer, file=sys.stderr)
print("Number of motifs \t=", str(d), file=sys.stderr)
print("Motif length \t\t=", str(m), file=sys.stderr)

# Device to run on
if(torch.cuda.is_available()):
    if(torch.cuda.device_count() > 1):
        device=torch.device('cuda:1')
    else:
        device=torch.device('cuda:0')
else:
    device=torch.device('cpu')
    torch.set_num_threads(24)

print("Device\t\t\t=", device, file=sys.stderr)


print(file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
print(file=sys.stderr)

exp_model_dir = model_output_dir
if not os.path.exists(exp_model_dir):
    os.makedirs(exp_model_dir)
    print("NOTE: \tCreated model directory", file=sys.stderr)
else:
    print("NOTE: \tDirectory exists, will overwrite previous files", file=sys.stderr)

print("NOTE: \tPreprocessing data", file=sys.stderr)

input_data = pd.read_csv(data, sep="\t")
input_seqs = np.asarray([su.onehot_encode_seq(x, m, True) 
                    for x in input_data["sequence"]])
n = len([q for q in input_data["sequence"][1]])
input_ind = np.asarray([[float(x) for x in y.split(",")] for y in input_data["ind"]])
cells = len(input_ind[0,:])

input_seqs = torch.from_numpy(input_seqs).unsqueeze(1)
input_ind = torch.from_numpy(input_ind)

input_seqs = input_seqs.to(device)
input_ind = input_ind.to(device)

train_data = TensorDataset(input_seqs, input_ind)
cross_vals = KFold(n_splits = val_factor, shuffle=True, random_state=None)
cross_vals_2 = KFold(n_splits = num_candidates, shuffle=True, random_state=None)

errtest_len = sm.get_errtest_len(epochs, num_errtest)


cell_type_names = pd.read_csv(celldata, sep="\t")
network_labels = cell_type_names["cell_type1"]
network_pal = sns.cubehelix_palette(network_labels.unique().size,
                                    light=.9, dark=.1, reverse=True,
                                    start=1, rot=-2)
network_lut = dict(zip(map(str, network_labels.unique()), network_pal))
network_colors = pd.Series(network_labels).map(network_lut)


print("NOTE: \tDone preprocessing data", file=sys.stderr)
print(file=sys.stderr)

print("---------------------------------------------------------------------", file=sys.stderr)
print("                             run start", file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)

all_val_losses = np.empty((num_calibrations,val_factor,errtest_len))
all_train_losses = np.empty((num_calibrations,val_factor,errtest_len))

fold_mean_minima = []

cal = 0
for i in range(num_calibrations):

    # Initialize model
    model = sm.ConvNet(optimizer, 
                    "MSE", 
                    epsilon_min, 
                    epsilon_max, 
                    sigma_motifs_min, 
                    sigma_motifs_max, 
                    sigma_net_min,
                    sigma_net_max, 
                    d, 
                    m, 
                    n, 
                    cells).double()

    model = model.to(device)

    model_lr = model.learning_rate

    # Save initial parameters for every calibration
    cal_intial_model_path = exp_model_dir + "/" + name + "_cal_" + str(cal+1) + "_initial_state_dict.pt"
    torch.save(model.state_dict(), cal_intial_model_path)

    # Perform K-fold cross-validation
    fold = 0
    for train_ind, val_ind in cross_vals.split(train_data):
        model.load_state_dict(torch.load(cal_intial_model_path))

        # Make this ridiculously high so that the model gets a better score than this
        best_val_loss = 1e100

        index_for_errtest = 0
        train_sampler = SubsetRandomSampler(train_ind)
        val_sampler = SubsetRandomSampler(val_ind)

        # To cycle through the data
        train_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, sampler=train_sampler, 
                                                num_workers = 0)
        val_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, sampler=val_sampler, 
                                                num_workers = 0)

        # Start out with loss being 0.0
        running_loss = 0.0
        val_running_loss = 0.0

        model_name = name + "_cal_" + str(cal+1) + "_fold_" + str(fold+1)
        best_model_save_path = exp_model_dir + "/" + model_name + "_state_dict.pt"
        # Progress bar of length num_learningsteps
        trange_name = "Run " + str(cal+1) + "_" + str(fold+1) + " of " + str(num_calibrations) + "_" + str(val_factor)
        t = trange(epochs, desc=trange_name, ncols=100) 
        for j in t:
            running_loss = 0.0
            val_running_loss = 0.0

            # Cycle through training set by batches
            for _, data in enumerate(train_loader, 0):
                # Get batch data
                batch_seqs, batch_ind = data

                # Zero the parameter gradients for each minibatch
                model.optimizer.zero_grad()

                # Forward + compare + backward + optimize
                outputs = model(batch_seqs)
                loss = model.criterion(outputs, batch_ind)
                loss.backward()
                model.optimizer.step()

                # Compute loss and add to running_loss
                running_loss += loss.item()

            # Cycle through validation set
            for _, data in enumerate(val_loader, 0):
                # Get validation data
                val_seqs, val_ind = data

                # Calculate outputs and loss, add to validation loss
                val_outputs = model(val_seqs)
                val_loss = model.criterion(val_outputs, val_ind)
                val_running_loss += val_loss.item()

            # When epoch is divisible by 'num_errtest', store loss and update progress bar
            if (j+1) % num_errtest == 0:      
                # Change the progress bar to show loss and validation loss
                t.set_description(trange_name + ' (loss=%g, val_loss=%g)' % (running_loss,val_running_loss))

                all_val_losses[cal,fold,index_for_errtest] = val_running_loss
                all_train_losses[cal,fold,index_for_errtest] = running_loss

                if(val_running_loss) < best_val_loss:
                    torch.save(model.state_dict(), best_model_save_path)
                    best_val_loss = val_running_loss

                running_loss = 0.0
                val_running_loss = 0.0
                index_for_errtest += 1

        fold += 1

    # Append mean of K-fold run minima to array
    #   -> to find out which cal performed best
    fold_mean_minima.append(np.mean(all_val_losses[cal,:,:].min(axis=1)))
    cal += 1    

# Get calibration with best mean minimum across candidates
best_model = fold_mean_minima.index(min(fold_mean_minima))
best_model_initial_path = exp_model_dir + "/" + name + "_cal_" + str(best_model+1) + "_initial_state_dict.pt"

print("Best model:", str(best_model+1), file=sys.stderr)

all_val_losses = np.empty((num_candidates,errtest_len))
all_train_losses = np.empty((num_candidates,errtest_len))

# Do 10x cross-validation, find motifs across 10 different subsets of data
# separately. Afterwards, compare which motifs were found consistently
new_splits = [(x,y) for (x,y) in cross_vals_2.split(train_data)]
train_inds2 = [x for (x,y) in new_splits]
val_inds2 = [y for (x,y) in new_splits]

candidate_train_val_ind_split_path = exp_model_dir + "/" + name + "_candidate_train_val_inds.p" 
pickle.dump(new_splits, open(candidate_train_val_ind_split_path, "wb" ) )

print("---------------------------------------------------------------------", file=sys.stderr)
print("                         training candidates", file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
can = 0
all_candidate_hdf5_paths = []

best_initial_model = sm.ConvNet(optimizer, 
                    "MSE", 
                    epsilon_min, 
                    epsilon_max, 
                    sigma_motifs_min, 
                    sigma_motifs_max, 
                    sigma_net_min,
                    sigma_net_max, 
                    d, 
                    m, 
                    n, 
                    cells).double()

# Load optimal initial parameters
best_initial_model.load_state_dict(torch.load(best_model_initial_path))

sigma_motifs = best_initial_model.sigma_motifs
sigma_net = best_initial_model.sigma_net
learning_rate = best_initial_model.learning_rate

for i in range(num_candidates):
    # Initialize candidate model with best initial parameters
    model = sm.BestInitialConvNet(optimizer, 
                    "MSE", 
                    learning_rate, 
                    sigma_motifs, 
                    sigma_net, 
                    d, 
                    m, 
                    n, 
                    cells).double()

    model = model.to(device)


    # Make this ridiculously high so that the model gets a better score than this
    best_val_loss = 1e100

    index_for_errtest = 0

    train_sampler = SubsetRandomSampler(train_inds2[i])
    val_sampler = SubsetRandomSampler(val_inds2[i])

    # To cycle through the data
    train_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, sampler=train_sampler, 
                                            num_workers = 0)
    val_loader = torch.utils.data.DataLoader(train_data, batch_size=batch_size, sampler=val_sampler, 
                                            num_workers = 0)

    # Start out with loss being 0.0
    running_loss = 0.0
    val_running_loss = 0.0

    trange_name = "Run_" + str(can+1) + " of " + str(num_candidates)

    model_name = name + "_candidate_" + str(can+1)
    candidate_save_path = exp_model_dir + "/" + model_name + "_state_dict.pt"

    # Progress bar of length num_learningsteps
    t = trange(epochs, desc=trange_name, ncols=100) 
    for j in t:
        running_loss = 0.0
        val_running_loss = 0.0

        # Cycle through training set by batches
        for _, data in enumerate(train_loader, 0):
            # Get batch data
            batch_seqs, batch_ind = data

            # Zero the parameter gradients for each minibatch
            model.optimizer.zero_grad()

            # Forward + compare + backward + optimize
            outputs = model(batch_seqs)
            loss = model.criterion(outputs, batch_ind)
            loss.backward()
            model.optimizer.step()

            # Compute loss and add to running_loss
            running_loss += loss.item()

        # Cycle through validation set
        for _, data in enumerate(val_loader, 0):
            # Get validation data
            val_seqs, val_ind = data

            # Calculate outputs and loss, add to validation loss
            val_outputs = model(val_seqs)
            val_loss = model.criterion(val_outputs, val_ind)
            val_running_loss += val_loss.item()

        # When epoch is divisible by 'num_errtest', store loss and update progress bar
        if (j+1) % num_errtest == 0:      
            # Change the progress bar to show loss and validation loss
            t.set_description(trange_name + ' (loss=%g, val_loss=%g)' % (running_loss,val_running_loss))

            all_val_losses[can,index_for_errtest] = val_running_loss
            all_train_losses[can,index_for_errtest] = running_loss

            if(val_running_loss) < best_val_loss:
                torch.save(model.state_dict(), candidate_save_path)
                best_val_loss = val_running_loss

            running_loss = 0.0
            val_running_loss = 0.0
            index_for_errtest += 1

    can += 1

candidate_minima = [x for x in all_val_losses.min(axis=1)]
best_candidate = candidate_minima.index(min(candidate_minima))
best_candidate_err = [x for x in all_val_losses[best_candidate,:]]
best_candidate_best_tp = best_candidate_err.index(min(best_candidate_err))

with open(exp_model_dir + "/Best_candidate_" + str(best_candidate+1) + ".txt", 'w') as the_file:
    the_file.write(str(best_candidate+1))

best_train_ind = train_inds2[best_candidate]
best_val_ind = val_inds2[best_candidate]

device = torch.device('cpu')

print("Best candidate:", str(best_candidate+1), file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
print("                         motif alignments", file=sys.stderr)
print("---------------------------------------------------------------------", file=sys.stderr)
best_candidate_best_model_path = exp_model_dir + "/" + name + "_candidate_" + str(best_candidate+1) + "_state_dict.pt"
bestmodel = sm.BestInitialConvNet(optimizer, 
                    "MSE", 
                    learning_rate, 
                    sigma_motifs, 
                    sigma_net, 
                    d, 
                    m, 
                    n, 
                    cells).double()
bestmodel = bestmodel.to(device)
bestmodel.load_state_dict(torch.load(best_candidate_best_model_path))
input_seqs = input_seqs.to(device)

# Aligning convolutional filters to train indices sequences, creating pwms
# based on which subsequences activate the different convolutional filters
best_motifs_pfm_dict, best_motifs_ppm_dict = su.align_conv_filters(bestmodel, input_seqs, input_data, m, best_train_ind)
meme_output_path = exp_model_dir + "/Motifs_" + name + "_Candidate_best_MEME.txt"

su.save_meme(best_motifs_ppm_dict, meme_output_path)

if(torch.cuda.is_available()):
    torch.cuda.empty_cache()

print("---------------------------------------------------------------------")


M = bestmodel.conv_1.weight.cpu().detach().numpy().squeeze()
w = bestmodel.fc.weight.cpu().detach().numpy()

w_df = pd.DataFrame(data=w, index=network_colors.index)
w_path = exp_model_dir + "/" + name + "_candidate_best_w_matrix.tsv.gz"
w_df.to_csv(w_path, sep='\t', compression='gzip')

clustermap_path = exp_model_dir + "/" + name + "_candidate_best_clustermap.png"
g = sns.clustermap(w_df, row_cluster=True, col_cluster=True, row_colors=network_colors)
for label in network_labels.unique():
    g.ax_col_dendrogram.bar(0, 0, color=network_lut[label],
                            label=label, linewidth=0)
g.ax_col_dendrogram.legend(loc='center left',bbox_to_anchor=(1.72,-2),frameon=True)
g.cax.set_position([1.24, .35, .02, .3])
g.savefig(clustermap_path)





# Losses plot ================================================================

plt.figure(figsize=(14,10))
for i in range(all_val_losses.shape[0]):
    plt.plot(np.transpose(all_val_losses[i,:]), label=str(i+1))
plt.legend(loc="upper left")
plt.title("Validation set losses for all " + str(all_val_losses.shape[0]) + " candidates")
plt.ylabel("Validation loss")
plt.xlabel("Epoch")
pltname = exp_model_dir + "/Plot_CandidateValLosses_" + name + ".png"
plt.savefig(pltname, bbox_inches="tight")
plt.close()


# Neural network weights plot ================================================

pltname = exp_model_dir + "/Plot_NeuralNetWeights_" + name + ".png"
plt.figure(figsize=(25,10))
ax = sns.heatmap(w, center=0)
plt.savefig(pltname, bbox_inches="tight")
plt.close()

# Convolutional filter plots =================================================

motif_dir = exp_model_dir + "/Motifs"
if not os.path.exists(motif_dir):
    os.makedirs(motif_dir)
    print("Created motif dir", file=sys.stderr)


# Plot filter weights
for i in range(w.shape[1]):
    motif_num = i
    logo_df = pd.DataFrame(data=M[i,:,:], 
                               index=range(M[motif_num,:,:].shape[0]),
                               columns=["A", "C", "G", "T"])
    crp_logo = logomaker.Logo(logo_df,
                              shade_below=.5,
                              fade_below=.5,
                              color_scheme='classic',
                              figsize=(4,3))
    plt.title("Motif "+str(motif_num))
    pltname = motif_dir + "/Plot_Motif_FilterWeight_" + name + "_Filter" + str(i) + ".png"
    plt.savefig(pltname, bbox_inches="tight")
    plt.close()

    
# Per-cell and per-gene correlations ========================================

val_actual = input_ind[best_val_ind].cpu().numpy()
val_prediction = bestmodel(input_seqs[best_val_ind]).detach().cpu().numpy()

randomly_permuted_val_seqs = np.asarray([su.onehot_encode_seq(x, m, True) 
                         for x in su.randomize_sequences(input_data["sequence"][best_val_ind])])
randomly_permuted_val_seqs = torch.from_numpy(randomly_permuted_val_seqs).unsqueeze(1)
randomly_permuted_val_seqs = randomly_permuted_val_seqs.to(device)
val_random_prediction = bestmodel(randomly_permuted_val_seqs).detach().cpu().numpy()

gene_corrs = []
cell_corrs = []
gene_corrs_rand = []
cell_corrs_rand = []

gene_corrs_spear = []
cell_corrs_spear = []
gene_corrs_spear_rand = []
cell_corrs_spear_rand = []


for gene_i in trange(val_actual.shape[0], desc="Calculating per-gene correlations"):
    if (np.std(val_actual[gene_i,:]) != 0) & (np.std(val_prediction[gene_i,:]) != 0):
        gene_corrs.append(np.corrcoef(val_actual[gene_i,:],val_prediction[gene_i,:])[0,1])
        gene_corrs_spear.append(spearmanr(val_actual[gene_i,:],val_prediction[gene_i,:])[0])
        gene_corrs_rand.append(np.corrcoef(val_actual[gene_i,:],val_random_prediction[gene_i,:])[0,1])
        gene_corrs_spear_rand.append(spearmanr(val_actual[gene_i,:],val_random_prediction[gene_i,:])[0])
    else:
        gene_corrs.append(0)
        gene_corrs_spear.append(0)
        gene_corrs_rand.append(0)
        gene_corrs_spear_rand.append(0)
    

for cell_i in trange(val_actual.shape[1], desc="Calculating per-cell correlations"):
    if (np.std(val_actual[:,cell_i]) != 0) & (np.std(val_prediction[:,cell_i]) != 0):
        cell_corrs.append(np.corrcoef(val_actual[:,cell_i],val_prediction[:,cell_i])[0,1])
        cell_corrs_spear.append(spearmanr(val_actual[:,cell_i],val_prediction[:,cell_i])[0])
        cell_corrs_rand.append(np.corrcoef(val_actual[:,cell_i],val_random_prediction[:,cell_i])[0,1])
        cell_corrs_spear_rand.append(spearmanr(val_actual[:,cell_i],val_random_prediction[:,cell_i])[0])
    else:
        cell_corrs.append(0)
        cell_corrs_spear.append(0)
        cell_corrs_rand.append(0)
        cell_corrs_spear_rand.append(0)

gene_corrs = np.asarray(gene_corrs)
gene_corrs[np.isnan(gene_corrs)] = 0
gene_corrs_rand = np.asarray(gene_corrs_rand)
gene_corrs_rand[np.isnan(gene_corrs_rand)] = 0
    
cell_corrs = np.asarray(cell_corrs)
cell_corrs[np.isnan(cell_corrs)] = 0
cell_corrs_rand = np.asarray(cell_corrs_rand)
cell_corrs_rand[np.isnan(cell_corrs_rand)] = 0

gene_corrs_spear = np.asarray(gene_corrs_spear)
gene_corrs_spear[np.isnan(gene_corrs_spear)] = 0
gene_corrs_spear_rand = np.asarray(gene_corrs_spear_rand)
gene_corrs_spear_rand[np.isnan(gene_corrs_spear_rand)] = 0
    
cell_corrs_spear = np.asarray(cell_corrs_spear)
cell_corrs_spear[np.isnan(cell_corrs_spear)] = 0
cell_corrs_spear_rand = np.asarray(cell_corrs_spear_rand)
cell_corrs_spear_rand[np.isnan(cell_corrs_spear_rand)] = 0

all_corrs = np.concatenate((gene_corrs, cell_corrs))
all_corrs_spear = np.concatenate((gene_corrs_spear, cell_corrs_spear))

all_corrs_rand = np.concatenate((gene_corrs_rand, cell_corrs_rand))
all_corrs_spear_rand = np.concatenate((gene_corrs_spear_rand, cell_corrs_spear_rand))

corrs_df = pd.DataFrame(data={"type": ["Per-gene"] * gene_corrs.shape[0] + ["Per-cell"] * cell_corrs.shape[0],
                              "pearson_correlation": all_corrs,
                              "spearman_correlation": all_corrs_spear})
corrs_rand_df = pd.DataFrame(data={"type": ["Per-gene"] * gene_corrs_rand.shape[0] + ["Per-cell"] * cell_corrs_rand.shape[0],
                                   "pearson_correlation": all_corrs_rand,
                                   "spearman_correlation": all_corrs_spear_rand})

plt.figure(figsize=(8,8))
ax = sns.violinplot(x="type", y="pearson_correlation", inner="quart",
               data=corrs_df)
ax.set(xlabel='', ylabel='Pearson correlation')
ax.set(ylim=(-1, 1))
plt.title("Per-gene and per-cell Pearson correlations between observed and predicted values")
pltname = exp_model_dir + "/Plot_PerGeneCellCorrelationPearson_" + name + ".png"
plt.savefig(pltname, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,8))
ax = sns.violinplot(x="type", y="spearman_correlation", inner="quart",
               data=corrs_df)
ax.set(xlabel='', ylabel='Spearman correlation')
ax.set(ylim=(-1, 1))
plt.title("Per-gene and per-cell Spearman correlations between observed and predicted values")
pltname = exp_model_dir + "/Plot_PerGeneCellCorrelationSpearman_" + name + ".png"
plt.savefig(pltname, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,8))
ax = sns.violinplot(x="type", y="pearson_correlation", inner="quart",
               data=corrs_rand_df)
ax.set(xlabel='', ylabel='Pearson correlation')
ax.set(ylim=(-1, 1))
plt.title("Per-gene and per-cell Pearson correlations between observed and predicted values using randomly permuted sequences")
pltname = exp_model_dir + "/Plot_PerGeneCellCorrelationPearsonRANDOM_" + name + ".png"
plt.savefig(pltname, bbox_inches="tight")
plt.close()

plt.figure(figsize=(8,8))
ax = sns.violinplot(x="type", y="spearman_correlation", inner="quart",
               data=corrs_rand_df)
ax.set(xlabel='', ylabel='Spearman correlation')
ax.set(ylim=(-1, 1))
plt.title("Per-gene and per-cell Spearman correlations between observed and predicted values using randomly permuted sequences'")
pltname = exp_model_dir + "/Plot_PerGeneCellCorrelationSpearmanRANDOM_" + name + ".png"
plt.savefig(pltname, bbox_inches="tight")
plt.close()


metrics_filename = exp_model_dir + "/" + name + "_metrics_for_best_candidate.txt"
avg_per_gene_spearmanr = np.mean(gene_corrs_spear)
avg_per_cell_spearmanr = np.mean(cell_corrs_spear)
dataset_sparsity = 1.0 - (np.count_nonzero(input_ind.detach().cpu().numpy())/float(input_ind.detach().cpu().numpy().size))
prediction_spearmanr = spearmanr(val_actual.flatten(), val_prediction.flatten())[0]
with open(metrics_filename, 'w') as curr_file:
    curr_file.write("Run metrics for " + name + "\n")
    curr_file.write("Dataset sparsity = " + str(dataset_sparsity) + "\n")
    curr_file.write("Prediction spearman R = " + str(prediction_spearmanr) + "\n")
    curr_file.write("Mean per-gene spearman R = " + str(avg_per_gene_spearmanr) + "\n")
    curr_file.write("Mean per-cell spearman R = " + str(avg_per_cell_spearmanr) + "\n")


# Aligning other candidates

other_candidate_numbers = [x for x in np.asarray(range(all_val_losses.shape[0]))[np.asarray(range(all_val_losses.shape[0])) != best_candidate]]
can_motifs_ppm_dicts = []
can_motifs_pfm_dicts = []
d = M.shape[0]
for can_nr in other_candidate_numbers:
    candidate_model_path = exp_model_dir + "/" + name + "_candidate_" + str(can_nr+1) + "_state_dict.pt"
    candidate_model = sm.BestInitialConvNet(optimizer, 
                    "MSE", 
                    learning_rate, 
                    sigma_motifs, 
                    sigma_net, 
                    d, 
                    m, 
                    n, 
                    cells).double()
    candidate_model = candidate_model.to(device)
    candidate_model.load_state_dict(torch.load(candidate_model_path))
    input_seqs = input_seqs.to(device)
    
       
    w_path = exp_model_dir + "/" + name + "_candidate_" + str(can_nr+1) + "_w_matrix.tsv.gz"
    w_df = pd.DataFrame(data=candidate_model.fc.weight.cpu().detach().numpy(), index=network_colors.index)
    w_df.to_csv(w_path, sep='\t', compression='gzip')

    g = sns.clustermap(w_df, row_cluster=True, col_cluster=True, row_colors=network_colors)
    for label in network_labels.unique():
        g.ax_col_dendrogram.bar(0, 0, color=network_lut[label],
                                label=label, linewidth=0)
    g.ax_col_dendrogram.legend(loc='center left',bbox_to_anchor=(1.72,-2),frameon=True)
    g.cax.set_position([1.24, .35, .02, .3])

    clustermap_path = exp_model_dir + "/" + name + "_candidate_" + str(can_nr+1) + "_clustermap.png"
    g.savefig(clustermap_path)

    # Aligning convolutional filters to train indices sequences, creating pwms
    # based on which subsequences activate the different convolutional filters
    can_motifs_pfm_dict, can_motifs_ppm_dict = su.align_conv_filters(candidate_model, input_seqs, input_data, m, train_inds2[can_nr])
    can_motifs_ppm_dicts.append(can_motifs_ppm_dict)
    can_motifs_pfm_dicts.append(can_motifs_pfm_dict)
    
all_candidates = [str(x+1) for x in np.asarray(range(all_val_losses.shape[0]))]    
all_candidates_with_best = [x for x in all_candidates]
all_candidates_with_best[best_candidate] = "best"

# Save w (final layer weights) and M (convolutional kernels) to hdf5 file for
# each candidate. 
all_model_HDF5_path = exp_model_dir + "/All_model_HDF5_" + name + "_M_w.h5"
fid = h5py.File(all_model_HDF5_path, 'w')
for idx, can in enumerate(all_candidates):
    candidate_state_dict = exp_model_dir + "/" + name + "_candidate_" + can + "_state_dict.pt"
    candidate_model = sm.BestInitialConvNet(optimizer, 
                    "MSE", 
                    learning_rate, 
                    sigma_motifs, 
                    sigma_net, 
                    d, 
                    m, 
                    n, 
                    cells).double()
    candidate_model.load_state_dict(torch.load(candidate_state_dict, map_location=torch.device('cpu')))
    candidate_model = candidate_model.eval()
    candidate_model = candidate_model.double()
    
    M = candidate_model.conv_1.weight.cpu().detach().numpy().squeeze()
    w = candidate_model.fc.weight.cpu().detach().numpy()
    
    fid["M_" + all_candidates_with_best[idx]] = M
    fid["w_" + all_candidates_with_best[idx]] = w
fid.close()

all_motifs_ppm_dict = {}
all_motifs_pfm_dict = {}
for i,can_nr in enumerate(other_candidate_numbers):
    for key in can_motifs_ppm_dicts[i].keys():
        all_motifs_ppm_dict[str(can_nr+1)+"_"+key] = can_motifs_ppm_dicts[i][key]
        all_motifs_pfm_dict[str(can_nr+1)+"_"+key] = can_motifs_pfm_dicts[i][key]
for key in best_motifs_ppm_dict.keys():
    all_motifs_ppm_dict["best_"+key] = best_motifs_ppm_dict[key]
    all_motifs_pfm_dict["best_"+key] = best_motifs_pfm_dict[key]
    
# Plot all motif logos (also from other candidates)
all_motif_dir = exp_model_dir + "/AllMotifs"
if not os.path.exists(all_motif_dir):
    os.makedirs(all_motif_dir)
    print("Created all_motif_dir", file=sys.stderr)
    
for key in all_motifs_ppm_dict.keys():
    motif_ic_mat = logomaker.transform_matrix(pd.DataFrame(all_motifs_ppm_dict[key], 
                                                          columns=['A','C','G','T']), 
                                              from_type='probability',
                                              to_type='information')
    crp_logo = logomaker.Logo(motif_ic_mat,
                              shade_below=.5,
                              fade_below=.5,
                              color_scheme='classic',
                              figsize=(4,3))
    plt.title("Motif "+str(key))
    pltname = all_motif_dir + "/Plot_Motif_Alignment_" + name + "_Filter" + str(key) + ".png"
    plt.savefig(pltname, bbox_inches="tight")
    plt.close()
    
    
# Save big MEME file with all motifs
    
all_meme_output_path = exp_model_dir + "/All_MEME_motifs_" + name + ".txt"
su.save_meme(all_motifs_ppm_dict, all_meme_output_path)    

all_motif_ppm_object_path = exp_model_dir + "/All_motifs_ppm_dict_" + name + ".p"
pickle.dump(all_motifs_ppm_dict, open(all_motif_ppm_object_path, "wb" ) )

all_motif_pfm_object_path = exp_model_dir + "/All_motifs_pfm_dict_" + name + ".p"
pickle.dump(all_motifs_pfm_dict, open(all_motif_pfm_object_path, "wb" ) )

# Cleaning up calibration files ===========================================

files_to_delete = []
# r=root, d=directories, f = files
for r, d, f in os.walk(exp_model_dir):
    for file in f:
        if 'pt' in file and ('fold' in file) or ('initial' in file):
            files_to_delete.append(os.path.join(r, file))
            # (we don't care about initial model states anymore)
if len(files_to_delete) > 0:
    for file in files_to_delete:
        print('rm ' + file)
        os.remove(file)
        
print("Done", file=sys.stderr)
