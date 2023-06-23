import random
import re
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
from sklearn.model_selection import KFold, train_test_split
import scanpy as sc
import scipy
import anndata as annd
from tqdm.auto import tqdm
from scipy.stats import pearsonr, spearmanr


def get_expression_hits_df(ad, 
                           all_loo_scores_df: pd.DataFrame,
                           motif_families: pd.DataFrame,
                           fdr_correct = True,
                           gene_symbols: str = 'Symbol', 
                           correlation = 'spearman',
                           use_cat = 'subfamily'):
    r"""
    Experimental.
    """
    assert correlation.lower() in ['pearson', 'spearman']
    all_expression_hits_df = []
    # cycle through motif sub-families
    for subfam in motif_families[use_cat].unique():
        curr_sub_df = motif_families[motif_families[use_cat] == subfam]
        cluster_filters = list(curr_sub_df['motif'].unique())
        cluster_genes = curr_sub_df['family_genes'][0]
        cluster_mean_weights = np.array(all_loo_scores_df.loc[cluster_filters].mean(0))
        cluster_genes = [x for x in cluster_genes if x in list(ad.var[gene_symbols])]
        
        if len(cluster_genes) == 0:
            print("Cluster " + subfam + " has no expressed genes")
            break
        all_corrs = []
        all_pvals = []
        all_mean_expression = []
        used_genes = []
        for j in cluster_genes:
            gene_idx = np.argwhere(np.array(ad.var[gene_symbols]) == j)[0,0]
            values = np.array(ad[:,gene_idx].X[:,0])
            if sum(values) == 0:
                pass
            else:                
                if correlation.lower() == 'pearson':
                    corr, p_val = pearsonr(values, cluster_mean_weights)
                else:
                    corr, p_val = spearmanr(values, cluster_mean_weights)
                all_corrs.append(corr)
                all_pvals.append(p_val)
                all_mean_expression.append(values.mean())
                used_genes.append(j)
        if fdr_correct:
            all_pvals = fdr(np.array(all_pvals))
        expression_hits_df = pd.DataFrame({'family': [subfam for x in range(len(used_genes))],
                                           'gene': used_genes, 
                                           'mean_expression': all_mean_expression,
                                           'correlation': all_corrs, 
                                           'pval': all_pvals})
        all_expression_hits_df.append(expression_hits_df)
    return pd.concat(all_expression_hits_df)


def onehot_seq(sequence, m=0, padding=False) -> np.ndarray:
    r"""
    Converts IUPAC sequence to one-hot-encoded numpy array with
      colums corresponding to ['A','C','G','T'].
    """
    import numpy as np
    import sys
    valid_keys = ['a', 'c', 'g', 't', 'u', 'n', 'r', 'y', 's', 'w', 'k', 'm']
    nucs = {'a': 0, 'c': 1, 'g': 2, 't': 3, 'u': 3}  # u allows for RNA seq
    if padding:
        assert m != 0, "If using padding, m should be bigger than 0"
        padding_mat = np.tile(0.25, (m - 1, 4))
    onehot = np.tile(.0, (len(sequence), 4))
    for i, char in enumerate(sequence.lower()):
        if char not in valid_keys:
            sys.exit("invalid char in sequence (choose from acgt and nryswkm)")
        elif char == 'n':  # if unknown char: equal p across ACGT
            onehot[i, :] = 0.25
        elif char == 'r':
            onehot[i, (0, 2)] = 0.5
        elif char == 'y':
            onehot[i, (1, 3)] = 0.5
        elif char == 's':
            onehot[i, (1, 2)] = 0.5
        elif char == 'w':
            onehot[i, (0, 3)] = 0.5
        elif char == 'k':
            onehot[i, (2, 3)] = 0.5
        elif char == 'm':
            onehot[i, (0, 1)] = 0.5
        else:
            onehot[i, nucs[char]] = 1
    if padding:
        onehot = np.concatenate((padding_mat, onehot, padding_mat))
    return onehot


def seq_list_to_conv(seq_list,
                     m=0,
                     padding=False,
                     expand=False,
                     expand_dim=1) -> np.ndarray:
    r"""
    Converts list of sequences to numpy tensor of one-hot encoded
      numpy arrays where dimension 1 corresponds to the sequence
      index, dimension 2 corresponds to the sequence position, 
      and dimension 3 corresponds to ['A','C','G','T'].
    """
    import numpy as np

    if expand:
        return np.expand_dims(
            np.stack([onehot_seq(x, m, padding) for x in seq_list]),
            expand_dim)
    else:
        return np.stack([onehot_seq(x, m, padding) for x in seq_list])


def pool_anndata(adata, 
                 neighbors=40, 
                 n_seed_cells=500, 
                 seed_selection=None,
                 mean=False, 
                 return_selection=False,
                 do_neighbors=True,
                 log_transform=True,
                 use_rep=None):
    r"""
    EXPERIMENTAL
    """
    import anndata as annd
    
    if do_neighbors:
        if use_rep is not None:
            sc.pp.neighbors(adata, n_pcs = 40, n_neighbors = neighbors,
                            use_rep=use_rep)
        else:
            sc.pp.neighbors(adata, n_pcs = 40, n_neighbors = neighbors)
    connectivities = adata.obsp['connectivities'].toarray()
    idx = connectivities.argsort()
    kn_neighbours=np.ones((connectivities.shape[0],neighbors), dtype=np.int64)
    for i in range(connectivities.shape[0]):
        kn_neighbours[i,:] = idx[i][::-1][0:neighbors]
    if seed_selection is None:
        seed_selection = np.random.choice(connectivities.shape[0], 
                                          n_seed_cells, 
                                          replace=False)
    pseudobulk_idx = kn_neighbours[seed_selection,:]
    if scipy.sparse.issparse(adata.X):
        X_data = np.copy(adata.X.toarray())
    else:
        X_data = np.copy(adata.X)
    pseudobulk_data = np.empty((n_seed_cells,adata.shape[1]))
    if mean:
        for idx,nbs in enumerate(tqdm(pseudobulk_idx)):
            pseudobulk_data[idx,:] = X_data[nbs,:].mean(0)
    else:
        for idx,nbs in enumerate(tqdm(pseudobulk_idx)):
            pseudobulk_data[idx,:] = X_data[nbs,:].sum(0)
    del X_data
    if log_transform:
        pseudobulk_data = np.log10(1 + pseudobulk_data)
    pseudo = pd.DataFrame(pseudobulk_data, 
                      index=["pseudobulk_" + str(x) for x in range(pseudobulk_data.shape[0])],
                      columns=adata.var_names.copy())
    if return_selection:
        return annd.AnnData(pseudo, var=adata.var.copy()), seed_selection, pseudobulk_idx
    else:
        return annd.AnnData(pseudo, var=adata.var.copy())


def pool_anndata_given_pseudobulk_idx(adata, 
                                      pseudobulk_idx,
                                      mean=False,
                                      log_transform=True):
    r"""
    EXPERIMENTAL
    """
    import anndata as annd
    if scipy.sparse.issparse(adata.X):
        X_data = np.copy(adata.X.toarray())
    else:
        X_data = np.copy(adata.X)
    n_seed_cells = pseudobulk_idx.shape[0]
    pseudobulk_data = np.empty((n_seed_cells,adata.shape[1]))
    if mean:
        for idx,nbs in enumerate(tqdm(pseudobulk_idx)):
            pseudobulk_data[idx,:] = X_data[nbs,:].mean(0)
    else:
        for idx,nbs in enumerate(tqdm(pseudobulk_idx)):
            pseudobulk_data[idx,:] = X_data[nbs,:].sum(0)
    del X_data
    if log_transform:
        pseudobulk_data = np.log10(1 + pseudobulk_data)
    pseudo = pd.DataFrame(pseudobulk_data, 
                      index=["pseudobulk_" + str(x) for x in range(pseudobulk_data.shape[0])],
                      columns=adata.var_names.copy())
    return annd.AnnData(pseudo, var=adata.var.copy())


class SeqDataset(Dataset):
    r"""
    Sequence dataset
    """
    def __init__(self, X_seqs, X_data, dtype='double'):
        """
        Args:
            X_seqs: array containing sequences for each gene (size [N x C x H x W])
            X_data: array containing values across samples for each gene (size [N x O]),
            dtype: data type for torch tensors. Choose from 'double' or 'float'
        """
        assert dtype in ["double", "float"
                         ], "invalid dtype (choose from 'double' or 'float')"
        self.X_seqs = torch.tensor(X_seqs)
        self.X_data = torch.tensor(X_data)
        if dtype == 'double':
            self.X_seqs = self.X_seqs.double()
            self.X_data = self.X_data.double()
        if dtype == 'float':
            self.X_seqs = self.X_seqs.float()
            self.X_data = self.X_data.float()
        self.output_size = self.X_data.shape[1]
        self.seq_length = self.X_seqs.shape[2]

    def __len__(self):
        return self.X_seqs.shape[0]

    def __getitem__(self, idx):
        return self.X_seqs[idx], self.X_data[idx]


def get_splits(ind_list: list, n_splits: int = 10) -> dict:
    all_splits = dict(
        zip(["test", "outer_train", "inner_train", "val"],
            [[] for x in range(n_splits)]))
    for idx, i in enumerate(KFold(n_splits=n_splits).split(ind_list)):
        i_outer_train = i[0]
        i_outer_test = i[1]
        i_inner_train, i_inner_val = train_test_split(i_outer_train,
                                                      test_size=0.2)
        all_splits["outer_train"].append(list(i_outer_train))
        all_splits["test"].append(list(i_outer_test))
        all_splits["inner_train"].append(list(i_inner_train))
        all_splits["val"].append(list(i_inner_val))
    return all_splits


def shan_ent(vect: np.ndarray) -> int:
    t_ent = 0
    vect = np.abs(vect)
    for x in vect:
        if x > 0:
            a = x / vect.sum()
            t_ent += (a) * np.log2(a)
    return -t_ent


def fdr(p_vals):
    r"""
    From https://stackoverflow.com/a/30748899
    """
    from scipy.stats import rankdata
    ranked_p_vals = rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_vals
    fdr[fdr > 1] = 1
    return fdr


def get_group_name(group_motifs, max_names=3):
    dot_dot_dot = False
    if len(group_motifs) > max_names:
        dot_dot_dot = True
    string_parts = group_motifs[0:min(max_names, len(group_motifs))]
    if dot_dot_dot:
        string_parts = string_parts + ["..."]
    return "/".join(string_parts)


def align_conv_filters(model, data_loader, device=torch.device("cuda")):
    """
    Aligns the convolutional filters of a given convolutional layer 
    to the given sequences.
    """
    # Motif analysis
    import numpy as np
    import torch
    from tqdm.auto import tqdm
    from torch.utils.data import DataLoader
    input_seqs = data_loader.dataset.X_seqs
    conv_layer = model.conv_1
    conv_layer = conv_layer.to(device)
    conv_layer.eval()
    n_seq = len(data_loader.dataset)
    seq_len = data_loader.dataset[0][0].shape[1]
    d = conv_layer.weight.shape[0]
    m = conv_layer.weight.shape[2]
    activations = torch.FloatTensor(size=(n_seq, d, seq_len-m+1, 1)) # assumes padding=0, strides=1
    curr_i = 0
    for idx, (seqs, vals) in enumerate(tqdm(data_loader)):
        curr_batch_size = len(seqs)
        seqs = seqs.to(device)
        with torch.no_grad():
            activations[curr_i:curr_i+curr_batch_size,:,:,:] = conv_layer(seqs).detach().cpu()
        curr_i += curr_batch_size
    activations = activations.squeeze().numpy()
    input_seqs = input_seqs.squeeze()
    motifs_pfm_dict = dict() # store pfms in this dict
    motifs_ppm_dict = dict() # store pwms in this dict
    # cycle through convolutional filters
    for filter_num in tqdm(range(d)):
        # select activations for filter. new array = nseq x length seq
        curr_activation = activations[:,filter_num,:] 
        # get those sequences that have positive values
        seq_has_pos_vals = np.argwhere(np.amax(curr_activation, axis=1) > 0)[:,0]
        # in the case that there is a minmum of 10 sequences that activate the filter
        if seq_has_pos_vals.shape[0] > 10: 
            # per sequence, get position of maximum activation
            per_seq_where_max_pos = np.argmax(curr_activation[seq_has_pos_vals], axis=1)
            curr_input_seqs = input_seqs[seq_has_pos_vals]
            curr_str_list = []
            # go through sequences and save to curr_str_list
            for i in range(seq_has_pos_vals.shape[0]): 
                # maximum activation
                curr_max = per_seq_where_max_pos[i] 
                # get subsequence that activated filter (max 1 per seq)
                curr_str_list.append(curr_input_seqs[i][curr_max:(curr_max+m)]) 
            # put them together in a numpy array
            sequence_array = np.stack(curr_str_list)
            # get sum per position
            sequence_array_summed = np.sum(sequence_array,axis=0) 
            # save pfm
            motifs_pfm_dict[str(filter_num)] = sequence_array_summed 
            # get counts per row
            row_sums = np.sum(sequence_array_summed, axis=1)
            # convert pfm to ppm
            sequence_array_summed = np.nan_to_num(sequence_array_summed / row_sums[:, np.newaxis])
            motifs_ppm_dict[str(filter_num)] = sequence_array_summed
    return motifs_pfm_dict, motifs_ppm_dict


def save_meme(motifs_ppm_dict, output_file="found_motifs.meme"):
    r"""
    Saves the found PPMs (given as dictionary) to a file that's
    compatible with MEME suite applications.
    """
    import pandas as pd
    meme_string = [
        "MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", ""
    ]
    for idx, key in enumerate(motifs_ppm_dict.keys()):
        curr_motif = pd.DataFrame(motifs_ppm_dict[key])
        s1 = "MOTIF " + str(key)
        s2 = "letter-probability matrix: alength= " + str(
            curr_motif.shape[1]) + " w= " + str(curr_motif.shape[0])
        s3 = curr_motif.to_csv(sep="\t", index=False, header=False)
        meme_string = meme_string + [s1, s2, s3]
    meme_string = "\n".join(meme_string)
    with open(output_file, 'w') as the_file:
        the_file.write(meme_string)
    print("wrote meme list")


def read_meme(meme_path: str) -> dict:
    r"""
    Reads PPMs from a MEME file. Buggy but seems to work for CIS-BP.
    """
    import numpy as np
    with open(meme_path, "r") as meme_file:
        meme_lines = meme_file.readlines()
    meme_lines = [x.strip() for x in meme_lines]
    MOTIF_lines = [x.startswith("MOTIF") for x in meme_lines]
    motif_ppm_dict = {}
    motif_len = int(meme_lines[np.argwhere(MOTIF_lines)[:, 0][0] +
                               1].split("w= ")[1])
    for i, line_i in enumerate(np.argwhere(MOTIF_lines)[:, 0]):
        curr_mot_name = meme_lines[line_i].split("MOTIF ")[1]
        curr_mot_lines = meme_lines[line_i + 2:line_i + 2 + motif_len]
        motif_ppm_dict[curr_mot_name] = np.array(
            [x.split("\t") for x in curr_mot_lines], dtype=float)
    return motif_ppm_dict


def to_z(mat, axis=0) -> np.ndarray:
    r"""
    Bit of an ugly bit of code
    """
    assert axis in [0,1], "axis should be 0 or 1"
    if isinstance(mat, pd.DataFrame):
        mat_index = mat.index.copy()
        mat_columns = mat.columns.copy()
        if axis == 1:
            curr_arr = np.nan_to_num(((mat.T - np.expand_dims(mat.mean(axis=0), 1)) / np.expand_dims(mat.std(axis=0), 1)).T)
            return pd.DataFrame(curr_arr, index = mat_index, columns = mat_columns)
        else:
            curr_arr = np.nan_to_num((mat - np.expand_dims(mat.mean(axis=1), 1)) / np.expand_dims(mat.std(axis=1), 1))
            return pd.DataFrame(curr_arr, index = mat_index, columns = mat_columns)
    else:
        if axis == 1:
            return np.nan_to_num(((mat.T - np.expand_dims(mat.mean(axis=0), 1)) / np.expand_dims(mat.std(axis=0), 1)).T)
        else:
            return np.nan_to_num((mat - np.expand_dims(mat.mean(axis=1), 1)) / np.expand_dims(mat.std(axis=1), 1))



def get_activations(model, data_loader, device=torch.device("cuda")):
    """
    Returns activations given model and data_loader
    """
    # Motif analysis
    import numpy as np
    import torch
    from tqdm.auto import tqdm
    conv_layer = model.conv_1
    conv_layer = conv_layer.to(device)
    conv_layer.eval()
    n_seq = len(data_loader.dataset)
    seq_len = train_loader.dataset[0][0].shape[1]
    d = conv_layer.weight.shape[0]
    m = conv_layer.weight.shape[2]
    activations = torch.FloatTensor(size=(n_seq, d, seq_len-m+1, 1)) # assumes padding=0, strides=1
    curr_i = 0
    for idx, (seqs, vals) in enumerate(tqdm(data_loader)):
        curr_batch_size = len(seqs)
        seqs = seqs.to(device)
        with torch.no_grad():
            activations[curr_i:curr_i+curr_batch_size,:,:,:] = conv_layer(seqs).detach().cpu()
        curr_i += curr_batch_size
    return activations


def create_alignment_df(tomtom_path, 
                        threshold = 0.05,
                        translate_ids = False,
                        db_meme_file = None,
                        id_translate_dict = None):
    r"""
    Experimental. 
    
    Create an alignment DataFrame from a TomTom alignment result.
    
    If field Target_ID in the dataframe contains IDs such as 
      'M5998_1.02', it might be useful to translate these using 
      translate_ids = True, in which case you need to provide
      db_meme_file or id_translate_dict. This is buggy; providing
      db_meme_file might only work for some CIS-BP databases
      provided by the MEME suite. 
    """
    if translate_ids:
        assert False in (db_meme_file == None, id_translate_dict == None), \
          "If translate_ids = True, provide db_meme_file or id_translate_dict"
    aln = pd.read_csv(tomtom_path, sep="\t", 
                      comment="#", engine="python", 
                      error_bad_lines=False)
    if db_meme_file != None:
        with open(db_meme_file, "r") as meme_file:
            db_lines = meme_file.readlines()
        all_ids = []
        all_names = []
        for x in db_lines:
            if x.startswith("MOTIF"):
                c_id, c_name = x.strip().split(" ")[1:]
                if c_name.startswith("("):
                    c_name = c_name.split("(")[1].split(")")[0]
                all_ids.append(c_id)
                all_names.append(c_name)
        id_translate_dict = dict(zip(all_ids, all_names))
    if translate_ids:
        aln["Target_code"] = aln["Target_ID"]
        aln["Target_ID"] = [id_translate_dict[x] for x in aln["Target_ID"]]
    aln = aln[aln["q-value"] <= threshold]
    return aln


def generate_alignment_graph(alignment_df,
                             cluster=True,
                             return_communities=True):
    r"""
    Experimental.
    Create an alignment graph from an alignment DataFrame.
    """
    from igraph import Graph
    query_set = set(alignment_df['Query_ID'])
    target_set = set(alignment_df['Target_ID'])
    n1 = len(query_set)
    n2 = len(target_set)
    g = Graph(n1 + n2, directed=False)
    g.vs["type"] = 0
    g.vs[n1:]["type"] = 1
    idx_mapping = dict(
        zip(
            sorted(query_set) + sorted(target_set),
            range(len(query_set) + len(target_set))))
    edges = list(
        zip([idx_mapping[x] for x in alignment_df["Query_ID"]],
            [idx_mapping[x] for x in alignment_df["Target_ID"]]))
    g.add_edges(edges)
    g.es["weight"] = list(-np.log10(alignment_df["q-value"]))
    g.vs["name"] = sorted(query_set) + sorted(target_set)
    if cluster:
        comms = g.community_walktrap().as_clustering()
    if return_communities:
        return g, comms
    else:
        return g


def plot_alignment_graph(alignment_graph, communities=None, plot_labels=False):
    r"""
    Experimental.
    Plot alignment graph. Optionally, add communities (clusters)
      and specify if you want to plot labels.
    """
    visual_style = {}
    visual_style["vertex_size"] = 4
    shape_dict = {0: 'rectangle', 1: 'circle'}
    color_dict = {0: 'lightblue', 1: 'salmon'}
    visual_style["vertex_shape"] = [
        shape_dict[x] for x in alignment_graph.vs["type"]
    ]
    visual_style["vertex_color"] = [
        color_dict[x] for x in alignment_graph.vs["type"]
    ]
    if plot_labels:
        visual_style["vertex_label"] = alignment_graph.vs["name"]
        visual_style["vertex_label_size"] = 8
        visual_style["vertex_label_angle"] = 1
    layout = alignment_graph.layout("kk")
    if communities != None:
        return igraph.plot(communities, mark_groups=True, **visual_style)
    else:
        return igraph.plot(alignment_graph, **visual_style)


def generate_motif_cluster_df(alignment_graph,
                              communities,
                              k_outer,
                              rep_cutoff: float = 0.4):
    r"""
    Experimental. 
    """
    groups = []
    for i in sorted(set(communities.membership)):
        groups.append(
            list(
                np.array(alignment_graph.vs["name"])[np.array(
                    communities.membership) == i]))
    curr_r = re.compile(
        '^[0-9]+_')  # For detecting convolutional filters VS database motifs
    group_reproducibilities = []
    for i in [list(filter(curr_r.match, x)) for x in groups]:
        group_reproducibilities.append(
            len(set([x.split("_")[0] for x in i])) / k_outer)
    group_filters = [
        i for i in [list(filter(curr_r.match, x)) for x in groups]
    ]
    curr_r = re.compile('^(?!^[0-9]+_)')
    group_motifs = [i for i in [list(filter(curr_r.match, x)) for x in groups]]
    groups_r = [
        x for i, x in enumerate(groups)
        if group_reproducibilities[i] >= rep_cutoff
    ]
    group_filters_r = [
        x for i, x in enumerate(group_filters)
        if group_reproducibilities[i] >= rep_cutoff
    ]
    group_motifs_r = [
        x for i, x in enumerate(group_motifs)
        if group_reproducibilities[i] >= rep_cutoff
    ]
    motif_cluster_df = pd.DataFrame({
        'group_name': [get_group_name(x) for x in group_motifs_r],
        'group':
        groups_r,
        'group_filters':
        group_filters_r,
        'group_motifs':
        group_motifs_r,
        'reproducibility':
        list(filter(lambda x: x >= rep_cutoff, group_reproducibilities)),
        'num_motifs': [len(x) for x in group_filters_r]
    })
    return motif_cluster_df
