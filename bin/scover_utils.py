def onehot_encode_seq(sequence, m=0, padding=False):
    """Converts a given IUPAC DNA sequence to a one-hot 
    encoded DNA sequence. 
    """
    import numpy as np
    import torch
    
    valid_keys = ['a','c','g','t','u','n','r','y','s','w','k','m']
    
    nucs = {'a':0,'c':1,'g':2,'t':3,'u':3}
    
    if padding:
        assert m != 0, "If using padding, m should be bigger than 0"
        padding_mat = np.tile(0.25,(m-1,4))
    onehot = np.tile(.0,(len(sequence),4))
    
    for i,char in enumerate(sequence.lower()):
        if char not in valid_keys:
            sys.exit("invalid char in sequence (choose from acgt and nryswkm)")
        elif char == 'n':
            onehot[i,:] = 0.25
        elif char == 'r':
            onehot[i,(0,2)] = 0.5
        elif char == 'y':
            onehot[i,(1,3)] = 0.5
        elif char == 's':
            onehot[i,(1,2)] = 0.5
        elif char == 'w':
            onehot[i,(0,3)] = 0.5
        elif char == 'k':
            onehot[i,(2,3)] = 0.5
        elif char == 'm':
            onehot[i,(0,1)] = 0.5
        else:
            onehot[i,nucs[char]] = 1
    
    if padding:
        onehot = np.concatenate((padding_mat, onehot, padding_mat))

    return onehot

def save_meme(motifs_ppm_dict, output_file="found_motifs.meme"):
    """Saves the found PPMs (given as dictionary) to a file that's
    compatible with MEME suite applications.
    """
    import pandas as pd

    meme_string = ["MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", ""]
    for idx,key in enumerate(motifs_ppm_dict.keys()):
        curr_motif = pd.DataFrame(motifs_ppm_dict[key])
        s1 = "MOTIF " + str(key)
        s2 = "letter-probability matrix: alength= " + str(curr_motif.shape[1]) + " w= " + str(curr_motif.shape[0])
        s3 = curr_motif.to_csv(sep="\t", index=False, header=False)
        meme_string = meme_string + [s1, s2, s3]
    
    meme_string = "\n".join(meme_string)
    
    with open(output_file, 'w') as the_file:
        the_file.write(meme_string)
        
    print("wrote meme list")

def align_conv_filters(model, input_seqs, m, train_ind):
    """Aligns the convolutional filters of a given scover model back
    to the given input sequences at the given indices. 
    """
    # Motif analysis
    import numpy as np
    import torch
    from tqdm import trange

    activation_seqs = input_seqs[train_ind]
    
    with torch.no_grad():
        model.eval()
        activations = model.conv_1(activation_seqs).cpu().detach().numpy().squeeze()
    
    n_seq = activation_seqs.shape[0]    
    activation_seqs = activation_seqs.squeeze()
    seq_len = activation_seqs.shape[1]
    d = activations.shape[1]
    
    motifs_pfm_dict = dict() # store pfms in this dict
    motifs_ppm_dict = dict() # store pwms in this dict
    
    # cycle through convolutional filters
    for filter_num in trange(d):
        
        # select activations for filter. new array = nseq x length seq
        curr_activation = activations[:,filter_num,:] 
        
        # get those sequences that have positive values
        seq_has_pos_vals = np.argwhere(np.amax(curr_activation, axis=1) > 0)[:,0]
        
        # in the case that there is a minmum of 10 sequences that activate the filter
        if seq_has_pos_vals.shape[0] > 10: 
            
            # per sequence, get position of maximum activation
            per_seq_where_max_pos = np.argmax(curr_activation[seq_has_pos_vals], axis=1)
            curr_activation_seqs = activation_seqs[seq_has_pos_vals]
            curr_str_list = []
            # go through sequences and save to curr_str_list
            for i in range(seq_has_pos_vals.shape[0]): 
                
                # maximum activation
                curr_max = per_seq_where_max_pos[i] 
                # get subsequence that activated filter (max 1 per seq)
                curr_str_list.append(curr_activation_seqs[i][curr_max:(curr_max+m)]) 
    
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

def randomize_sequences(sequences):
    """Randomly permutes a set of DNA sequences.
    """
    import random
    shuffled_seqs = []
    for seq in sequences:
        shuffled_seqs.append(''.join(random.sample(seq, len(seq))))
    return shuffled_seqs
