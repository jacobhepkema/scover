# Model files

## General guidelines
Each of the trained models contains 10 models: one for each of the data folds.

## Use of model files
To load the models, adapt the code from the third notebook in the `example_notebooks` directory, for instance, the following code loads all ten models and looks at test set prediction. This assumes that you already have loaded the dataset etc; see the notebook for the full run-through. I've highlighted the part of the code that concerns loading the model parameters.

```python
from scipy.stats import pearsonr

corrs = []
for curr_outer in tqdm(range(k_outer)):
    curr_test_idx = idx_sets['test'][curr_outer]
    test_data = SeqDataset(X_seqs[curr_test_idx], X_data[curr_test_idx], dtype="float")
    test_loader = DataLoader(test_data, batch_size=128, num_workers=0, shuffle=False)
    best_config = all_best_configs[curr_outer]
    sn = SeqNet(seq_length=seq_length, 
                     output_size=output_size, 
                     learning_rate=best_config["learning_rate"],
                     motif_length=12,
                     num_motifs=600,
                     sigma_motifs=best_config["sigma_motifs"],
                     sigma_net=best_config["sigma_net"],
                     use_elu=True)

    ################## Loading the model ##############################

    state_dict_path = f"models/{exp_name}_trained_outer_" + str(
        curr_outer) + ".p"      # Adjust this to reflect the model of choice

    sn.load_state_dict(torch.load(state_dict_path))

    ###################################################################

  
    if torch.cuda.is_available():
        sn = sn.cuda()
    sn = sn.eval()
    
    curr_predictions = np.zeros((len(test_data), output_size))
    curr_ground_truth = np.zeros((len(test_data), output_size))
    curr_start_idx = 0
    for idx, (X,Y) in enumerate(tqdm(test_loader)):
        if torch.cuda.is_available():
            X = X.to(torch.device('cuda'))
        curr_pred = sn(X)
        curr_pred_len = curr_pred.shape[0]
        curr_predictions[curr_start_idx:curr_start_idx+curr_pred_len,:] = curr_pred.detach().cpu().numpy()
        curr_ground_truth[curr_start_idx:curr_start_idx+curr_pred_len,:] = Y.detach().cpu().numpy()
        curr_start_idx += curr_pred_len
    corrs.append(pearsonr(curr_ground_truth.flatten(), curr_predictions.flatten())[0])
```
