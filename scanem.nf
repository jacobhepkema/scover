params.name = "example"
save_dir = file("./output/$params.name/")

params.data = "data/smalldataset.tsv"
data = file(params.data)

params.celldata = "add"
if(params.celldata == "add"){
    exit 1, "--celldata is a required argument"
}

celldata = file(params.celldata)

// Delete param & cals folders if they already exist                            
if( save_dir.exists() ) {                                                       
    println "Deleting existing save_dir"                                        
    del_result_save = save_dir.deleteDir()                                      
    println del_result_save ? "Deleted" : "Cannot delete!"                      
} 

// What do we align to?
params.tomtom = "resources/Mus_musculus.meme"
tomtom = file(params.tomtom)

// Threshhold q-value for TOMTOM alignment
params.tom_thresh = 0.05

// Run parameters
params.epochs = 24
params.num_candidates = 10
params.num_calibrations = 30
params.motif_length = 12
params.motif_amount = 300
params.val_factor = 10
params.batch_size = 128
params.seed = 42
params.opt = "SGD"

// Ranges for hyperparameter search, see README to see what these do
params.sigma_motifs_min = 1e-7
params.sigma_motifs_max = 1e-3
params.sigma_net_min = 1e-5
params.sigma_net_max = 1e-2
params.epsilon_min = 5e-4
params.epsilon_max = 5e-2


println """\
=========================================================================
=========================================================================

  scanem  v0.2 (Jun 9 2020) 

=========================================================================

  run name             : ${params.name}
  data path            : ${params.data}
  cell label data path : ${params.celldata}
  motif length         : ${params.motif_length}
  amount of motifs     : ${params.motif_amount}
  epochs               : ${params.epochs}
  batch size           : ${params.batch_size}
  K in K-fold CV       : ${params.val_factor}
  number of cal        : ${params.num_calibrations} 
  number of candidates : ${params.num_candidates}
  tomtom db file       : ${params.tomtom}  
  random seed          : ${params.seed}     

=========================================================================
=========================================================================

         """
         .stripIndent()

process scanem {
    publishDir "$save_dir", pattern: "*.pt", mode: "link"  
    publishDir "$save_dir", pattern: "AllMotifs/*.png", mode: "link"
    publishDir "$save_dir", pattern: "Motifs/*", mode: "link"  
    publishDir "$save_dir", pattern: "*png", mode: "link"
    publishDir "$save_dir", pattern: "All_MEME_*", mode: "link"
    publishDir "$save_dir", pattern: "*.p", mode: "link"
    publishDir "$save_dir", pattern: "Best_*", mode: "link"
    publishDir "$save_dir", pattern: "*.tsv.gz", mode: "link"
    publishDir "$save_dir", pattern: "*for_best_candidate.txt", mode: "link"
    publishDir "$save_dir", pattern: "All_model_HDF5_*", mode: "link"

    output:
    file "*.tsv.gz" 
    file "*.p" optional true
    file "Best_*" optional true
    file "*.pt" 
    file "AllMotifs/*.png"
    file "Motifs/*"
    file "*.png"
    file "All_MEME_*" into all_motifs_tmp_ch 
    file "*metrics_for_best_candidate.txt"
    file "All_model_HDF5_*"

    script:
    """
    scanem.py \
    $data \
    $celldata \
    -name $params.name \
    -epochs $params.epochs \
    -num_candidates $params.num_candidates \
    -c $params.num_calibrations \
    -v $params.val_factor \
    -batch_size $params.batch_size \
    -m $params.motif_length \
    -d $params.motif_amount \
    -seed $params.seed \
    -opt $params.opt \
    -sigma_motifs_min $params.sigma_motifs_min \
    -sigma_motifs_max $params.sigma_motifs_max \
    -sigma_net_min $params.sigma_net_min \
    -sigma_net_max $params.sigma_net_max \
    -epsilon_min $params.epsilon_min \
    -epsilon_max $params.epsilon_max
   """
}

// Feed into two different tomtom alignments 
all_motifs_tmp_ch
    .into { all_motifs_ch; all_motifs_tom_ch }

// Align to database
process tomtom {
    publishDir "$save_dir", mode: "link"

    input:
    file meme from all_motifs_tom_ch 

    output:
    file "all_tomtom/tomtom.tsv" into all_motifs_tom_tsv_ch
    file "all_tomtom/tomtom.html" optional true
    file "all_tomtom/tomtom.xml" optional true

    script:
    """
    /opt/bin/tomtom \
    -thresh $params.tom_thresh \
    -o all_tomtom \
    $meme $tomtom 
    """
}

// Align motifs to itself, output text only to tsv file
process tomtom_allmotifs {
    publishDir "$save_dir", mode: "link"
        
    input:
    file all_meme from all_motifs_ch
    
    output:
    file "motifmotif_alignment.tsv" into all_motifs_tsv_ch
    
    script:
    """
    /opt/bin/tomtom -norc -text \
    -thresh $params.tom_thresh \
    -o all_motifs \
    $all_meme $all_meme \
    > motifmotif_alignment.tsv
    """
}

// Run downstream analysis for network outputs
process motif_analysis {
    publishDir "$save_dir", mode: "link"

    input:
    file tom_out from all_motifs_tom_tsv_ch
    file mot_out from all_motifs_tsv_ch

    output:
    file "*.png"
    file "*.tsv"
    file "*.html" 
    file "*.RDS"

    script:
    """
    motif_graph_analysis.r \
    -d $tomtom \
    -a $tom_out \
    -m $mot_out \
    -t $params.tom_thresh 
    """
}

