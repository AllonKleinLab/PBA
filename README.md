Population balance analysis
=======================

Population balance analysis (PBA) relates the observed states of a system to its steady-state dynamics by applying the law of population balance (ref 1). The dynamics are represented as a potential field, _V_, over the space of gene expression. Knowledge of V allows one to calculate possible dynamical trajectories and their properties. PBA requires enough cells to be sampled to approximate a continuum of gene expression states; and some prior estimates of the relative rates of proliferation and loss in different gene expression states (_R_). The potential _V_ is the solution to the partial differential equation,

<p align="center">
<img src="https://github.com/AllonKleinLab/PBA/blob/master/aux_files/diff_drift_eq.png" width=170 />
</p>

where _c_ is the density of cells in gene expression space, and _R_ is the net rate of cell division minus cell loss. In practice, this equation cannot be numerically solved, but for a set of observed transcriptional states _x<ub>i</sub>_ sampled from _c_, PBA calculates an approximate solution for _V_. We prove mathematically in (ref 1) that this approximation converges in the high sampling limit. 

## Usage ##

### Inputs ###

1. **Expression matrix**. This (.npy or .csv) file should contain a matrix of single-cell gene expression values and is used to generate a k-nearest neighbor (knn) graph adjacency matrix. Rows represent cells and columns represent genes. 
2. **Edge list** [alternative to 1.]. Instead of an expression matrix (input 1.), users can upload a list of edges representing a (knn) graph over sampled gene expession states. The file should contain an edge in the format "_i,j_" on each line (0-based numbering). Users can generate, visualize and then export knn graph edge lists in our companion software _SPRING_ (ref 2), available as a [webserver](https://kleintools.hms.harvard.edu/tools/spring.html) or a [standalone program](https://github.com/AllonKleinLab/SPRING/). 
3. **Source/sink vector**. This (.npy or .csv) file should contain a vector of source/sink terms representing the relative rates of proliferation and loss at each sampled gene expession state. Note that uniformly changing R by a scalar factor _f_ is equivalent to changing the diffusion rate (level of stochasticity) by _1/f_.
4. **Lineage-specific sink matrix** [optional]. If provided, this matrix can be used to define terminal lineages and compute the fate probabilities of sampled gene expession state. This (.npy or .csv) file should contain a matrix with one column for each lineage and one row for each cell. The _i,j_ entry represents the flux of cells from gene expression state _i_ into lineage _j_. 
5. **Diffusion constant (_D_)** [defaut: 1.0]. This parameter controls the level of stochasticity in the model. It is redundant: setting _D = x_ is equivalent to multiplying inputs 3., 4. by the factor _1/x_.

We provide example datasets from (ref 1) in `example_datasets/`. 

### Run PBA ###

PBA applies a sequence of calculations (see below). Each one can be run as a separate script, which relies on the output from the previous script. To execute all the steps at once, run `PBA_pipeline.py` as follows: 

        Inputs: Expression matrix (Input 1. above) or edge list (Input 2. above), 
                A gobal source/sink vector (Input 3. above) 
                And optionally a lineage specific sink matrix (Input 4. above)
        Output: If no edge list was inputted: an edge list (a file called "edge_list.csv" saved to the same directory as the expression matrix)
                The pseudoinverse of the knn graph Laplacian (a matrix called Linv.npy saved to the same directory as the edge list or expression matrix)
                The potential (an array called V.npy saved to the same directory as the edge list or expression matrix)
                If lineage specific fates were inputted: a fate probability matrix where rows are cells and columns are fates (an array called B.npy saved to the same diectory as the edge list or expression matrix)

        Usage: python PBA_pipeline.py -X <path_to_expression_matrix>            (required if no edge list is supplied; .csv or .npy)
                                      -E <minimum_mean expression>              (default = 0.05; used to filter genes for knn graph)
                                      -V <minimum_CV>                           (default = 2; used to filter genes for knn graph)
                                      -N <Normalize>                            (default = False; used to normalize expression data for knn graph)
                                      -p <PCA dimension>                        (default = 50; used to compute distance matrix for knn graph)
                                      -k <number of nearest neighbors>          (default = 10; used to compute edge list for knn graph)
                                      -e <path_to_edge_list>                    (required if no expression matrix is supplied)
                                      -R <path_to_sources_sinks_vector>         (required; .npy or .csv)
                                      -S <path_to_lineage_specific_sink_matrix> (optional, needed to compute fate probabilities; .csv or .npy)
                                      -D                                        (default = 1.0; controls the level of stochasticity in the model)


Alternatively, it is possible to run each step separately, as follows: 

_(Note: Some of the following functions differ from the pseudocode of (ref 1) in that they take as input an edge list as opposed to an expression matrix. This configuration improves computational efficiency and allows generation of edge lists from outside software.)_ 

1. **Compute a knn graph from an expression matrix**. Use `compute_knn_graph.py`. This script, together with graph visualizaion tools, is also available in our companion software [_SPRING_ (ref 2)](https://github.com/AllonKleinLab/SPRING/tree/master)


        Inputs: Expression matrix (see Input 1. above)
        Output: Edge list (a file called "edge_list.csv" saved to the same directory as the expression matrix)

        Usage: python compute_knn_graph.py -X <path_to_expression_matrix>   (required; .csv or .npy)
                                           -E <minimum_mean expression>     (default = 0.05; used to filter genes)
                                           -V <minimum_CV>                  (default = 2; used to filter genes)
                                           -N <Normalize>                   (default = False; used to normalize expression data for knn graph)
                                           -p <PCA dimension>               (default = 50; used to compute distance matrix)
                                           -k <number of nearest neighbors> (default = 10; used to compute edge list)

2. **Compute the pseudo-inverse Laplacian**. Use `compute_Linv.py`. 


        Inputs: Edge list (see Input 2. above)
        Output: Pseudoinverse of the graph Laplacian (a matrix called Linv.npy saved to the same directory as the edge list)

        Usage: python compute_Linv.py -e <path_to_edge_list> (required)

3. **Compute the potential (V)**. Use `compute_potential.py`. 


        Inputs: Source/sink vector (see Input 3. abobe) and the pseudoinverse of the graph Laplacian 
        Output: The potential (an array called V.npy saved to the same directory as the pseudoinverse Laplacian)

        Usage: python compute_potential.py  -R <path_to_sources_sinks_vector>    (required; .npy .csv)
                                            -L <path_to_pseudoinverse_laplacian> (default: "Linv.npy" in the same directory as R; .npy)
        
4. **Compute fate probabilities**. Use `compute_fate_probabilities.py`. 


        Inputs: Lineage-specific sink matrix (see Input 4. above), the potential function V and a knn graph edge list (see Input 2. above)
        Output: Fate probability matrix where rows are cells and columns are fates (an array called B.npy saved to the same diectory as the potential V)

        Usage: python compute_fate_probabilities.py -S <path_to_lineage_specific_sink_matrix> (required; .csv or .npy)
                                                    -V <path_to_potential_vector>             (default: "V.npy" in same directory as S; .npy or .csv)
                                                    -e <path_to_edge_list>                    (default: "edge_list.csv" in same directory as S)
                                                    -D <diffusion_constant>                   (default = 1.0; controls the level of stochasticity in the model)\n

5. **Compute mean first passage times**. Use `compute_mean_first_passage_times.py`. 


        Inputs: The potential function V and a knn graph edge list (see Input 2. above)
        Output: Mean first passage time matrix T where T[i,j] = MFPT(i -> j) (an array called T.npy saved to the same diectory as the potential V)

        Usage: python compute_fate_probabilities.py -V <path_to_potential_vector>             (required; .npy or .csv)
                                                    -e <path_to_edge_list>                    (default: "edge_list.csv" in same directory as V)
                                                    -D <diffusion_constant>                   (default = 1.0; controls the level of stochasticity in the model)\n


## Testing ##


0. Unzip the example_datasets folder (can be done in terminal with the following command)



		unzip example_datasets.zip

1. To reproduce the PBA results used for Figure 1 of (ref 1), run the following command
 


        python PBA_pipeline.py -X example_datasets/bifurcation/X.npy -R example_datasets/bifurcation/R.npy -S example_datasets/bifurcation/S.npy

2. To reproduce the PBA results used for Figures 2-4 of (ref 1), run the following command



        python PBA_pipeline.py -X example_datasets/hematopoiesis/X.npy -R example_datasets/hematopoiesis/R.npy -S   example_datasets/hematopoiesis/S.npy
 
To visualize the results, run `python view_testing_results.py`

## Support ##

Contact calebsw@gmail.com

## References ##

Ref 1 (https://www.pnas.org/content/115/10/E2467)

Ref 2 (https://www.nature.com/articles/nature25741)

