Population balance analysis
=======================

Population balance analysis (PBA) is a method for reconstructing dynamics from static-snapshots of single-cell gene expression. The inputs to PBA are:

1. A sample of single-cell gene expression profiles
2. A vector (_R_) with estimates of the relative rates of proliferation and loss at each sampled gene expression state
3. A scalar diffusion constant (_D_) that reflects the level of stochasticity in the dynamics*

Under mild assumptions, there is a unique dynamical system that could have produced these data in steady-state. PBA uncovers the dynamical system by applying the law of population balance, which simply states that the flux of cells into and out of a small region of gene expression space must balance. More specifically, PBA models gene expression dynamics as a diffusion-drift process with diffusion rate (D) and inhomogenous boundary conditions (_R_). From the observed cell density (_c_), PBA calculates the unique potential field (_V_) satisfying the steady-state diffusion-drift equation:

<p align="center">
<img src="https://github.com/AllonKleinLab/PBA/blob/master/diff_drift_eq.png" width=270 />
</p>


The final outputs of PBA are:

1. A potential landscape _V_ (encoded by the values it takes at each sampled gene expession state)
2. Transition probabilities between the sampled gene expession states (derived from _V_ and _D_)
3. If terminal fates are given: the fate probabilities of each sampled gene expession state

*_Note: The output of PBA does not change when if R and D are scaled by a common factor. So in practice D is redudant to R._ 

## Usage ##

### Input ###

Input1. **Expression matrix**. This (.npy or .csv) file should contain a matrix of single-cell gene expression values and is used to generate a k-nearest neighbor (knn) graph adjacency matrix. Rows represent cells and columns represent genes. 
Input2. **Edge list** [alternative to 1.]. Instead of an expression matrix (input 1.), users can upload a list of edges representing a (knn) graph over sampled gene expession states. The file should contain an edge in the format "_i,j_" on each line. Users can generate, visualize and then export knn graph edge lists in our companion software _SPRING_, available as a [webserver](https://kleintools.hms.harvard.edu/tools/spring.html) or a [standalone program](https://github.com/AllonKleinLab/SPRING/). 
Input3. **Source/sink vector**. This (.npy or .csv) file should contain a vector of source/sink terms representing the relative rates of proliferation and loss at each sampled gene expession state. Note that uniformly changing R by a scalar factor _f_ is equivalent to changing the diffusion rate (level of stochasticity) by _1/f_.
Input4. **Lineage-specific sink matrix** [optional]. If provided, this matrix can be used to define terminal lineages and compute the fate probabilities of sampled gene expession state. This (.npy or .csv) file should contain a matrix with one column for each cell and one row for each terminal lineage. The _i,j_ entry represents the flux of cells into lineage _i_ from gene expression state _j_. 

We provide example data in `example_datasets/`. 

### Run PBA ###

PBA applies a sequence of transformations to the data (see below). Each transformation can be run as a separate script, but relies on the output from the previous script. To execute all the steps at once, run `PBA_pipeline.py.` All argument flags listed below can be used with PBA_pipeline.py. 

1. **Compute a knn graph from an expression matrix**. Use `compute_knn_graph.py`. This script, together with graph visualizaion tools, is also avaiable in our companion software [_SPRING_](https://github.com/AllonKleinLab/SPRING/tree/master)


        Input: blah
        Output: blah
        Usage: blah blah

  Then direct your browser to `http://localhost:8000`.



        Input: blah
        Output: blah
        Usage: blah blah

  Then direct your browser to `http://localhost:8000`.


       hello

## Testing ##

To test CoMEt, run the following commands:

    cd test
    python test.py

The tests are successful if the last line of the text printed to the terminal is `"PASS"`.

## Reference ##

_(in submission)_

When using `runHotNet2.py`, you may also optionally provide any or all of the parameters listed
below. If one of these parameters is not provided, it will be set to the default value shown below.

        =============================================================================================================
        | PARAMETER NAME          | DEFAULT            | DESCRIPTION                                                |
        =============================================================================================================
        |-r/--runname             | None               |Name of run / disease.                                      |
        -------------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size       | 2                  |Minimum size connected components that should be returned.  |
        -------------------------------------------------------------------------------------------------------------
        |-c/--num_cores           | 1                  |Number of cores to use for running permutation tests in     |
        |                         |                    |parallel. If -1, all available cores will be used.          |
        -------------------------------------------------------------------------------------------------------------
        |-dp/--delta_permutations | 100                |Number of permutations to be used for delta parameter       |
        |                         |                    |selection.                                                  |
        -------------------------------------------------------------------------------------------------------------
        |-sp                      | 100                |Number of permutations to be used for statistical           |
        |--significance_permutatio|                    |significance testing.                                       |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_directory    | hotnet_output      |Output directory. Files results.json, components.txt, and   |
        |                         |                    |significance.txt will be generated in subdirectories for    |
        |                         |                    |each delta.                                                 |
        -------------------------------------------------------------------------------------------------------------

For simple runs on classic HotNet, use the `runClassicHotNet.py` Python script.  The following
parameters are required:

        ========================================================================================
        | PARAMETER NAME          | DESCRIPTION                                                |
        ========================================================================================
        |-mf/--infmat_file        |Path to .mat file containing influence matrix               |
        ----------------------------------------------------------------------------------------
        |-if/--infmat_index_file  |Path to tab-separated file containing an index in the first |
        |                         |column and the name of the gene represented at that index   |
        |                         |in the second column of each line.                          |
        ----------------------------------------------------------------------------------------
        |-hf/--heat_file          |Path to heat file containing gene names and scores. This    |
        |                         |can either be a JSON file created by generateHeat.py        |
        |                         |(described below), in which case the file name must end in  |
        |                         |.json, or a  tab-separated file containing a gene name in   |
        |                         |the first column and the heat score for that gene in the    |
        |                         |second  column of each line.                                |
        ----------------------------------------------------------------------------------------
        |-pnp                     |Path to influence matrices for permuted networks. Include   |
        |--permuted_networks_path |##NUM## in the path to be replaced with the iteration       |
        |                         |number                                                      |
        ----------------------------------------------------------------------------------------

Running with only the parameters specified above will create a 'hotnet_output' directory in your
current working directory that contains 5 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
classic HotNet algorithm.  The contents of the directories are identical to those described above
for simple runs of HotNet2 algorithm using `runHotNet2.py`. Similarily, the `runClassicHotNet.py`
script accepts the same optional parameters as the `runHotNet2.py` script.
