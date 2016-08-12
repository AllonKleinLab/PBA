Population balance analysis
=======================

Population balance analysis (PBA) is a method for reconstructing dynamics from static-snapshots of single-cell gene expression. The inputs to PBA are:

1. A sample of single-cell gene expression profiles (X)
2. A vector (R) with estimates of the relative rates of proliferation and loss at each sampled cell
3. A scalar diffusion constant (D) that reflects the level of stochasticity in the dynamics*

Under mild assumptions, there is a unique dynamical system that could have produced these data in steady-state. PBA uncovers the dynamical system by applying the law of population balance, which simply states that the flux of cells into and out of a small region of gene expression space must balance. More specifically, PBA models gene expression dynamics as a diffusion-drift process with diffusion rate (D) and inhomogenous boundary conditions (R). From the observed cell density (c; approximated by X), PBA calculates the unique potential field (V) satisfying the steady-state diffusion-drift equation:

<p align="center">
<img src="https://github.com/AllonKleinLab/PBA/blob/master/diff_drift_eq.png" width=270 />
</p>


The final outputs of PBA are:

1. A potential landscape V (encoded by the values it takes at each sampled cell)
2. Transition probabilities between the sampled cells (derived from V and D)
3. If terminal fates are given: the fate probabilities of each sampled cell

*_Note: The output of PBA does not change when if R and D are scaled by a common factor. So in practice D is redudant to R._ 

## Requirements ##

PBA requires the following Python modules.

2. SciPy
3. NumPy
4. sciki-learn (optional speed-up)

## Usage ##

### Input ###

The input data for PBA consists of a:

1. **Expression matrix (X)**. This *.npy file should contain a matrix of single-cell gene expression values. Rows represent cells and columns represent genes. X is used to generate a k-nearest neighbor (knn) graph adjacency matrix
2. **Edge list (A)** [alternative to 1.]. Instead of an expression matrix (input 1.), users can upload a list of edges representing a (knn) graph over sampled cells. The file should contain an edge in the format "i,j" on each line. Users can generate, visualize and then export knn graph edge lists in our companion software _SPRING_, available as a [webserver](https://kleintools.hms.harvard.edu/tools/spring.html) or a standalone program. 
2. **Source/sink vector (R)**. This *.npy file should contain a vector of source/sink terms representing the relative rates of proliferation and loss at each sampled cell. Note that uniformly changing R by a scalar factor f is equivalent to changing the diffusion rate (level of stochasticity) by 1/f.
3. **Lineage-specific sink matrix (S)** [optional]. If provided, this matrix can be used to define terminal lineages and compute the fate probabilities of sampled cells. This *.npy file should contain a matrix with one column for each cell and one row for each terminal lineage. S<sub>i,j</sub> represents the flux of cells into lineage i from cell j. 

We provide example data in `example_datasets/`. 

### Run CoMEt ###

We provide two pipelines for performing CoMEt:

1. **Run COMEt MCMC algorithm on real data and create output website**. Use the `run_comet_simple.py` script to run the Markov chain Monte Carlo (MCMC) algorithm on the given mutation matrix. `run_comet_simple.py` outputs a [JSON](http://json.org/) file that stores the parameters of the run, a tab-separated file that lists the collections identified by CoMEt (sorted descending by sampling frequency), and a website that can be used to visualize the results.
2. **Run CoMEt MCMC algorithm on real data, assess the significance against permuted data, and create output website**. Use the `run_comet_full.py` script to perform CoMEt with the same output as the `run_comet_simple.py` but with significant test. This pipeline computes the collections with statistical significance and identifies the consensus modules. The output of this pipeline contains a [JSON](http://json.org/) file that stores the parameters of the run, a tab-separated file that lists the collections identified by CoMEt (sorted descending by sampling frequency), and a website that can be used to visualize the results.

To view the results website, download the required Javascript files (see Requirements above) and start a Python web server:

        cd OUTPUT_DIRECTORY # the output directory you provided to run_comet_simple.py or run_comet_full.py
        bower install
        python -m SimpleHTTPServer 8000

  Then direct your browser to `http://localhost:8000`.

### Compute weights exhaustively ###

We also provide the script `run_exhaustive.py` as a simple way to compute the weight &phi;(M) for *all* gene sets *M* in a given dataset (using the same input format as above). The output of `run_exhaustive.py` is a tab-separated file that lists the weight &phi;(M) for all gene sets in the dataset (sorted ascending by &phi;(M)).

## Support ##

Please visit [our Google Group](https://groups.google.com/forum/#!forum/dendrix) to post questions and view discussions from other users, or contact us through [our research group's website](http://compbio.cs.brown.edu).

## Testing ##

To test CoMEt, run the following commands:

    cd test
    python test.py

The tests are successful if the last line of the text printed to the terminal is `"PASS"`.

## Reference ##

[Mark D.M. Leiserson](http://maxleiserson.com)\*, [Hsin-Ta Wu](http://cs.brown.edu/~bournewu/)\*, [Fabio Vandin](http://www.imada.sdu.dk/~vandinfa/), [Benjamin J. Raphael](http://compbio.cs.brown.edu). CoMEt: A Statistical Approach to Identify Combinations of Mutually Exclusive Alterations in Cancer. In *Proceedings of the 19th Annual Conference on Research in Computational Molecular Biology (RECOMB)* 2015. [Extended abstract](http://link.springer.com/chapter/10.1007%2F978-3-319-16706-0_19#page-1) and [preprint](http://arxiv.org/abs/1503.08224).
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
