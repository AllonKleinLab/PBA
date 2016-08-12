Population balance analysis
=======================

Population balance analysis (PBA) is a method for reconstructing dynamics from static-snapshots of single-cell gene expression. The inputs to PBA are:

1. A set of single-cell gene expression profiles (X)
2. Prior estimates (R) of the relative rates of proliferation and loss in different gene expression states
3. A diffusion constant (D) that reflects stochasticity in the dynamics*

Under mild assumptions, there is a unique dynamical system that could have produced these data in steady-state. PBA uncovers the dynamical system by applying the law of population balance, which simply states that the flux of cells into and out of a small region of gene expression space must balance.  

*Note that multiplying R and D by a common factor changes the time-scale but nothing else, so these parameters are redundant in practice. 

HotNet2 vs. Classic HotNet
------------------------
This distribution contains two related algorithms: the original HotNet algorithm "classic HotNet",
and an updated version "HotNet2". For details on the two algorithms, please refer to the publications
listed at the end of this README.

In brief, HotNet2 uses a new heat diffusion kernel analogous to random walk with restart that better
captures the local topology of the interaction network as compared to the general heat diffusion
process used by classic HotNet. HotNet2 also uses an asymmetric influence score and different
permutation testing and parameter selection procedures. Although classic HotNet is included for
completeness, we recommend using HotNet2.


Requirements
------------------------

* Linux/Unix
* [Python 2.7](http://python.org/)
* [NumPy 1.6.2](http://www.numpy.org/)
* [SciPy 0.10.1](http://www.scipy.org/)
* [NetworkX 1.7](http://networkx.github.io/)
* [h5py 2.4.0](http://www.h5py.org/)
* Fortran or C compiler (optional but recommended for performance)

HotNet2 will likely work with additional versions of Python, h5py, NetworkX, NumPy, and SciPy, but
alternative configurations have not been tested.

Support
------------------------
For support using HotNet, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).

Setup
------------------------

### Compilation
For best performance, install a Fortran or C complier and run one of the following commands
(or some appropriate variation of them) prior to running HotNet for the first time:

With a Fortran compiler:

    python hotnet2/setup_fortran.py build_src build_ext --inplace

With a C compiler:

    python hotnet2/setup_c.py build_src build_ext --inplace

If you are unable to perform these steps, the code will transparently fall back to a pure Python
implementation.

### Influence matrix creation

For each gene-gene interaction network you want to use with HotNet2, you must perform a one-time
step to generate the corresponding influence matrix. Use the provided `makeRequiredPPRFiles.py`
script to create the real and permuted personalized pagerank influence matrices. We have included
configuration files for creating the required matrices for the HPRD and iRefIndex PPI networks,
which you can run as follows:

    python makeRequiredPPRFiles.py @influence_matrices/hprd.config

    python makeRequiredPPRFiles.py @influence_matrices/irefindex.config

This will create the following files in the `influence_matrices/hprd` / `influence_matrices/irefindex`
directories:

* `{hprd/iref}_index_genes`: Gene-index file for the largest component in the given network.
* `{hprd/iref}_edge_list`: Edge list file for the largest component in the given network.
* `{hprd/iref}_ppr_{beta}.h5`: Personalized page-rank influence matrix in HDF5 format.
* `permuted`: directory containing 100 subdirectories with the above files for permuted matrices

For other networks, and for creating influence matrices for use with the classic HotNet algorithm,
see the "Advanced use" section below.

Note that this step will take a long time. Fortunately, though, you only need to do it once per interaction network you wish to use. We provide 1000 permuted networks for the HINT+HI2012,
iRefIndex 9, and Multinet interaction networks for download:

* [HINT+HI2012](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/hint+hi2012.tar) (~500Mb)
* [iRefIndex 9](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/iref9.tar) (~1.2Gb)
* [Multinet](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/multinet.tar) (~1.4Gb)

You can use these permuted networks instead of generating your own. However, you will
still need to generate the influence matrices for each of the permuted networks. See
the "Advanced Use" section below.

Simple runs
------------------------
Once you have performed the influence matrix creation step described above, you can use the
`runHotNet2.py` script to get started running HotNet2 quickly and easily. You must provide the
following parameters:

        ========================================================================================
        | PARAMETER NAME          | DESCRIPTION                                                |
        ========================================================================================
        |-mf/--infmat_file        |Path to HDF5 (.h5) file containing influence matrix. NumPy  |
        |                         | (.np) and MATLAB files also supported.                     |
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
current working directory that contains 4 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
HotNet2 algorithm. The output files are:

* `components.txt`: Lists subnetworks identified as significantly altered, one per line. Genes
  within each subnetwork are separated by tabs.
* `significance.txt`: For k from 2 - 10, lists the number of subnetworks of size >= k found in the
  real data, the expected number of subnetworks of size >= k based on permuted data, and the p-value
  for the observed number of subnetworks.
* `results.json`: Contains all of the above information plus the parameters used for the run in
  JSON format to faciliate further automated processing

The `runHotNet2.py` script can also be used to create a web visualization of the output subnetworks.
To do so, include the `--edge_file` parameter, and, optionally, other visualization-related parameters:

        =============================================================================================================
        | PARAMETER NAME          | DEFAULT   | DESCRIPTION                                                         |
        =============================================================================================================
        |-ef/--edge_file          | None               |Path to TSV file listing edges of the interaction network,  |
        |                         |                    |where each row contains the indices of two genes that are   |
        |                         |                    |connected in the network. This is used to create subnetwork |
        |                         |                    |visualizations; if not provided, visualizations will not be |
        |                         |                    |made.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-dsf                     | None               |Path to a tab-separated file containing a gene name in the  |
        |--display_score_file     |                    |first column and the display score for that gene in the     |
        |                         |                    |second column of each line.                                 |
        -------------------------------------------------------------------------------------------------------------
        |-nn/--network_name       | Network            |Display name for the interaction network.                   |
        -------------------------------------------------------------------------------------------------------------

This will result in a a `viz` subdirectory of the output directory. To view the visualizations,
navigate to the `viz` directory and run `python -m SimpleHTTPServer`, then visit `http://localhost:8000`
in a browser.

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
