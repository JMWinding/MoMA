# README

This git repository is the artifact of the Sigcomm 2023 paper "Towards Practical and Scalable Molecular Networks".

[![DOI](https://zenodo.org/badge/679340304.svg)](https://zenodo.org/badge/latestdoi/679340304)

Our code as well as data are both available, while documentation is still under construction.

## How to download this artifact?

Use `git clone` to download this repository.

## How to install this artifact?

This project is composed of the experiment part and the data processing part.

### Experiment

For experiment part, please prepare two Arduinos (one for TXs and one for RX) and install the Arduino IDE. The following two header files required to run the code, which can be found in the Arduino IDE library (In the Arduino IDE, click the `Tools > Manage Libraries...`).

* `Ticker.h`
* `OneWire.h`

And of course, you need to prepare pumps, tubes, salt (and/or baking soda) and the [EC reader](https://www.dfrobot.com/product-1123.html).

### Data processing

For data processing part, please install MATLAB R2019b or any later version. The following toolboxes are **_highly recommended_**.

* `Communication Toolbox`
    + This toolbox is used in `code_algo/get_gold_code.m` to generate Gold code of certain parameters, and `code_algo/get_gold_code2.m` to generate MoMA code (which is Gold code modulated with Manchester code). 
    + However, these two functions are only used in simulation for the early stage of the research. Once we go to experiment, we only need to generate the codes of required length and store in a TXT file before hand. Thus, during transmission the TX Arduino will pick the code for each TX from this file. Since MoMA assumes that the RX knows the pre-assigned code to each TX, the decoding will proceed depending on this assignment. 

* `Parallel Computing Toolbox`
    + This toolbox is used in `loop_emulates_txrx_all.m` and `loop_emulates_txrx_noncoherent.m` to parallel the decoding of multiple traces, thus accelerate the whole result generating procedure. Note that this toolbox is NOT used to speed up the decoding of ONE trace. Instead, each CPU core is used to process different traces concurrently. 
    + You can get rid of this toolbox by making the following modifications to the code, **_which will definitely and extremely increase the time to generate the results_**.
        - In `main_emulates_txrx_allloop.m`, comment the following code from line 5 to 7.
        ```
        if isempty(gcp("nocreate"))
            parpool(feature("numcores"));
        end
        ```
        This code creates the maximum number of workers, which is the same as the number of the CPU cores. But you can also specify a suitable number you want, e.g. 8 workers with the following code
        ```
        if isempty(gcp("nocreate"))
            parpool(8);
        end
        ```
        - In `loop_emulates_txrx_all.m` and `loop_emulates_txrx_noncoherent.m`, change the code in line 108 from
        ```
        parfor kk = 1:totalruns
        ```
        to
        ```
        for kk = 1:totalruns
        ```

## How to run this artifact?

The experiment is conducted on a testbed controlled by two Arduinos. If you would like to set up the testbed yourself, please refer to [Setting up the testbed](/documentation/testbed.md). We also provide our experimental data in this repository, which can be directly applied for the data processing part.

The data processing codes are MATLAB scripts and functions.

## How to compare this artifact’s outputs to outputs described in the paper?

We provide two MATLAB scripts that can directly generate the figures in the paper.

[`main_emulates_txrx_allloop.m`](/main_emulates_txrx_allloop.m) processes the experiment data and save the results as `.mat` files in corresponding folders. The script contains multiple loops, each defining a variable that specifies either the experiment data or the algorithm. Since not all combinations are reported in the paper, we add multiple `if` conditions and mark in comment their corresponding figure indices. Please feel free to modify the `if` conditions if you would like to see other results. The script then calls [`loop_emulates_txrx_all.m`](/loop_emulates_txrx_all.m) and [`loop_emulates_txrx_noncoherent.m`](/loop_emulates_txrx_noncoherent.m) to process the specified data in the specified way.

[`code_figure/generate_all_fig.m`](/code_figure/generate_all_fig.m) reads the results in `.mat` files and generate the figures. It calls multiple scripts, each generating one figure, and corresponding indices in the paper are marked in comment. And you can run each of them separately as well.

Simply running these two scripts in sequence will generate all the reported figures, but `main_emulates_txrx_allloop.m` may take a long time.

## How to reuse this artifact?

[`main_emulate_txrx.m`](/main_emulate_txrx.m) is a cleaner version of `main_emulates_txrx_allloop.m`, which is used to debug the algorithms with a single experiment trace.

For more details about the scripts, please refer to [Processing the data](/documentation/data_process.md).
