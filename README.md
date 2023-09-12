# README

This git repository is the artifact of the Sigcomm 2023 paper "Towards Practical and Scalable Molecular Networks".

[![DOI](https://zenodo.org/badge/679340304.svg)](https://zenodo.org/badge/latestdoi/679340304)

Our code as well as data are both available, while documentation is still under construction.

## How to download this artifact?

Use `git clone` to download this repository.

## How to install this artifact?

This project is composed of the experiment part and the data processing part.

For experiment part, please prepare two Arduinos (one for TXs and one for RX) and install the Arduino IDE. The required `Ticker.h` and `OneWire.h` can be downloaded from the Arduino IDE library.

For data processing part, please install MATLAB R2019b or any later version. Please install any toolbox that is reported missing when you run the scripts.

## How to run this artifact?

The experiment is conducted on a testbed controlled by two Arduinos. If you would like to set up the testbed yourself, please refer to [Setting up the testbed](/documentation/testbed.md). We also provide our experimental data in this repository, which can be directly applied for the data processing part.

The data processing codes are MATLAB scripts and functions.

## How to compare this artifact’s outputs to outputs described in the paper?

We provide two MATLAB scripts that can directly generate the figures in the paper.

[`main_emulates_txrx_allloop.m`](/main_emulates_txrx_allloop.m) processes the experiment data and save the results as `.mat` files in corresponding folders. The script contains multiple loops, each defining a variable that specifies either the experiment data or the algorithm. Since not all combinations are reported in the paper, we add multiple `if` conditions and mark in comment their corresponding figure indices. Please feel free to modify the `if` conditions if you would like to see other results.

[`code_figure/generate_all_fig.m`](/code_figure/generate_all_fig.m) reads the results in `.mat` files and generate the figures. It calls multiple scripts, each generating one figure, and corresponding indices in the paper are marked in comment. And you can run each of them separately as well.

Simply running these two scripts in sequence will generate all the reported figures, but `main_emulates_txrx_allloop.m` may take a long time. 

For more details about  the scripts, please refer to [Processing the data](/documentation/data_process.md).
