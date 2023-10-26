# Processing the data

## Variables

As said, `main_emulates_txrx_allloop.m` processes the data in the way that are reported in the paper. However, there are multiple loops in the file that specifies which files to process and how to process them.

### Data related variables

| Variable  | Meaning                                                                                 | Potential value                                                                                   |
|-----------|-----------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
| `datanote0` | The single-molecule channel in experiment or the multi-molecule channel to be emulated. | n-by-1 string array. Only the combinations of all fork or line channels are meaningful.           |
| `nMo`       | The number of molecule in the experiment data or the emulated data.                     | Any positive integer when $datanote0$ is a single string; or equals to the length of $datanote0$. |
| `code`      | The code used by the TXs.                                                               | Check [dataset](/dataset) and [documentation](/documentation/testbed.md)                            |
| `T`         | The chip interval in milliseconds.                                                      | Check [dataset](/dataset) and [documentation](/documentation/testbed.md)                            |
| `pumpstr`   | The active TX pins on TX Arduino.                                                       | Check [dataset](/dataset) and [documentation](/documentation/testbed.md)                            |
| `Lp`        | The length of the preamble as multiples of CDMA code length.                            | Check [dataset](/dataset) and [documentation](/documentation/testbed.md)                            |

### Decoder

`main_emualtes_tx_allloop.m` calls two different scripts which actually process the raw data and output the results.

* [`loop_emulates_txrx_noncoherent.m`](/loop_emulates_txrx_noncoherent.m) uses the [noncoherent decoder](/code_algo/decode_mmo_noncoherent_MMoNTx.m) (constant thresholder), which serves as a baseline. The only result assumes that the RX knows the arriving time of each TX.
* [`loop_emulates_txrx_all.m`](/loop_emulates_txrx_all.m) uses the [proposed decoder](/code_algo/decode_mmo_coherent_MMoNTxSW11loop.m). All results in the paper (except one bar) are processed by this decoder.

### Decoder related variables

| Variable | Meaning                                                                                             | Potential value                                     |
|----------|-----------------------------------------------------------------------------------------------------|-----------------------------------------------------|
| `cenote`   | Preset adaptive filtering loss weights for channel estimation.                                      | Check [GetCEWeights()](/code_algo/GetCEWeights.m).   |
| `pdnote`   | Preset adaptive filtering loss weights for packet detection.                                        | Check [GetCEWeights()](/code_algo/GetCEWeights.m).   |
| `T2`       | $T/T2$ is the oversampling rate.                                                                    | No lower than the actual Arduino sampling interval. |
| `Lp2`      | Conceptual preamble length, which treats some data bits as preamble.                                | No lower than $Lp$.                                 |
| `algoPD`   | Packet detection algorithm.                                                                         | "gt", "gt1" or "sc".      |
| `algoCE`   | Channel estimation algorithm.                                                                       | "gt"  (only if `algoPD` is "gt") or "af0".     |
| `mode_pd`  | Packet detection only. The decoder stops as long as all TXs are detected, whether correct or wrong. | True (only if `algoPD` is "sc") or false.                                      |
| `debug_pd` | Debugging packet detection. It corrects false packet detection and records all metrics.             | True (only if `algoPD` is "sc") or false.            |

The `algoPD` has three options

* "gt" means ground truth, which assumes the RX knows exactly when each packet arrvies. The ground truth is based on the peak of the CIR of each TX.
* "gt1" also assumes the RX knows exactly when each packet arrvies, but with a few bits lag. This mode is meant to test the decoding performance after correct packet detection, because we observe that the packet detection is usually a few bits behind the peak of the CIR and such lag is larger at higher data rate (i.e. higher ISI).
* "sc" performs packet detection normally.

The `algoCE` has two options

* "gt" mean ground truth. But this is not the real ground truth, which is estimated from the whole packet. It is only legitimate when `algoPD` is "gt".
* "af0" performs channel estimation normally, where the CIR is estimated from the preamble only. It is legitimate with any `algoPD` option.

## Results format

If you run `main_emulates_txrx_allloop.m` yourself, the results can be found in folder [`mat_temp`](/mat_temp).

The results are saved in the same naming manner for directories and `.mat` files.

### First layer directory

```
mat[datanotes](PD(debug))
```

* `datanotes` correpsonds to the `datanote0` in [Data related variables](#data-related-variables), where the n-by-1 string array is converted to a string connected by dashes.
* `PD` exists only if `mode_pd` is true, and `PDdebug` exists only if `debug_pd` is true.

### Second layer directory

```
ce[cenote](-pd[pdnote])
```

* `cenote` corresponds to the `cenote` in [Algorithm related variables](#algorithm-related-variables).
* `(-pd[pdnote])` exists only if `pdnote` is different from `cenote`.

### .mat file

```
emulates_[T]ms_[pumpstr]_[Lp(-Lp2)]_[code]_[nMo]_[algoPD]-[algoCE].mat
```

* These variables exactly matches their description in [Data related variables](#data-related-variables) and  [Algorithm related variables](#algorithm-related-variables).
* `(-Lp2)` exists only if `Lp2` is different from `Lp`.

## MATLAB functions

Here we only list the general purpose of each MATLAB functions in [folder](/code_algo). For more details, please refer to the comments in each function.

* `MatlabIndex2BinarySeq.m` converts an integer to a binary sequence. It is used in Viterbi decoder to reduce the storage for Markov states of transmitted bit sequence.
* `ComputeViterbiLogProb.m` computes the probability of a Viterbi states based on its observation.
* `DecodeSequence_ViterbiForward4.m` performs the forward precedure in Viterbi algorithm, which computes the probability of each trellis.
* `DecodeSequence_ViterbiBackward.m` performs the backward procedure in Viterbi algorithm, which returns the bit sqeuence of highest probability.
* `decode_mmo_noncoherent_MMoNTx.m` decodes the packets with a simple thresholder.
* `decode_mmo_coherent_MMoNTxSW11loop.m` decodes the packet with Viterbi algorithm, which is proposed in the paper.
* `GeneratePreambleChips.m` generates the sequence of preamble chips based on the CDMA code and preamble bits. The sequence is of `1`s and `-1`s.
* `GeneratePreambleBits.m` generates the sequence of preamble bits. The sequence is of `1`s and `-1`s.
* `GenerateOversampleChips.m` upsample the transmitted chips by a specified rate. The original sequence is of `1`s and `-1`s, but it is padded with `0`s.
* `GenerateDataChips.m` generates the sequence of data chips based on the CDMA code and data bits. The sequence is of `1`s and `-1`s.
* `GenerateCodeChips.m` generates the CDMA code chips based on the basic code and code type. The sequence is of `1`s and `-1`s.
* `ToPos.m` convert a seqeuence of `1`s and `-1`s to a sequence of `1`s and `0`s.
* `GeneratePreambleChipsPos.m`  generates the sequence of preamble chips based on the CDMA code and preamble bits. The sequence is of `1`s and `0`s.
* `GenerateDataChipsPos.m` generates the sequence of data chips based on the CDMA code and data bits. The sequence is of `1`s and `0`s.
* `EstimateChannelMMoSW.m` estimate the channel of the specified TXs.
* `emulates_construct_rxIn.m` load the data from the dataset as the input of decode like `decode_mmo_coherent_MMoNTxSW11loop.m`.
* `GetCEWeights.m` generates the preset adaptive filtering loss weights for channel estimation.
* `get_gold_code2.m` generates MoMA code.
* `get_gold_code.m` generates Gold code.
* `get_decode_matrix.m` generates the matrix used for LS channel estimation. It is the same matrix that is also used in our adaptive filtering algorithm.
* `GeneratePreambleDetectionChips.m` generates the sequence for preamble detection. 
* `read_txrx.m` load the TX and RX data from the dataset.
* `PreambleDetectionSW8.m` detects the location of all packets in the current sliding window.
* `RebuildKnownPacket.m` rebuilds the received signal of TXs with their decoded chips and estimated CIR.
* `RemoveKnownPacket.m` removes the rebuilt signal from the whole RX signal.
* `sim_mc_cir.m` simulates a CIR based on 1-D model.
* `sim_mc_cir3.m` simulates a CIR based on 3-D model.
* `sim_tx.m` simulates the TX signal.
* `sim_mmo_tx.m` simulates the TX signal on multiple molecules.