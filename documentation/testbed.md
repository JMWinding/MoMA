# Setting up the testbed

## Channel

Please set up arbitrary molecular channels for your own purpose. Here, we provide the two topologies used in our paper, i.e. the line and the fork.

![testbed_line](/documentation/figures/testbed_line.jpg "Line channel testbed.")

![topo_line](/documentation/figures/topo_line.jpg "Line channel topology.")

![testbed_fork](/documentation/figures/testbed_fork.jpg "Fork channel testbed.")

![topo_fork](/documentation/figures/topo_fork.jpg "Fork channel topology.")

Three different types of containers are prepared:

* Pure water tank: the background flow pump draws water from this container.
* Saline solution tank: the TX pumps draws information solution from this tank (or tanks).
* Waste tank: solution goes to this tank after passing by the RX EC meter.

## Circuit diagram

Here is a overview of the testbed circuit.

![circuits](/documentation/circuits/circuit.png "Arduino circuits.")

### Transmitter(s)

Each pump serves as one TX. It is turned on and off under the control of a transistor (NTE373) whose base is connected to the output pin of the TX Arduino.

### Receiver

The RX Arduino setting is much more easier, which only requires connecting the [EC meter](https://www.dfrobot.com/product-1123.html) and the SD card. The RX Arduino runs the script under folder [`Arduino/Rx`](/Arduino/Rx).

### Synchronizing TX and RX

The TX Arduino and the RX Arduino are synchronized through "Serial1", which requires connection of TX Arduino Pin18 to RX Arduino Pin19, and TX Arduino Pin19 to RX Arduino Pin18.

## Coding

We tested several different coding, each running a different Arduino code and requiring a separate `.txt` file in the TX Arduino SD card. The example scripts are for Gold code (length 7), Goldman code (length 14), and OOC code (length 14). Following are the explanations for each scheme.

### Goldman (MoMA CDMA)

This CDMA scheme uses MoMA code (naming Goldman, which is the abbrevation for gold plus Manchester). For a given length of gold code, it modulate every chip with Manchester code, which doubles the code length but at the same time balance all the codes.

To send a bit `1`, the TX sends the code with OOK. To send a bit `0`, the TX sends the complementary code.

The preamble is constructed by repeating __each chip of the code__ $Lp$ times.

![goldman](/documentation/figures/goldman.jpg "Goldman scheme.")

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_goldman`, where `n` represents the number of active TXs. Besides, `goldman.txt` is required in the TX SD card.

### Goldman0

This CDMA scheme is almost identical to [Goldman](#goldman-moma-cdma), except that the TX sends __nothing__ for a bit `0`.

![goldman0](/documentation/figures/goldman0.jpg "Goldman0 scheme.")

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_goldman0`, where `n` represents the number of active TXs. Besides, `goldman.txt` is required in the TX SD card.

### OOC

This CDMA scheme uses OOC code [[1]](#1).

To send a bit `1`, the TX sends the code with OOK. To send a bit `0`, the TX sends the complementary code.

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_ooc`, where `n` represents the number of active TXs. Besides, `ooc.txt` is required in the TX SD card.

![ooc](/documentation/figures/ooc.jpg "OOC scheme.")

### OOC0

This CDMA scheme is almost identical to [OOC](#ooc), except that the TX sends __nothing__ for a bit `0`.

![ooc0](/documentation/figures/ooc0.jpg "OOC0 scheme.")

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_ooc0`, where `n` represents the number of active TXs. Besides, `ooc.txt` is required in the TX SD card.

### Gold

This CDMA scheme uses Gold code. However, only balanced codes (with the number of `1`s and `0`s differing by at most one) are used.

The preamble is constructed by repeating __each chip of the code__ $Lp$ times.

![gold](/documentation/figures/gold.jpg "Gold scheme.")

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_gold`, where `n` repres

### Plain0

This is __NOT__ a CDMA scheme. Instead, it is only the basic OOK without multiplexing purpose. To send a bit `1`, the TX sends one pulse in the beginning of the symbol. To send a bit `0`, the TX sends nothing.

However, each TX is still assigned a Goldman code to construct preamble.

![plain0](/documentation/figures/plain0.jpg "Plain0 scheme.")

For this scheme, the TX Arduino runs the scripts in folders naming `Arduino/[n]Tx_plain0`, where `n` represents the number of active TXs. Besides, `goldman.txt` is required in the TX SD card.

## Data collection

The experiment is conducted as the following steps.

1. Check the connectivity of the circuit.
2. Refill the containers and turn on the power supplies.
3. Start the receiver code `RX.ino`, until the RX Arduino terminal shows `initialization done`.
4. Start the transmitter code. If everything is right, the TX Arduino terminal will show the information of loading pumps. Then, the TX and RX Arduino terminal will concurrently write new files for each experiment trace in SD cards.
5. When enough data is collected or any of the container depletes, please stop the testbed (by turning of the power supplies) and retrieve the data from SD cards.

## Data format

### TX data

Take [/dataset/data5/125ms_2-3-4-5_16_goldman/tx/00.TXT](/dataset/data5/125ms_2-3-4-5_16_goldman/tx/00.TXT) as an example of TX data.

#### Experiment setup

The first line shows the basic setup of the experiments.

```
125 100 4 7 16
```

* `125` is the chip interval in milliseconds. But coding `plain0` is an exception, where the chip interval is double the value (e.g. `875ms_2_2_plain0` has 1.75s symbol interval).
* `100` is the number of data bits in each packet.
* `4` is the number of active TX.
* `7` is the length of the gold code. It is either used directly as CDMA or to generate the MoMA codebook, depending on the coding scheme explained [in this part](#coding).
* `16` is the length of the preamble as multiples of the code length.

#### Transmitter setup

The next several lines, which is equal to the number of active TX, records the information of each TX.

```
(2, 20, 0, 10.00), (code 0)
(3, 258, 0, 10.00), (code 4)
(4, 496, 3, 10.00), (code 3)
(5, 734, 5, 10.00), (code 1)
```

For example, in `(2, 20, 0, 10.00), (code 0)`,

* `2` denotes that this pump is controled by Arduino pin 2. Note that the TX Arduino pin 2 to 5 actually match the TX 1 to 4 in both [line and fork channels](#channel), because pin 1 __cannot__ serve as the pump control.
* `20` is the large TX packet offset in the unit of chips. 
* `0` is the small TX packet offset in the unit of chips. The sum of the two offsets is the actual TX packet offset. Dividing it into two parts is for the purpose of studying the influence of transmission delay across molecules, but this part is currently in the discussion section of the paper and no results reported.
* `10.00` was reserved for the concentration of molecular solution. However, it is a useless variable for now, since all solutions of the same molecule are of equal concentration.
* `code 0` denotes that this TX uses the first code in the coding file, which is explained [here](#coding). The index starts from 0.

#### Transmitter behavior

The rest of the file after the line `START` records the behavior of each TX, and there are three different lines.

```
38332, 0, prem 0
80340, 0, bit 0
38332, 0, 0
```

* The first value is the timestamp in milliseconds.
* The second value denotes the index of the TX. `0` refers to the first pump recorded in the file, which is `(2, 20, 0, 10.00), (code 0)`.
* The third value makes a difference for these three lines: `prem 0` (`prem 1`) marks the transmission of bit 0 (1) in preamble; `bit 0` (`bit 1`) marks the transmission of bit 0 (1) in data; a single digit `0` (`1`) marks the transmission of chip 0 (1). Typically, `prem` and `bit` are only conceptual representations showing the state transitions, while they are always followed by a sequence of `0`s and `1`s showing the actual behavior of transmitters.

### RX data

Take [/dataset/data5/125ms_2-3-4-5_16_goldman/rx/00.TXT](/dataset/data5/125ms_2-3-4-5_16_goldman/rx/00.TXT) as an example of RX data.

```
24087	21
24730	22
24735	21
24740	20
24745	19
24750	20
24755	19
24760	18
24765	18
24770	18
24775	18
...
```

RX data are composed of two columns:

* The first column is the timestamp of RX in milliseconds. Since the TX Arduino and the RX Arduino are synchronized, the first line of RX data is recorded exactly as the start of the first TX packet. Because the TX code and the RX code start at different times, the absolute value are different but the offset remains constant.
* The second column is the reading of the analog pin A1, i.e. the EC reader. Note that the actual value of EC requires a linear conversion from the measurement. But since such conversion is LINEAR, it does not really matter we use either the original reading or the actual EC.


## Our dataset

Our dataset can be found in folder [dataset](/dataset).

```
├── dataset
    ├── data1
    ├── data2
    ├── data3
    ├── data4
    └── data5
        ├─125ms_2-3-5_16_goldman
        ├─125ms_3-4_16_goldman
        ├─125ms_2-3_16_goldman
        ├─125ms_2-3-4-5_16_goldman
        └─125ms_2-3-4_16_goldman
            ├─rx
            │   ├─00.TXT               
            │   └─...              
            └─tx    
                ├─00.TXT           
                └─...
```

### Second layer

In the second layer, data is organized in five folders by the correpsonding testbed and type of molecules. The details are described in the following table. By default, none of the transmitters share the same CDMA code. However, `data2` is the only exception, where all transmitters use the same code.

| Folder | Channel |  Molecule | Comment                            |
|:------:|:-------:|:---------:|------------------------------------|
|  data1 |   line  |   $NaCl$  |                                    |
|  data2 |   line  |   $NaCl$  | Repeat usage of the same CDMA code |
|  data3 |   line  |   $NaCl$  |                                    |
|  data4 |   fork  | $NaHCO_3$ |                                    |
|  data5 |   fork  | $NaHCO_3$ |                                    |

### Third layer

In the third layer, data is organized by different test cases. The folders are named in the way as

```
[T]ms_[pumpstr]_[Lp]_[code]
```

Take `125ms_2-3-4-5_16_goldman` as an example.

| Variable | Value   | Meaning                                                  |
|----------|---------|----------------------------------------------------------|
| T        | 125     | Chip interval in milliseconds.                           |
| pumpstr  | 2-3-4-5 | Indices of active pumps (TXs) separated by '-'.          |
| Lp       | 16      | Length of preamble as multiples of code length.          |
| code     | goldman | Coding type. More details are explained [in this part](#coding). |

### Fourth layer

The transmitter and the receiver data are stored in folder `tx` and `rx` separately, while each trace are in the form of `.TXT` naming as a two-digit interger.

## Reference

[1] Chu, Wensong, and Charles J. Colbourn. "Optimal (n, 4, 2)-OOC of small orders." Discrete Mathematics 279.1-3 (2004): 163-172.