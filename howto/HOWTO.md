# Pysca toolbox HOWTO

This is a walkthrough showing how to put pysca in action with the provided example tracesets. It does not go deep under the hood; if you like that, the best way is to dive into the code starting from the scripts used here.

## Environment setup

Pysca needs python 2.7 with numpy and matplotlib, and jupyter if you like to work with notebook format in a browser. Here is how to create a minimal environment for pysca with Anaconda python distribution in Linux and activate it:

    elbrus:pysca ilya$ conda create --name py27min python=2.7 numpy matplotlib jupyter
    [... all the conda printouts will be here ..]
    elbrus:pysca ilya$ source activate py27min
    (py27min) elbrus:pysca ilya$

Note that Anaconda provides numpy built against Intel MKL (Math Kernel Library), which is essential for performance. If you use other python dstribution, ensure that you use numpy-MKL.

## Traceset conversion
 
Convert an example set of power traces obtained from a software AES running on an ATmega microcontroller.

    (py27min) elbrus:pysca ilya$ python trs2npz.py traces\swaes_atmega_power
    Number of traces: 2000
    Samples per trace: 2800
    Samples datatype: int8
    Data bytes: 16
    Trace block size:
    Header size:
    Preallocating arrays
    Populating arrays
    Saving file
    Done

## CPA and LRA attacks on SW AES

For this example we will use the commad line script. Execute the script that performs the attacks, recovering a single key byte.

    (py27min) elbrus:pysca ilya$ attackaessbox.py
    ---
    Attack parameters
    Intermediate function   : sBoxOut
    CPA leakage function    : leakageModelHW
    LRA basis functions     : basisModelSingleBits
    Encryption              : True
    S-box number            : 2
    Known key               : 0x2b7e151628aed2a6abf7158809cf4f3c
    Known roundkey          : 0x2b7e151628aed2a6abf7158809cf4f3c
    ---
    Loading traces/swaes_atmega_power.npz
    Number of traces loaded : 100
    Trace length            : 700
    Loading time            : 0.03 s
    ---
    Attack
    ConditionalAverager: initialized for 256 values and trace length 700
    ---
    Results after 20 traces
    CPA
    Winning candidate: 0xd3, peak magnitude 0.827157
    Correct candidate: 0x15, peak magnitude 0.748879, rank 22
    LRA
    Winning candidate: 0x7d, peak magnitude 0.951581
    Correct candidate: 0x15, peak magnitude 0.884216, rank 50
    ---
    [...]
    Results after 100 traces
    CPA
    Winning candidate: 0x15, peak magnitude 0.481228
    Correct candidate: 0x15, peak magnitude 0.481228, rank 1
    LRA
    Winning candidate: 0x15, peak magnitude 0.512743
    Correct candidate: 0x15, peak magnitude 0.512743, rank 1
    ---
    Cumulative timing
    24.56 s
    ---
    Plotting...

Observe the result visualization. The plots show results of CPA (correlation traces) and LRA (R2 traces and matrix of basis fucntion coefficients characterising the leakage function) for the maximum amonut of traces, and evolution of the correct key candidate rank with the increasing amount of traces.

<img src="howto-script-aes-result.png" width="640">

The parameters of the attack can be adjusted in the configuration section of the script.

## CPA and LRA attack on HW DES

For this example, we will use another (convenient) way to work: a Jupyter noteboook in a browser. Launch the notebook server:

    (py27min) elbrus:pysca ilya$ jupyter notebook
    [I 13:15:38.823 NotebookApp] Serving notebooks from local directory: /Users/ilya/pysca
    [I 13:15:38.823 NotebookApp] 0 active kernels 
    [I 13:15:38.824 NotebookApp] The Jupyter Notebook is running at: http://localhost:8888/?token=a75f5aa53be646d4a96bedc760728c9baea3a3a72b9111be
    [I 13:15:38.824 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).

You will see the browser poping up with the directory contents, with two .ipynb notebooks.

<img src="howto-notebooks.png" width="640">

Open the notebook with the attack on DES round XOR.

<img src="howto-notebook-des.png" width="640">

The code in the notebook is arranged in cells. Here is the cell with attack settings:

<img src="howto-notebook-des-settings.png" width="640">

The example traceset comes in npz format. Execute cells one-by-one and get an attack result plot in-line.

<img src="howto-notebook-des-result.png" width="640">

This is the instantaneous attack result. To see the evolution, you can proceed with the following cells.

This is it so far for the basics.
