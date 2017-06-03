# Pysca toolbox

This toolbox was started in 2014 to experiment with efficient differential power analysis (DPA) techniques from the paper "Behind the Scene of Side Channel Attacks" by Victor Lomné, Emmanuel Prouff, and Thomas Roche (https://eprint.iacr.org/2013/794).

To clone this repo with the included example traces you will need Git-LFS. Without Git-LFS, only pointers to traces will be cloned.

## Why
The toolbox was designed with the following in mind:
* state-of-the-art DPA techniques
* performance
* visualization of metrics for security evaluations purpose (and not just attack)
* simplicity and flexibility through use of a language suitable for scientific computing

In terms of these points, Pysca (still) outperforms some commercial tooling. Pysca is nowadays mostly superseded by https://github.com/Riscure/Jlsca.

## What
Pysca implements:
* non-profiled linear-regression analysis (LRA) with configurable basis functions
* classical correlation power analysis (CPA)
* significant speed-up of the above by conditional averaging
* targets: AES (S-box out) and DES (round in XOR round out, round out, S-box out)
* visualization of results

## How

For usage basics refer to the [HOWTO](howto/HOWTO.md).

For a deeper dive into leakage modelling using linear regression, clone the tutorial from a separate repo into the subfolder:

    git clone https://github.com/ikizhvatov/leakage-modelling-tutorial.git

## Details
Pysca works on traces stored in npz (numpy zipped) format. Example tracesets are included in the repo using git-lfs. The conversion script from Riscure Inspector trs format is included. The trs reader was originally implemented by Erik van den Brink.

Under the hood, the most interesting technical tricks in pysca are perhaps:
* fast computation of correlation (see https://github.com/ikizhvatov/efficient-columnwise-correlation for a dedicated study)
* conditional averaging implementation for DES (because of all the bit permutations, it requires splitting the leakage function into two stages)

Author: Ilya Kizhvatov<br>
Version: 1.0, 2017-05-14
