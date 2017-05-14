# Pysca toolbox

This toolbox was started in 2014 to experiment with efficient differential power analysis (DPA) techniques from the paper "Behind the Scene of Side Channel Attacks" by Victor Lomné, Emmanuel Prouff, and Thomas Roche (https://eprint.iacr.org/2013/794).

The toolox can do:
* non-profiled linear-regression analysis (LRA)
* classical correlation power analysis (CPA)
* speed-up of the above by conditional averaging
* AES (S-box out) and DES (round in XOR round out, round out, S-box out)
* visualization of results

It works on traces converted from trs format into numpy npz. The conversion script is included. The trs reader was originally implemented by Erik van den Brink.

Pysca runs under python 2.7, numpy 1.12.1, matplotlib 2.0.2. See project wiki for usage examples. The example traces are included in the repo using git-lfs.

Pysca is mostly superseded by https://github.com/Riscure/Jlsca.

Author: Ilya Kizhvatov

Version: 1.0, 2017-05-14
