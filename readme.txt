This is pysca toolbox by Ilya.
Version 0.3, 2015-10-22

The toolbox can do LRA and CPA on AES and DES.
It works traces converted from Inspector trs format.
The conversion script is included.

Can do, unlike Inspector:
- LRA on DES
- user-defined "leakage models" for LRA on AES and DES
- CPA faster than Inspector multi-threaded (thanks conditional averaging)
- flexible visualization
- hackable code

Cannot do so far:
- inner round intermediates
- work with trs directly
- GUI