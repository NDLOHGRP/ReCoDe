## :warning: Deprecated
This repository is deprecated in favor of [pyReCoDe](https://github.com/NDLOHGRP/pyReCoDe)   
Please check out pyReCoDe to see the latest developments.


# ReCoDe
A Reduced Compressed Description for Direct Electron Microscopy Data

The ReCoDe package can be used to reduce and compress sparse electron microscopy data produced by direct electron detectors. The current limited release supports only DEFLATE compression and has been tested on CentOS 7.

# Installation
The package requires the libz and OpenMP libraries to be present in your PATH.<br>
To build from source use the following command: <br>
<b>make recode</b><br>
This will build the executable in "bin"

# Sample Data
ReCoDe requires both the image frames and a dark reference dataset as input. Sample image and dark data files are provided in "data/Frame_0_12-13-59.436.bin" and "data/Dark_Frame_12-23-00.232.bin", respectively.

# Usage
To reduce and compress use the following command:<br>
<b>./recode -rc &lt;input image filename&gt; &lt;dark image filename&gt; &lt;input parameter filename&gt; &lt;output directory&gt;</b><br>
E.g.: ./recode -rc data/Frame_0_12-13-59.436.bin data/Dark_Frame_12-23-00.232.bin config/recode_params.txt scratch/<br>
  
To decompress and expand use the following command:<br>
<b>./recode -de &lt;recode filename&gt; &lt;output directory&gt;</b><br>
E.g.: ./recode -de scratch/Frame_0_12-13-59.436.bin.rc1 scratch/<br>
  
In the first example the "config/recode_params.txt" file specifies the reduction compression parameters, such as reduction level, number of OpenMP threads to be used etc.
