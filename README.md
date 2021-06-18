# disruptr
Framework to facilitate in-silico repression or network disruption experiments

Disruptr provides functions to facilitate the computation of network potential and change in network potential from user-provided mRNA expression data. 
Comprehensive methods can be found here (https://www.biorxiv.org/content/10.1101/854695v3), or by perusing function specific documentation. 
A few of the functions in this package rely on another package that I wrote which can be found at https://github.com/DavisWeaver/crosstalkr

## Features
 
The package ships with all the processed data from the paper above to allow users to mess around with the different functionality quickly and easily. 
Users can use these analytic tools in their own studies by inputting an expression matrix to the `compute_np()` function. 
Some of the slower operations can be sped up via parallelization by setting the `ncores` variable >1. 
The package also ships with some methods for downloading and processing PPI databases such as biogrid and stringdb

## Future work

Future work will include methods to process and incorporate GDSC data into these types of analysis, expanded metrics of cell state to allow different kinds of in-silico repression to be run, and methods to facilitate filtering larger networks into functionally relevant subnetworks.

## Support
 
Please feel free to reach out to me at dtw43@case.edu or create an issue here on github if something isn't working as expected. 
