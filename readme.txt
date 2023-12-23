The codes here are used for generating the corresponding figures for the manuscript:

Mechanisms of Dominant Electrophysiological Features of Four Subtypes of Layer 1 Interneurons

John Hongyu Meng 1, Benjamin Schuman 2, Bernardo Rudy 2,3, Xiao-Jing Wang (2023) Journal of Neuronscience.
https://doi.org/10.1523/JNEUROSCI.1876-22.2023

Previous bioAxiv versions:
https://www.biorxiv.org/content/10.1101/2022.08.23.505010v2.abstract

The files may not be complete in the review page on the ModelDB page. If happens, download the whole .zip file and everything should be there. 

For models, which are under ./Models folder:

The models for 4 subtypes of L1 interneurons are saved under corresponding folders. 
*  The <modelName>.py files under each subfolder are the main models
*  The corresponding plot script may be required to generate corresponding figures.

The codes are tested under Python 3.9. To run the code, the Brian2 package is required.
https://brian2.readthedocs.io/en/stable/index.html

The OU process coding of Brian2 is based on 
https://brian2.readthedocs.io/en/stable/examples/advanced.Ornstein_Uhlenbeck.html?highlight=ou%20process#example-ornstein-uhlenbeck
For sufficient reference:
http://var.scholarpedia.org/article/Stochastic_dynamical_systems


For the Dynamic gain, the example analysis code for the VIP cell model is saved under the ./DG folder
*  To generate figures, 1st to run the VIPlongsingle.py script, then run DG_main.py script.

For tools, which are under ./Tools folder
The codes are tested under Matlab R2021b. Data from one alpha7 cell is included as an example. The original data was collected by pClamp.
Run plotv_L1.m to generate the sample analysis. 

Editted by John Meng, 02/28/2023. Reach me at john.meng@nyu.edu
