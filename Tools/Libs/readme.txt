08242021
This folder contains all the functions I need to do the analysis of the pClamp data.
The goals are the following:
1st, generate the indexes to classify PV vs. SST. What Shlomo needs are
    a. Slope of the f-I curve.
    b. Rise rate/ decay rate. 
    c. Adaptation index. 
    d. (Max-Min)/Max of ISI. 
    e. Loop over I_inj/firing rate.
Ideally I would like to generate an analyzed file, so I don't need to run the code everytime. It will takes a while.
So I would prefer sth. like a state{} function that collect everything, and save it for latter use.

2nd, analyze the data for L1 interneurons.
    a. Reset
    b. Threshold.
    c. Standard for irregularity, ramping, and initial bursting.
    d. Plotting all those stuff.
    e. generate a pathway depended thing.


Let me think about the work flow. Usually, comes from 3 parts

Loading. Guess I will do it by hand. 

Loop through sweeps.

1a, 1c, 1d, 2c, 2d.

Loop through spikes. 

1b, 2a, 2b, 

Saving

2e

Plotting

2d.
