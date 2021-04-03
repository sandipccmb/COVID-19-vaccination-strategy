# COVID-19-vaccination-strategy
Code for the article:
PRAGMATIC COVID-19 VACCINATION STRATEGY FOR INDIA: A MATHEMATICAL MODELLING BASED ANALYSIS
Sandip Mandal; Nimalan Arinaminpathy; Balram Bhargava; Samiran Panda
 
This repository contains codes and data used to simulate and analyze COVID-19 vaccination strategies for India under various scenarios.
1.	Impact of vaccination among priority groups, on symptomatic cases and on deaths 
2.	Deaths averted by different strategies for an infection-preventing vaccine and for a symptomatic disease preventing vaccine 

Note: The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 
 
OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). However as Matlab and Octave are available for most operating systems, codes should run on Mac OSX and Linux as well.

Installation guide
MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2018b.

Codes & their functionality

Setup_model.m
This model has six state variables which are further divided into three age groups, two risk group (co-morbid or not) and two groups according to their vaccination status. All the variables and parameter values are assigned in this script. 

Make_model2.m
This is a function which specifies the transmission model in matrix form, capturing the linear and non-linear parts separately.

goveqs_basis3.m
Using the matrix formulation constructed in Make_model2, this computes the time derivatives for each state variable given values for those states.

linspecer.m
This function is used to plot multiple lines with distinguishable and nice looking colors.

jbfill.m
This routine will shade the area of a 2-D plot between two user defined vectors. 
(Ref. John Bockstege (2021). Shade area between two curves (https://www.mathworks.com/matlabcentral/fileexchange/13188-shade-area-between-two-curves), MATLAB Central File Exchange. Retrieved April 2, 2021.)

get_init.m
This function assigns the initial conditions for simulation, depending on vaccination coverage, and existing seroprevalence. 

Get_aggregators.m
Function input ind_sets is a cell array, each element a vector showing the set of indices that need to be aggregated over.

Get_address.m
We make use of Matlab ‘structures’ (the equivalent of ‘lists’ in R) to assist in book-keeping on the indices corresponding to different state variables. This function assists in constructing those indices. For example, the index named i.A.v0.r1.a3 is the index corresponding to the asymptomatic, non-vaccinated, with co-morbidity and among elderly population. 

Alloc_parameters.m
For a given parameter vector x, this function allocates the corresponding values to each of the parameters being sampled.



Instructions for use

In the above titled article, there are two figures (Figure 3 and Figure 4) in the main text and four figures (figure S1-S4) in the supplementary document. To generate these figures and output table (Table 1) select the appropriate script and run it in MATLAB. Save the data file as instructed at the bottom of each script and then run the file named as ‘Figure X.m’ to get the figure. Parameters should be adjustment for infection preventing vaccine and symptomatic disease preventing vaccine. 

Setup_model.m
Model variables and parameter values are assigned in this script. To setup the model for infection preventing vaccine choose the parameter value p.c1 = 0.6 (when vaccine efficacy is 60%) and 0.9 (if efficacy = 90%) and p.c3 = 0. For disease preventing vaccine p.c3 = 0.6 (when vaccine efficacy is 60%) and 0.9 (if efficacy = 90%) and p.c1 = 0;

Simulate.m
After setting up the model, run this code to get data for the plot “Figure 3” of this article. This figure requires six data set. First choose the priority groups whom to vaccinate. Select scenario = 1 (when priority group is only HCW+FW); scenario = 2 (when priority group is only HCW+FW + Co-morbid (<60y)); scenario = 3 (when priority group is only HCW+FW + Co-morbid (<60y)+ elderly);   

For susceptibility reducing (infection preventing) vaccine save scenario1 as ‘simulate1_sus.mat”. Similarly for the scenario 2 and 3 save the data as ‘simulate2_sus.mat” and ‘simulate3_sus.mat” respectively. Likewise, for symptomatic disease preventing vaccine save three scenarios as ‘simulate1_sev.mat’, ‘simulate2_sev.mat’ and ‘simulate3_sev.mat’.
  
Figure3.m

Illustration of vaccine impact of each of the priority groups can be obtained by running this code. 

Table1.m

Load appropriate data (such as ‘simulate1_sus.mat’ to find the impact of vaccination among HCW+EW only) and select ‘target_pop’ accordingly. This code gives percentage reduction in peak symptomatic incidence, percentage reduction in cumulative mortality and number needed to vaccinate to avert one death.     
 
Plot_results.m

First set up the model either for ‘Infection preventing vaccine’ or ‘Disease preventing vaccine’. Then choose R0 value (for the current article we have chosen 1.25, 2, 2.5) and save the data as Out_1p25.m, Out_2.m and Out_2p5.m corresponding to each R0 value. 

FigureS1.m

This gives the plot of deaths averted by vaccination of different priority based scenarios and uniform scenario using infection preventing vaccine of efficacy 60%. Required data files are Out_1p25.m, Out_2.m and Out_2p5.m (corresponding to this scenario).

FigureS2.m

This gives the plot of deaths averted by vaccination of different priority based scenarios and uniform scenario using symptomatic disease preventing vaccine of efficacy 60%. Required data files are Out_1p25.m, Out_2.m and Out_2p5.m (corresponding to this scenario).

FigureS3.m

This gives the plot of deaths averted by vaccination of different priority based scenarios and uniform scenario using infection preventing vaccine of efficacy 90%. Required data files are Out_1p25.m, Out_2.m and Out_2p5.m (corresponding to this scenario).

FigureS4.m

This gives the plot of deaths averted by vaccination of different priority based scenarios and uniform scenario using symptomatic disease preventing vaccine of efficacy 90%. Required data files are Out_1p25.m, Out_2.m and Out_2p5.m (corresponding to this scenario).

Figure4_toppanel.m

Optimal priority strategy for an infection preventing vaccine for three different R0 values are plotted with this code. Optimal strategy for each R0 value has been chosen from FigureS1.m.

Figure4_bottompanel.m

Optimal priority strategy for a symptomatic disease preventing vaccine for three different R0 values are plotted with this code. Optimal strategy for each R0 value has been chosen from FigureS2.m.   
   


MAT files:

Out_1p25.mat
Out_2.mat
Out_2p5.mat
Simulate1_sus.mat
Simulate2_sus.mat
Simulate3_sus.mat
Simulate1_sev.mat
Simulate2_sev.mat
Simulate3_sev.mat




