# TPD-Kinetic-Model
Matlab scripts for modeling targeted protein degradation. Please see model details in "A quantitative approach for defining degradability landscape of protein degraders".  


# System Requirements
- OS Requirements: Supported for macOS, Windows, Linux. The scripts have tested on macOS Ventura 13.4.1.
- MATLAB Requirements: Supported for MATLAB 2023a and later versions. Require the Curve Fitting Toolbox.


# Installation Guide
Download full directories from GitHub to local. (Estimated installation time: <5 min)


# Demo Instructions
Purpose,output,and usage for each directory of scripts:
- BRD: BRD system, reproduce Figure 2A/2B. Adapt this to simulate time course degradation response of your target of interest. Update measured degradation time course data file BRD_tc.csv and parameters in BRD_main.m with your target of interest. (Estimated run time: <10min)
- KINOME: KINOME system, generate results for Figure 3. Adapt this to fit degradation rate from single dose degradation level measurement of your target of interest. Update measured degradation level data file Foretinib_data.csv or TAE684_data.csv and parameters in KINOME_main.m with your target of interest. (Estimated run time: <10min)
- LANDSCAPE: BTK system, reproduce Figure 4A/4B. Adapt this to simulate degradability landscape of your target of interest. Update measured degradation level data file BTK_data.csv and parameters in LANDSCAPE_main.m with your target of interest.(Estimated run time: <30min)
- WEE1: WEE1 system, reproduce Figure 5A/5B/5E. Adapt this to simulate functional inhibition, inhibition landscape and hook effect of your target of interest. Update measured degradation and acitivity level data file WEE1_data.csv and parameters in WEE1_main.m with your target of interest.(Estimated run time: <30min)

Directory structure:\
Each directory hosts three Matlab scripts: 
- *_main.m: the main program (load input data, specify parameters, and run simulations)
- *_model.m: ODE model, called by *_main.m   
- events_time.m: controls ODE time steps, called by *_main.m

Run guide:\
Run *_main.m from within each directory. 

# License
This code is released under the MIT License. 

# Contact
Please E-mail weidu.pku@gmail.com for inquiries. 
