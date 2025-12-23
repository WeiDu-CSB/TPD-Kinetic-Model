# TPD-Kinetic-Model
Matlab scripts for modeling targeted protein degradation. See model details in "A quantitative approach for defining degradability landscape of protein degraders".  

- BRD: BRD system, reproduce Figure 2A/2B. Adapt this to simulate time course degradation response of your target of interest. 
- KINOME: KINOME system, generate results for Figure 3. Adapt this to fit degradation rate from single dose degradation level measurement of your target of interest. 
- LANDSCAPE: BTK system, reproduce Figure 4A/4B. Adapt this to simulate degradability landscape of your target of interest. 
- WEE1: WEE1 system, reproduce Figure 5A/5B/5E. Adapt this to simulate functional inhibition, inhibition landscape and hook effect of your target of interest. 

Each directory hosts three Matlab scripts: 
- *_main.m: the main program (load input data, specify parameters, and run simulations)
- *_model.m: ODE model, called by *_main.m   
- events_time.m: controls ODE time steps, called by *_main.m
Run *_main.m from within each directory. 

E-mail weidu.pku@gmail.com for inquiries. 
