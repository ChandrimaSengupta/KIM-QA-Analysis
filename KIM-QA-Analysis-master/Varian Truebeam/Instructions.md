# Instructions

## Requirement: MATLAB, MCR 2018b 64 bit. 
MCR has been used for the purpose when the computers at different sites do not have MATLAB installed. In that case, install MCR on a computer that has MATLAB and create a standalone executable by writing 'mcc -e ***.mlapp' in the command prompt. This will create a standalone executable file. This executable file then can be used in any computer that has MCR istalled. 

## Instructions to run: 
1. Download static localisation codes (App_Static_loc.mlapp and Staticloc.m), dynamic localisation codes (App_Dynamic_loc.mlapp and Dynamicloc.m) and treatment interruption codes (App_treat_int.mlapp and TreatmentInt.m).
2. If MATLAB has been installed, Open the MATLAB app (.mlapp file) and run the code. The code looks for:
   a) 'Select KIM Traj folder' - The folder that contains KIM log files.
   b) 'Select Robot Traj file' - The robot trace file.
   c) 'Select coord file' - The patient coordinate file. 
   d) 'Select param file' - The parameters (a, b, c) that shift the KIM traces to match with the robot trace in time in the interval between a and b at a step of c. 
   e) Click on 'analyse'.
 
 ## Output
 This will produce an output file with the information mean, standard deviation and percentiles in LR, SI and AP directions for all the data points in the KIM input trace and a figure with KIM trace and robot trace on top of each other.
   
