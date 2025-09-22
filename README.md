<h1>ObjFind-M </h1>

ObjFind-M is a Python framework for uncovering the metabolic objectives of mammalian cells directly from experimental data.
Instead of assuming what cells optimize (like growth or productivity), ObjFind-M uses fluxomic and metabolomic datasets to infer the actual reaction-level goals that cells appear to follow.

<img width="2666" height="1033" alt="New ObjFindCHO Methodology (13)" src="https://github.com/user-attachments/assets/d1b8a8cf-35bb-4ca5-b4be-9b30268c253c" />

We provide two ready-to-use scripts:

•	ObjFind-M CHOmpact.py → runs on a small, simplified CHO cell model (“CHOmpact”), which is easy to interpret and useful for testing.

•	ObjFind-M iCHO2441.py → runs on a large genome-scale model (iCHO2441) to capture the full complexity of CHO cell metabolism.


What does it do?

-Reads in published fluxomics and metabolomics data.  
-Fits these data to a metabolic network model.  
-Infers a “fitness function” (a weighted set of reactions) that best explains how the cells distribute their metabolism.  
-Outputs ranked reaction coefficients and predicted fluxes to Excel files.  

This lets you see which reactions (like ATP synthase, TCA cycle steps, or lipid shuttles) are most central to cell objectives under different conditions.

How to run

1.	Make sure you have Python 3 and the following libraries installed:
o	COBRApy
o	Pyomo
o	pandas
o	NumPy
o	A solver such as Ipopt.
2.	Adjust the CONFIG section at the top of either script to point to your model file, exchange data, intracellular flux data, and mapping file.
3.	Run "ObjFind-M CHOmpact.py" or "ObjFind-M iCHO2441.py"
4.	Results are saved to Excel files (coefficients and fluxes).
   
Repository contents
•	ObjFind-M CHOmpact.py – example run using the reduced CHOmpact model.  
•	ObjFind-M iCHO2441.py – example run using the iCHO2441 genome-scale model.  
•	Mapping and data files (Excel sheets) – link experimental data to model reactions.  


Citation
If you use this code, please cite:
Morrissey et al. (2025). Inferring the Metabolic Objectives of Mammalian Cell Culture through Inverse Modeling of Metabolomics and Fluxomics.

