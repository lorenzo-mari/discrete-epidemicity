The scripts in this folder replicate the full set of reactivity/epidemicity analyses reported in the paper and reproduce all figures included within the main text and the supporting information (SI). Note that, in some cases, results are pre-calculated in order to save lengthy computational times. Users who whish to repeat all calculations should set the flag_comp variable equal to 1 in the relevant routines. 

Specifically:
- plot_figures.m is the master script through which the individual routines (included in the scripts folder) that generate the figures are called; 
- scripts/recurrent.m generates Figure 1 in the main text (COVID-19 outbreak under recurrent imports of infectious individuals);
- scripts/all_diseases.m generates Figure 2 in the main text (epidemicity analysis for 15 respiratory viruses) and Figure S10 in the SI (same as Figure 2, but using the L2 norm); 
- scripts/ticks_control.m generates Figure 3 in the main text (epidemicity analysis for seven tickborne zoonoses) and Figure S11 in the SI (same as Figure 3, but using different control targets);
- scripts/parametric.m generates Figure 4 in the main text (parametric epidemicity analysis for emerging infectious diseases);
- scripts/COVID_19.m generates Figure S1 in the SI (COVID-19 case study);
- scripts/measles.m generates Figure S2 in the SI (measles case study);
- scripts/measles_Delta.m generates Figure S3 in the SI (measles case study, sensitivity analysis to different discretization timesteps);
- scripts/sensitivity.m generates Figures S4--S9 in the SI (respiratory viruses, sensitivity analysis to parameter variations);
- scripts/ticks_control_rnd.m generates Figure S12 in the SI (tickborne zoonoses, sensitivity analysis to parameter variability);
- scripts/invasive_num.m generates Figure S13 in the SI (reactivity analysis for eight invasive plant species).

All scripts have been developed and tested with MATLAB R2024b and require the Statistics and Machine Learning Toolbox. They are fully compatible with the online version of MATLAB, for which a free access option exists. Users looking for an open-source alternative can try with GNU Octave (https://octave.org/).
