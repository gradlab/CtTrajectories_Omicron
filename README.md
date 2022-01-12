## Code and data associated with "Viral dynamics and duration of PCR positivity of the SARS-CoV-2 Omicron variant"

Stephen M. Kissler`*`, James Hay`*`, Joseph R. Fauver`*`, Christina Mack`*`, Caroline G. Tai, Radhika M. Samant, Deverick J. Anderson, Gaurav Khullar, Matthew MacKay, Miral Patel, Shannan Kelly, April Manhertz, Isaac Eiter, Daisy Salgado, Tim Baker, Ben Howard, Joel T. Dudley, Christopher E. Mason, David D. Ho, Nathan D. Grubaugh`*`, Yonatan H. Grad`*`

`*` denotes equal contribution

__proportion_over_time.R__ generates the time course of RT-qPCR positivity depicted in Figure 1 (requires data/ct_dat_subset_figure1.RData).

__run_analysis.R__ calls the remaining files to generate and save the MCMC viral trajectory fits depicted in Figure 2 (requires data/ct_dat_refined.RData).