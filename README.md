# MHC-Shake: MHCflurry benchmark recreation and application on antigen presenting prediction for TCGA NSCLC Patients

* ECBM 4060: Genomic Science Final Project
* Group Members: Harry Lee, Jason Mohabir, Lillian Wang, Shuangyi Xu

## Introduction

Major histocompatibility complex (MHC) is a set of surface proteins that plays an essential role in the adaptive immune system. MHC presents the pathogen-derived antigen onto the cell surface, and the presented peptide will be later recognized by T cells to initiate an immune response. The peptide is loaded onto MHC, therefore accurately predict the binding affinity for MHC with the antigen peptide has found wide application in autoimmunity, vaccine design, and cancer immunotherapy. 

A recently published allele-specific class I MHC binding prediction tool: MHCflurry 1.2.0 developed by O’Donnell et al. showed great potential for becoming a reliable resource in predicting affinity of MHC and antigen binding. The original intention that leads to the creation of MHCflurry is to facilitate cancer vaccine development, and it is currently in the process of improving affinity prediction and expand the use to non-allele specific. 

In their published paper, the benchmark results showed MHCflurry outperformance over the standard predictors NetMHC 4.0 and NetMHCpan 3.0 in terms of accuracy and speed overall. Since the publication, MHCflurry has been updated to improve the affinity prediction. In the currently newest version MHCflurry 1.4.0, the approach has been switched from allele specific to pan allele to accommodate the undersampled diversity of HLA genes. 

The first objective of this paper is to recreate the benchmark results and compare with the benchmark results using the updated MHCflurry 1.4.0. The second objective is to perform MHCflurry affinity analysis on TCGA NSCLC patients by creating pseudo cancer antigen based on the mutation data. By investigating the relation between mutation load burden and maximum MHC binding affinity, we expect to provide some guidance for patients when choosing a particular cancer immunotherapy. For simplicity and demonstration, we choose to use the most common and well characterized HLA allele A*01:01.

## Methods and Results

Please see detailed discussion of our work linked here: 

## Acknowledgments

* Thank you to Wei-Yi Cheng for your project guidance and support this semester. 
* Thank you to Timothy O'Donnell for meeting with our group and providing us with useful tips and discussion. 

## References

T. J. O’Donnell, A. Rubinsteyn, M. Bonsack, A. B. Riemer, U. Laserson, and J. Hammerbacher, "MHCflurry: Open-Source Class I MHC Binding Affinity Prediction," Cell Systems, 2018. Available at: https://www.cell.com/cell-systems/fulltext/S2405-4712(18)30232-1.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
