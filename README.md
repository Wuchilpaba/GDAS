# GDAS
Developed By Python and deployed with Docker

Supported running glycoproteomics analysis in GlycReSoft or O-Pair with images, other analysis process supported by Ray images 

GDAS's core novelty lies in a streamlined, multi-step workflow that strategically reduces computational load: it first employs an ultrafast open search (e.g., MSFragger-Glyco) on mass spectrometry data to rapidly screen for significantly regulated glycoproteins, statistically reducing the vast proteome database to a manageable subset. This resource-conserving step enables subsequent, in-depth, and targeted analysis for both N- and O-glycosylation using specialized tools (e.g., GlycReSoft and O-Pair). Furthermore, a unique Final Analysis Module utilizes an advanced statistical and machine learning pipeline (incorporating Bootstrap/Bayesian methods, XGBoost, and Random Forest) to integrate quantitative results and generate a robust, comprehensive glycosylation score for confident identification of disease-specific glycoforms. We validate GDAS using published Alzheimer's disease data, demonstrating its power to pinpoint biologically relevant glycosylation changes in targeted proteins.  


# Citation:
Kong, A. T., Leprevost, F. V., Avtonomov, D. M., Mellacheruvu, D., & Nesvizhskii, A. I. (2017). MSFragger: ultrafast and comprehensive peptide identification in mass spectrometry–based proteomics. Nature Methods, 14(5), 513-520.

Yu, F., Teo, G. C., Kong, A. T., Haynes, S. E., Avtonomov, D. M., Geiszler, D. J., & Nesvizhskii, A. I. (2020). Identification of modified peptides using localization-aware open search. Nature Communications, 11(1), 1-9.

Polasky, D. A., Yu, F., Teo, G. C., & Nesvizhskii, A. I. (2020). Fast and Comprehensive N-and O-glycoproteomics analysis with MSFragger-Glyco. Nature Methods, 17(11), 1125-1132.

Yu, F., Haynes, S. E., Teo, G. C., Avtonomov, D. M., Polasky, D. A., & Nesvizhskii, A. I. (2020). Fast Quantitative Analysis of timsTOF PASEF Data with MSFragger and IonQuant. Molecular & Cellular Proteomics, 19(9), 1575-1585.

Polasky, D., Geiszler, D., Yu, F., Li, K., Teo, G. C., & Nesvizhskii, A. I., (2023). MSFragger-Labile: A flexible method to improve labile PTM analysis in proteomics. Molecular & Cellular Proteomics, 22(5), 100538.

Yu, F., Teo, G. C., Kong, A. T., Fröhlich, K., Li, G. X., Demichev, V., & Nesvizhskii, A. I., (2023). Analysis of DIA proteomics data using MSFragger-DIA and FragPipe computational platform. Nature Communications, 14, 4154.

Yu, F., Yamei Deng, & Nesvizhskii, A. I. (2025). MSFragger-DDA+ Enhances Peptide Identification Sensitivity with Full Isolation Window Search. Nature Communications, 16, 3329.

Klein, J., Carvalho, L., & Zaia, J. (2018). Application of Network Smoothing to Glycan LC-MS Profiling. Bioinformatics.

Lei Lu, Nicholas M. Riley, Michael R. Shortreed, Carolyn R. Bertozzi & Lloyd M. Smith (2020). O-Pair Search with MetaMorpheus for O-glycopeptide characterization. Nature Communications, 17, 1133-1138.

M. Bern, Yong J. Kil, Christopher H. Becker (2012). Byonic: Advanced Peptide and Protein Identification Software. Current Protocols in Bioinformatics, 40.
