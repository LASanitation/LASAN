# LASAN
Project to develop indices of ecological condition for the city of Los Angeles.

Source data for species observations are available in this directory: https://drive.google.com/drive/folders/1F_LaKPq1fOZKxmS6985-bZXL_umzn35z?usp=sharing

Source data for environmental layers are linked and described here: https://docs.google.com/spreadsheets/d/1JPwPqgNaD17vz2PMCRRdepr3YSNCH3t6Nkw56Bj3-Po/edit?usp=sharing

Generating bias rasters to account for geographic clustering in our observational data.
Script to generate bias raster for non-native species observations in Los Angeles: https://github.com/LASanitation/LASAN/blob/main/LABiasRasterV1.R
Script to generate bias raster for native species observations in Los Angeles: https://github.com/LASanitation/LASAN/blob/main/LABiasRasterV2.R

Evaluating species distribution models generated using Maxent.
Determine the sensitivity of non-native species to environmental gradients: https://github.com/LASanitation/LASAN/blob/main/LAIndicatorTaxaV1.R
Determine the sensitivity of native species to environmental gradients: https://github.com/LASanitation/LASAN/blob/main/LAIndicatorTaxaV2.R

Select indicator species based on their relatively high sensitivity to environmental gradients.
Generate species distribution models of indicator species, collapse their prediction maps in order to generate ecological indices for Los Angeles.
Generate the Los Angeles Ecological Index (LAEI): https://github.com/LASanitation/LASAN/blob/main/LAEcologyIndexGeneratorV1.R
Generate the Native Los Angeles Ecological Index (nLAEI): https://github.com/LASanitation/LASAN/blob/main/LAEcologyIndexGeneratorV2.R

Use random forest to evaluate the reliability of ecological indices across the extent of Los Angeles.
Evaluate the LAEI: https://github.com/LASanitation/LASAN/blob/main/LAEcologyIndexEvaluatorV1.R
Evaluate the nLAEI: https://github.com/LASanitation/LASAN/blob/main/LAEcologyIndexEvaluatorV2.R
