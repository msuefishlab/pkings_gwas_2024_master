*Paramormyrops kingsleyae* Phenotyping
==================

**Data Analysis**

This directory contains source `MATLAB` code for phenotyping EODs as well as source data and output files for inspection.

Original EOD recordings are located in `input_data/02_Phenotyping/eods`, PDF file plots of EOD landmarks are located in `output_data/02_EOD_Plots` along with the measurements of these landmarks in `output_data/02_EOD_Plots/measurement_data.csv`

The analysis can be regenerated issuing the following commands to `MATLAB` in this directory:

```matlab
EOD_Analysis
Plot_EODs
```

`EOD_Analysis.m` performs the following operations:

> For each file in `input_data/02_Phenotyping/eods` directory:
>  1. Extract "short" EOD recordings (recordings that encompass only a single EOD)
>  2. Find the peaks of the EOD waveform
>  3. Flip so that P1 is always positive (if necessary)
>  4. Normalize P-T-P amplitude = 1V
>  5. Center waveforms so that T0=P1
>  6. Record these EODs to an array with metadata
>  7. Average the EODs to a second array with metadata

`Plot_EODs.m` performs the following operations:

>For each file in `input_data/02_Phenotyping/eods` directory:
>  1. Create a PDF figure.  The top row represents the centered, normalized EODs and the bottom row represents the averaged EOD.
>  2. Perform landmark analysis on the averaged EOD (using `standard_eod_measurement.m`).  Plot these points on the averaged EOD
>  3. Write a PDF file to `output_data/02_EOD_Plots`.
>  3. Write a row to `output_data/02_EOD_Plots` containing the landmark analysis results.

**Phenotype Calls**

Phenotypes were categorized based on manual inspection of PDF plots above for TP/BP/Wobbles (and TypeI/TypeII Magnostipes).  These assignments are summarized in `input_data/01_Terra/data_model/fish_data_2023.txt`
