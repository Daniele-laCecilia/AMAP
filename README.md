# AMAP
Probabilistic indicator for threshold-based risk assessment

The AMAP indicator quantifies the impact of uncertain parameters on the probability of a target quantity to exceed user-defined thresholds. The indicator can be readily employed in probabilistic risk assessment.

We applied AMAP in la Cecilia et al., (2020) to quantify the sensitivity of soil and aquifer contamination following herbicide glyphosate (GLP), and its transformation product aminomethylphosphonic acid (AMPA), dispersal to soil hydraulic parameters (k, b, psi_s, phi, Sr_l, Sr_g).
In the example provided here, the target quantity is GLP concentration in the aquifer.

# Curiosity:
"AMA" stands for Aronne Monica Alberto, who are the authors of the AMA family of indicators.
"AMAP" builds on the AMA family and adds the Probability framework.

# How to:
This function computes AMA sensitivity indices and principal Sobol Indices (Si) and AMAP starting from a set of Monte Carlo realizations. It also gives conditional moments with respect to one parameter.

To run the example:
If you use R, copy the .R scripts from the folder "R" and paste them in the folder "data"; then, run the main script "AMAP_analysis.R";
If you use Matlab, copy the .m scripts from the folder "Matlab" and paste them in the folder "data"; then, run the main script "AMAP_analysis.m".
In the main scripts (.R or .m), the Users can:
* Decide the indicators to plot by setting "var_out"
* Switch between an analysis of the parameters sensitivity to the target quantities at a fixed point in time or over time by setting the variable dotime=1 or dotime=0, respectively.

We invite the Users to test AMAP on their data, if suitable.
* Copy and paste the scripts in your working directory containing your data;
* In the INPUT section of the main scripts, the Users define:
  * The indicators to plot by setting "var_out";
  * The value of "dotime" to switch between an analysis at a fixed point in time or over time;
  * The .csv file name of the matrix "Sampling_points", containing Monte Carlo realizations of the parameters [n simulations x n parameters]. The first row contains the name of the uncertain parameters. The column delimiter is set to semicolomn;
  * The .csv file names (fixed time and over time) of the matrix "Output_Mat", containing full model evaluations [n simulations x n time points]. Should you have only one file, either fixed time or over time, please fill both entries with the same file name and then set "dotime" accordingly. The first row contains the name of the output variables. The column delimiter is set to semicolomn;
  * The number of classes for "nclass" used to compute conditional statistics;
  * The threshold value for "thr".
* The script automatically get "name_var", the name of each output variable (text string);
* In Matlab, the color palette shall be edited if the number of uncertain parameters is different than 6;
* Run.

# Plot
The plot generated is a barplot, which can have 1 column if a fixed time is used in the analysis or a number of columns equal to the number of time points in case an over time analysis is performed.
The higher the value, the higher the contribution of an uncertain parameter with respect to the selected target variable.
In analysis over time, the overall and single contribution of uncertain parameters can change. That is, at a particular time, the uncertainty of the parameters can be more or less significant than at other times.

# Requirements
The scripts were tested in R 3.5.0 and Matlab 2019a.

# Getting help
Should you encounter a clear bug, please describe it at https://github.com/Daniele-laCecilia/AMAP/issues.
Should you have questions/feedbacks about the analysis/function, you may reach us at our email addresses reported below.
Should you will to run a full uncertainty quantification analysis on your model, feel free to get in contact with us.

# Community engagement
Different indices may differently apportion the model parameters uncertainty. If this is the case for you using the AMA indices, please share your comparison with us via email. We will compile all the contributions and eventually deliver a final review with all the sources acknowledged.

# Authors
Daniele, la Cecilia (The University of Sydney, now at eawag, @ daniele.lacecilia at sydney.edu.eu)
Giovanni M., Porta (Politecnico di Milano @ giovanni.porta at polimi.it)

# Acknowledgements
Please feel free to use these functions and properly acknowledge this contribution by citing:
la Cecilia et al., (2020), should you use the AMAP function;
Dell'Oca et al., (2017), should you use the AMA family (AMAE, AMAV, AMAskew, AMAkurt).

Data and script were used in the journal article:
Daniele la Cecilia, Giovanni M. Porta, Fiona H.M. Tang, Monica Riva, Federico Maggi
Probabilistic indicators for soil and groundwater contamination risk assessment
Ecological Indicators
Volume 115
2020
106424
https://doi.org/10.1016/j.ecolind.2020.106424
http://www.sciencedirect.com/science/article/pii/S1470160X20303617

The AMA family of moment-based sensitivity indices was introduced in:
Aronne Dell’Oca, Monica Riva, Alberto Guadagnini
Moment-based metrics for global sensitivity analysis of hydrological systems
Hydrology and Earth System Sciences
Volume 21
2017
6219-6234
https://doi.org/10.5194/hess-21-6219-2017
https://www.hydrol-earth-syst-sci.net/21/6219/2017/

# Disclaimer
The concentration values are the result of a numerical model and were not validated against field data.
