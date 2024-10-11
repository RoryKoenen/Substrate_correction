Topic: -Biochemistry -Enzymology -Thrombosis and hemostasis.

The script is for the use of a continuous monitoring assay, in which a pro-enzyme is activated by an enzyme (complex). The activated proenzyme converts a substrate that is present in the solution. The resulting curve ideally is a parabola, but it may be distorted by the consumption of the substrate (depending on its kinetic parameters). The script can be used to correct the substrate conversion under conditions where there is an inhibitor of the proenzyme activating enzyme is present.
An example is tissue factor / factor VIIa that activates factor X, the factor Xa generated is measured by substate CS(11)65. An inhibitor (TFPI) can be titrated to the enzyme complex.
Reference: Peraramelli S, et al., Inhibition of tissue factor:factor VIIa-catalyzed factor IX and factor X activation by TFPI and TFPI constructs. J Thromb Haemost. 2014; 12: 1826-37. 10.1111/jth.12713.

This repository contains an R script that allows the correction of enzymatic substrate conversion for substrate consumption, recorded in time. The initial 5 minutes of the curves are then fitted with a parabola, of which the derivative is taken. This represents the amount of enzyme generated (activated) as a function of time. The slope can be used to calculate percentage of inhibition by the inhibitor.
The script returns the corrected curves, the fitted curves, first derivative, the fitting parameters, and the percentage of enzyme activity by increasing concentrations of inhibitor. All files are saved in time-labelled .csv formats in folders labelled as date and can be readily imported in spreadsheet software (provided that the locale is US or UK, ie. periods as decimal separators). 

The data can be imported directly from the clipboard. Important considerations are:

There is a script for Mac OSX and for Windows. This is because the clipboard is addressed differently in these environments.
The first column of the data always contains the substrate conversion curve on which all corrections are based. This generally means that the enzymatic substrate conversion is recorded in the absence of an inhibitor but else under the same conditions as all other recorded curves. Make sure that this curve is straight and has no dips or lags.
The dimensions are minutes for time and OD for the optical densities.
The locale is British-American, meaning periods to separate decimal numbers.
The data to be imported should be "rectangular", meaning only the named columns with the respective readings. No notes, dates, figures, filenames etc.
The titles or annotations are always in the first row of the data-array. These will be converted into syntactically correct titles for R.
The script corrects for empty cells, but be aware that all readings will be truncated to the length of the shortest measurement. This is because the used R data.frame class can only contain vectors with the same length.
The repository contains sample data that contains a control curve without inhibitor (is the reference for correction) and 2 curves that have different concentrations of inhibitor.
