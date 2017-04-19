#
#							the MERpsychophysics (an R based software)
#	Version:0
#	Author: Alessandro Moscatelli <moskante@gmail.com>
#	Depends: R(>= 2.15.0), lme4, mnormt, Matrix, lattice
#	Date: 20.10.2012

#	Description: The software provides functions to estimate the Point of Subjective Equality (PSE) and the Just Noticeable Difference (JND) within the GLMM framework. Fuctions require as input a GLMM fitted with glmer{lme4} (a "mer" object. See package lme4 for further details). By default it only allows a probit link function (corresponding to the cumulative Gaussian distribution in the psychometric function). The function MERplot() plot a mer object as psychometric functions. The function MERboot() estimate the PSE, the JND and their confidence interval by means of bootstrap method. The function delta.psy.probit estimate the PSE, the JND and their confidence interval by means of delta method. This function can be used on glm object (lme4 = F) or on a mer object (lme4 = T) having a probit link function. The other two functions resample.mer() and suff.GLMM() are called by MERboot() for resampling and ordering the output. The function MERsimulate() simulate a dataset from a psychophysical experiment.
source("MERsimulate.r")
source("kombo.r")
source("xplode.mer.r")
source("MERdelta.probit.r")
source("MERtreatment.r")
source("pseMer.R")

#for the single-subject analysis (psychometric function)
source("delta.psy.probit.r")
source("psych.function.r")

cat("The MERpsychophysics proto-package\n",
	"Author: Alessandro Moscatelli\n",
	"\n",
	"List of Functions:\n", ls(),"\n",
	"\n",
	"The MERpsychophysics is a free R-based software \n",
	"The MERpsychophysics is distributed as is, and comes with ABSOLUTELY NO WARRANTY.\n",
	"No liability is accepted for any damage or loss resulting from the use of these routines.\n",
	"\n",
	"This software is distributed under the terms of the GNU General Public License, either Version 2, June 1991 or Version 3, June 2007.\n",
	"The terms of version 2 of the license are in a file called LICENSE.txt which you should have received with this software.\n",
	"Copies of both versions 2 and 3 of the license can be found at http://www.R-project.org/Licenses/.\n",
	"\n",
	"If you use the MERpsychophysics software please cite the following article:\n",
	"	Moscatelli, A; Mezzetti, M; Lacquaniti, F. (2012). Modelling \n",
	"	Psychophysical Data at the Population-Level: The Generalized Linear \n",
	"	Mixed Model. Journal of Vision, 12(11):26\n",
	"See also citation() and citation('lme4') to cite R and the lme4 package\n")

library(lme4)
library(mnormt)
#function pseMer also requires boot and beepr