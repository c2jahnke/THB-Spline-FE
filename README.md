# THB-Spline-FE
1D finite element analysis with truncated hierarchical B-splines (THB-splines)

###### DISCLAIMER ######
November 2017, Jonathan Jahnke

NO WARRANTIES OR GARANTIES

THIS CODE WAS DEVELOPPED AS A MASTER'S THESIS PROJECT

TO USE IT, ALL FOLDERS AND SUBFOLDERS HAVE TO BE ADDED TO THE MATLAB PATH
########################

# CONTENT
2DPlots:
scripts to generate plots of 2D basis functions

AdditionalFunctions:
contains additional functions, e.g for HB and THB refinement, for the generation of Gauss points and weights and for error analysis.

CurvePlot:
algorithms to generate B-Spline curves

NurbsBookAlgorithms:
contains algorithms from the NURBS book by Piegl and Tiller for basic evaluation of B-splines and B-spline curves
 
BSplClasses:
contains bSplBas class and inherited classes bSplBasFun (one basis function), hbSplBas (hierarchical) and a multilevel class hbSplBasML. 
Further, thbSplBas inherits from hbSplBas and thbSplBasML from hbSplBasML.
boundCond is a class to save Dirichlet or Neumann boundary conditions.
PoissonSolv is a class to solve a Poisson model problem using multi basis such as hbSplBasML and thbSplBasMl.

Poisson:
Example file to solve a Poisson problem for different source terms

