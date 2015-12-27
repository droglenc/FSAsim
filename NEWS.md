# FSAsim 0.0.4 ongoing
* Compiling under R 3.2.3 and using roxygen2 5.0.1
* DESCRIPTION: Updated.
* Added Imports for `manipulate`, `relax`, `graphics`, and `stats`.  Added Suggests for `tcltk` (ultimately want to get rid of this and `relax`).
* `growthModelSim()`: Modified.  Fixed links to functions in `FSA` that had changed names.
* `leslieSim()`: Modified.  Major overhaul including changing to use `manipulate` rather than `relax`, moved `removals` to an argument rather than a slider bar, chagned `type=` to `sim=`, and dramatically changed the plot for `sim="montecarlo"`.
* `mrClosed1Sim()`: Modified.  Changed to use `manipulate` rather than `relax`.  Removed the `type=` argument and put it as a dynamic choice.  Added the `sim=` argument which allows simulations of the `distribution` or for assessing `assumptions` violations.  Changed look of the histogram (more like `hist()` from `FSA`).  Force more bins in the histogram.  Modified the legend look.  Changed default number of `rsmpls=`.  Streamlined code some.
* `vbComp()`: Modified.  Reorganized code and deleted bad fishR links in documentation.
* `vbDataGen()`: Modified.  Removed bad link to `gReshape()` in `FSA`, which no longer exists.

# FSAsim 0.0.3 May15
* DESCRIPTION: Updated.
* Added importFrom for `tidyr` (for `spread()` in `vbDataGen()`).
* `growthModelSim()`: Modified.  Changed links to `FSA` in the help file.
* `vbDataGen()`: A complete rebuild.  See help file.

# FSAsim 0.0.2 Mar15
* `growthModelSim()`: Added.
* `leslieSim()`: Modified.  Updated for changes to `depletion()` in `FSA`.
* `simAgeBias()`: Modified.  Updated for changes to `ageBias()` in `FSA`.
* `srSim()`: Added.  Still needs work to add functionality for "Shepherd" and "SailaLorda" methods.

# FSAsim 0.0.1 May14
* Compiling under R 3.1.0, using roxygen2 4.0.0, and using github as a repository.
* `catchCurveSim()`: Added.
* `cohortSim()`: Added.
* `lengthWeightSim()`: Added.
* `leslieSim()`: Added.
* `lwModelSim()`: Added.
* `sample4ALK()`: Added.
* `simAgeBias()`: Added.
* `simAges()`: Added.
* `simLenFromAge()`: Added.
* `simLenSelect()`: Added.
* `mrClosed1Sim()`: Added.
* `srCobWeb()`: Added.
* `vbComp()`: Added.
* `VBGMlit`: Added.
* `vbDataGen()`: Added.
