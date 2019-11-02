# FSAsim 0.0.6 ongoing
* Removed dependency on `relax` and `tcltk`.
* `srSim()`: Modified. Replaced dependency on `relax` with `manipulate`. Removed `formula` and `data` arguments and, thus, the functionality to plot with data as this was used to find reasonable starting values for the model fitting which is better accomplished with `srStartsDP()`.
* `srStartsDP()`:  Modified. Replaced dependency on `relax` with `manipulate`.
* `vbStartsDP()`:  Modified. Replaced dependency on `relax` with `manipulate`.


# FSAsim 0.0.5 10-Jan-17
* Added `::` to functions in many functions.
* `.onAttach()`: Removed.
* `mrClosed1Sim()`: Modified. Fixed bug in plot when an infinite population was estimated. Tried to streamline code some (did not seem to make faster).
* `sample4ALK()`: Modified. Fixed misuse of `df` where should have been `data`.
* `simLenSelectP()`: Modified. Changed plot with `show=TRUE` and actual probabilities so that they actually represent the beta distribution over the range of values. Removed an unnecessary catch for `max.height=`. Minor spelling corrections.
* `srStartsDP()`: Added from `FSA`.
* `vbDataGen()`: Modified. Changed ouput variable `ageFrac` to `ageFracG` and added `ageFracY`. Fixed a bug related to repeating the lengths-at-capture (appeared only when the sampling season was the whole year and a fish happened to be sampled on the first day). Streamlined some code and updated the documentation.
* `vbStartsDP()`: Added from `FSA`. Fixed a bug related to call to `vbStarts()` and how `type=` worked with `vbStarts()` (thanks to Jake Lowe for pointing this out).

# FSAsim 0.0.4 1-Feb-16
* Compiling under R 3.2.3 and using roxygen2 5.0.1
* DESCRIPTION: Updated.
* Added Imports for `manipulate`, `relax`, `graphics`, and `stats`. Added Suggests for `tcltk` (ultimately want to get rid of this and `relax`).
* `catchCurveSim()`: Modified. Major overhaul including changing to use `manipulate` rather than `relax`, removed `Zsteady=` and `NoSteady=` in favor of a checkbox (but must be the same for both Z and No), removed `ZDeltaAge=` and `NoDeltaAge=` in favor of one `deltaAge=` to be used for both Z and No. Moved legend to bottomleft. Changed lines such that including randomness in Z or No is considered an assumption violation.
* `cohortSim()`: Modified. Changed to use `manipulate` rather than `relax`.
* `growthModelSim()`: Modified. Changed to use `manipulate` rather than `relax`. Removed the ability to plot actual data (see `vbStarts`) for this functionality. Changed parameter limits for most models. Fixed links to functions in `FSA` that had changed names.
* `lengthWeightSim()`: Modified. Changed to use `manipulate` rather than `relax`.
* `leslieSim()`: Modified. Major overhaul including changing to use `manipulate` rather than `relax`, moved `removals` to an argument rather than a slider bar, chagned `type=` to `sim=`, and dramatically changed the plot for `sim="montecarlo"`.
* `lwModelSim()`: Deleted. Same as `lengthWeightSim()`.
* `mrClosed1Sim()`: Modified. Changed to use `manipulate` rather than `relax`. Removed the `type=` argument and put it as a dynamic choice. Added the `sim=` argument which allows simulations of the `distribution` or for assessing `assumptions` violations. Changed look of the histogram (more like `hist()` from `FSA`). Force more bins in the histogram. Modified the legend look. Changed default number of `rsmpls=`. Streamlined code some.
* `vbComp()`: Modified. Reorganized code and deleted bad fishR links in documentation.
* `vbDataGen()`: Modified. Removed bad link to `gReshape()` in `FSA`, which no longer exists.

# FSAsim 0.0.3 May15
* DESCRIPTION: Updated.
* Added importFrom for `tidyr` (for `spread()` in `vbDataGen()`).
* `growthModelSim()`: Modified. Changed links to `FSA` in the help file.
* `vbDataGen()`: A complete rebuild. See help file.

# FSAsim 0.0.2 Mar15
* `growthModelSim()`: Added.
* `leslieSim()`: Modified. Updated for changes to `depletion()` in `FSA`.
* `simAgeBias()`: Modified. Updated for changes to `ageBias()` in `FSA`.
* `srSim()`: Added. Still needs work to add functionality for "Shepherd" and "SailaLorda" methods.

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
