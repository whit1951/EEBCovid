# Lessons from movement ecology for the return to work: modeling contacts and the spread of COVID-19

This repository is associated with the publication "Lessons from movement ecology
for the return to work: modeling contacts and the spread of COVID-19," currently
[published as a pre-print on medR$\chi$iv](https://www.medrxiv.org/content/10.1101/2020.05.27.20114728v2).

There are four subdirectories in this repository:

## [data](data/)

This folder contains a single file of annonymized office and lab locations for members of the
[University of Minnesota's Department of Ecology, Evolution, and Behavior](https://cbs.umn.edu/academics/departments/eeb).
These data were used to construct networks of shared spaces for the network model
portion of the aforementioned manuscript.

## [code](code/)

This folder contains several scripts to run analyses and generate figures associated
with the publication. In particular, see [EEBNetwork.Rmd](code/EEBNetwork.Rmd) for a
walkthrough of the network generation.

Code associated wtih analyses run for the supporting information exploring the
sensitivity of results in the main text are in the further sub-directory
[Supplemenatary_Analyses](code/Supplementary_Analyses/).

## [results](results/)

Many of the scripts in the [code](code/) directory will save their outputs for
convenience in this [results](results/) folder...

## [figures](figures/)

as well as save some figures in the [figures](figures/) folder. Images saved here
include those similar to the final figures used in the text. Note that most of the
main text figures have their own sub-directory which additionally contains code to
generate those figures, e.g. [Figure 7 from the main text](figures/disease_dynamics/),
showing a summary of disease dynamics in different networks.
