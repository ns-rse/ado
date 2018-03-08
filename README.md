# Overview

This repository contains two sets of [Stata](https://www.stata.com/) files...

* `genass` a Stata programme I wrote to facilitate analysing large numbers of [SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism).
* `clayton` a copy of [David Clayton](https://en.wikipedia.org/wiki/David_Clayton) `genassoc` package which used to be available at [http://www-gene.cimr.cam.ac.uk/clayton/software/stata/](http://www-gene.cimr.cam.ac.uk/clayton/software/stata/), but no longer is.  These are a dependency of `genass`.


I do not claim any copyright or authority over David Clayton's files they are provided here for convenience  as `genass` uses them and over the years I've been contacted by a few people who wished to use `genass` but could not access David's `genassoc` package.

## Installation

You can clone this repository locally using...

    git clone https://github.com/ns-ctru/ado

You then have to ensure that Stata knows where these files are and this is achieved by modifying your `adopath` (see `man adopath` in Stata if you do not know how to do this).
