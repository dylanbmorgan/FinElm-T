# FinElm-T
Coursework for PX912 as part of HetSys CDT at the University of Warwick, to perform finite element analysis on a shape, 'T'.

This has been written (pretty much) entirely in Julia so it's pretty quick. Functionality for hyperthreading has also been included to make things even faster, so on most systems, this code is likely to be limited by memory.

## Installation 
To install, clone this repository somewhere on your filesystem.

Ensure you have an up to date version of [Julia](https://julialang.org/downloads/) installed, in addition to the following dependencies:
- [Plots](https://docs.juliaplots.org/stable/install/)
- [PyCall](https://docs.juliahub.com/PyCall/GkzkC/1.92.0/)
- [matplotlib](https://matplotlib.org/stable/users/installing/index.html)
- [LatexStrings](https://docs.juliahub.com/LaTeXStrings/H4HGh/1.2.0/)
- [PrettyTables](https://ronisbr.github.io/PrettyTables.jl/stable/)

## Usage
This has been written with the intention of being using through the REPL. Navigate to `/src` and open the REPL.

Load the code using `include("Main.jl")`. The code can be executed simply by executing `finelm()`, where it runs with sane defaults, however there are also several optional arguments:

1. `time`: Boolean option on whether to print timings. Ensure to use with `@time finelm(...)`. Default = false
2. `points`: Integer option for the number of nodes and elements to use. Default = 2000
3. `scalfac`: Float64 option for the scaling factor when visualising the deformation. Default = 1e11
4. `interactive`: Integer in range 1:3:
	- `1`: Save generated figures without opening interactively. 
	- `2`: View figures interactively without saving. 
	- `3`: Both save and view figures interactively. 
5. `writefile`: Booleon option for whether to write the program output to a file 
6. `filename`: String option for file path of output file. Must point to a file, not directory. 
7. `figpath`: String option for file path of figures to be saved to. Must point to a directory, not file. 

***
Copyright © 2022, Dylan Morgan
