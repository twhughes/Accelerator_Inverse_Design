# DLA Optimization Software

## Contents
This matlab program is used to generate optimized dielectric laser accelerator structures.  
The code produces two separate structures:
- optimized for maximum acceleration gradient
- optimized for maximum acceleration gradient divided by the maximum electric field amplitude (either in material or design region).

## Paper information
This code was used in work related to this paper:

[Method for Computationally Efficient Design of Dielectric Laser Accelerator Structures](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-13-15414 "Method for Computationally Efficient Design of Dielectric Laser Accelerator Structures")

## How to run simulation
Just run the script "optimization.m".  The FDFD code and sparse matrix LU-factoring code is located in "dependencies/".
The entire directory needs to be added to the path to work correctly

## Citing our work

If you would like to use this code, please cite the paper:


        @article{Hughes:17,
            author = {Tyler Hughes and Georgios Veronis and Kent P. Wootton and R. Joel England and Shanhui Fan},
            journal = {Opt. Express},
            keywords = {Gratings; Optical devices; Subwavelength structures},
            number = {13},
            pages = {15414--15427},
            publisher = {OSA},
            title = {Method for computationally efficient design of dielectric laser accelerator structures},
            volume = {25},
            month = {Jun},
            year = {2017},
            url = {http://www.opticsexpress.org/abstract.cfm?URI=oe-25-13-15414},
            doi = {10.1364/OE.25.015414},
        }

