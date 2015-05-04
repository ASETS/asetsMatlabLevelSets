## [ASETS](http://www.advancedsegmentationtools.org) LevelSet repository 

### License:  
BSD (see license.md)

### Developers:
- Martin Rajchl (@mrajchl), Imperial College London (UK)
- John SH. Baxter (@jshbaxter), Robarts Research Institute (CAN)
- Jing Yuan, Robarts Research Institute (CAN)

### Features: 
- Fast time-implicit levelsets in 2D/3D
    - Binary 
    - Coupled
    - Multi-phase
- Implemented in multiple languages
   - Matlab/mex/C
   - Matlab/CUDA
- Tutorials
    - T01 Time-implicit levelset propagation
    - T02 Binary levelset image segmentation
- Application examples for (medical) image segmentation:
    - 3D multi-phase levelset segmentation


### Overview of folder structure:   
*./*: Compile scripts, readme, license and todo list  
*./applications*: Contains examples of typical applications in image segmentation and analysis  
*./data*: Example data to run the applications  
*./lib*: Is created by compile.m and contains the compiled C/mex files  
*./maxflow*: asetsMaxFlow optimization code in C/mex and Matlab 
*./tutorials*: Contains available tutorials   

### Compile/Installation instructions:  
To compile the C/mex code run:
```matlab
compile.m
```
which creates the folder *./lib*. For testing purposes run any script in *./tests*.   

### Tests:  
- Matlab 2014a, 64-bit Linux (Ubuntu 12.04 LTS)  
- Matlab 2015a, 64-bit Windows 7
- Matlab 2012a, 32-bit WinXP

