Contents under \EighteenStoryStructure:

The 18-story example is based on a one-third scale 18-story steel frame structure tested at E-defense in Japan.  Model updating is performed entirely through simulation models, and then with the experimental data. Detailed descriptions of the structure and model updating studies can be found in the companion monograph, ["Formulation and application of SMU - an open-source MATLAB package for structural model updating."](https://github.com/ywang-structures/Structural-Model-Updating/blob/master/Formulation%20and%20application%20of%20SMU%20%E2%80%93%20an%20open-source%20MATLAB%20package%20for%20structural%20model%20updating.pdf)

1. 18-Story Steel Structure - Simulated Data

A study is first performed entirely in simulation to demonstrate the model updating procedure of the 18-DOF shear-frame model. The following table shows the masses and the initial/nominal inter-story stiffness values. To construct an “actual” structure in the simulation study, the table also shows the "actual" inter-story stiffness values. Modal properties of the “actual” structure are directly used as “experimental” modal properties. It is assumed that only seven out of the eighteen DOFs (#1, #4, #7, #9, #12, #15, and #18) are instrumented with sensors and only the lowest four modes are “measured” and available for model updating.

<img src="https://github.com/ywang-structures/Structural-Model-Updating/blob/master/Examples/EighteenStoryStructure/Figures/18DOF_1.png?raw=true" width="700" height="550" />



2. 18-Story Steel Structure - Experimental Data

Model updating is also performed using the experimental data from shaking table tests of the one-third scale 18-story steel frame structure. Floor masses are accurately weighed, while all 18 inter-story stiffness values require updating. The structure is mounted on a shake table which provides a uniaxial base excitation for the structure. To make the problem more challenging, it is assumed that sensor data from only six DOFs (#3, #6, #9, #12, #15, and #18) are available. Modal properties of the 18-story structure are extracted from the experimental data. The first and the second extracted modes are used for model updating.


<img src="https://github.com/ywang-structures/Structural-Model-Updating/blob/master/Examples/EighteenStoryStructure/Figures/18DOF_2.png?raw=false" width="250" height="300" />
Photo credit: The National Research Institute for Earth Science (NIED), Japan
