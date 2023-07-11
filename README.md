
# Hybrid DODE

A hybrid modelling framework for the estimation of dynamic origin-destination flows

This repository includes the codes for the hybrid dynamic origin-destination flow estimation (hybrid-DODE) method proposed in [LINK of the paper]

The novel hybrid DODE framework introduces synchronous modelling of region-level and centroid-level traffic dynamics into DODE prolems. 

The region-level traffic flows are described by the MFD, while the centroid-level traffic flows are represented by the linear mapping of origin-destination flows onto link counts.
### Concept

The hybrid DODE problem presented in this study introduces the hybrid modelling of traffic dynamics. This is accomplished by incorporating a centroid-level linear approximation and a region-level traffic model which operates simultaneously. 

<p align="center">
<img src="BilvFml.png " width="40%" height="40%">
</p>

Figure shows a conceptual diagram of the hybrid DODE framework. We have hybrid DODE as the upper level and traffic assignment as the lower-level. Note that there are three traffic models in this framework. The traffic assignment level is a simulation based model which operates closer to real traffic networks, while there are two other traffic models being used in the OD estimation level. These two models are analytical models based on the hierarchy of operation. The centroid-level is modelled as a linear traffic approximation model which maps the OD flows onto link counts via an assignment matrix.
### Content

This repositiory contains following functions that could be used to execute the hybrid-DODE in any urban network 

- **Main_script.m** : Main script that executes the Hybrid DODE.
- **run_simulator_CN.m** : script that executes the traffic assignment and retrieves centroid-level data from the simulator.
- **run_simulator_RG.m** : script that executes the traffic assignmnet and retrieves region-level data.
- **build_ODestimation_combined3.m** : CasADI based script that builds the hybride DODE optimisation problem.
- **solve_ODestimation_combined3.m** : CasADI based script that solves the hybrid DODE optimisation problem.
- **purturb_OD.m** : script that perturbs ground-truth OD matrices.
- **CalibrateMFD.m** : The script to calibrate MFD parameters according to traffic conditions.
- **build_calibration2.m** : CasADI based script that builds the calibration problem for the regional traffic model.
- **solve_calibration2.m** : CasADI based script that solves the calibrations problem.
### Documentation

An implementation of the proposed Hybrid DODE framework is documented and given in :
[Link to Paper after publication] 
### Software Requirements

The hybrid DODE framework requires three main tools that could be listed as:

- a programming language
- traffic assignment software
- non-linear, non-convex optimization solver

This repository and its associated scripts are designed to work with the software tools listed below.

| Software     | Version Required |
| :---         |    :----:        |  
| MATLAB       |    2020b        |
| AIMSUN       | 8.4              |
| CasADi       | 5.5             |


### Acknowledgment

This research is funded by iMOVE CRC and supported by the Cooperative Research Centres program, an Australian Government initiative. This research is also partly funded by the Australian Research Council (ARC) through Dr. Mehmet Yildirimoglu's Discovery Early Career Researcher Award (DECRA; DE220101320).
### License

[MIT](https://choosealicense.com/licenses/mit/)


### Main contact

For more information on the hybrid DODE, please contact : 

- Dr Mehmet Yildirimoglu (m.Yildirimoglu@uq.edu.au) 
- Dr Sakitha Kumarage (s.kumarage@uq.edu.au) 
