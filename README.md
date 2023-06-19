
# Hybrid DODE

A hybrid modelling framework for the estimation of dynamic origin-destination flows

This repository keeps the code base of hybrid dynamic origin-destination flow estimation (hybrid-DODE) method proposed in [LINK of the paper]

The novel hybrid DODE framework introduces synchronous modelling of region-level and centroid-level traffic dynamics into DODE prolems. 

The region-level traffic flows are described by the MFD, while the centroid-level traffic flows are represented by the linear mapping of origin-destination flows onto link counts.
### Concept

The hybrid DODE problem presented in this study introduces the hybrid modelling of traffic dynamics. This is accomplished by incorporating a centroid-level linear approximation and a region-level traffic model which operates simultaneously. 

<p align="center">
<img src="BilvFml.png " width="50%" height="50%">
</p>

Figure shows a conceptual diagram of the hybrid DODE framework. We have hybrid DODE as the upper level and traffic assignment as the lower-level. Note that there are three traffic models in this framework. The traffic assignment level is a simulation based model which operates closer to real traffic networks, while there are two other traffic models being used in the OD estimation level. These two models are analytical models based on the hierarchy of operation. The centroid-level is modelled as a linear traffic approximation model which maps the OD flows onto link counts via an assignment matrix.
### Content

This repositiory contains following functions that could be used to execute the hybrid-DODE in any urban network 

- **Main_script.m** : Main script that execute the Hybrid DODE
- **run_simulator_CN.m** : script that execute the traffic assignmnet and retrieve centroid-level data
- **run_simulator_RG.m** : script that execute the traffic assignmnet and retrieve region-level data
- **build_ODestimation_combined3.m** : CasADI based script that develop the hybride DODE optimisation problem
- **solve_ODestimation_combined3.m** : CasADI based script that solve the hybride DODE optimisation problem
- **purturb_OD.m** : script that purthurb ground-truth OD matrices
- **CalibrateMFD.m** : The script to calibrate MFD parameters acrroding to traffic conditions
- **build_calibration2.m** : CasADI based script that develop the calibraition of the regional parameters.
- **solve_calibration2.m** : CasADI based script that solvesthe calibraition optimization problem built by  build_calibration2.m.
### Documnetation

An implementation of the proposed Hybrid DODE framework is documented and given in :
[Link to Paper after publicatoion] 
### Software Requirements

The hybrid DODE framework requires three main software categories that could be listed as:

- a programming language
- traffis assignment software 
- non-linear, non-convex optimization tool

This repository and associated scripts are developed to execute along with following software tools.

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
