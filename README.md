# Calorimeter simulation for NNbar/HIBEAM experiment


<h3>How to Build</h3>

<h5>Within the directory that contains nnbar-calo-sim, run </h5>
  
$ mkdir nnbar-calo-sim-build <br>
$ cd nnbar-calo-sim-build  <br>
$ cmake -DGeant4_DIR=path_to_Geant4_installation/lib[64]/Geant4-10.0.0/ ../nnbar-calo-sim <br>
$ make <br>              

<h5>You made need to run </h5>

$ source path_to_Geant4_installation/bin/geant4.sh 
