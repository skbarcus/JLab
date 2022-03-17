The file target_data_SBS4.npy is a numpy array of labels indicating whether an event was an elastic event(1) or a non-elastic event (0).

The file training_data_SBS4.npy is the corresponding numpy array of the detector data before any normalization is applied for each event (normalization is done in the notebook). Each event has 6 variables.

Event Data Structure:
Element     Variable
0           HCal Cluster Energy [GeV]
1           HCal Cluster Time [ns]
2           HCal Cluster X-Position (vertical) [m]
3           HCal Cluster Y-Position (horizontal) [m]
4           BigBite Preshower Energy [GeV]
5           BigBite Shower Energy [GeV]
