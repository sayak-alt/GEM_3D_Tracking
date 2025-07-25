Scripts to do cosmic muon tracking through three layers of GEM chambers (using APV+MPD based readout)
Pedcal.C calculates the pedestal from the individual APV cards
Sorting.C calculates the hit and charge information from the individiaul GEMs
GEMs_residual_calculation.C sortes the events passing through all the GEMs
GEMs_residual_x and y.C calculates the GEM residuals in X & Y direction (offset correction)
muon_tracker.C does the 3D track reconstruction of the muon events
