function [epsilon] = strain(PhysicalChemicalDatasetCombo,Porosity,taubulk)
%chemloss_factor function
% by R.E. Breunig
% last update: 2020.11.11
% This function is based on the Hayes 2019 ( DOI: 10.1126/sciadv.aao0834)
% equation for weathered porosity. The form of this function here in
% particular is capable of calculating  the strain, epsilon, of a
% sample from the porosity and tau_bulk of a sample. Epsilon is a unitless
% value. Negative values of epsilon indicate compression while positive
% values indicate augmentation.
% Strain relies on comparison of weathered material to protolith and
% normalizing change in bulk density of a sample against immobile element
% depletion/augmentation. The choice of immobile element and protolith are
% implicit in this function via taubulk. 

%
%
%Summary of inputs:
%Physical&ChemicalDataset
%   The input for Physical&ChemicalDataset is a struct that must contain columns
%   'AvgDepthBelowGroundSurface_m' and a column of porosity.
%Porosity
%   The input for porosity should be the string title of the column storing
%   porosity (in decimal form) from the Physical&ChemicalDataset struct.
%taubulk
%taubulk should be a vertical vector input of taubulk calculated with respect to an
%immobile element and protolith type fit for the landscape.



epsilon= ((taubulk+1)./(1-PhysicalChemicalDatasetCombo.(Porosity)))-1;
if abs(epsilon)>2 
    warning{'Epsilon value exceeds absolute value of 2-- check inputs'}
else 
    epsilon;

end

