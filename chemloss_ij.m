function [tau_ij] = chemloss_ij(Geochem,MobileElementJ,ImmobileElementI,ProtolithType,XofLowermost)
%chemloss_ij function
% by R.E. Breunig
% last update: 2020.11.11
% This function is a take from Brimhall & Dietrich 1987(DOI:
% 10.1016/0016-7037(87)90070-6) and detailed description behind the general
% theory of chemical depletion can be found there. The form of the function here in
% particular is capable of calculating chemical depletion, tau, of a mobile
% element(j) in respect to an immobile element (i). Tau_i;j is a unitless
% value where negative values indicate depletion and positive values
% indicate augmentation.
%Common immobile elements in literature for calculating tau_i;j are Zr, Ti,
% & Nb. Choice of immobile should be specific to site.
%
%
%Summary of inputs:
%Geochem
%   The input for Geochem is a struct that must contain columns
%   'AvgDepthBelowGroundSurface_m', 'HoleID', and user choice of
%   ImmobileElementI and MobileElementJ.
%MobileElementJ
%   The input of MobileElementJ should be the title of the column
%   corresponding to the mobile element of interest in the form of a
%   string. (ex. 'SiO2_pct')
%ImmobileElementI
%   The input of ImmobileElementI should be the title of the column
%   corresponding to the immobile element of interest in the form of a
%   string. (ex. 'TiO2_ppm')
%ProtolithType
%   Chemical depletion and the calculation of tau_i;j relies on the
%   comparison of regolith (i.e. weathered material) to protolith (i.e.
%   unweathered material). This function has 4 valid string inputs for
%   ProtolithType: 'TotalLowermost', 'AverageLowermostPerHole',
%   'LowermostPerHole', and 'XLowermostAveraged'.
%       'TotalLowermost' assumes that the single deepest well in the sampled
%       landscape is representative of protolithic material for all holes and
%       depths.
%       'AverageLowermostPerHole' assumes that the depth to which each well
%       is drilled to is significant and represents a lower boundary of
%       interest. This calculation averages the lowermost geochemical data
%       i and j of each borehole to represent a single protolith across the
%       full landscape.
%       'LowermostPerHole' assumes that each hole has a unique parent
%       lithology. This calculation finds the deepest interval of each
%       borehole, which represents the protolith for it's specific hole. 
%       'XLowermostAveraged' assumes that the deepest intervals across the
%       landscape represent protolithic material for all holes and depths,
%       but this calculation in particular might be best suited for
%       lithology with slight heterogeneity. 
%XofLowermost
%   XofLowermost input is only required if ProtolithType is
%   'XLowermostAveraged'. XofLowermost must be a positive integer and
%   indicates the number of depth intervals averaged together to form the
%   protolith. ex. XofLowermost=3 indicates the 3 lowest samples in the
%   landscape averaged together represent the protolithic material.




if contains(ProtolithType,'TotalLowermost')
    [~,indexmaxd]= max(Geochem.('AvgDepthBelowGroundSurface_m')); %find 
    iprotolith=(Geochem.(ImmobileElementI)(indexmaxd));
    jprotolith=(Geochem.(MobileElementJ)(indexmaxd));
    tau_ij=(Geochem.(MobileElementJ)/jprotolith).*(iprotolith./Geochem.(ImmobileElementI))-1;
    if tau_ij >= 5
         warning('Value is high inidicating high augmentation. Review inputs.')
    else tau_ij;
    end 
    
    

    
elseif contains(ProtolithType,'AverageLowermostPerHole')
    HoleIDList=unique(Geochem.('HoleID')); %HoleIDList gives all unique names 
    for i=1:numel(HoleIDList) %%this loop finds the row of the deepest interval @ each hole
     rowsofinterest=find(contains(Geochem.('HoleID'),(HoleIDList(i))));
     [~,index_inrowsofinterest]=max(Geochem.('AvgDepthBelowGroundSurface_m')(rowsofinterest));
     ProtolithRow_inGeochem=rowsofinterest(index_inrowsofinterest);
     DeepestPerHole(i,:)=ProtolithRow_inGeochem; 
    end
    for i=1:numel(DeepestPerHole) %%this loop finds the deepest immobile and mobile chemistry measurements for each hole
    deepestiperhole(i,:)=Geochem.(ImmobileElementI)(DeepestPerHole(i));
    deepestjperhole(i,:)=Geochem.(MobileElementJ)(DeepestPerHole(i));
    end
    iprotolith=mean(deepestiperhole); %%immobile element chemistry averaged across the deepest interval of each hole
    jprotolith=mean(deepestjperhole);%%mobile element chemistry averaged across the deepest interval of each hole
    tau_ij=(Geochem.(MobileElementJ)/jprotolith).*(iprotolith./Geochem.(ImmobileElementI))-1;
    if tau_ij >= 5
         warning('Value is high inidicating high augmentation. Review inputs.')
    else tau_ij;
    end 
    
    
    
    
elseif contains(ProtolithType,'LowermostPerHole')
   HoleIDList=unique(Geochem.('HoleID')); %HoleIDList gives all unique names 
   for i=1:numel(HoleIDList) 
       rowsofinterest=find(contains(Geochem.('HoleID'),(HoleIDList(i)))); %%finds rows corresponding to particular hole
      [~,index_inrowsofinterest]=max(Geochem.('AvgDepthBelowGroundSurface_m')(rowsofinterest)); 
       ProtolithRow_inGeochem=rowsofinterest(index_inrowsofinterest); %%finds index of deepest interval of particular hole within full Geochem struct
       iprotolith=Geochem.(ImmobileElementI)(ProtolithRow_inGeochem);%%immobile element concentration of protolith at particular hole
       jprotolith=Geochem.(MobileElementJ)(ProtolithRow_inGeochem); %%mobile element concentration of protolith at particular hole
       tau_ij(min(rowsofinterest):max(rowsofinterest),:)=(Geochem.(MobileElementJ)(min(rowsofinterest):max(rowsofinterest))/jprotolith).*(iprotolith./Geochem.(ImmobileElementI)(min(rowsofinterest):max(rowsofinterest)))-1; 
       %%^for the hole of interest in each rep of the loop, this function
       %%calculates tau_ij where the protolith is the lowest depth in the
       %%hole.
   end   
   if tau_ij >= 5
         warning('Value is high inidicating high augmentation. Review inputs.')
         else tau_ij;
   end

   
 
elseif contains(ProtolithType,'XLowermostAveraged')
    [~,indexmaxd]= maxk(Geochem.('AvgDepthBelowGroundSurface_m'),XofLowermost);
    iprotolith=mean(Geochem.(ImmobileElementI)(indexmaxd));
    jprotolith=mean(Geochem.(MobileElementJ)(indexmaxd));
    tau_ij=(Geochem.(MobileElementJ)/jprotolith).*(iprotolith./Geochem.(ImmobileElementI))-1;
   if tau_ij >= 5
         warning('Value is high inidicating high augmentation. Review inputs.')
         else tau_ij;
   end

else 
    warning('Input valid string for ProtolithType')
end

end

