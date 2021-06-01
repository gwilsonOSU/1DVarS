% Parametric Beach Modeling
% 
% 1D Barred Beach Model:
%   make1DBeachEngine   - given inputs, compute full 1D beach
%   makeBackgroundProfileComposite - make combo linear-exp h0
%   makeBarredBeach     - create 1D barred profile from given h0
%   findKGammaComposite - helpful function for makeBarredBeach
%
% 2D Barred Beaches:
%   make2DBeachEngine   - given inputs, compute full 2D beach for a grid
%   sortInputs          - sort shore, bar and deep inputs by increasing y
%   smoothInputs        - cubic spine smooth of shore, bar and deep data
%   designGrid          - default design grid to span the input data
%   findShoreNormalSmooth - finds the shore normal for any offshore x,y
%   findBarAndDeep      - find intersection of transect with bar and deep
%   findDepthOnTransect - find actual depth at x,y on transect
%
% Examples
%   example1DCase       - compute and compare transect to Duck survey
%   testCase1D.mat      - data for example1DCase
%   example2DCase       - compute and compare model to Duck CRAB survey.
%   testCase2D.mat      - data for example2DCase
%
% Docs
%   README              - tips on using this package
%   Ruessink2003.pdf    - original paper on parametric sand bar forms
%   HLEV14.pdf          - 2014 Coastal Engineering paper for 1D version
%   HLEV16.pdf          - 2016 Coastal Engineering paper for 2D version
