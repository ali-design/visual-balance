%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author: Ali Jahanian: Ali.Jahanian.HCI@gmail.com                    %%%                                                                   %%%
%%% http://people.csail.mit.edu/jahanian/index.html                     %%%
%%% Version: 2.0                                                        %%%
%%% Date: 05/01/2014                                                    %%%
%%% For sharing with others, please first ask permission from           %%% 
%%% Ali Jahanian.                                                       %%%
%%% For your publication, please cite these two refernces:              %%%
%%% 1) Jahanian, A., Vishwanathan, S. V. N., and Allebach, J. P.,       %%%
%%% Learning visual balance from large-scale datasets of aesthetically  %%%
%%% highly rated images. In Proc. IS&T/SPIE Electronic Imaging,         %%%
%%% International Society for Optics and Photonics, (2015).             %%%
%%% 2) Jahanian A., Quantifying Aesthetics of Visual Design Applied to  %%%
%%% Automatic Design. Ph.D. Dissertation, Purdue University,            %%%
%%% West Lafayette, IN, December 2014.                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mixture] = process_GMM_Layout(dataset)

if isempty(dataset)
    fprintf('Loading the dataset ... \n');
    dataset = load('dataset');
    dataset = dataset.dataset;
end

K = 5;
ConditionNumber = 1e5;

fprintf('Start processing ... \n');
tic
[mixture] = GaussianMixture(dataset, K, ConditionNumber);
toc
% plotGMM(mixture);