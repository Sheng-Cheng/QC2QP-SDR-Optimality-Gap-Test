% QC2QP_SDR_optimalityGap_test_example.m
% This is an examplary script that shows users how to call the function 
% QC2QP_SDR_optimalityGap_test to conduct an optimality gap test
% test on a quadratic program with two quadratic constraints (QC2QP). 

% Feel free to remove the comment on trial example 1-10 to see the output
% of these worked out examples. More data ara available in 
% 'QC2QP_dataset.m'. You are also welcome to try your own problem data.

% Sheng Cheng (cheng@terpmail.umd.edu), Nov. 2018.


clear all;clc;
% Trial example 1: nonconvex QC2QP with no optimality gap
    M0 = [ 0    -2     0;
          -2    -1    -2;
           0    -2     1];
   
M1 = [-2     3     2;
       3    -3     1;
       2     1     1];
   
M2 = [ 4    -1     5;
      -1     4     5;
       5     5     3];  
% =========================

% Trial example 2: nonconvex QC2QP with no optimality gap
% M0 = [ 0    -2     0;
%       -2    -1    -2;
%        0    -2     1];
%    
% M1 = [-2     3     2;
%        3     3     1;
%        2     1    -2];
%    
% M2 = [ 4    -1     5;
%       -1     4     5;
%        5     5     1];
% =========================

% Trial example 3: dual infeasible
% M0 = [ 0    -2     0;
%       -2    -1    -2;
%        0    -2     1];
%    
% M1 = [-2     3     2;
%        3    -3     1;
%        2     1     1];
%    
% M2 = [ 4    -1     5;
%       -1     4    -5;
%        5    -5    -1];
% =========================

% Trial exmaple 4: semidefinite relaxation infeasible
% M0 = [ 0     5     2;
%        5     5     3;
%        2     3     1];
%    
% M1 = [ 2    -2     0;
%       -2     5    -3;
%        0    -3     5];
%  
% M2 = [-4     2     2;
%        2    -5    -1;
%        2    -1    -1];
% =========================

% Trial example 5: nonconvex QC2QP with no optimality gap (n = 3)
% M0 = [14    11     7    11;
%       11    13     9     8;
%        7     9     9    16;
%       11     8    16     3];
% 
% M1 = [1     6    13     9;
%       6    14     9     8;
%      13     9    10    16;
%       9     8    16    13];
% 
% M2 = [12     9    10     2;
%        9     7    10    15;
%       10    10    19     9;
%        2    15     9    14];
% ========================= 
   
% Trial example 6: nonconvex QC2QP with no optimality gap (n = 4)   
% M0 = [13    11    14    13     6;
%       11    13     9     9     8;
%       14     9    18     5     4;
%       13     9     5    12    14;
%        6     8     4    14    11];
%    
% M1 = [15     9     9    12    15;
%        9     7     5     3    13;
%        9     5     7    11     8;
%       12     3    11    18     7;
%       15    13     8     7    16];
% 
% M2 = [20    12     3    11    11;
%       12    20    15     9    15;
%        3    15     2    12    15;
%       11     9    12    10    11;
%       11    15    15    11    19];
% ========================= 

% Trial example 7: nonconvex QC2QP with no optimality gap (n = 5)   
% M0 = [9     2    14     2     8    18;
%       2    16     4    12     9     1;
%      14     4    18     7     5     4;
%       2    12     7    13     8    11;
%       8     9     5     8     9     4;
%      18     1     4    11     4    19];
% 
% M1 = [7    12    10    13    10    10;
%      12    17     8     5     4    14;
%      10     8     1     9    12     4;
%      13     5     9     9    12     7;
%      10     4    12    12     4    13;
%      10    14     4     7    13     1];
% 
% M2 = [9     2    15     9     4    10;
%       2    14    11    14    10    10;
%      15    11    18    18     4     4;
%       9    14    18     4     4    11;
%       4    10     4     4     6     2;
%      10    10     4    11     2    10];
% ========================= 

% Trial example 8: nonconvex QC2QP with an optimality gap (n = 3)
% M0 = [7     9     8     8;
%       9    14    11     6;
%       8    11    15    15;
%       8     6    15    14];
%  
% M1 = [6    13     8     9;
%      13     9     5    14;
%       8     5    17     8;
%       9    14     8    15];
%  
% M2 = [16     2    13     7;
%        2    13     7    12;
%       13     7     8    14;
%        7    12    14    12];
% =========================

% Trial example 9: nonconvex QC2QP with an optimality gap (n = 4)
% M0 = [7    13    12    10    14;
%      13    16    14     8    12;
%      12    14    20     8     8;
%      10     8     8    19    15;
%      14    12     8    15    11];
% 
% M1 = [ 7    10     8    13     5;
%       10     6     4    12    13;
%        8     4    10    14     5;
%       13    12    14     5     6;
%        5    13     5     6    19];
% 
% M2 = [8    11    11     1     9;
%      11    12     7    13     6;
%      11     7     4     6    11;
%       1    13     6    13     6;
%       9     6    11     6    15];
% ========================= 

% Trial example 10: nonconvex QC2QP with an optimality gap (n = 5)
% M0 = [3    13     7    13    12    10;
%      13    18     7    13    16    11;
%       7     7    18     7     9    17;
%      13    13     7     9     5    10;
%      12    16     9     5     6    14;
%      10    11    17    10    14    11];
% 
% M1 = [14     4    17     6    15    10;
%        4    11    14    16    16     8;
%       17    14    17     7    10     3;
%        6    16     7     7     8    17;
%       15    16    10     8    14    11;
%       10     8     3    17    11     5];
% 
% M2 = [4    15    12    12     9     8;
%      15     7     4     5     2     8;
%      12     4    12     8    14     9;
%      12     5     8    12     6    11;
%       9     2    14     6    16    10;
%       8     8     9    11    10    20];
% ========================= 
% More trial examples:
% load('QC2QP_dataset.mat');
% ========================= 

epsilon2 = 1e-5;

status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2,1);