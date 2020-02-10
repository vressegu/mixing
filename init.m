

clear;
close all;
clc;
dbstop if error;

current_tests=[ pwd '/current_tests' ];
fct = genpath([ pwd '/functions' ]);
mains = [ pwd '/mains'];
addpath(pwd)
addpath(fct)
addpath current_tests mains;
clear current_tests mains fct
home;