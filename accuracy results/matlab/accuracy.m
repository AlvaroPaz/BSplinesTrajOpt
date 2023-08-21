clc; clear all; close all;

old_analytic_a_c = importdata('old_analytic/a_analytic/cData.txt');
old_bfgs_a_c = importdata('old_bfgs/a_bfgs/cData.txt');

sum(sum(abs(old_analytic_a_c-old_bfgs_a_c)))





old_analytic_a_c = importdata('cData.txt');
old_bfgs_a_c = importdata('cData_.txt');

sum(sum(abs(old_analytic_a_c-old_bfgs_a_c)))
