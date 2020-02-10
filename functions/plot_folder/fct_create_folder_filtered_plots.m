function fct_create_folder_filtered_plots(model)
% Create folders to save plots and files
%

folder_simu = model.folder.folder_simu;
model.sigma.k_c = inf;

if ~ (exist([folder_simu '/files'],'dir')==7)
    mkdir([folder_simu '/files']);
end
if ~ (exist([folder_simu '/one_realization'],'dir')==7)
    mkdir([folder_simu '/one_realization']);
end
if ~ (exist([folder_simu '/Spectrum'],'dir')==7)
    mkdir([folder_simu '/Spectrum']);
end
