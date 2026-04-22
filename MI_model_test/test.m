this_dir = fileparts(mfilename('fullpath'));   % 当前脚本所在文件夹
project_root = fileparts(this_dir);            % 上一级：project_code
addpath(genpath(fullfile(project_root,'HOIs_tools')));

X = table2array(data);

I = gcmi_cc(X(:,1),X(:,2));

