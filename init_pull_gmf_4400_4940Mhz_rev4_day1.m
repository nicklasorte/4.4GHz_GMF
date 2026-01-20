close all;
close all force;
clc;
app=NaN(1);  %%%%%%%%%This is to allow for Matlab Application integration.
format shortG
top_start_clock=clock;
folder1='C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\4.4GHz';
cd(folder1)
addpath(folder1)
pause(0.1)
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\GMF_functions')
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Basic_Functions')
addpath('C:\Users\nlasorte\OneDrive - National Telecommunications and Information Administration\MATLAB2024\Census_Functions')
pause(0.1)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%Split GMF into Service Types and Agency
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%Step 1: Pull GMF and do light cleaning.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GMF inputs:
cd(folder1)
pause(0.1)
gmf_MinMHz=4400;
gmf_MaxMHz=4940;
rev_num=01202026;
[gmf_table]=pull_gmf_excel_onedrive_rev2(app,gmf_MinMHz,gmf_MaxMHz,rev_num);
gmf_header=gmf_table.Properties.VariableNames;
cell_gmf=table2cell(gmf_table);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%gmf_header IDX
col_service_idx=find(matches(gmf_header,'Service'));
col_agency_idx=find(matches(gmf_header,'Agency'))
city_col_idx=find(matches(gmf_header,'XAL'))
col_gmf_ser_idx=find(matches(gmf_header,'SER'))
col_nts_idx=find(matches(gmf_header,'NTS'));
col_latdd_idx=find(matches(gmf_header,'XLatDD'));
col_londd_idx=find(matches(gmf_header,'XLonDD'));
col_freq1_idx=find(matches(gmf_header,'FRQMHz'));
col_ems_idx=find(matches(gmf_header,'EMS'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cut
%%%%%%%%%%%%%Filter out CAN and NG
ng_idx=find(contains(cell_gmf(:,col_gmf_ser_idx),'NG'));  %%%%%%Non-Government
cell_gmf(ng_idx,:)=[];
can_idx=find(contains(cell_gmf(:,col_gmf_ser_idx),'CAN'));   %%%%Cananda
cell_gmf(can_idx,:)=[];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Simplifiy the Agency Name in the GMF
[cell_gmf]=simplify_gmf_agency_name_rev1(app,cell_gmf,gmf_header);
%%%%%%%%%%Unique Agency Types
unique(cell_gmf(:,col_agency_idx))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Unique other strings in the GMF and remove N/A  or NaN
header_string='XAG';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='XAH';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='STC';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='PWR';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='RAG';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='RAH';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='RAL';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='RSC';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='Service';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

header_string='EMS';
[cell_gmf]=unique_and_remove_nan_rev1(app,gmf_header,header_string,cell_gmf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Clean up the Tx EUT name and fill in holes
[cell_gmf]=clean_gmf_tx_eut_rev1(app,gmf_header,cell_gmf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Clean up the Rx EUT name and fill in holes
[cell_gmf]=clean_gmf_rx_eut_rev1(app,gmf_header,cell_gmf);

%%%%%%%%%Clean up the Lat/Lon, if missing --> NAN
[cell_gmf]=clean_gmf_tx_latlon_rev1(app,gmf_header,cell_gmf);


size(cell_gmf)

horzcat(gmf_header',cell_gmf(1,:)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sort the GMF into Service Types and Agency
service_type='Aeronautical Mobile'
print_excel_service_agency_rev1(app,cell_gmf,service_type,gmf_header,rev_num)

service_type='Maritime Mobile'
print_excel_service_agency_rev1(app,cell_gmf,service_type,gmf_header,rev_num)

service_type='Land Mobile'
print_excel_service_agency_rev1(app,cell_gmf,service_type,gmf_header,rev_num)

service_type='Fixed'
print_excel_service_agency_rev1(app,cell_gmf,service_type,gmf_header,rev_num)







end_clock=clock;
total_clock=end_clock-top_start_clock;
total_seconds=total_clock(6)+total_clock(5)*60+total_clock(4)*3600+total_clock(3)*86400;
total_mins=total_seconds/60;
total_hours=total_mins/60;
if total_hours>1
    strcat('Total Hours:',num2str(total_hours))
elseif total_mins>1
    strcat('Total Minutes:',num2str(total_mins))
else
    strcat('Total Seconds:',num2str(total_seconds))
end

'DONE'
cd(folder1)
