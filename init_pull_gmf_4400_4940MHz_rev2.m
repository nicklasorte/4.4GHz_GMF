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


% %%%%%%%%Pull in James p2p pairing
% excel_filename='4400 to 4940 Fixed P2P Assignments James.xlsx'
% mat_filename_str='James_p2p_pairing_rev1'
% tf_repull=0
% [cell_p2p_pair]=load_full_excel_rev1(app,mat_filename_str,excel_filename,tf_repull);
% size(cell_p2p_pair)
% cell_p2p_pair(1:3,:)'
% james_header=cell_p2p_pair(1,:)'
% 
% 'Need to fix the missing'
% pause;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%james_header IDX
% gmf_num_col_idx=find(matches(james_header,'Applicant_Serial_Number'))
% pair_col_idx=find(matches(james_header,'Paired_With_Serial_Number'))
% tx_lat1_col_idx=find(matches(james_header,'Tx_Antenna_Latitude_Decimal_Degrees'))
% tx_lon1_col_idx=find(matches(james_header,'Tx_Antenna_Longitude_Decimal_Degrees'))
% 
% cell_p2p_pair_data=cell_p2p_pair(:,[gmf_num_col_idx,pair_col_idx,tx_lat1_col_idx,tx_lon1_col_idx])
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 'Next step, plot the p2p links (red lines) with the 50 states outline, use James since it has pairings'
% 
% 
% 'start here'
% pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GMF inputs:
cd(folder1)
pause(0.1)
gmf_MinMHz=4400;
gmf_MaxMHz=4940;
rev_num=07182025;
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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%%%%%%%%%%If there is more than 1 entry for the Receiver location, break it into separate rows.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%Check: RLA	RLG
%%%%[cell_gmf]=expand_gmf_rx_loc_rev1(app,gmf_header,cell_gmf);

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



%%%%%%%%Cut the 'Fixed' from the list
mobile_row_idx=find(contains(cell_gmf(:,col_service_idx),'Mobile'));
cell_gmf_mobile=cell_gmf(mobile_row_idx,:);
size(cell_gmf_mobile)
table_mobile_gmf=cell2table(cell_gmf_mobile);
table_mobile_gmf.Properties.VariableNames=gmf_header;
%writetable(table_mobile_gmf,strcat('Mobile.xlsx'));



%%%%%%%%%%%%%%%%%%%%Find the different notes
%%%%%%%S296--Not to preclude assignment of this frequency to other agencies at specific locations.
cell_temp_gmf=cell_gmf_mobile;
temp_cell_row_notes=cell_temp_gmf(:,col_nts_idx);
substring = 'S296';
logical_indices = cellfun(@(x) ischar(x) && contains(x, substring), temp_cell_row_notes);
note_idx=find(logical_indices);
note_table=cell2table(cell_temp_gmf(note_idx,:));
note_table.Properties.VariableNames=gmf_header;
%writetable(note_table,strcat('MOBILE_Note_Table_',substring,'.xlsx'));

%%%%%%%%%%%%%%%%%%%%Find the different notes
%%%%%%%P032 NIB
cell_temp_gmf=cell_gmf_mobile;
temp_cell_row_notes=cell_temp_gmf(:,col_nts_idx);
substring = 'P032';
logical_indices = cellfun(@(x) ischar(x) && contains(x, substring), temp_cell_row_notes);
note_idx=find(logical_indices);
note_table=cell2table(cell_temp_gmf(note_idx,:));
note_table.Properties.VariableNames=gmf_header
%writetable(note_table,strcat('MOBILE_Note_Table_',substring,'.xlsx'));



%%%%Split, Unique, Combine
uni_mobile_type=unique(cell_gmf_mobile(:,col_service_idx))
split_cell=cellfun(@(x) strsplit(x, ','), uni_mobile_type, 'UniformOutput', false);
expand_split_cell=horzcat(split_cell{:});
uni_mobile=unique(expand_split_cell)'


%%%%%%%%%%Split GMF into these sub types: uni_mobile
%%%Remove the "Fixed"
temp_fix_idx=find(contains(uni_mobile,'Fixed'))
uni_mobile{temp_fix_idx}=[];
uni_mobile=uni_mobile(~cellfun('isempty',uni_mobile))
num_mob_type=length(uni_mobile)
cell_split_mobile_gmf=cell(num_mob_type,2);

for i=1:1:num_mob_type
    uni_mobile{i}
    temp_idx=find(contains(cell_gmf_mobile(:,col_service_idx),uni_mobile{i}));
    sub_cell_gmf_mobile=cell_gmf_mobile(temp_idx,:);
    cell_split_mobile_gmf{i,1}=uni_mobile{i};
    cell_split_mobile_gmf{i,2}=sub_cell_gmf_mobile;

    table_sub_mobile_gmf=cell2table(sub_cell_gmf_mobile);
    table_sub_mobile_gmf.Properties.VariableNames=gmf_header;
    writetable(table_sub_mobile_gmf,strcat('Sub_',uni_mobile{i},'.xlsx'));

    %%%%Check
    %unique(sub_cell_gmf_mobile(:,col_service_idx))

end
cell_split_mobile_gmf

size(cell_split_mobile_gmf)

tf_mobile_tree=0
if tf_mobile_tree==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Make the square tree map

    [num_rows,~]=size(cell_split_mobile_gmf)
    for i=1:1:num_rows
        cell_square_gmf=cell_split_mobile_gmf{i,2};
        str_service_type=cell_split_mobile_gmf{i,1}
        string_title=strcat(str_service_type,':','4400-4940MHz')
        filename1=strcat(str_service_type,'_4400-4940MHz.png')
        treemap_gmf_rev1(app,cell_square_gmf,col_agency_idx,string_title,filename1)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%'Next Step: Make the square for the sub mobile with the frequency just between 4800-4940MHz'
    [num_rows,~]=size(cell_split_mobile_gmf)
    for i=1:1:num_rows
        cell_square_gmf=cell_split_mobile_gmf{i,2};
        str_service_type=cell_split_mobile_gmf{i,1}

        %%%Do it for 4.8-4.94GHz
        array_freq=cell2mat(cell_square_gmf(:,col_freq1_idx));
        above_idx=find(array_freq>=4800);
        size(cell_square_gmf)
        cell_square_gmf=cell_square_gmf(above_idx,:);
        size(cell_square_gmf)
        string_title=strcat(str_service_type,':','4800-4940MHz')
        filename1=strcat(str_service_type,'_4800-4940MHz.png')
        treemap_gmf_rev1(app,cell_square_gmf,col_agency_idx,string_title,filename1)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end






% % 'Next Step:Try to find the channels used, Bandwidth and Frequency, and Histogram those for the different sub mobile gmf'
% % pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Separate the Fixed
%%%%%%%%Cut the 'Fixed' from the list
fixed_row_idx=find(matches(cell_gmf(:,col_service_idx),'Fixed'));
size(fixed_row_idx)
non_fixed_row_idx=find(matches(cell_gmf(:,col_service_idx),'Fixed')==0);
size(non_fixed_row_idx)

cell_gmf_non_fixed=cell_gmf(non_fixed_row_idx,:);
cell_gmf_fixed=cell_gmf(fixed_row_idx,:);
size(cell_gmf)
size(cell_gmf_fixed)  %%%%4215 fixed
size(cell_gmf_non_fixed)

%%%%%'Separate the fixed that have a radius of operation. --> These are transportables'
col_txrad_idx=find(matches(gmf_header,'TxRad'));
col_rxrad_idx=find(matches(gmf_header,'RxRad'));
cell_gmf_fixed_radius=cell_gmf_fixed(:,[col_txrad_idx,col_rxrad_idx]);
non_empty_idx=find(~cellfun(@isempty,cell_gmf_fixed_radius));
[row_indices, col_indices]=ind2sub(size(cell_gmf_fixed_radius), non_empty_idx);
uni_non_emt_idx=unique(row_indices);
cell_gmf_fixed_transportable=cell_gmf_fixed(uni_non_emt_idx,:);
size(cell_gmf_fixed_transportable)
ft_gmf=cell2table(cell_gmf_fixed_transportable);
ft_gmf.Properties.VariableNames=gmf_header;
%writetable(ft_gmf,strcat('Fixed_Transportable.xlsx'));


%%%%%%%%%%%%%%%%%%%%Find the different notes
%%%%%%%S296--Not to preclude assignment of this frequency to other agencies at specific locations.
cell_temp_gmf=cell_gmf_fixed_transportable;
temp_cell_row_notes=cell_temp_gmf(:,col_nts_idx);
substring = 'S296';
logical_indices = cellfun(@(x) ischar(x) && contains(x, substring), temp_cell_row_notes);
note_idx=find(logical_indices);
note_table=cell2table(cell_temp_gmf(note_idx,:));
note_table.Properties.VariableNames=gmf_header;
%writetable(note_table,strcat('Transportable_Note_Table_',substring,'.xlsx'));



%%%%%%%Need to find the idx of the fixed p2p
[num_fix_rows,~]=size(cell_gmf_fixed)
full_fixed_idx=[1:1:num_fix_rows]';
fixed_p2p_idx=setxor(full_fixed_idx,uni_non_emt_idx);
% % size(fixed_p2p_idx)
% % size(full_fixed_idx)
% % size(uni_non_emt_idx)
pre_cell_gmf_fixed_p2p=cell_gmf_fixed(fixed_p2p_idx,:);


%%%Cut those with no Lat/Lon
pre_cell_gmf_fixed_p2p_latlon=pre_cell_gmf_fixed_p2p(:,[col_latdd_idx,col_londd_idx]);
array_isnan=any(cell2mat(cellfun(@isnan,pre_cell_gmf_fixed_p2p_latlon, 'UniformOutput', false)),2);
non_nan_idx=find(array_isnan==0);
cell_gmf_fixed_p2p=pre_cell_gmf_fixed_p2p(non_nan_idx,:);
table_p2p=cell2table(cell_gmf_fixed_p2p);
table_p2p.Properties.VariableNames=gmf_header;
%writetable(table_p2p,strcat('Fixed_P2P.xlsx'));



%%%%%%%%%%See how this list compares with James List
%%%%%%Do it with the GMF serial number
cell_gmf_fixed_p2p_num=cell_gmf_fixed_p2p(:,1);
for i=1:1:length(cell_gmf_fixed_p2p_num)
    temp_string=cell_gmf_fixed_p2p_num{i};
    cell_gmf_fixed_p2p_num{i}=temp_string(find(~isspace(temp_string)));  %%%%%%%%%%Remove the White Spaces
end

% cell_p2p_pair_num=cell_p2p_pair([2:end],1);
% for i=1:1:length(cell_p2p_pair_num)
%     temp_string=cell_p2p_pair_num{i};
%     cell_p2p_pair_num{i}=temp_string(find(~isspace(temp_string)));  %%%%%%%%%%Remove the White Spaces
% end

size(cell_gmf_fixed_p2p_num)
% size(cell_p2p_pair_num)


%%%%%Find the max bandwidth of each fixed p2p EMS
cell_fixed_p2p_ems=cell_gmf_fixed_p2p(:,col_ems_idx);
num_ems=length(cell_fixed_p2p_ems);
array_ems=NaN(num_ems,2);
array_channel_bw=vertcat(1.25,2.5,5,10,20,30,40);
for i=1:1:num_ems
    temp_ems_str=cell_fixed_p2p_ems{i};
    if contains(temp_ems_str,',')
        temp_cell_ems=strsplit(temp_ems_str,',')';
    else
        temp_cell_ems=cell(1,1);
        temp_cell_ems{i}=temp_ems_str;
    end

    %%%%Convert to num
    temp_cell_ems=temp_cell_ems(~cellfun('isempty',temp_cell_ems));
    num_cell=length(temp_cell_ems);
    temp_array_ems=NaN(num_cell,1);
    for j=1:1:length(temp_cell_ems)
        temp_str=temp_cell_ems{j};
        temp_array_ems(j)=convert_str_ems2num_mhz_rev1(app,temp_str);
    end
    array_ems(i,1)=max(temp_array_ems);
    %%%%%%Snap to 40,30,20,10,5,2.5,1.25MHz channels
    temp_idx=find(array_ems(i,1)<=array_channel_bw);
    if isempty(temp_idx)
        %%%%%%44mhz seems to keep popping up
        %%%array_ems(i,1)
        array_ems(i,2)=array_channel_bw(end);
    else
        array_ems(i,2)=array_channel_bw(min(temp_idx));
    end
end

edges=horzcat(0,1.5,3,7.5,15,25,35,45);
[hist_num,edges] = histcounts(array_ems(:,2),edges);

%%'Find the Histogram bandwidth of the fixed p2p channels'
figure;
bar(array_channel_bw, hist_num);
grid on;
xlabel('P2P Channel Size');
ylabel('Counts');
filename1=strcat('P2P_Channel_Size_Hist.png');
pause(0.1)
saveas(gcf,char(filename1))
pause(0.1)






%%%%%%%%%%%%%%%%%%%%%%%%%%%Try Mapping the Links to PEA/County using
%%%%%%%%%%%%%%%%%%%%%%%%%%%(inti_graph_short_links_rev1.m) old code as a
%%%%%%%%%%%%%%%%%%%%%%%%%%%starting point.
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Find the total number in each PEA, all P2P links
tic;
load('cell_pea_census_data.mat','cell_pea_census_data') %%%%%PEA Name, PEA Num, PEA {Lat/Lon}, PEA Pop, PEA Centroid, Census {Geo ID},Census{Population},Census{NLCD}, Census Centroid
toc;

array_p2p_latlon=cell2mat(cell_gmf_fixed_p2p(:,[col_latdd_idx,col_londd_idx]));
[num_pea,~]=size(cell_pea_census_data)
cell_pea_inside_idx=cell(num_pea,1);
tic;
for pea_idx=1:1:num_pea
    temp_pea_bound=cell_pea_census_data{pea_idx,3};
    [inside_idx]=find_points_inside_contour_two_step(app,temp_pea_bound,array_p2p_latlon);
    cell_pea_inside_idx{pea_idx,1}=inside_idx;
end
toc;  %%%%0.5 seconds

cellsz =cell2mat(cellfun(@size,cell_pea_inside_idx,'uni',false));
pea_count_data=cellsz(:,1);
other_label='4GHz_P2P_links1_Count'
sum(pea_count_data)
%pea_heatmap_graph_rev2(app,pea_count_data,other_label)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Do it for 4.8-4.94GHz
array_p2p_freq=cell2mat(cell_gmf_fixed_p2p(:,col_freq1_idx));
above_idx=find(array_p2p_freq>=4800);
above_array_p2p_latlon=array_p2p_latlon(above_idx,:);
[num_pea,~]=size(cell_pea_census_data)
cell_pea_inside_idx=cell(num_pea,1);
tic;
for pea_idx=1:1:num_pea
    temp_pea_bound=cell_pea_census_data{pea_idx,3};
    [inside_idx]=find_points_inside_contour_two_step(app,temp_pea_bound,above_array_p2p_latlon);
    cell_pea_inside_idx{pea_idx,1}=inside_idx;
end
toc;  %%%%0.5 seconds

cellsz =cell2mat(cellfun(@size,cell_pea_inside_idx,'uni',false));
pea_count_data=cellsz(:,1);
other_label='4GHz_P2P_links_Count_Above4800MHz'
sum(pea_count_data)
%pea_heatmap_graph_rev2(app,pea_count_data,other_label)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Make the square tree map for number of links per agency'

[uni_agency,~,c_idx] = unique(cell_gmf_fixed_p2p(:,col_agency_idx));
[count_ag,edges]=histcounts(c_idx,unique(c_idx),'Normalization', 'probability');
[sort_count_ag,sort_idx]=sort(count_ag,'descend');
agency_rectangles = treemap_topleft(sort_count_ag);
color_set=flipud(plasma(length(uni_agency)));

%%%Need to find Abbrevations for Agencies
temp_idx=find(contains(uni_agency,'Energy'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='DOE';
end
temp_idx=find(contains(uni_agency,'Air Force'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='AF';
end
temp_idx=find(contains(uni_agency,'Coast Guard'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='CG';
end
temp_idx=find(contains(uni_agency,'Commerce'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Com';
end
temp_idx=find(contains(uni_agency,'Justice'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='J';
end
temp_idx=find(contains(uni_agency,'Marine Corps'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='MC';
end
temp_idx=find(contains(uni_agency,'Interior'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Int';
end
temp_idx=find(contains(uni_agency,'Agriculture'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Ag';
end
sort_uni_agency=uni_agency(sort_idx)

%%%%%%%%%%%%%%%%%%%%%%%Labels with Percentages
num_ag=length(sort_count_ag);
cell_ag_per=cell(num_ag,1);
for i=1:1:num_ag
    temp_str=sort_uni_agency{i} + "\n" + strcat(num2str(round(sort_count_ag(i)*100)),'%');
    cell_ag_per{i}=compose(temp_str);
end

close all;
f1=figure;
hold on;
plotRectangles(agency_rectangles,cell_ag_per,color_set)
outline(agency_rectangles)
title('P2P:4400-4940MHz')
f1.Position
f1.Position = [100 100 600 600];
pause(1)
filename1=strcat('TreeMap_P2P_All.png');
pause(0.1)
saveas(gcf,char(filename1))
pause(0.1)
close(f1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Just for Above 4800Mhz
[uni_agency,~,c_idx] = unique(cell_gmf_fixed_p2p(above_idx,col_agency_idx));
[count_ag,edges]=histcounts(c_idx,unique(c_idx),'Normalization', 'probability');
[sort_count_ag,sort_idx]=sort(count_ag,'descend');
agency_rectangles = treemap_topleft(sort_count_ag);
color_set=flipud(plasma(length(uni_agency)));


%%%Need to find Abbrevations for Agencies
temp_idx=find(contains(uni_agency,'Energy'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='DOE';
end
temp_idx=find(contains(uni_agency,'Air Force'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='AF';
end
temp_idx=find(contains(uni_agency,'Coast Guard'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='CG';
end
temp_idx=find(contains(uni_agency,'Commerce'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Com';
end
temp_idx=find(contains(uni_agency,'Justice'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='J';
end
temp_idx=find(contains(uni_agency,'Marine Corps'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='MC';
end
temp_idx=find(contains(uni_agency,'Interior'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Int';
end
temp_idx=find(contains(uni_agency,'Agriculture'));
if ~isempty(temp_idx)
    uni_agency{temp_idx}='Ag';
end
sort_uni_agency=uni_agency(sort_idx)

%%%%%%%%%%%%%%%%%%%%%%%Labels with Percentages
num_ag=length(sort_count_ag);
cell_ag_per=cell(num_ag,1);
for i=1:1:num_ag
    temp_str=sort_uni_agency{i} + "\n" + strcat(num2str(round(sort_count_ag(i)*100)),'%');
    cell_ag_per{i}=compose(temp_str);
end

close all;
f1=figure;
hold on;
plotRectangles(agency_rectangles,cell_ag_per,color_set)
outline(agency_rectangles)
title('P2P:4800-4940MHz')
f1.Position
f1.Position = [100 100 600 600];
pause(1)
filename1=strcat('TreeMap_P2P_Above4800MHz.png');
pause(0.1)
saveas(gcf,char(filename1))
pause(0.1)
close(f1)









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