function date_str = fct_data_2011(day)
% Write the date with format yyyymmdd from the day-th of year 2011
%

year = 2011;
ref_idx_day = datenum(['1-January-' num2str(year)])-1;
date_idx = ref_idx_day + day;
date_str = datestr(date_idx, 30);
date_str(9:end) = [];