function [output] = CleanMyTable(TableInput,Max_NaN)

% For any data given as a table, this function will clean the data under 2
% cases:
% 1) remove the columns that has more than Max_NaN.
% 2) remove the rows with NaN's.

% Loading the data:
data_clean = TableInput;
data_array = table2array(data_clean);
% Removing columns with more than Max_NaN:
len.columns = length(data_array(1,:));  % Saves the column length
len.rows = length(data_array(:,1));     % Saves the row length
counter = zeros(1,len.columns);
for n = 1:len.columns
    counter(len.columns - n + 1) = sum(isnan(data_array(:,...
        len.columns - n + 1)));
    if counter(len.columns - n + 1)>Max_NaN
        data_clean(:,len.columns - n + 1) = [];
    end
end
% Removing rows with NaN's:
data_clean = rmmissing(data_clean);
% Outputting the clean data:
[output] = data_clean;
end