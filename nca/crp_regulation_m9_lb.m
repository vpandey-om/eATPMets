% Define the column names
CRP_col = 'CRP';
ratio_cols = {'ratio_t0', 'ratio_t5', 'ratio_t30', 'ratio_t60', 'ratio_t120', 'ratio_t180'};

% Path to the Excel file
table_path = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/crp_net_data_lb.xlsx';

% Read the table from the specified path
data_table = readtable(table_path);

% Extract the specified ratio columns as a matrix
ratio_matrix = table2array(data_table(:, ratio_cols));

% Create a mask for rows where all values are zero
all_zeros_mask = all(ratio_matrix == 0, 2);

% Filter out rows where all ratio values are zero
filtered_table = data_table(~all_zeros_mask, :);

% Extract the ratio columns from the filtered table
filtered_ratio_matrix = table2array(filtered_table(:, ratio_cols));

% Calculate log2 ratios for non-zero elements only
% Create a mask for non-zero values
nonzero_mask = filtered_ratio_matrix ~= 0;

% Initialize log2ratio_matrix with zeros
log2ratio_matrix = zeros(size(filtered_ratio_matrix));

% Calculate log2 for non-zero elements
log2ratio_matrix(nonzero_mask) = log2(filtered_ratio_matrix(nonzero_mask));

% Replace the original ratio columns in the filtered table with the log2 values
filtered_table(:, ratio_cols) = array2table(log2ratio_matrix, 'VariableNames', ratio_cols);

% Display the filtered data table
disp('Filtered Data Table with All Columns:');
disp(filtered_table);

% Save the filtered table to an Excel file
output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/filtered_data_table_with_all_columns.xlsx';
writetable(filtered_table, output_file);
disp(['Filtered data table with all columns saved to: ', output_file]);

load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')

% Identify CRP-regulated genes by intersecting with the model's gene list
crp_regulated_genes = intersect(filtered_table{:,10}, iML1515.genes);

% Initialize a cell array to store stacked results
stacked_results = {};

% Iterate through each CRP-regulated gene
for i = 1:length(crp_regulated_genes)
    % Get the index of the current gene in the model's gene list
    gene = crp_regulated_genes{i};
    [~, gind] = ismember(gene, iML1515.genes);
    
    % Find reactions associated with this gene
    [rind] = find(iML1515.rxnGeneMat(:, gind));
    
    % Extract relevant reaction information
    rxns = iML1515.rxns(rind);
    rxnNames = iML1515.rxnNames(rind);
    subSystems = iML1515.subSystems(rind);
    
    % Append the results for this gene to the stacked_results array
    for j = 1:length(rxns)
        stacked_results = [stacked_results; {gene, rxns{j}, rxnNames{j}, subSystems{j}}];
    end
end

% Convert the stacked results to a table
stacked_results_table = cell2table(stacked_results, ...
    'VariableNames', {'Gene', 'ReactionID', 'ReactionName', 'Subsystem'});

% Display the stacked table
disp('Stacked Results Table:');
disp(stacked_results_table);

% % Save the stacked table to an Excel file
% output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/stacked_results_table_lb.xlsx';
% writetable(stacked_results_table, output_file);
% disp(['Stacked results saved to: ', output_file]);

%%

% Add reaction information to the filtered table
filtered_genes = filtered_table{:, 10}; % Assuming column 10 contains the genes
[is_match, idx] = ismember(filtered_genes, stacked_results_table.Gene);

% Create columns for reaction information in the filtered table
reaction_id = cell(height(filtered_table), 1);
reaction_name = cell(height(filtered_table), 1);
subsystem = cell(height(filtered_table), 1);

% Populate reaction information for matching genes
reaction_id(is_match) = stacked_results_table.ReactionID(idx(is_match));
reaction_name(is_match) = stacked_results_table.ReactionName(idx(is_match));
subsystem(is_match) = stacked_results_table.Subsystem(idx(is_match));

% Add reaction information to the filtered table
filtered_table.ReactionID = reaction_id;
filtered_table.ReactionName = reaction_name;
filtered_table.Subsystem = subsystem;

% Save the updated filtered table to an Excel file
output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/filtered_table_with_reactions_lb.xlsx';
writetable(filtered_table, output_file);

disp(['Updated filtered table saved to: ', output_file]);


%%
% Define the column names
CRP_col = 'CRP';
ratio_cols = {'ratio_t0', 'ratio_t5', 'ratio_t30', 'ratio_t60', 'ratio_t120', 'ratio_t180'};
% Path to the Excel file
table_path = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/crp_net_data_m9.xlsx';

% Read the table from the specified path
data_table = readtable(table_path);

% Extract the specified ratio columns as a matrix
ratio_matrix = table2array(data_table(:, ratio_cols));

% Create a mask for rows where all values are zero
all_zeros_mask = all(ratio_matrix == 0, 2);

% Filter out rows where all ratio values are zero
filtered_table = data_table(~all_zeros_mask, :);

% Extract the ratio columns from the filtered table
filtered_ratio_matrix = table2array(filtered_table(:, ratio_cols));

% Calculate log2 ratios for non-zero elements only
% Create a mask for non-zero values
nonzero_mask = filtered_ratio_matrix ~= 0;

% Initialize log2ratio_matrix with zeros
log2ratio_matrix = zeros(size(filtered_ratio_matrix));

% Calculate log2 for non-zero elements
log2ratio_matrix(nonzero_mask) = log2(filtered_ratio_matrix(nonzero_mask));

% Replace the original ratio columns in the filtered table with the log2 values
filtered_table(:, ratio_cols) = array2table(log2ratio_matrix, 'VariableNames', ratio_cols);

% Display the filtered data table
disp('Filtered Data Table with All Columns:');
disp(filtered_table);

% Save the filtered table to an Excel file
output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/filtered_data_table_with_all_columns.xlsx';
writetable(filtered_table, output_file);
disp(['Filtered data table with all columns saved to: ', output_file]);

load('/Users/vikash/Documents/MATLAB/eatp_metabolism/models/iML1515.mat')

% Identify CRP-regulated genes by intersecting with the model's gene list
crp_regulated_genes = intersect(filtered_table{:,10}, iML1515.genes);

% Initialize a cell array to store stacked results
stacked_results = {};

% Iterate through each CRP-regulated gene
for i = 1:length(crp_regulated_genes)
    % Get the index of the current gene in the model's gene list
    gene = crp_regulated_genes{i};
    [~, gind] = ismember(gene, iML1515.genes);
    
    % Find reactions associated with this gene
    [rind] = find(iML1515.rxnGeneMat(:, gind));
    
    % Extract relevant reaction information
    rxns = iML1515.rxns(rind);
    rxnNames = iML1515.rxnNames(rind);
    subSystems = iML1515.subSystems(rind);
    
    % Append the results for this gene to the stacked_results array
    for j = 1:length(rxns)
        stacked_results = [stacked_results; {gene, rxns{j}, rxnNames{j}, subSystems{j}}];
    end
end

% Convert the stacked results to a table
stacked_results_table = cell2table(stacked_results, ...
    'VariableNames', {'Gene', 'ReactionID', 'ReactionName', 'Subsystem'});

% Display the stacked table
disp('Stacked Results Table:');
disp(stacked_results_table);

% % Save the stacked table to an Excel file
% output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/stacked_results_table_lb.xlsx';
% writetable(stacked_results_table, output_file);
% disp(['Stacked results saved to: ', output_file]);

%%

% Add reaction information to the filtered table
filtered_genes = filtered_table{:, 10}; % Assuming column 10 contains the genes
[is_match, idx] = ismember(filtered_genes, stacked_results_table.Gene);

% Create columns for reaction information in the filtered table
reaction_id = cell(height(filtered_table), 1);
reaction_name = cell(height(filtered_table), 1);
subsystem = cell(height(filtered_table), 1);

% Populate reaction information for matching genes
reaction_id(is_match) = stacked_results_table.ReactionID(idx(is_match));
reaction_name(is_match) = stacked_results_table.ReactionName(idx(is_match));
subsystem(is_match) = stacked_results_table.Subsystem(idx(is_match));

% Add reaction information to the filtered table
filtered_table.ReactionID = reaction_id;
filtered_table.ReactionName = reaction_name;
filtered_table.Subsystem = subsystem;

% Save the updated filtered table to an Excel file
output_file = '/Users/vikash/Documents/MATLAB/eatp_metabolism/revision/filtered_table_with_reactions_m9.xlsx';
writetable(filtered_table, output_file);

disp(['Updated filtered table saved to: ', output_file]);





