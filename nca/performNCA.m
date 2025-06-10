function [Ae, Se, Ae_log2, Se_log2] = performNCA(table_path, CRP_col, ratio_cols)
    % Read the table from the specified path
    data_table = readtable(table_path);
    
    % Extract the 'CRP' column as a matrix
    CRP_matrix = data_table.(CRP_col);
    
    % Extract the specified ratio columns as a matrix
    ratio_matrix = table2array(data_table(:, ratio_cols));
    
    % Create a mask for non-zero values
    nonzero_mask = ratio_matrix ~= 0;

    % Initialize log2ratio_matrix with zeros
    log2ratio_matrix = zeros(size(ratio_matrix));

    % Calculate log2 for non-zero elements only
    log2ratio_matrix(nonzero_mask) = log2(ratio_matrix(nonzero_mask));

    % Perform NCA on the raw ratio matrix
    [Ae, Se] = FastNCA(ratio_matrix, CRP_matrix);

    % Perform NCA on the log2 transformed matrix
    [Ae_log2, Se_log2] = FastNCA(log2ratio_matrix, CRP_matrix);

   
end
