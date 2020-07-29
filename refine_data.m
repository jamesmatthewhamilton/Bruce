load("BRCA_patient.mat")
load("consolidated_data.mat")
load("consolidated_features.mat")
load("consolidated_patient_info.mat")

feature_idx = 41;

data_CDE_ID = CelltoArray(data(3:end, 2));
patient_info_column_3 = CelltoArray(patient_info(:, 3));

data_CDE_ID = CutArrayElements(data_CDE_ID, 12);
patient_info_column_3 = CutArrayElements(patient_info_column_3, 12);
data_CDE_ID = DisfigureDuplicates(data_CDE_ID);
patient_info_column_3 = DisfigureDuplicates(patient_info_column_3);

matches_rev = CountMatches(patient_info_column_3, data_CDE_ID);
fprintf("Matches between patient_info_column_1 and data_CDE_ID: %i\n", matches_rev)
matches= CountMatches(data_CDE_ID, patient_info_column_3);
fprintf("Matches between data_CDE_ID and patient_info_column_1: %i\n", matches_rev)

%data_pathologic_stage = CelltoArray(data(3:end, 41));
[x, y] = ConcateGivenKey(data_CDE_ID, patient_info_column_3, matches, ...
                         all_data, data(3:end, feature_idx));
                                  
x = x';
writecell(all_features, "feature_lables.csv");
writematrix(x, "x.csv");
writecell(y, strcat("y_f", int2str(feature_idx), ".csv"));

function [matches, array_1] = CountMatches(array_1, array_2)
    matches = 0;
    for ii = 1 : length(array_1)
        for jj = 1 : length(array_2)
            %fprintf("Comparing: %s to %s\n", array_1(ii), array_2(jj));
            if (array_1(ii) == array_2(jj))
               matches = matches + 1;
               break
            end
            if (jj == length(array_2))
                fprintf("No match at %s.\n", array_1(ii))
            end
        end
    end
end

function [matches, array_1] = CountMatchesOLD(array_1, array_2, cut_at)
    matches = 0;
    found_at_pre = zeros(size(array_2)*2);
    found_at = zeros(size(array_2)*2);
    next = 1;
    for ii = 1 : length(array_1)
        for jj = 1 : length(array_2)
            %fprintf("Comparing: %s to %s\n", array_1(ii), array_2(jj));
            char_1 = char(array_1(ii));
            char_2 = char(array_2(jj));
            if (strcmp(char_1(1:cut_at), char_2(1:cut_at)))
               for kk = 1 : next
                   if (found_at(kk) == jj)
                        fprintf("Warning duplicate data at key1:%d and key2:%d element %s\n", ii, jj, array_2(jj))
                        fprintf("    -> Previously found at key1:%d \n", found_at_pre(kk))
                        fprintf("    -> DISFIGURING ELEMENT AT key1:%d\n", ii)
                        array_1(ii) = "NULL" + array_1(ii);
                        fprintf("    -> %s\n", array_1(ii))
                   end
               end
               found_at_pre(next) = ii;
               found_at(next) = jj;
               matches = matches + 1;
               %fprintf("^^^ Match Found!\n");
               next = next + 1;
               break
            end
            if (jj == length(array_2))
                fprintf("No match at %s.\n", array_1(ii))
            end
        end
    end
end


function array_1 = DisfigureDuplicates(array_1)
    num_disfigs = 0;
    for ii = length(array_1) : -1 : 1
        for jj = 1 : length(array_1)
            if (array_1(ii) == array_1(jj) && ii ~= jj)
                fprintf("Warning duplicate data at index %d and %d element %s\n", ii, jj)
                fprintf("    -> %s \n", array_1(jj))
                fprintf("    -> DISFIGURING ELEMENT AT :%d\n", ii)
                array_1(ii) = "NULL" + array_1(ii);
                fprintf("    -> %s\n", array_1(ii))
                num_disfigs = num_disfigs + 1;
                break
            end
        end
    end
    fprintf("Total Disfigurations: %i\n", num_disfigs)
end

function [x, y] = ConcateGivenKey(key_array_1, key_array_2, matches, x_in, y_in) % 1222, 1099, 1091, 1222, 1099
    [r,c] = size(x_in);
    y = cell(matches, 1);
    x = zeros(matches, r);
    next = 1;
    for ii = 1 : length(key_array_1)
        for jj = 1 : length(key_array_2)
            %fprintf("Comparing: %s to %s\n", key_array_1(ii), key_array_2(jj));
            if (key_array_1(ii) == key_array_2(jj))
               fprintf("^^^ Match found at key1:%d and key2:%d... (%d of %d)\n", ii, jj, next, matches);
               y(next,1) = y_in(ii);
               x(next,:) = x_in(:,jj);
               next = next + 1;
               break
            end
            if (jj == length(key_array_2))
                fprintf("No match at %s.\n", key_array_1(ii))
            end
        end
    end
end

function array_string = CelltoArrayLengthXCheck(cell, cut_at)
    matrix = cell2mat(cell);
    array_string = strings(1, length(matrix));
    for ii = 1 : length(matrix)
        cell_contents = matrix(ii,:);
        if (cell_contents(end-3:end) == "-01A")
            %fprintf("%s\n", cell_contents(end-3:end))
            array_string(ii) = lower(convertCharsToStrings(cell_contents(1:cut_at)));
        else
            array_string(ii) = "DUPLICATE" + lower(convertCharsToStrings(cell_contents)); 
        end
    end
end

function array_string = CelltoArray(cell)
    matrix = cell2mat(cell);
    array_string = strings(1,length(matrix));
    for ii = 1 : length(matrix)
        cell_contents = matrix(ii,:);
        array_string(ii) = lower(convertCharsToStrings(cell_contents));
    end
end

function array = CutArrayElements(array, cut_off)
    for ii = 1 : length(array)
        array_char = char(array(ii));
        array(ii) = string(array_char(1:cut_off));
    end
end