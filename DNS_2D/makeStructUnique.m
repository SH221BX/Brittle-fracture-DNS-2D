function uniqueStructArray = makeStructUnique(structArray)
    % Find the number of structs in the array
    numStructs = numel(structArray);
    % Initialize a logical index vector to mark unique structs
    isUnique = true(1, numStructs);
    % Get all field names
    fieldNames = fieldnames(structArray);
    
    % Iterate through each struct
    for i = 1:numStructs
        if ~isUnique(i)
            % If this entry is already marked as not unique, skip it
            continue;
        end
        
        % Iterate through each other struct for comparison
        for j = i+1:numStructs
            % Compare each field of the current and other struct
            areEqual = true;
            for k = 1:numel(fieldNames)
                field = fieldNames{k};
                % Check if the fields are equal
                if ~isequal(structArray(i).(field), structArray(j).(field))
                    areEqual = false;
                    break; % If any field is not equal, the structs are not the same
                end
            end
            
            % If all fields are equal, mark the other struct as not unique
            if areEqual
                isUnique(j) = false;
            end
        end
    end
    
    % Select the structs that are marked as unique
    uniqueStructArray = structArray(isUnique);
end

