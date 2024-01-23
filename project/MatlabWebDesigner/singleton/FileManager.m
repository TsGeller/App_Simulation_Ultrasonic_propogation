classdef FileManager < handle
    properties (Access = private)
        FileID
    end
    
    methods (Access = private)
        function obj = FileManager()
            % Ouvrir le fichier pour l'écriture
            obj.FileID = fopen('dataSimulation.json', 'w');
        end
    end
    
    methods (Static)
        function instance = getInstance()
            persistent uniqueInstance
            if isempty(uniqueInstance) || ~isvalid(uniqueInstance)
                uniqueInstance = FileManager();
            end
            instance = uniqueInstance;
        end
    end
    
    methods
        function writeToFile(obj, data)
            json_str =jsonencode(data);
            fprintf(obj.FileID, '%s\n', json_str);
        end
        
        function closeFile(obj)
            fclose(obj.FileID);
        end
    end
end
