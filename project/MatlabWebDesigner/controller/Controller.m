classdef Controller <handle
    properties
        simulation
        view
        fileManager
    end
    methods(Access = public)
        function setUpDataConfig(obj)
            obj.view.defaultParam(obj.simulation)
        end
        function controller = Controller(appParam)
            controller.view = View(appParam);
            controller.simulation = Simulation();            
        end

        function model = makeSimulation(model)
            model.view.disableButtonWhilesimulation();
            model.view.showMessageStartSimulation();
            try
            model.simulation.makeSimulation();
            model.view.enableButtonWhilesimulation();
            model.view.showMessageEndSimulation();
            model.exportData();
            catch exeption
                model.view.enableButtonWhilesimulation();
                model.view.showMessageUI("Error: "+exeption.message);
            end
            
        end

        function obj =openCTFile(obj)
            try
                obj.simulation.openCTFile();
                obj.view.checkRequirementFiles(obj.simulation.filenameCT,obj.simulation.filenameIMR);
            catch MException
                obj.view.showMessage(MException);
            end
            obj.view.displaypathOfFileCT(obj.simulation);
        end

        function obj = openIRMFile(obj)
            try
                obj.simulation.openIRMFile();
                obj.view.checkRequirementFiles(obj.simulation.filenameCT,obj.simulation.filenameIMR);
            catch MException
                obj.view.showMessage(MException);
            end
            obj.view.displayPathOfFileIrm(obj.simulation);
        end
        function validate = saveData(obj, app)
            validate = true;

            try
                validateInput = @(x) ~isempty(regexp(x, '^\d*\.?\d*$', 'once')); % Vérifie que x contient uniquement des chiffres et des décimales

                if validateInput(app.MaxwaterthresholdHUEditField.Value) && ...
                        validateInput(app.MinbonethresholdHUEditField.Value) && ...
                        validateInput(app.MaxbonethresholdHUEditField.Value) && ...
                        validateInput(app.WaveFrequencyMHzEditField.Value) && ...
                        validateInput(app.WaveSpeedmsEditField.Value) && ...
                        validateInput(app.PointsperwavelengthEditField.Value) && ...
                        validateInput(app.DensitywaterKgm3EditField.Value) && ...
                        validateInput(app.SpeedsofttissuemsEditField.Value) && ...
                        validateInput(app.SofttissuedensityKgm3EditField.Value) && ...
                        validateInput(app.SpeedcorticalbonemsEditField.Value) && ...                        
                        validateInput(app.TransducercurvatureradiusmEditField.Value) && ...
                        validateInput(app.EffectivetransducerdiametermEditField.Value) && ...
                        validateInput(app.PointperperiodEditField.Value) && ...
                        validateInput(app.FocalDistancemmEditField.Value) && ...
                        validateInput(app.TransducersurfaceacousticpressurePaEditField.Value)
                        

                    obj.simulation.animalName = app.NameAnimaldefaultMacaqueEditField.Value;
                    obj.simulation.CT_th1 = str2double(app.MaxwaterthresholdHUEditField.Value);
                    obj.simulation.CT_th2 = str2double(app.MinbonethresholdHUEditField.Value);
                    obj.simulation.CT_th3 = str2double(app.MaxbonethresholdHUEditField.Value);
                    obj.simulation.wave_frequency = str2double(app.WaveFrequencyMHzEditField.Value);
                    obj.simulation.wave_speed = str2double(app.WaveSpeedmsEditField.Value);
                    obj.simulation.ppw = str2double(app.PointsperwavelengthEditField.Value);
                    obj.simulation.dwater = str2double(app.DensitywaterKgm3EditField.Value);
                    obj.simulation.stissue = str2double(app.SpeedsofttissuemsEditField.Value);
                    obj.simulation.dtissue = str2double(app.SofttissuedensityKgm3EditField.Value);
                    obj.simulation.sbone = str2double(app.SpeedcorticalbonemsEditField.Value);
                    obj.simulation.PPP = str2double(app.PointperperiodEditField.Value);
                    obj.simulation.source_roc = str2double(app.TransducercurvatureradiusmEditField.Value);
                    obj.simulation.source_diameter = str2double(app.EffectivetransducerdiametermEditField.Value);
                    obj.simulation.source_amp = str2double(app.TransducersurfaceacousticpressurePaEditField.Value);
                    obj.simulation.focol_distance = str2double(app.FocalDistancemmEditField.Value);
                else
                    error('One or more entries are not valid numbers.');
                end
            catch exception
                disp(exception);
                validate = false;
            end
        end
        function updateStatus(obj,boolean)
            if(boolean)
                obj.view.saveData();
            else
                obj.view.unsaveData();
            end
        end
        function exportData(obj)
            obj.fileManager = FileManager.getInstance();
            disp('file manager')
            disp(obj.fileManager)
            data_str = obj.simulation.exportData();
            obj.fileManager.writeToFile(data_str);
            %obj.fileManager.closeFile();
            obj.view.exportData();          
        end
        function loadConfig(obj)
            try
                pathConfig =obj.simulation.loadConfig();
                obj.view.showPathConfig(pathConfig);
                obj.setUpDataConfig();                
            catch
                error = 'None config loaded';
                obj.view.showPathConfig(error);
                pause(2);
                error = '';
                obj.view.showPathConfig(error);
            end            
            
        end
    end
end