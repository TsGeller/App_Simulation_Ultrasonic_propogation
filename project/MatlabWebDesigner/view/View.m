classdef View
    properties
        app
    end
    methods
        function obj = View(appParam)
            obj.app = appParam;
        end
        function defaultParam(obj,simulation)
            obj.app.NameAnimaldefaultMacaqueEditField.Value = simulation.animalName;
            obj.app.MaxwaterthresholdHUEditField.Value = num2str(simulation.CT_th1);
            obj.app.MinbonethresholdHUEditField.Value = num2str(simulation.CT_th2);
            obj.app.MaxbonethresholdHUEditField.Value = num2str(simulation.CT_th3);
            obj.app.WaveSpeedmsEditField.Value = num2str(simulation.wave_speed);
            obj.app.WaveFrequencyMHzEditField.Value = num2str(simulation.wave_frequency);
            obj.app.PointsperwavelengthEditField.Value = num2str(simulation.ppw);
            obj.app.DensitywaterKgm3EditField.Value  = num2str(simulation.dwater);
            obj.app.SpeedsofttissuemsEditField.Value  = num2str(simulation.stissue);
            obj.app.SofttissuedensityKgm3EditField.Value = num2str(simulation.dtissue);
            obj.app.SpeedcorticalbonemsEditField.Value = num2str(simulation.sbone);            
            obj.app.TransducercurvatureradiusmEditField.Value = num2str(simulation.source_roc);
            obj.app.EffectivetransducerdiametermEditField.Value = num2str(simulation.source_diameter);
            obj.app.TransducersurfaceacousticpressurePaEditField.Value = num2str(simulation.source_amp);
            obj.app.PointperperiodEditField.Value = num2str(simulation.PPP);
            obj.app.FocalDistancemmEditField.Value = num2str(simulation.focol_distance);
            return ;
        end

        function showMessage(obj,exception)
            fprintf('Une erreur s est produite : %s\n', exception.message);
        end
        function disableButtonWhilesimulation(obj)
            obj.app.MRIButton.Enable = 'off';
            obj.app.CTButton.Enable = 'off';
            obj.app.SimulateButton.Enable = 'off';
        end
        function enableButtonWhilesimulation(obj)
            obj.app.MRIButton.Enable = 'on';
            obj.app.CTButton.Enable = 'on';
            obj.app.SimulateButton.Enable = 'on';
        end

        function displaypathOfFileCT(obj, simulation)
            if(simulation.filepath1CT ~= 0)
                parts = strsplit(simulation.filepath1CT, filesep); % Utilisez filesep pour assurer la portabilité entre les systèmes d'exploitation

                % Récupérez le dossier parent en utilisant l'avant-dernière partie
                if length(parts) >= 2
                    dirPart = parts{end - 1};
                else
                    dirPart = '';  % Si le chemin n'a pas d'élément avant-dernier (dossier parent)
                end                
                obj.app.pathoffileCTLabel.Text = "/"+dirPart+"/";
                obj.app.pathoffileCTLabel.Visible = 'on';
            else
                obj.app.pathoffileCTLabel.Text = 'No files selected';
                obj.app.pathoffileCTLabel.Visible = 'on';
            end
        end

        function displayPathOfFileIrm(obj,simulation)
            if(simulation.filepathIMR ~=0)
                parts = strsplit(simulation.filepathIMR, filesep); % Utilisez filesep pour assurer la portabilité entre les systèmes d'exploitation

                % Récupérez le dossier parent en utilisant l'avant-dernière partie
                if length(parts) >= 2
                    dirPart = parts{end - 1};
                else
                    dirPart = '';  % Si le chemin n'a pas d'élément avant-dernier (dossier parent)
                end
                
                obj.app.pathoffileimrLabel.Text = "/"+dirPart+"/";
                obj.app.pathoffileimrLabel.Visible = 'on';
            else
                obj.app.pathoffileimrLabel.Text = 'No files selected';
                obj.app.pathoffileimrLabel.Visible = 'on';

            end
        end
        function checkRequirementFiles(obj,filepathCt,filepathIrm)
            if ~isempty(filepathCt) && ~isempty(filepathIrm)
                obj.app.SimulateButton.Enable = 'on';
            end
        end
        function showMessageStartSimulation(obj)
            obj.app.statusSimulationLabel.Text = 'Loading simulation...';
            obj.app.statusSimulationLabel.Visible = 'on';
        end
        function showMessageEndSimulation(obj)
            obj.app.statusSimulationLabel.Text = 'Simulation terminated.';
        end
        function showMessageUI(obj,string)
             obj.app.statusSimulationLabel.Text = string;
            obj.app.statusSimulationLabel.Visible = 'on';
        end
        
        function saveData(obj)
        obj.app.message_statusLabel.Visible= 'on';
        obj.app.message_statusLabel.Text = 'Data saved';
        pause(2)
        obj.app.message_statusLabel.Visible= 'off';
        end
        function unsaveData(obj)
            obj.app.message_statusLabel.Visible= 'on';
            obj.app.message_statusLabel.Text = 'Data Unsaved';
            pause(2)
            obj.app.message_statusLabel.Visible= 'off';
        end
        function exportData(obj)
            obj.app.status_saveLabel.Visible = 'on';
            obj.app.status_saveLabel.Text = 'Data Exported !';
            pause(2);
            obj.app.status_saveLabel.Visible = 'off';
        end
        function showPathConfig(obj,string)
            obj.app.path_configjsonLabel.Visible = 'on';
            obj.app.path_configjsonLabel.Text = string;
        end
    end
end

