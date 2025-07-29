
for i=1:length(averaged_eods)

measurement_data(i)=standard_eod_measurement(averaged_eods(i).wave,averaged_eods(i).sampRate,averaged_eods(i).sample_name);

disp("Creating Figure.  Ready for Plotting!")
    
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'units','normalized')

%% Plot All Normalized EODs
subplot(2,2,1)
        hold on;
            [~,n_waves]=size(normalized_eods(i).wave);
            for j=1:n_waves

                if isempty(normalized_eods(i).wave(:,j))
                    continue
                end
                time=linspace(0,1000*length(normalized_eods(i).wave(:,j))/normalized_eods(i).sampRate,length(normalized_eods(i).wave(:,j)));
                [~,ivmax]=max(normalized_eods(i).wave(:,j));
                time = time - (time(ivmax));
                plot(time,normalized_eods(i).wave(:,j)); %,'DisplayName',normalized_eods(i).info_text{j});
            end

            title(['N= ', num2str(n_waves) , ' EODs from ',string(normalized_eods(i).sample_name)],'Interpreter', 'none')
            hold off;

%% Plot All Normalized EODs, Focusing on P0    
subplot(2,2,2)
    hold on;
        [~,n_waves]=size(normalized_eods(i).wave);
        for j=1:n_waves

            if isempty(normalized_eods(i).wave(:,j))
                continue
            end
            time=linspace(0,1000*length(normalized_eods(i).wave(:,j))/normalized_eods(i).sampRate,length(normalized_eods(i).wave(:,j)));
            [~,ivmax]=max(normalized_eods(i).wave(:,j));
            time = time - (time(ivmax));
            plot(time,normalized_eods(i).wave(:,j)) %,'DisplayName',normalized_eods(i).info_text{j});
        end

        title(['N= ', num2str(n_waves) , ' EODs from ',string(normalized_eods(i).sample_name)],'Interpreter', 'none')
        axis( [ time(measurement_data(i).iZC1)-1 (time(measurement_data(i).iZC1)+1) -4e-3 3e-3])
        hold off;

%% Plot the Averaged EOD with Landmarks
subplot(2,2,3)
    [~,ivmax]=max(averaged_eods(i).wave);
    time = time - (time(ivmax));
    plot(time,averaged_eods(i).wave);
    hold on;
    plot(measurement_data(i).tT1,measurement_data(i).vT1,'k+'); %begin and end time landmarks are always present
    text(measurement_data(i).tT1+.05,measurement_data(i).vT1,'\leftarrow EOD start', 'HorizontalAlignment','left')
    plot(measurement_data(i).tT2,measurement_data(i).vT2,'k+');
    text(measurement_data(i).tT2+.05,measurement_data(i).vT2,'\leftarrow EOD end', 'HorizontalAlignment','left')
    plot(measurement_data(i).tP1,measurement_data(i).vP1,'r+');					%phase1 landmarks (P1, S1, S2) always present
    text(measurement_data(i).tP1+.05,measurement_data(i).vP1,'\leftarrow P1', 'HorizontalAlignment','left')
    plot(measurement_data(i).tS1,measurement_data(i).vS1,'r+');
    text(measurement_data(i).tS1+.05,measurement_data(i).vS1,'\leftarrow Inflection 1', 'HorizontalAlignment','left')
    plot(measurement_data(i).tS2,measurement_data(i).vS2,'r+');
    text(measurement_data(i).tS2+.05,measurement_data(i).vS2,'\leftarrow Inflection 2', 'HorizontalAlignment','left')

    plot(measurement_data(i).tP0,measurement_data(i).vP0,'r+');
    text(measurement_data(i).tP0+.05,measurement_data(i).vP0,'\leftarrow P0', 'HorizontalAlignment','left','Rotation', -30)

    plot(measurement_data(i).tZC1,measurement_data(i).vZC1,'k+');
    text(measurement_data(i).tZC1+.05,measurement_data(i).vZC1,'\leftarrow Zero Crossing 1', 'HorizontalAlignment','left','Rotation', -30)
    plot(measurement_data(i).tZC2,measurement_data(i).vZC2,'k+');
    text(measurement_data(i).tZC2+.05,measurement_data(i).vZC2,'\leftarrow Zero Crossing 2', 'HorizontalAlignment','left', 'Rotation', -30)
    plot(measurement_data(i).tP2,measurement_data(i).vP2,'r+');
    text(measurement_data(i).tP2+.05,measurement_data(i).vP2,'\leftarrow P2', 'HorizontalAlignment','left')
    
    title(averaged_eods(i).sample_name,'Interpreter', 'none')
        
%% Plot Averaged EOD, with Focus on P0
subplot(2,2,4);
    plot(time,averaged_eods(i).wave,'k')
    width=fix(averaged_eods(i).sampRate*.0005);
    hold on;
%     plot(time(istart:izc1),wave(istart:izc1),'r')
    x=time(measurement_data(i).iZC1-width:measurement_data(i).iZC1);
    y=averaged_eods(i).wave(measurement_data(i).iZC1-width:measurement_data(i).iZC1);
    areashade(x,y,0,'r','h')
    areashade(x,y,0,'g','l')
    mark=plot(measurement_data(i).tP0,measurement_data(i).vP0,'ok');
    text(measurement_data(i).tP0+.025,measurement_data(i).vP0,['\leftarrow P0.  aP0=', num2str(measurement_data(i).aP0)], 'HorizontalAlignment','left','rotation', -25)

    set(mark,'MarkerSize',15);
    axis( [ time(measurement_data(i).iZC1)-1 (time(measurement_data(i).iZC1)+1) -4e-3 3e-3])
    hold off;
    
%% Export PDF
pdfprinfig=gcf;

file_name = strcat(averaged_eods(i).sample_name, '.png'); % Ensure it's a single string
file_path = fullfile('EOD_Plots', file_name);
print(pdfprinfig, '-dpng', file_path); 
close(pdfprinfig)
disp("Closed Figure.  Ready for Next!")

end

%% Write Measurement File

colnames=flip(fieldnames(measurement_data));
measurement_data = orderfields(measurement_data, colnames);
colnames=fieldnames(measurement_data)';
ds=squeeze(struct2cell(measurement_data))';
ds_f=vertcat(colnames,ds);
fid = fopen('EOD_Plots/measurement_data.csv','wt');

if fid>0
 for k=1:size(ds_f,1)
     if k==1
        fprintf(fid,'%s,%s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',ds_f{k,:});
     else
        fprintf(fid,'%s,%f,%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',ds_f{k,:});
     end
 end
  fclose(fid);
end

%% Clean Up
clearvars -except normalized_eods averaged_eods subjectlist measurement_data 