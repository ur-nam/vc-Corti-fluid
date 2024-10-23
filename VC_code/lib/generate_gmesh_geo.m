function generate_gmesh_geo(xx)

    np = length(xx);

    fileID = fopen('./hinput/gmesh_geo.geo','w');
    fprintf(fileID,'//+\nSetFactory("OpenCASCADE")\n');


    id = 1:np;
    fprintf(fileID,'//+\nPoint(%d) = {%.3f, %.3f, 1.0};\n',transpose([id(:),xx]));

    conn = [[id(1:end-1);id(2:end)],[id(end);id(1)]];
    nl = length(conn);
    id = 1:nl;
    fprintf(fileID,'//+\nLine(%d) = {%d, %d};\n',[id;conn(1,:);conn(2,:)]);
    fprintf(fileID,'//+\nCurve Loop(1) = {');
    fprintf(fileID,'%d, ',id);
    fprintf(fileID,'\b\b};\n');
    fprintf(fileID,'//+\nPlane Surface(1) = {1};\n');

    fclose(fileID);    

end