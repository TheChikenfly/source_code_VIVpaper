% This function writes every example relevant values in the printParams struct into a .tex file

function texGenerator( printParams, exampleName, texFolderPath )

  %Name of the file
  nameExampleFileVals = ['valsExample' exampleName ];

  % open the tex file
  exampleTexFile = fopen( [ texFolderPath  nameExampleFileVals '.tex' ] ,'w') ;

  % Read print params list
  paramNames = fieldnames( printParams ) 

  % Print each param
  for paramIndex = 1:length( paramNames )
    nameParam  = paramNames{ paramIndex } ;
    valueParam = getfield(printParams, nameParam) ;
    fprintf(exampleTexFile, [ '\\newcommand{\\' exampleName nameParam '}{' num2str(valueParam) '}' '\n' ] ) ;
  end

  % close file
  fclose(exampleTexFile);

end