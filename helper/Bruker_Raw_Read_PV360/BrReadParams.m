function arr = BrReadParams(params, path)
% Reads specified parameters from a Bruker Parameter file (p.ex. acqp)
% Opens the file path or, if not given, opens file selection dialogue 
% to select a file, reads it and returns the values of the parameters 
% contained in the string array params.
% For the params-array, the following possibilities exist:
% - A single parameter name: In this case function returns the value of
% this parameter, with the correct type.
% - An arrays of parameters: The function returns a structure, the tag
% names of which are the parameter names.
% The parameter names in the params-array may also contain array indices ()
% or structure tag-numbers (p.ex. 'TPQQ.1').
% If a parameter is an array, the names of the returned array are named
% value1, ....

% Written by Rolf Pohmann, Max Planck Institute for Biological Cybernetics,
% Tübingen, Germany

if nargin < 1
    disp('Reads Bruker Parameter files and returns parameter values.');
    disp('Usage: result = BrReadParams(ParameterArray, Filename');
    arr = -1;
    return;
end
if ~isequal(class(params), 'char') &  ~isequal(class(params), 'cell')
    disp('Reads Bruker Parameter files and returns parameter values.');
    disp('Usage: result = BrReadParams(ParameterArray, Filename');
    arr = -1;
    return;
end
if isequal(class(params), 'char')
    params = {params};
end
if nargin < 2
    [f,p] = uigetfile([{'acqp;method;imnd;subject' ...
            'Bruker parameter files';'*.*','All Files'}], ... 
            'Select parameter-file for reading')
    if f == 0
        arr = -1;
        return;
    end
    path = strcat(p,f);
end
%% Open file 
file = fopen(path,'r');
if file == -1
    arr = -2;
    return
end
line = 0;
paramfound = zeros(numel(params),1);
if numel(params)>1
    for cnt = 1: numel(params)
        arr.(params{cnt})={};
    end
end
%% Read file line by line
while line ~= -1
  line = fgets(file);
  if line ~= -1
    line = strtrim(line);
    %Does the line contain a parameter name?
    if (strfind(line, '##$'))==1
        [pname, pvalue] = strtok(line(4:end), '=');
        %Does this parameter appear in the params-array?
        nfoundparams = 0;
          for cnt=1:numel(params)

            if strcmpi(pname,strtok(params{cnt},'([.'))   %Parameter found
               nfoundparams = nfoundparams + 1;
               %disp([pname,pvalue]);
               if nfoundparams > 1
                   res = saveres;
               else
                   IsStruct = 0;
                   %What type of parameter is it?
                   pvalue = pvalue(2:end);  %Remove =
                   pvalue = strtrim(pvalue);
                   if pvalue(1)~='('    %It's not an array: number or String
                       num = find((pvalue>57 | pvalue<43 | pvalue==44) ...
                           & pvalue~=101 & pvalue~=69);
                       if numel(num)==0    %It's a number
                           %disp([pname,': ',pvalue]);
                           res = str2double(pvalue);
                       else %end of parameter is a single number
                           res = pvalue;                       
                       end  %end of parameter is a single string
                   else     %end of parameter is a scalar
                       %parameter is an array or a structure
                       %Is it an array?
                       if numel(find((pvalue(2:end-1)>57 | pvalue(2:end-1)<43) ...
                               & pvalue(2:end-1)~=32 & pvalue(2:end-1)~=44 ) ) == 0 ... %for Jobs value
                               & pvalue(end)==')'
                           NumDim = str2num(pvalue(2:end-1));
                           if NumDim==10240
                               NumDim=32;
                           end
                           if numel(NumDim) == 1 & (strcmpi(pname,'ACQ_jobs')==0) %for PV360
                               NumDim = [1,NumDim];
                           end
                           Ind = 0;
                           r = '';
                           res = '';
                           
                           while Ind < prod(NumDim)
                               l = strtrim(fgets(file));
                               %Determine type
                               if numel(res)==0
                                    if numel(find((l>57 | l<43) & l~=32 & l~=44)) == 0  %only numbers
                                        IsNum = 1;
                                        res = zeros(NumDim);
                                    elseif l(1)=='<'    %A single string
                                        IsNum = -1;
                                        res = '';
                                    elseif l(1)=='('    %An array of structures
                                        IsNum = -2;
                                        res = {prod(NumDim)};
                                        IsStruct = 1;
                                    else                %Any other kind of array
                                        IsNum = 0;
                                        res = cell(NumDim);
                                    end
                                    res = res';
                               end
                               while numel(l) > 0
                                   if IsNum==0
                                       [res{Ind+1}, l] = strtok(l);
                                       Ind = Ind+1;
                                   elseif IsNum == -1
                                       if Ind==0 & l(end)=='>'
                                           res = l(2:end-1)';
                                           l='';Ind =prod(NumDim); 
                                       elseif Ind==0 & l(end)~='>'
                                           res = l(2:end)';
                                       elseif Ind~=0 & l(end)=='>'
                                           res = [res,l(1:end-1)];
                                           l='';Ind =prod(NumDim);
                                       elseif Ind~=0 & l(end)~='>'
                                           res = [res,l(1:end)];
                                       end
                                       Ind = Ind+1;
                                    elseif IsNum==-2
                                        if numel(r)>0
                                            l = [r,l];
                                        end
                                       [r, l] = strtok(l(strfind(l,'(')+1:end),')');
                                       if numel(l)>0
                                           if Ind == 0
                                               res = BrConvertToStruc(r); 
                                           else
                                               res(Ind+1) = BrConvertToStruc(r);
                                           end
                                           r='';
                                           Ind = Ind+1;
                                       else
                                           r = ['(',r];
                                       end
                                    elseif IsNum==1
                                       [r, l] = strtok(l);
                                       res(Ind+1) = str2num(r);
                                       Ind = Ind+1;
                                   end

                               end
                        
                           end
                           res = res';
                       else
                           IsStruct=1;
                           while (pvalue(end) ~= ')')
                               l = strtrim(fgets(file));
                               pvalue = [pvalue,l];
                           end
                       end  %end of parameter is an array else ...
                   
               end  %end of parameter is array or structure
            end
            saveres = res;
            nameappend = '';
            %If structure: Check whether the params{cnt} contains a .
            if IsStruct==1 & (strcmpi(pname,'ACQ_jobs')==0) %for PV360
                if numel(strfind(params{cnt},'.'))~=0
                    [n,tag] = strtok(params{cnt},'.');
                    field = ['value',tag(2:end)];
                    FieldExist=0;
                    fields=fieldnames(res);
                    for fieldcnt=1:numel(fields)
                        if strcmp(fields{fieldcnt},field)
                            FieldExist=1;
                        end
                    end
                    if FieldExist==1
                        if isnumeric(res(1).(field))
                            structres =  zeros(numel(res),1);
                            for fieldcnt=1:numel(res)
                                structres(fieldcnt)=res(fieldcnt).(field);
                            end
                        else
                            structres = cell(numel(res),1);
                            for fieldcnt=1:numel(res)
                                structres{fieldcnt}=res(fieldcnt).(field);
                            end
                        end
                        res = structres;
                        nameappend=['__',tag(2:end)];

                    end
                else
                    pvalue = pvalue(2:end-1);
                    ind = 1;
                    res = struct;
                    while numel(pvalue)>0
                        [n,pvalue] = strtok(pvalue,',');
                        [num, isnum] = str2num(n);
                        if isnum == 0
                            res.(['val',num2str(ind)]) = n;
                        else
                            res.(['val',num2str(ind)]) = num;
                        end
                        ind = ind+1;
                    end
                end
            end
            %Return result
            if (numel(res)==1)
              if (numel(params)==1)
                  arr = res;
                  fclose(file);
                  return
              else
                  arr.(params{cnt})=res;
              end
            else
              % If the params-Entry contains an index: Only return
              % this element
              [p,ind] = strtok(params{cnt},'([');
              if numel(ind) > 1
                ind = strtok(ind(2:end),'])');
                index = str2num(ind);
                if index<-1
                    index=-1;
                end
                if index>Ind-1
                    index=Ind-1;
                end
                res = res(index+1);
                params{cnt} = [strtok(params{cnt},'(['),'_',num2str(index)];
              end
              %Add data to output structure of return value
              if (numel(params)==1)
                arr = res;
                fclose(file);
                return
              else
                arr.([strtok(params{cnt},'[(.'),nameappend]) = res;
              end
                
            end
            paramfound(cnt)=1;   
            end  %end of parameter found
           
          end      %end of loop through parameters
      end          %end of line contains parameter
  end  %end of next line read
end  %end of loop through file  
%Are there any not-found parameters?
if numel(params)==1
    arr = [];
else
    for cnt=1:numel(params)
        if paramfound(cnt)==0
            arr.(strtok(params{cnt},'[(.'))=[];
        end
    end
end
fclose(file);
return

function struc = BrConvertToStruc(String)

cnt = 1;
while numel(String)>0
    [item, String] = strtok(String,',');
    item = strtrim(item);
    if item(1) == '<'         %String
        item = item(2:end-1);
    elseif  numel(find((item>57 | item<43 | item==44) ...
                       & item~=101 & item~=69))==0;     %It's a number
        item = str2num(item);
    end
    struc.(['value',num2str(cnt)]) = item;
    cnt = cnt+1;
end    

return


