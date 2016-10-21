function pgfplot_tbl(tbl_data,N,tbl_name,format,options)
% pgfplot_tbl(tbl_data,N,tbl_name,format,options)
% pgfplot_tbl creates table datafiles that can be read by pgfplots
% The data is given as a struct where each field in the struct is
% considered a column in the table. The field name of the struct is used to
% access the data of a particular column in pgfplots. The fields in the
% struct have to be vectors
%
% Input:
%         tbl_data: struct containing the data
%                N: number of rows to write (not counting the header)
%                   empty argument uses the length of the first field
%         tbl_name: filename (including path) of the output file
% Optional input arguments:
%           format: struct with strings. The fields of the struct have the
%           same name as the fields in 'tbl_data'. Each field in 'format'
%           is a fprintf-compatible format-string. 
%           format.(field) = '%f' (default for all fields)
%
%         options: struct with options (defaults are given)
%                  options.no_header = 0, write header
%                  options.row_sep   = '', separator used at the end of a line 
% 
%
% See also fprintf

%Get column headers
col_names = fieldnames(tbl_data);

if(~exist('format','var') || isempty(format))
  format  = default_formatting(col_names); 
end
%Load default options for unspecified fields
if(~exist('options','var'))
   options = struct();
end
options = default_options(options);
[fID msg] = fopen(tbl_name,'w');
assert(isempty(msg),msg);
if(~options.no_header)
  write_header(fID,col_names,options);
end

if(isempty(N) || N <= 0)
  N = length(tbl_data.(col_names{1}));
end

write_body(fID,tbl_data,col_names,format,N,options)

fclose(fID);

disp(['Wrote pfgplot data to ' tbl_name ]);

end

function format = default_formatting(col_names)
  format = struct();
  for i=1:length(col_names)
      format.(col_names{i}) = '%d';
  end

end

function options = default_options(options)
  no_header = 0;
  row_sep = '';
  if(~isfield(options,'no_header'))
    options.no_header = no_header;
  end
  if(~isfield(options,'row_sep'))
    options.row_sep = row_sep;
  end
end

function write_header(fID,col_names,options)
  for i=1:length(col_names)
    fprintf(fID,'%s ',col_names{i});
  end
  fprintf(fID,'%s \n',options.row_sep);
end

function write_body(fID,tbl_data,col_names,format,N,options)
for i=1:N
  for j=1:length(col_names)
    len = length(tbl_data.(col_names{j}));
    assert(len >= N,[col_names{j} ...
    ' is out of bounds.' 'length(' col_names{j} ') = ' num2str(len) ...
    ', N = ' num2str(N) ]);
    fprintf(fID,[format.(col_names{j}) ' '],tbl_data.(col_names{j})(i));
  end
    fprintf(fID,'%s \n',options.row_sep);
end

end

