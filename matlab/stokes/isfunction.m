function [ result ] = isfunction( function_string )
%% Returns whether function_string is a known function.
% From
% http://commoditymodels.files.wordpress.com/2010/03/matlab-parallel-computing-toolbox-and-interpolating-futures-curves.pdf
% Will allow backwards compatibility of scripts with Parallel Computing
% Toolbox commands on older versions of Matlab, or versions without the
% toolbox add-on installed.
result = 0; %#ok<NASGU>
% If it's a keyword, I'd say 'yes'.  Although some might dispute this.
if iskeyword(function_string)
    result = 1;
else
    fh = str2func(function_string);    
f = functions(fh);   
% If it's a function, functions() will return the file it's in.
if (~isempty(f.file))
        result = 1;
    else
        result = 0;
    end
end
end