function varargout = get_outputs(fn, ixsOutputs)
% see http://stackoverflow.com/questions/7921133/anonymous-functions-calling-functions-with-multiple-output-forms
output_cell = cell(1,max(ixsOutputs));
[output_cell{:}] = (fn());
varargout = output_cell(ixsOutputs);