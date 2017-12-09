function [ y ] = printc( string )
%PRINTC Print centered string
	SLASH80 = '///////////////////////////////////////////////////////////////////////////\n';
	comSize = get(0, 'CommandWindowSize');
	WIDTH = comSize(1);
	ml = length(string);
	fprintf(SLASH80);
	WIDTH=80;
	fprintf(blanks(floor((WIDTH-ml) / 2)));
	fprintf('%s\n', string);
	fprintf(SLASH80);
	y = 0;
end