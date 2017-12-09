classdef Element
	properties
		name = ''
		peak1 = 0;
		peak2 = 0;
		peak3 = 0;
		peak4 = 0;
		peak5 = 0;
		background = 0;
	end
	
	methods
		function e = Element(name, peak1, peak2, peak3, peak4, peak5, back)
			e.name = name;
			e.peak1 = peak1; e.peak2 = peak2;
			e.peak3 = peak3; e.peak4=peak4; e.peak5 = peak5;
			e.background = 0;
		end
	end
end