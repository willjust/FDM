%% Spectroscopy, NE412
e(1) = Element('Galium-72', 834.0, 2201.7, 629.9, 2507.8, 894.2, 0); 
e(2) = Element('Silver-110m',657.7,884.2,937.4,1383.9, 763.8, 0);
e(3) = Element('Iron-59', 5.9, 6.5, 1099.2,1291.6, 192.3, 0);
e(4) = Element('Arsenic-76',  559.1, 657, 1216, 1212.7, 1228.5, 0); 
e(5) = Element('Gold-198', 411.8, 70.8, 675.9, 68.9, 80.2, 0);
e(6) = Element('Nickel', 511, 810.8, 863.9, -5, -5, 0);

e(7) = Element('Sm-155', 104.3, 41.51, 245.7, 141.4, 78.5, 0);
e(8) = Element('U', 99.5, 103.7, 106.4, 228.1, 277.9, 0);
e(9) = Element('Mg', 1368.6, 2754.1, -5, -5, -5,0)
e(10) = Element('Sc', 889.3, 1120.5, 142.5, -5, -5,0)
e(11) = Element('', , , , , ,0)

sample = [834.2 630 74.22 2205.3 601.1 894.4;
		  657.8 884.8 1873.2 764.0 1385.2 706.8;
		  74.49 411.8 85.23 1099.4 1292.3 1461.6;
		  0 0 0 0 0 0; %we did not measure sample 4
		  559.2 657.1 1216.1 1229.1 42.56 111.6;
		  411.9 70.29 80.63 250.3 675.9 42.44;
		  810.9 411.9 511.0 74.35 1461.7 0];

montana = [41.51 103.5 1369.5 889.4  564.2 69.11 1120.9 2761 487.2 1597.6 329 1099.7 145.6 312 603 1292 815.8 280.5 411.7 159.5 396.8 209.1 320.2 2246.9 868];

N_sample = size(sample,1);
N_elem = size(e,2);
percent = zeros(N_sample, N_elem);
error = 1;

% Calculate the match with each element
for i=1:N_sample
	for j=1:N_elem
		for k = 1:6
			if(abs(e(j).peak1 - sample(i,k)) < error)
				percent(i,j) = percent(i,j) + 10;
			end
			if(abs(e(j).peak2 - sample(i,k)) < error)
				percent(i,j) = percent(i,j) + 10;
			end
			if(abs(e(j).peak3 - sample(i,k)) < error)
				percent(i,j) = percent(i,j) + 10;
			end
			if(abs(e(j).peak4 - sample(i,k)) < error)
				percent(i,j) = percent(i,j) + 10;
			end
			if(abs(e(j).peak5 - sample(i,k)) < error)
				percent(i,j) = percent(i,j) + 10;
			end
		end
	end
end

monSoilElem = zeros(1,N_elem);

for j=1:N_elem
	for k = 1:25
		if(abs(e(j).peak1 - montana(k)) < error)
			monSoilElem(j) = monSoilElem(j) + 10;
		end
		if(abs(e(j).peak2 - montana(k)) < error)
			monSoilElem(j) = monSoilElem(j) + 10;
		end
		if(abs(e(j).peak3 - montana(k)) < error)
			monSoilElem(j) = monSoilElem(j) + 10;
		end
		if(abs(e(j).peak4 - montana(k)) < error)
			monSoilElem(j) = monSoilElem(j) + 10;
		end
		if(abs(e(j).peak5 - montana(k)) < error)
			monSoilElem(j) = monSoilElem(j) + 10;
		end
	end
end