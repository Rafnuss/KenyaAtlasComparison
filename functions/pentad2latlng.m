function [lat0, lng0] = pentad2latlng(p)

	lat0 = str2double(p(1:2)) + str2double(p(3:4)) / 60;
	lng0 = str2double(p(6:7)) + str2double(p(8:9)) / 60;
	if (p(5) == "_" || p(5) == "a")
		lat0 = -lat0;
    end
	if (p(5) == "b" || p(5) == "a") 
		lng0 = -lng0;
    end

end