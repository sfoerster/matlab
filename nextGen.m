
function [out,cbin] = nextGen(dist,bin,cbin)
	
	
	out = 0;
	cbin = 1;
	for n=1:numel(dist)
		
		temp = dist(n)*cbin;
		cbin = compbinom(cbin,bin);
		
		out = paddadd(out,temp);
		%shis = paddmatadd(shis,out);
	end
	
	
	
end


function out = paddadd(shortvec,longvec)

	out = zeros(size(longvec));
	
	out(1:numel(shortvec)) = shortvec;
	
	out = out + longvec;	
	
end

function cbin = compbinom(cbin,bin)

	if isempty(cbin)
		cbin = bin;
	else
		cbin = conv(cbin,bin);
	end	
	
end
