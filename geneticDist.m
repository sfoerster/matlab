
function [dist,shis,ex] = geneticDist(initial,p,trialsPerGen,g)
	
	bin = binom(p,trialsPerGen);
	
	shis = initial;
	cbin = [];
	dist = initial;
	for n=2:g
		disp(['Processing Generation: ' num2str(n)])
		[dist,cbin] = nextGen(dist,bin,cbin);
		shis = paddmatadd(shis,dist);
	end
	
	ex = sum(repmat(0:size(shis,2)-1,size(shis,1),1) .* shis,2);
	
	
	
end

function out = paddmatadd(shis,dist)

	[m,n] = size(shis);
	
	temp = zeros(m,numel(dist));
	temp(1:m,1:n) = shis;
	
	out = [temp; dist];	
	
end