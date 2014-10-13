
function bin = binom(p,n)
	
	bin = [];
	for k=0:n
		bin(k+1) = nchoosek(n,k)*(p^k)*((1-p)^(n-k));
	end
	
end