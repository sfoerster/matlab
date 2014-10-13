function den = getDen(n)
	
	if n==1
		den = 1;
		return;
	end
	
	prev = getDen(n-1);
	
	den = [prev n*[1 prev]];
	
end