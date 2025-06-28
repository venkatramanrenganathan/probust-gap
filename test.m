a = 0.95;
e = 0.1;
if(e <= a)
    LHS = (a+e)/(1-e)
    RHS = (1.75/(1.01-a)) + e
else
    fprintf('Wrong Choice');
end
