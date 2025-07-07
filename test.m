barB = 0.9;
gap = 0.2;
LHS = (barB+gap)/(1-gap)
RHS = exp(7*gap)*barB
if(LHS <= RHS)
    fprintf('Correct Choice \n ');
else
    fprintf('Wrong Choice \n ');
end
