function LHS = BDF2_LHS(u,um1,um2)
    LHS = u - (4/3)*um1 + (1/3)*um2;
end