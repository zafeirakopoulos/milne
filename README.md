# milne
Isolating real roots of bivariate polynomials using Milne's volume function.

The Milne package uses SLV. 

How to use
-----------

 
    MYPATH:="{path to the code}/":
    SLV_PATH:=cat(MYPATH,"SLV/"):
    read cat(MYPATH, "Milne.mpl"):
    read cat(SLV_PATH, "SLV.mpl"):

    h1:=2*x^2+3*x*y+y^2-3:
    h2:=4*x^2+7*y^2+9*y:
    Milne:-solveMilne(h1,h2,1);


Results to 

                           2            2    
                  h1 := 2 x  + 3 x y + y  - 3

                             2      2      
                    h2 := 4 x  + 7 y  + 9 y

    [[[    -3  -5  1]  [-3  1  -5  1]]                                   ]
    [[[-2, --, --, -], [--, -, --, -]], 0.030, 0.001, 0.001, 0.003, 0.035]
    [[[    4   4   2]  [4   2  4   2]]                                   ]

Except if your computer is 15 years old (as the code is essentially) and the times are 
closer to those (from the original readme file).

    [[[-5  -3  -5  -3]  [-3  1  -5  -3]]                                   ]
    [[[--, --, --, --], [--, -, --, --]], 0.120, 0.016, 0.084, 0.072, 0.292]
    [[[4   8   4   8 ]  [8   2  4   8 ]]                                   ]

The output [ [a_1,b_1,c_1,d_1], [a_2,b_2,c_2,d_2],...,[a_n,b_n,c_n,d_n], t_1, t_2, t3, t4, t5] (soon to be in better format) 
means that there is a (unique up to multiplicity) real root in the box [a_i,b_i]x[c_i,d_i]. 
The times are for computing the volume function, the initial box containing all real roots, a polynomial remainder sequence for the volume function, and the 
subdivision for the real root isolation respectively. 



 
