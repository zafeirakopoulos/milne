Milne := module()

export
    solveMilne,
    computeVF,
    computeInitialBox,
    computePRS,
    subdiv2,
    per,
    computeVF_multires;
  
local 
    MaxBoundRoots,
    primpartPRS;
    
#############################################3
    # Primitive-Part Sequence.
    primpartPRS := proc (f1, f2, x)
        #
        local g1, g2, g3, m, q, PRS;
        #
        g1 := sort( collect( f1, x ), x);
        g2 := sort( collect( f2, x ), x);
        g3 := sort( collect( primpart( prem( g1, g2, x, 'm', 'q' ), x ), x ), x );
        if ( signum( lcoeff( m ) ) > 0) then g3 := -g3; fi;
        if ( g3 = 0 ) then
            return [g1, g2];
        end if;
        PRS := [g1, g2, g3];
        while true do
            g1 := g2;
            g2 := g3;
            g3 := sort( collect( primpart( prem( g1, g2, x, 'm', 'q' ), x ), x ), x );
            if ( signum( lcoeff( m ) ) > 0) then g3 := -g3; fi;
            if ( g3 = 0 ) then
                return PRS;
            else
                PRS := [op(PRS), g3];
            end if;
        end do;

        return PRS;
    end proc:

############################################
per := proc(L)

local tmp, v, i;

if (nops(L)=0) then
	return 0;
end if;
#print(L);
v:=0:
tmp:=op(1,L);
for i from 2 to nops(L) do
if ((op(i,L)<>0) and (tmp*op(i,L)>=0)) then
v:=v+1:
end if;
tmp:=op(i,L):
end do;
#print(v);
return v;
end proc:


############################################
solveMilne := proc (h1,h2, sel)

local vf, bbox, M, sols, time1, time2, time3, time4, time5  ;

time1:=time();
vf:=computeVF(h1,h2):
time2:=time();
bbox:=computeInitialBox(h1,h2):
time3:=time();
M:=computePRS(op(1,vf),op(2,vf), sel):
time4:=time();
sols:=subdiv2(M, bbox):
time5:=time();
return [sols, time2-time1, time3-time2, time4-time3, time5-time4, time5-time1];
end proc:


MaxBoundRoots:=proc(p)
    local _c,_r,tmp;
    _c := [coeffs(p)];
    tmp :=1;
    for _r in _c do;
        tmp:=max(tmp, abs(op(1,_c)/_r));
    end do;
    return tmp;
end proc;

############################################
computeVF_multires:=proc(h1,h2)

local S0, S1, upol:


upol:=u+(x-t)*(y-s):
S0:=Determinant(Matrix(mbezout([h1,h2,upol], [x,y]))):
S0:=simplify(S0/u^ldegree(S0,u)):
S1:=diff(S0,u):
return [S0,S1]:

end proc:

############################################
computeVF:=proc(h1,h2)

local S0, S1, upol:

upol:=u+(x-t)*(y-s):
S0:= resultant(resultant(h1,upol,y),resultant(h2,upol,y),x):
S0:=simplify(S0/u^ldegree(S0,u)):
S1:=diff(S0,u):
return [S0,S1]:
end proc:


#########################
computeInitialBox:=proc(h1,h2)

local b1, b2 :

b1:=MaxBoundRoots(resultant(h1,h2,y),x):
b2:=MaxBoundRoots(resultant(h1,h2,x),y):
return [-(ceil(abs(b1))),ceil(abs(b1)+1), -(ceil(abs(b2))),ceil(abs(b2))+1]:
end proc:


#########################
computePRS:=proc(S0, S1, selector)

if (selector = 1)  then
	return subs({u=0},PRS:-StHa(S0,S1,u)):
elif (selector = 2)  then
	return subs({u=0},PRS:-subresPRS(S0,S1,u)):
elif (selector = 3)  then
	return subs({u=0},PRS:-primpartPRS(S0,S1,u)):
end if;
end proc:


################################
subdiv2:=proc(M, bbox)

local  lob, box, E, sols,  loops, M1, M2, M3, M4:


sols:=[];

lob:=[bbox]:
#print("In subdivision");
loops:=0: 
#while ((nops(lob)>0) and (max_loops>0)) do
while (nops(lob)>0)  do
loops:=loops+1:
box:=op(1,lob):

M1:=subs({t=op(2,box),s=op(3,box)}, M);
M2:=subs({t=op(1,box),s=op(4,box)}, M);
M3:=subs({t=op(2,box),s=op(4,box)}, M);
M4:=subs({t=op(1,box),s=op(3,box)}, M);

if ((op(nops(M1),M1)=0) or (op(nops(M2),M2)=0) or (op(nops(M3),M3)=0) or (op(nops(M4),M4)=0)) then
print("FAIL! Last element was 0...");
return [];
end  if;

#print("M1",M1);
#print("M2",M2);
#print("M3",M3);
#print("M4",M4);


E:=-0.5*(per(M1)+per(M2)-per(M3)-per(M4)):
#print("E",E);
if (E<>0) then
if (E=1) then
#printf("* Single root in: [%G,%G] x [%G,%G]\n", op(1,box), op(2,box), op(3,box), op(4,box)):
sols:=[op(sols),[op(1,box), op(2,box), op(3,box), op(4,box)]]:
else
lob:=[op(lob), [op(1,box)  , op(1,box)+((op(2,box)-op(1,box))/2)	, op(3,box)+((op(4,box)-op(3,box))/2) ,op(4,box)]]:
lob:=[op(lob), [op(1,box)+((op(2,box)-op(1,box))/2), op(2,box)     , op(3,box)+((op(4,box)-op(3,box))/2) ,op(4,box)]]:
lob:=[op(lob), [op(1,box)  , op(1,box)+((op(2,box)-op(1,box))/2)   , op(3,box) , op(3,box)+((op(4,box)-op(3,box))/2)]]:
lob:=[op(lob), [op(1,box)+((op(2,box)-op(1,box))/2), op(2,box)     , op(3,box) , op(3,box)+((op(4,box)-op(3,box))/2)]]:
end if:
end if:
lob:=subsop(1=NULL, lob):
end do:

#print("loops", loops);
return sols:
end proc:

end module: # end Milne
