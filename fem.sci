//                           Equipe
//                     AGOSSOU Hermann G1B
//                    AMANI Jean-Yves G1B


//définition du second membre
function y = f(x)
    y = cos(x)
endfunction


// introduction du vecteur X
x = [0:0.1:1] 


//tracé du second membre
//plot2d(x,f(x))


// définition de la solution théorique
function y= u(x)
    y = -1 + (1-cos(1))*x + cos(x)
endfunction

//tracé de cette solution
plot2d(x,u(x))  //ca marche



//introduction de la dimension et du vecteur maillage
n = 5
disp("")
disp("dimension n = ")
disp(n)
disp("")

//definition du pas du maillage
p = 1/(n+1)

//
X = linspace(0,1,n+2)




//définition des fonction phi(i,)
function y = phi(i,x,X)
    y = ((x-X(i))/(X(i+1)-X(i))).*bool2s(X(i)<x & x<=X(i+1)) + ((x-X(i+2))/(X(i+1)-X(i+2))).*bool2s(X(i+1)<x & x<=X(i+2))
endfunction

//tracé de la fonction phi(i=3 , ) cas i = 3
//x0 = [0:0.00001:1]
//i=3
//plot2d(x0,phi(i,x0,X))  //ca marche



// définition de la matrice de rigidité A
T = zeros(n,n)
T0 = T(1,:)
T0(1)= 2*(n+1)
T0(2) = -(n+1)
A = toeplitz(T0)

//affichage de deux elements de A
disp("affichage de la matriice de rigidité A")
disp("A = ")
disp(A)
disp("")


// definition de la fonction phi(i,x,X).cos(x))
function y = phiCos(n,j,x,X)
    y = phi(j,x,X).*cos(x)
endfunction
//affichage de la fonction phiCos(n=5,i=3)
//n = 5
//j = 2
//x0 =0:0.001:1
//plot2d(x0,phiCos(n,j,x0,X)) 



// initialisation du second membre, le vecteur B
B = linspace(1,n,n)
for j =1:n
    B0 = integrate('phiCos(n,j,x,X)','x',0,[1])
    B(j) = B0(1)
end

//affichage du vecteur B
disp("")
disp("affichage du vecteur B")
disp("B = ")
disp(B) 
disp("")


//résolution du système linéaire AU=B
sp = sparse(A)
Un = umfpack(sp,"\",B') //Un est un vecteur contement les valeurs de la solution approchée "un" aux points du maillage

//affichage des composantes de la solution du approchée
//disp(Un) 


//tracé de la solution approchée "un" avec le vecteur Un
Un = Un'
//X1 = linspace(p,1-p,n) //un vecteur des points du maillage sauf 0 et 1 donc de x1 à xn
//plot2d(X1,Un) 


// définition de la fonction un
X = linspace(0,1,n+2)
function y = un(x)
    s = 0
    for j = 1:n
        s = s + Un(j)*phi(j,x,X)
    end
    y = s
endfunction

//tracé de la solution approchée "un" avec la fonction un
//x2 = 0:0.00001:1
//plot2d(x2,un(x2)) 
    





////Approximation de l'erreur en norme L2

I = 0 //initialisation de la valeur de l'intégrale de (u - un)^2 
for i = 0:n
    I = I + (   u(i/(n+1)) - un(i/(n+1))   )^2
end
I = I/(n+1)

// calcul de la norme L2 de (u - un)
Err_L2 = sqrt(I)
disp("")
disp("Erreur en norme L2")
disp("Err_L2 = ")
disp(Err_L2)
disp("")




////Approximation de l'erreur en norme H1
    
      //calcul du vecteur des évaluations de UnPrime aux points du maillage
UnPrime         = zeros(1,n+2)
UnPrime(1)      = + (un(p))/p
UnPrime(n+2)    = - (un(1-p))/p
for j = 2:(n+1)
    UnPrime(j)  =(  un((j+1)*p) )/p - ( un(j*p) )/p
end

      //calcul de la valeur de l'intégrale de (uPrime - unPrime)^2 
I = 0 //initialisation de la valeur de l'intégrale de (uPrime - UnPrime)^2 
for i = 0:n
    I = I + (   u(i/(n+1)) - UnPrime(i+1)   )^2
end
I = I/(n+1)

      //calcul de la norme H1 de (u - un)
Err_H1 = sqrt(Err_L2^2 + I)
disp("")
disp("Erreur en norme H1")
disp("Err_H1 = ")
disp(Err_H1)
disp("")

// tracé de l'erreur = fct(n) en échelle logarithmique ie fct(n)=log(Err_H1(10^n))
plot2d(["ln"],[5 10 20 40 80 160 320 640], [0.2796968 0.2665998 0.2609445 0.2014898  0.1231447 0.1088801 0.1100329   0.1104319 ]) 

//Commentaire
//On n'obtient pas une pente de -1 ou inférieur à -1 donc ce n'est pas conforme à la théorie
