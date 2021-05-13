/* g is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
\\ q est l'ordre du groupe (F128)
\\ t est le nombre d'erreurs que l'on peut corriger (35)
\\ ie le code est de distance delta=2*t+1 =71
\\ m est le message

\\ Souci quand le générateur s'appelait x,
\\ g plus commode avec des polynômes.
g = ffgen(q,'a);

\\ subst(x,y,z): in expression x, replace the variable y by the expression z.
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);


\\ On cherche les erreurs dans le message, ce sont elles qui composent l'information secrète !
\\ pour les trouver, on applique la méthode de décodage des codes BCH utilisant le syndrome comme vu en cours.
\\ Une fois que l'on connait E le polynôme localisateur de l'erreur,  ses racines donneront les indices auxquels sont situées les erreurs.



\\ fonction qui construit le polynome syndrome d'un message
syndrome(msg, a, b) = {
	    delta = 2*t -1;
	    \\ le polynome syndrome est de degré <= delta
	    return (sum(l=0,delta , subst(msg, 'x, a^(b+l)) * 'x^l));
	    }


\\ Les polynomes E et R sont solution du problème de Padé: S*E=R mod (X^(2t)).
\\ La fonction de pari/gp bestapprPade donne la solution sous forme S=R/E :
\\*bestapprPade(x, {B}): returns a rational function approximation to x. This 
\\*function applies to series, polmods, and rational functions of course. 
\\ Si E(a^-i)=0, il y a une erreur en position i, et la valeur de l'erreur est déterminée en calculant :
\\ R(a^-i)/E(a^-i).


find_secret(msg)={
	for(k = 0, q, 
	      a = ffprimroot(g);
	      \\ On va parcourir les racines primitives a jusqu'à en trouver une qui convient,
	      \\ c'est à dire, qui pour j parcourant 2t indices consécutifs, a^(b+j) annule tout mot du code.
	      \\ Pour cela, il faut boucler sur la valeur de b (première puissance pour laquelle a est racine), qui est inconnue.
	      for(b = 0, q-1,
	      	    \\ on va remplir une liste initialement vide avec les erreurs détectée pour la valeur a
	      	    err = List();
		    pade = bestapprPade(Mod(syndrome(msg,a,b),x^(2*t)));
		    \\ on récupère nominateurs et dénominateurs
		    \\ le dénominateur E est le polynome localisateur recherché
		    E = denominator(pade);
		    R = numerator(pade);
		    for(i = 0, q-2,
		    	  V = subst(E,'x,a^(-i));
			  \\ si a^(-i) est racine de E (ie V==0), il y a une erreur en position i
			  if(V==0, val = subst((R/deriv(E))*(x^(b-1)), 'x , a^(-i));
			  new = fqx2int(val,g);
			  \\ on complète la liste des erreurs 
			  listput(err, new)));
	      \\ on suppose que le message fait au moins 5 caractères
	      \\ par exemple "salut"
	      if(#err >5, return (Strchr(Vec(err))));
	      );
	);
}


\\ On récupère le message sous forme de polynome
mot_recu = int2fqx(m, g);
\\ On retrouve les erreurs
secret=find_secret(mot_recu);
print(secret);

