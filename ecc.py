# Sorbonne Université LU3IN024 2024-2025
# TME 5 : Cryptographie à base de courbes elliptiques
#
# Etudiant.e 1 : KOWAL 28606347
# Etudiant.e 2 : ALVES 21105811

from math import sqrt
import matplotlib.pyplot as plt
from random import randint
import random

# Fonctions utiles

def exp(a, N, p):
    """Renvoie a**N % p par exponentiation rapide."""
    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L
    res = 1
    for Ni in binaire(N):
        res = (res * res) % p
        if (Ni == 1):
            res = (res * a) % p
    return res


def factor(n):
    """ Return the list of couples (p, a_p) where p is a prime divisor of n and
    a_p is the p-adic valuation of n. """
    def factor_gen(n):
        j = 2
        while n > 1:
            for i in range(j, int(sqrt(n)) + 1):
                if n % i == 0:
                    n //= i
                    j = i
                    yield i
                    break
            else:
                if n > 1:
                    yield n
                    break

    factors_with_multiplicity = list(factor_gen(n))
    factors_set = set(factors_with_multiplicity)

    return [(p, factors_with_multiplicity.count(p)) for p in factors_set]


def inv_mod(x, p):
    """Renvoie l'inverse de x modulo p."""
    return exp(x, p-2, p)


def racine_carree(a, p):
    """Renvoie une racine carrée de a mod p si p = 3 mod 4."""
    assert p % 4 == 3, "erreur: p != 3 mod 4"

    return exp(a, (p + 1) // 4, p)


# Fonctions demandées dans le TME

def est_elliptique(E):
    """
    Renvoie True si la courbe E est elliptique et False sinon.

    E : un triplet (p, a, b) représentant la courbe d'équation
    y^2 = x^3 + ax + b sur F_p, p > 3
    """
    p,a,b=E
    return (-(4*a**3+27*b**2)%p!=0)


def point_sur_courbe(P, E):
    """Renvoie True si le point P appartient à la courbe E et False sinon."""
    if ()==P:
        return True
    x,y=P
    p,a,b=E

    e1= y**2 % p
    e2= (x**3+a*x+b)%p
    return e1==e2


def symbole_legendre(a, p):
    """Renvoie le symbole de Legendre de a mod p."""
    N=(p-1)//2
    e=exp(a,N,p)

    if e==1:
        return 1

    if e==p-1:
        return -1 % p

    return 0


def cardinal(E):
    """Renvoie le cardinal du groupe de points de la courbe E."""
    p,a,b=E
    res = 1

    for x in range(p):
        v = (x**3 + a*x + b) % p
        sl=symbole_legendre(v,p)

        if sl==1:
            res+=2

        if sl==0:
            res+=1
            
    return res

def liste_points(E):
    """Renvoie la liste des points de la courbe elliptique E."""
    p, a, b = E
    
    res = [()]
    
    assert p % 4 == 3, "erreur: p n'est pas congru à 3 mod 4."

    for x in range(p):
        rhs = x**3 + a*x + b
        if symbole_legendre(rhs,p) == 1:
            racine = racine_carree(rhs,p)
            res.append((x,racine))
            res.append((x,(-racine)%p))
        elif symbole_legendre(rhs,p) == 0:
            res.append((x,0))
    return res


def cardinaux_courbes(p):
    """
    Renvoie la distribution des cardinaux des courbes elliptiques définies sur F_p.

    Renvoie un dictionnaire D où D[i] contient le nombre de courbes elliptiques
    de cardinal i sur F_p.
    """
    borne_min = int(p + 1 - 2 * sqrt(p))
    borne_max = int(p + 1 + 2 * sqrt(p))
    D = {c: 0 for c in range(borne_min, borne_max + 1)}
    
    for a in range(p):
        for b in range(p):
            delta = -16*(4*(a**3) + 27*(b**2))
            if delta % p == 0:
                continue
            E = (p,a,b)
            c = cardinal(E)
            D[c] = D.get(c,0) + 1
    D = {k: v for k, v in D.items() if v > 0}
    return D


def dessine_graphe(p):
    """Dessine le graphe de répartition des cardinaux des courbes elliptiques définies sur F_p."""
    bound = int(2 * sqrt(p))
    C = [c for c in range(p + 1 - bound, p + 1 + bound + 1)]
    D = cardinaux_courbes(p)

    plt.bar(C, [D[c] for c in C], color='b')
    plt.show()


def moins(P, p):
    """Retourne l'opposé du point P mod p."""
    if P == ():
        return ()
    return (P[0], (-P[1])%p)

def est_egal(P1, P2, p):
    """Teste l'égalité de deux points mod p."""
    if P1 == () and P2 == ():
        return True
    if P1 == () and P2 != ():
        return False
    if P1 != () and P2 == ():
        return False
    return (P1[0]%p == P2[0]%p) and (P1[1]%p == P2[1]%p)

def est_zero(P):
    """Teste si un point est égal au point à l'infini."""
    return P == ()



def addition(P1, P2, E):
    """Renvoie P1 + P2 sur la courbe E."""
    p, a, b = E
    if P1 == () and P2 == ():
        return ()
    if P1 == () and P2 != ():
        return P2
    if P1 != () and P2 == ():
        return P1
    if P1 == moins(P2,p):
        return ()
    if P1 == P2:
        ld = ((3*(P1[0]**2)+a)*inv_mod((2*P1[1]),p)) % p
        x_res = (ld**2 - 2*P1[0])%p
        y_res = (ld*(P1[0]-x_res)-P1[1])%p
        return (x_res,y_res)
    else:
        ld = ((P2[1]-P1[1])*inv_mod(P2[0]-P1[0],p))%p
        x_res = (ld**2 - P1[0] - P2[0])%p
        y_res = (ld*(P1[0]-x_res)-P1[1])%p
        return (x_res,y_res)
    


def multiplication_scalaire(k, P, E):
    """Renvoie la multiplication scalaire k*P sur la courbe E."""   
    if P == ():
        return ()
    if k<0:
        return moins(multiplication_scalaire(-k,P,E),E[0])
    p, _, _ = E
    P = (P[0]%p, P[1]%p)
    R = ()
    Q = P
    while k>0:
        if k%2 == 1:
            R = addition(R,Q,E)
        Q = addition(Q, Q, E)
        k = k // 2
    return R


def generate_divisors(factors):
    """Génère les diviseurs à partir d'une factorisation [(p1,e1), (p2,e2), ...]"""
    def backtrack(i, current):
        if i == len(factors):
            yield current
            return
        p, e = factors[i]
        for exp in range(e + 1):
            yield from backtrack(i + 1, current * (p ** exp))
    return sorted(set(backtrack(0, 1)))

def ordre(N, factors_N, P, E):
    """Renvoie l'ordre du point P dans les points de la courbe E mod p. 
    N est le nombre de points de E sur Fp.
    factors_N est la factorisation de N en produit de facteurs premiers."""
    if P == ():
        return 1
    
    for d in generate_divisors(factors_N):
        if multiplication_scalaire(d,P,E) == ():
            return d



def point_aleatoire_naif(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""
    p, a, b = E
    while True:
        x = random.randrange(p)
        y = random.randrange(p)
        print("Point courant :",(x,y))
        if (y * y) % p == (x**3 + a * x + b) % p:
            return (x, y)




def point_aleatoire(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""
    p, a, b = E
    while True:
        x = random.randrange(p)
        rhs = (x**3 + a * x + b) % p       
        try:
            y = racine_carree(rhs, p)
            return (x, y)
        except ValueError:
            continue


def point_ordre(E, N, factors_N, n):
    """Renvoie un point aléatoire d'ordre N sur la courbe E.
    Ne vérifie pas que n divise N."""
    while True:
        point = point_aleatoire(E)
        if ordre(N,factors_N,point,E) == n:
            return point

def keygen_DH(P, E, n):
    """Génère une clé publique et une clé privée pour un échange Diffie-Hellman.
    P est un point d'ordre n sur la courbe E.
    """
    sec = random.randrange(1,n)
    pub = multiplication_scalaire(sec,P,E)
    
    return (sec, pub)


def echange_DH(sec_A, pub_B, E):
    """Renvoie la clé commune à l'issue d'un échange Diffie-Hellman.
    sec_A est l'entier secret d'Alice et pub_b est l'entier public de Bob."""

    return multiplication_scalaire(sec_A,pub_B,E)
