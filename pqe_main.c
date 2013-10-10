/**********************************************************************/
/* PROGRAMMATION QUADRATIQUE EN NOMBRES ENTIERS                       */
/* Relaxation aggrégée (article Djerdjour, Mathur, Salkin, ORL, 1987) */
/* D. Quadri et E. Soutif                                             */
/**********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/***********************/
/* Types et constantes */
/***********************/

#include "pqe_types.h"

/* def du type vecteur_x nécessaire pour transférer les données  à y_ik*/
typedef struct {double val_x; bool modifier;} vecteur_x;

/* def du type vecteur_y nécessaire pour récupérer les données relatives à y_ik*/
typedef struct {int val_y ; bool modifier;} vecteur_y;

typedef vecteur_x elt_x;
typedef vecteur_y elt_y;

typedef st_tri_coeff elt_coeff;
typedef elt_coeff * coeff_1D;
typedef elt_x * eltx_1D;
typedef elt_y ** elty_2D;

typedef struct _contrainte {
    int    ind_cont;
    bool   triee;
    elt_1D tab;
} contrainte;

typedef contrainte * contrainte_1D;

#include "pqe_const.h"
const double epsilon = 0.000001;
const double precision_double = 0.00001; 
const double infini = 1E+8;
const int    diminue_pas = 5;


typedef struct _cellule {
    double           ecart;
    int              ind_cont;
    struct _cellule *suivant;
} cellule;

typedef cellule *p_cellule;

/**********************/     
/* Variables globales */
/**********************/     

#include "pqe_var.h"


/*************************************/
/* Prototypes des fonctions externes */
/*************************************/

#include "pqe_entete_cplex.h"
#include "pqe_prototypes.h"



void echanger (elt_1D tab, int i, int j)
{
    /* Echange les elts de tab situés aux indices i et j */
    elt temp;
    
    temp   = tab[i];
    tab[i] = tab[j];
    tab[j] = temp;
    
} /* echanger */


void pqe_main_tri_rapide (int g, int d, elt_1D tab)
{
    /* tri de la partie g..d du tableau */
    int i = g;
    int j = d+1;
    int pivot = tab[g].cle;
    
    
    if (g < d) {
        while (i < (j-1)) {
            /* de g a i (inclus) tous les elements sont plus petits ou = pivot
            de j (inclus) a d tous les elements sont plus grands que pivot */
            i = i + 1;
            j = j - 1;
            while ((i < j) && (tab[i].cle <= pivot)) i++;
            while ((j > i) && (tab[j].cle > pivot))  j--;
            if (i < j) echanger (tab, i, j);
            else {
                /* cas i=j, retablir l'invariant */
                if (tab[i].cle <= pivot) j++;
                else i--;
            }
        }
        echanger (tab, g, i); /* pivot est maintenant a la position i */
        if ((i-1) > g) pqe_main_tri_rapide (g, i-1, tab);
        if ((i+1) < d) pqe_main_tri_rapide (i+1, d, tab);
    }
} /* pqe_main_tri_rapide */

void echanger_coeff (coeff_1D tab, int i, int j)
{
    /* Echange les elts de tab situés aux indices i et j */
    elt_coeff temp;
    
    temp   = tab[i];
    tab[i] = tab[j];
    tab[j] = temp;
    
} /* echanger */

void tri_rapide_coeff (int g, int d, coeff_1D tab)
{
    /* tri de la partie g..d du tableau */
    int i = g;
    int j = d+1;
    int pivot = tab[g].cle;
    
    
    if (g < d) {
        while (i < (j-1)) {
            /* de g a i (inclus) tous les elements sont plus petits ou = pivot
            de j (inclus) a d tous les elements sont plus grands que pivot */
            i = i + 1;
            j = j - 1;
            while ((i < j) && (tab[i].cle <= pivot)) i++;
            while ((j > i) && (tab[j].cle > pivot))  j--;
            if (i < j) echanger_coeff (tab, i, j);
            else {
                /* cas i=j, retablir l'invariant */
                if (tab[i].cle <= pivot) j++;
                else i--;
            }
        }
        echanger_coeff (tab, g, i); /* pivot est maintenant a la position i */
        if ((i-1) > g) tri_rapide_coeff (g, i-1, tab);
        if ((i+1) < d) tri_rapide_coeff (i+1, d, tab);
    }
}/* tri_rapide_coeff */




/******************************************************************/

void mettre (int inf, int sup, double *capac, double *res, double *p, 
double *w, double *solution, st_tri_kp *tri)
{
    int i, j;
    double  sous_tot = 0.0;
    
    for (i = inf; i <= sup; i++)
    {
        j = tri[i].ind;
        sous_tot += w[j];
        solution[j] = 1.0;
        *res += p[j];
        print("*res = %f j = %d, p[%d] = %f\n", *res, j, j, p[j]);
    }
    *capac -= sous_tot;
    
} /* mettre */

/******************************************************************/

bool rentre (int inf, int sup, double capac, double *w, st_tri_kp *tri)
{
    int i, j;
    double sous_tot = 0.0;
    
    if (inf > sup) return VRAI; /* Necessaire et rien n'empeche de mettre du vide */
    
    for (i = inf; i <= sup; i++)
    {
        j = tri[i].ind;
        sous_tot += w[j];
    }
    return (sous_tot <= capac);
    
} /* rentre */

/******************************************************************/

void couper (int i, double capac, double *p_res, double *p, double *w, 
             double *solution, st_tri_kp *tri)
{
    int j;
    
    j            = tri[i].ind;
    solution[j]  = capac/ w[j];
    *p_res      += p[j] * solution[j];
    
} /* couper */

/******************************************************************/

void echange (int A, int B, st_tri_kp *tri)
{
    st_tri_kp  echange;
    
    if (A != B)
    {
        echange = tri[A];
        tri[A]  = tri[B];
        tri[B]  = echange;
    }
} /* echange */

/******************************************************************/

void met_a_zero (int inf, int sup, double *solution, st_tri_kp *tri)
{
    int i, j;
    
    for (i = inf; i <=sup; i++)
    {
        j = tri[i].ind;
        solution[j] = 0.0;
    }
} /* met_a_zero */

/******************************************************************/
/* Résolution du sac à dos linéaire version générale              */
/* c'est à dire dès le départ la fction objectif est linéaire     */
/******************************************************************/

void reso_kp (int nbvar, double *p, double *w, double cap, 
              int LIM_INF, int LIM_SUP, double *p_resu, double *solution, 
              int *indice_critique, st_tri_kp *tri)
{
    /* ------------------------------------------------------------------------------*/
    /* nbvar : nb de variables (n) p : coeff fonction obj,                           */
    /* w : coeff contrainte en y_ik                                                  */
    /* de capacite, cap : capacite, *p_resu : valeur objectif, solution : sol.       */
    /* optimale en x                                                                 */
    /* On suppose qu'avant le premier appel à reso_kp, solution a été initialisé à 0 */
    /* ----------------------------------------------------------------------------- */
    
    /* ---------------------------- */
    /* Methode de G. PLATEAU - 1987 */
    /* ---------------------------- */
    
    int     i, j;
    double  CLEF, cap_courante = cap;
    
    
    if (LIM_INF == LIM_SUP)
    {
        if (rentre(LIM_INF, LIM_SUP, cap, w, tri))
        mettre(LIM_INF, LIM_SUP, &cap_courante, p_resu, p, w, solution, tri);
        else {
            couper(LIM_INF, cap, p_resu, p, w, solution, tri);
            *indice_critique = LIM_INF;
        }
    }
    
    if (LIM_INF < LIM_SUP) {
        i = LIM_INF;
        j = LIM_SUP;
        CLEF = tri[LIM_INF].cle;
        while (i < j) {
            i++;
            while ((tri[i].cle > CLEF) && (i <= LIM_SUP)) i++;
            while (tri[j].cle < CLEF) j--;
            if (i < j) echange(i, j, tri);
        }
        
        echange(LIM_INF, j, tri);
        
        /* A ce niveau, tout ce qui est a gauche de j est de cle superieure */
        /* et tout ce qui est a droite est de cle inferieure                */
        
        cap_courante = cap;
        
        if (rentre(LIM_INF, j-1, cap_courante, w, tri)) {
            mettre(LIM_INF, j-1, &cap_courante, p_resu, p, w, solution, tri);
            if (rentre(j, j, cap_courante, w, tri)) {
                mettre(j,j, &cap_courante, p_resu, p, w, solution, tri);
                reso_kp(nbvar, p, w, cap_courante, j+1, LIM_SUP, p_resu, solution,
                indice_critique, tri);
            }
            else {
                couper(j, cap_courante, p_resu, p, w, solution, tri);
                *indice_critique = j;
            }
        }
        else {
            met_a_zero(j, LIM_SUP, solution, tri);
            reso_kp(nbvar, p, w, cap_courante, LIM_INF, j-1, p_resu, solution,
            indice_critique, tri);
        }
    }
} /* reso_kp */

/******************************************************************/
/* On résout maintenant le sac à dos linéaire en 0-1              */
/* La fction objectif est linéaire par morceaux en y_ik           */
/******************************************************************/

void pqe_main_reso_kp_y(int n, int_1D u, double_1D A, double b, double_2D s,
double_2D y, int *p_ietoile, int *p_pietoile, int_1D p)
{
    
    /* --------------------------------------------------------------------------------- */
    /* n        : nb de variables (n)                                                    */
    /* u        : borne des var x_i, ce sont elles qui permettent de calculer            */
    /*                               le nbre de var y_ik                                 */
    /* A        : coefficient dans la contrainte agrégée ie W inclu                      */
    /*                                                    (autant de coeff que de y_ik)  */
    /* b        : membre de gauche de la contrainte agrégée ie W inclu                   */
    /* s        : s_ik coefficients de la fction objectif LINEARISEE par morceaux        */
    /* y        : y_ik variables de la fction objectif LINEARISEE par morceaux           */
    /* p_ietoile: devient i*, indice critique, après lui la var est fractionnaire        */
    /* ppietoile: devient k*, indice k, correspond à i* mais en 2ème indice              */
    /* p        : tableau des indices k correpondant aux p_i et à p_i*, besoin pr derivee*/
    /*                                                                                   */
    /*                                                                                   */
    /* --------------------------------------------------------------------------------- */
    
    double_1D  solution, weight, obj;
    /* obj      : coefficients de la fction objectif QUADRATIQUE                         */
    double     valopt,max_ratio;
    int        somme_ui, i, j, k, cpt, ind_critique, atrier, indice;
    st_tri_kp *tri;
    bool       stop;
    
    
    print("Entree dans reso_kp_y \n");
    
    /* Nombre de variables y_ik : somme_ui */
    for (i = 1, somme_ui = 0; i <= n; i++) somme_ui = somme_ui + u[i];
    
    /* vecteur solution du sac-à-dos linéaire */
    pqe_init_1D((void **) &solution, sizeof(double), somme_ui);
    for (i=1; i <= somme_ui; i++) solution[i] = 0.0;
    
    
    /* coefficients de la contrainte de capacité */
    pqe_init_1D((void **) &weight, sizeof(double), somme_ui);
    cpt = 1;
    for (i = 1; i <= n; i++)
    for (j = 1 ; j <= u[i]; j++) weight[cpt++] = A[i];
    
    pqe_init_1D((void **) &tri, sizeof(st_tri_kp), somme_ui);
    
    /* coefficients de la fonction objectif */
    pqe_init_1D((void **) &obj, sizeof(double), somme_ui);
    atrier = 0;
    cpt = 0;
    for (i = 1; i <= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
        cpt++; /* indice mis à plat de la variable yik */
        obj[cpt]     = s[i][k];
        if (s[i][k] < 0)  y[i][k] = 0.0;
        else {
            atrier++;
            tri[atrier].cle = (obj[cpt] / weight[cpt]);
            tri[atrier].ind = cpt;
        }
    }
    
    print(" avant appel reso_kp atrier = %d\n", atrier);
    /*
    for (i = 1; i <= somme_ui; i++) printf(" obj[%d] = %f", i, obj[i]); printf("\n");
    for (i = 1; i <= somme_ui; i++) printf(" weight[%d] = %f", i, weight[i]); printf("\n");
    */
    
    /* résolution du sac-à-dos linéaire */
    valopt = 0.0;
    ind_critique = atrier + 1;
    reso_kp (atrier, obj, weight, b, 1, atrier, &valopt, solution,
    &ind_critique, tri);
    print("ind_critique = %d\n", ind_critique);
    
    /* Recupération du vecteur solution en y */
    cpt = 1;
    for (i = 1; i <= n; i++)
    for (k = 1; k <= u[i]; k++) {
        y[i][k] = solution[cpt];
        print("i = %d    k = %d   cpt = %d  y[%d][%d] = %lf  solution[%d] = %f\n",
        i, k, cpt, i,k, y[i][k], cpt, solution[cpt]);
        cpt++;
    }
    
    /* détermination de ietoile */
    
    
    if (ind_critique <= atrier) {
        /* Cas général : une variable a été coupée   -> fractionnaire */
        ind_critique = tri[ind_critique].ind;
        /* Calcul de l'indice critique (i*, pi*) */
        i = 1;
        while ((ind_critique - u[i]) >= 1) {
            ind_critique -= u[i];
            i++;
        }
        
        *p_ietoile = i;
        *p_pietoile = ind_critique - 1;
        /* et non plus ind_critique seult car c'est le "k"*/
        /* qui suit les variables entières                */
        
        /* analyse de la sol en y pour la détermination des pi */
        for (i = 1; i <= n; i++)
        if (i != *p_ietoile) {
            p[i] = 0;
            k = 1;
            stop = FAUX;
            while ((k <= u[i]) && !stop )
            if (y[i][k] == 0) {
                print("stop = vrai : %d %d \n", i,k);
                stop = VRAI;
            }
            else p[i] = k++;
        }
        
        p[*p_ietoile] = *p_pietoile;
        
    }
    else {
        /* Aucune varaiable n'est coupée -> indice critique :
        indice de la variable ayant le plus grand ratio obj/poids
        avec obj < 0 */
        
        *p_ietoile = 1;
        max_ratio = -infini;
        for (i = 1; i <= n; i++)
        for (k = 1; k <= u[i]; k++) {
            if ((s[i][k] < 0) && (s[i][k] / A[i] > max_ratio) ) {
                *p_ietoile = i;
                max_ratio = s[i][k] / A[i];
            }
        }
        
        /* Détermination des pi : analyse de la sol en y */
        for (i = 1; i <= n; i++) {
            p[i] = 0;
            k = 1;
            stop = FAUX;
            while ((k <= u[i]) && !stop )
            if (y[i][k] == 0) {
                print("stop = vrai : %d %d \n", i,k);
                stop = VRAI;
            }
            else p[i] = k++;
        }
        
        *p_pietoile = p[*p_ietoile];
    }
    
    print("ietoile = %d; pietoile = %d\n", *p_ietoile, *p_pietoile);
    
    
    /*    if (affichage) {      * /
    printf("cpt = %d somme_ui = %d, atrier = %d\n", cpt, somme_ui, atrier);
    for (i = 1; i <= atrier; i++) printf("obj[%d] = %f   ", i, obj[i]); print("\n");
    for (i = 1; i <= atrier; i++) printf("poids[%d] = %f   ", i, poids[i]); print("\n");
    printf("\n ind_critique = %d\n", ind_critique);
    / *   }        * /
    
    / * Calcul de l'indice critique (i*, pi*) * /
    i = 1;
    while ((ind_critique - u[i]) >= 1) {
    ind_critique -= u[i];
    i++;
}
    
    *p_ietoile = i;
    *p_pietoile = ind_critique - 1;
    / * et non plus ind_critique seult car c'est le "k"* /
    / * qui suit les variables entières                * /
    print("ietoile = %d; pietoile = %d\n", *p_ietoile, *p_pietoile);
    
    / * Recupération du vecteur solution en y * /
    cpt = 1;
    for (i = 1; i <= n; i++)
    for (k = 1; k <= u[i]; k++) {
    
    / * ATTENTION VERIFIER * /
    if (s[i][k] > 0) y[i][k] = solution[cpt];
    else  y[i][k] = 0.0;
    printf("i = %d    k = %d   cpt = %d  y[%d][%d] = %lf  solution[%d] = %f\n",
    i, k, cpt, i,k, y[i][k], cpt, solution[cpt]);
    cpt++;
}
    
    / * analyse de la sol en y pour la détermination des pi * /
    for (i = 1; i <= n; i++)
    if (i != *p_ietoile) {
    p[i] = 0;
    k = 1;
    stop = FAUX;
    while ((k <= u[i]) && !stop )
    if (y[i][k] == 0) {
    print("stop = vrai : %d %d \n", i,k);
    stop = VRAI;
}
    else p[i] = k++;
}
    
    p[*p_ietoile] = *p_pietoile;
    */
    
    
    pqe_init_free_1D(solution);
    pqe_init_free_1D(weight);
    pqe_init_free_1D(tri);
    pqe_init_free_1D(obj);
    
    
    print(" valopt = %lf  ietoile = %d ", valopt, *p_ietoile);
}   /* pqe_main_reso_kp_y */

/****************************************************************************************/

double calcul_derivee(int n, double_2D s, int ietoile, int pietoile, double_1D b,
double_2D a, int_1D p, double_1D A, double B, int j)
{
    double resu, somme1, somme2, num, denom;
    int i;
    
    for (i = 1, somme1 = 0.0; i <= n; i++) somme1 += p[i]* a[j][i];
    for (i = 1, somme2 = 0.0; i <= n; i++) somme2 += p[i]* A[i];
    
    num   = (b[j] - somme1) * A[ietoile] - (B - somme2) * a[j][ietoile];
    denom = (A[ietoile] * A[ietoile])  ;
    resu  = s[ietoile][pietoile+1] * ( num / denom );
    
    return (resu);
} /* calcul_derivee */

/*****************************************************************************************/

double valeur_obj(int n, double_2D s, int ietoile, int pietoile, double B,
double_1D A, int_1D p, double_2D y)
{
    /*     double somme2;    */
    double resu = 0.0;
    int i, k;
    
    /* for (i = 1, somme2 = 0.0; i <= n; i++) somme2 += p[i]* A[i]; */
    
    for (i = 1; i <= n; i++)
    for (k = 1; k <= p[i]; k++) resu += s[i][k];
    
    print("resu = %f ietoile = %d pietoile = %d\n", resu, ietoile, pietoile);
    
    resu += s[ietoile][pietoile+1] * y[ietoile][pietoile+1];
    /* ... * ( (B-somme2) / A[ietoile] ); */
    print("resu = %f\n", resu);
    
    return(resu);
} /* valeur_obj */

/******************************************************************/

double norme(double *d, int m)
{
    int    j;
    double somme =0.0;
    
    for (j = 1; j <= m; j++) somme += d[j] * d[j];
    return (sqrt(somme));
} /* norme */


/******************************************************************/

double min2(double x, double y)
{
    if (x < y) return(x);
    else       return(y);
} /* min2 */

/******************************************************************/

double min3(double x, double y, double z) {
    if (x < y) {
        return(min2(x,z));
    }
    else return(min2(y,z));
} /* min3 */

/******************************************************************/

double agreg_col(int n, int m, double_2D a, int i, double_1D w)
{
    /* fonction g de l'article de Djerdjour : retourne g_i(w) */
    int    j;
    double somme = 0.0;
    
    for (j = 1; j <= m; j++) somme += w[j] * a[j][i]; 
    
    return(somme);
} /* agreg_col */

/******************************************************************/

double agreg_b(int m, double_1D b, double_1D w)
{
    int    j;
    double somme = 0.0;
    
    for (j = 1; j <= m; j++) somme += w[j] * b[j]; 
    
    return(somme);
} /* agreg_b */

/******************************************************************/

void pqe_main_calcule_w_etoile(int n, int m, double_1D w0, double_2D s, double_2D a,
                               double_1D b, int_1D u, double_1D w_etoile, double *p_borne1,
                               int *p_nbiter, double *p_temps, double_2D y_cont)
{
    const double eps2 = 0.000001;
    const double eps3 = 0.000001;

    double_1D w, A, d;
    double_2D y, theta;
    int_1D    p;
    double    derivee_j, B, alpha, alpha1, alpha2, alpha21, alpha22; 
    double    alpha3, alpha31, alpha32,borne_sup, ancienne_borne = -infini;
    double    somme_d, somme_w, b_d, b_w, ratio, denominateur, x, z;
    int       k, i, j, ietoile, pietoile, iteration, palier_alpha;
    bool      arret = FAUX, tousnuls, alpha_D = FAUX;
    char      c;
    

    print("n = %d m = %d", n, m);
    for (j = 1; j <= m; j++) print(" w0[%d] = %lf",j, w0[j]);
    for (i = 1; i <= n; i++)
	for (k = 1; k <= u[i]; k++) print(" s[%d][%d] = %lf ", i, k, s[i][k]);







    pqe_init_1D((void **)  &w,     sizeof(double), m);
    pqe_init_1D((void **)  &d,     sizeof(double), m);
    pqe_init_1D((void **)  &A,     sizeof(double), n);
    pqe_init_1D((void **)  &p,     sizeof(int),    n);
    pqe_init_2D((void ***) &theta, sizeof(double), 4, n);
    y = (double **) malloc(sizeof(double *) * (1+n));
    for (i = 0; i < 1 + n; i++) y[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
    
    pqe_common_init_temps(p_temps);
    
    print("\n n = %d   m = %d\n", n, m);
    /* Initialisation */
    
    *p_borne1 = infini;
    for (j = 1; j <= m; j++) w[j] = w0[j];
    /*alpha = 2.0;*/
    /* alpha = norme (w0, n); */
    alpha = norme (w0, m)/4.0;
    iteration = palier_alpha = 0;
    *p_nbiter = 0;
    
    while (! arret) {
        
        (*p_nbiter)++;
        
        /* Aggregation des contraintes */
        for (i = 1; i <= n; i++) {
            A[i] = 0.0;
            for (j = 1; j <= m; j++) A[i] += w[j] * a[j][i];
            print("A[%d] = %lf ",i, A[i]);
        }
        
        B    = 0.0;
        for (j = 1; j <= m; j++) B += w[j] * b[j];
        print(" B = %lf\n", B);
        
        /* Resolution du pb de sac-à-dos linéaire continu */
        pqe_main_reso_kp_y(n, u, A, B, s, y, &ietoile, &pietoile, p);
        
        borne_sup = valeur_obj(n, s, ietoile, pietoile, B, A, p, y);
        *p_borne1 = borne_sup < (*p_borne1) ? borne_sup : (*p_borne1);
        
        print("\n Valeur courante de la borne supérieure : %f alpha = %f", borne_sup, alpha);
        /*
        for (i = 1; i <= n; i++) print("p[%d] = %d ",i, p[i]); print("\n");
        */
        
        /* Calcul de d */
        
        for (j = 1; j <= m; j++) {
            /* calcul de la dérivée partielle dZ/dWj */
            derivee_j = calcul_derivee(n, s, ietoile, pietoile, b, a, p, A, B, j);
            
	    /* Perturbation de la dérivee pour éviter les erreurs d'arrondis */
	    x = (double) rand() / (double) RAND_MAX;
	    z = (9.0*x+1)*eps3;
	    if ((derivee_j >= 0.0) && (derivee_j < eps3)) derivee_j = z;
	    if ((derivee_j > -eps3) && (derivee_j < 0.0)) derivee_j = -z;
            if (w[j] > precision_double) d[j] = -derivee_j;
            else                         d[j] = (derivee_j < 0.0 ? -derivee_j : 0.0);
            /*---------------------------------------------------------------------------*/
            /*																			  */
            /*  Changement ici pour afficher les w_j pour observer ce qu'il se           */
            /*  passe quand m=5 et n varie                                               */
            /*---------------------------------------------------------------------------
            print("w[%d] = %lf  derivee_%d = %lf d[%d] = %lf  \n",
            j,w[j],j, derivee_j,j, d[j]);
            */
	    if ((borne_sup > 82.8) && (borne_sup < 83.2))
                print("w[%d] = %lf  derivee_%d = %lf d[%d] = %.308lf  \n", j,w[j],j, derivee_j,j, d[j]);
        }
        
      /* scanf("%c", &c);  */
        print("fin calcul derivees\n");
        
        tousnuls = VRAI;
        for (j = 1; (j <= m) && (tousnuls); j++)

            if (d[j] != 0.0) { tousnuls = FAUX; print(" c'est d[%d] qui est non nul (%lf) ", j, d[j]); }
        
        if (tousnuls) arret = VRAI; 
	print(" tousnuls = %d", tousnuls);
        
        /* Calcul de alpha */
        
        if (alpha_D) {
            /* Calcul de aplha1 */
            alpha1 = INFINI;
            for (j = 1; j <= m; j++) {
                if ((d[j] < 0.0) && ((-(w[j]/d[j])) < alpha1)) alpha1 = (-(w[j]/d[j])); 
            }
            
            /* Calcul de aplha2 */
            alpha2 = INFINI;
            
            somme_d = 0.0;
            for (i = 1; i <= n; i++) somme_d += p[i] * agreg_col(n, m, a, i, d);
            
            somme_w = 0.0; 
            for (i = 1; i <= n; i++) somme_d += p[i] * agreg_col(n, m, a, i, w);
            
            b_d = agreg_b(m, b, d);
            b_w = agreg_b(m, b, w);
            
            /* Calcul de alpha21 */
            denominateur = somme_d - b_d;
            if (denominateur > 0.0 ) alpha21 = (b_w - somme_w) / denominateur;
            else                     alpha21 = INFINI;
            
            /* Calcul de alpha22 */
            denominateur = b_d - somme_d - agreg_col(n, m, a, ietoile, d);
            if (denominateur > 0.0) {
                alpha22 = (somme_w + agreg_col(n, m, a, ietoile, w) - b_w) / denominateur;
            }
            else alpha22 = INFINI;
            
            alpha2 = min2(alpha21, alpha22);
            
            /* Calcul de aplha3 */
            
            /* Calcul des theta_k_i */
            for (i = 1; i <= n; i++)
            if (i != ietoile) {
                theta[1][i] = s[i][p[i]] * agreg_col(n, m, a, ietoile, w) 
                              - s[ietoile][pietoile + 1] * agreg_col(n, m, a, i, w);
                theta[2][i] = s[ietoile][pietoile + 1] * agreg_col(n, m, a, i, d) 
                              - s[i][p[i]] * agreg_col(n, m, a, ietoile, d);
                theta[3][i] = s[ietoile][pietoile + 1] * agreg_col(n, m, a, i, w) 
                              - s[i][p[i] + 1] * agreg_col(n, m, a, ietoile, w);
                theta[4][i] = s[i][p[i] + 1] * agreg_col(n, m, a, ietoile, d) 
                              - s[ietoile][pietoile + 1] * agreg_col(n, m, a, i, d);
            }
            
            /* Calcul de alpha31 */
            alpha31 = INFINI;
            for (i = 1; i <= n; i++)
            if (i != ietoile) {
                if ((theta[2][i] > 0.0) && (theta[1][i] > 0.0)) { 
                    ratio = theta[1][i] / theta[2][i];
                    if (ratio < alpha31) { 
                        alpha31 = ratio;
                        print("\n i = %d ratio = %.2f", i, ratio);
                    }
                }
            }
            
            /* Calcul de alpha32 */
            alpha32 = INFINI;
            for (i = 1; i <= n; i++)
            if (i != ietoile) {
                if ((theta[4][i] > 0.0) && (theta[3][i] > 0.0)) { 
                    ratio = theta[3][i] / theta[4][i];
                    if (ratio < alpha32) alpha32 = ratio;
                }
            }
            
            alpha3 = min2(alpha31, alpha32);
            
            alpha = min3(alpha1, alpha2, alpha3);
            
            /* print("\n alpha = %f alpha1 = %.2f alpha21 = %.2f alpha22 = %.2f alpha31 = %.2f alpha32 = %.2f alpha2 = %.2f alpha3 = %.2f", 
            alpha, alpha1, alpha21, alpha22, alpha31, alpha32, alpha2, alpha3); */
            print("\n alpha = %f  borne_sup = %.2f", alpha, borne_sup);
        }
        else {
            
            if (ancienne_borne <= borne_sup + eps2) {
                /* Pas d'amélioration de la borne */
                if (iteration >= diminue_pas) {
                    alpha /= 2.0;
                    iteration = 0 ;
                }
                else {
                    iteration++;
                }
            }
            else {
                /* Amélioration de la borne */
                iteration = 0;
            }
            
            if (palier_alpha == 100) {
                alpha /= 2.0;
                iteration = 0 ;
                palier_alpha = 0;
            }
            
            palier_alpha++;
            print("\n alpha = %f  borne_sup = %.2f", alpha, borne_sup);
        }
        
        /* Calcul du nouveau w */
        
        for (j = 1; j <= m; j++) {
            if (d[j] != 0.0) {w[j] = w[j] + alpha * d[j] / norme(d, m);
            print("w[%d] modifie  = %lf  d[%d] = %.100lf norme(d) = %lf \n", j,w[j],j, d[j], norme(d, m));
	    }
            if (w[j] < 0.0) w[j] = 0.0;
        }
        
        if (alpha < epsilon) arret = VRAI;
        ancienne_borne = borne_sup;
        /*scanf("%c", &c);*/
        
        print(" fin boucle while !arret \n");
	print(" à la fin,  alpha=%f arret = %d", alpha, arret);
    }
    
    for (j = 1; j <= m; j++) { w_etoile[j] = w[j]; print("wetoile[%d] = %f\n", j, w_etoile[j]); }
    /* affichage de la solution en y pour Borne 1, pour voir si admissible pr MKP*/
    for (i=1;i<=n;i++){
        for(k=1;k<=u[i];k++){
            y_cont[i][k] = y[i][k];
            /*    print("y[%d][%d] = %lf\n", i, k, y[i][k]); */
        }
    }
    
    pqe_common_calcule_temps(p_temps);
    
    pqe_init_free_1D(w);
    pqe_init_free_1D(d);
    pqe_init_free_1D(A);
    pqe_init_free_1D(p);
    pqe_init_free_2D((void *) y, n);
    
    
    print("\n borne1 = %f (%.2f s)\n", borne_sup, *p_temps);
    fflush(stdin);
    /* scanf("%c", &c);
    scanf("%c", &c); */
} /* pqe_main_calcule_w_etoile */

/******************************************************************/

void pqe_main_linea_fction_quadra (int n, int_1D u, double_1D c, double_1D d, double_2D vec_s)
{
    static  int appel = 0;
    int       i, k;
    double_2D vec_f ;  
    /* double   vec_f[1000][1000]; */
    /* double f[1000][1000]; */
    double truc;
    /* double_1D g; */
    
    print("\n n = %d appel = %d\n", n, appel);
    pqe_init_2D((void ***) &vec_f, sizeof(double), n, pqe_main_max_ui(n,u));
    /*vec_f = (double **) malloc(sizeof(double *) * (1+n)); 
    for (i = 0; i < 1 + n; i++) {
    vec_f[i] = (double *) malloc(sizeof(double) * (1 + u[i])); 
    print("c[%d] = %lf   d[%d] = %lf u[%d] = %d\n", i, c[i], i, d[i], i, u[i]);
}
    */
    /*
    print(" \n avant essai ecriture dans f[10][15]\n");
    f[10][15] = 915.00;
    print("\nf[10][15] = %f\n", f[10][15]);
    */
    /*
    pqe_init_1D((void **) &g, sizeof(double), n);
    y = (double **) malloc(sizeof(double *) * (1+n));
    for (i = 0; i < 1 + n; i++) y[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
    */
    
    for (i = 1; i <= n; i++) {
        print("c[%d] = %lf   d[%d] = %lf u[%d] = %d\n", i, c[i], i, d[i], i, u[i]);
        vec_f[i][0] = 0.0;
        for (k = 1; k <= u[i]; k++) {
            /* f[i][k] = 0.0; */
            /* g[i] = 0; */
            print(" i= %d   k = %d\n", i, k);
            truc = (c[i]*(double) k - d[i]*(double) k*(double) k);
            print(" truc = %f \n", truc);
            print("\nf[10][15] = %f\n", vec_f[10][15]);
            /* print("\n truc = %f f[%d][%d] = %f  f[%d][%d] = %f \n", truc,i, k , f[i][k], i, k-1, f[i][k-1]);*/
            vec_f[i][k] = truc;
            print("\n apres affectation \n");
            /* f[i][k] = (c[i]*k - d[i]*k*k); */
            print("f[%d][%d] = %lf\n ", i, k, vec_f[i][k]);
            fflush(stdout);
            vec_s[i][k] = vec_f[i][k] - vec_f[i][k-1];
            /* g[i] = (g[i] + s[i][k]*y[i][k]); */
            print("s[%d][%d] = %lf \n", i,k, vec_s[i][k]);
        }
        
    }
    
    
    print("\n avant liberation de f\n");
    pqe_init_free_2D((void *) vec_f, n); 
    print("\n apres liberation de f\n");
    
    
    appel++;
} /*  pqe_main_linea_fction_quadra */



/******************************************************************/

int somme_ui(int n, int_1D u)
{
    int i, resu;
    
    /* Nombre de variables y_ik : somme_ui */
    for (i=1, resu= 0; i <=n; i++) resu = resu + u[i];
    
    return(resu);
} /* somme_ui */

/******************************************************************/

int pqe_main_max_ui(int n, int_1D u)
{
    
    int i, resu = 0;
    
    for (i=1; i <=n; i++) if (u[i] > resu) resu = u[i];
    
    return(resu);
} /* max_ui */

/******************************************************************/

void genere_exemple_1(int *p_n, int *p_m, double_1D *p_c, double_1D *p_d,
double_2D *p_a, double_1D *p_b, int_1D *p_u,
double_1D *p_w0)
{
    int i;
    
    *p_n = 2;
    *p_m = 2;
    
    /* initialisation des vecteurs */
    
    
    pqe_init_1D((void **)  p_w0, sizeof(double), *p_m);
    pqe_init_2D((void ***) p_a,  sizeof(double), *p_m, *p_n);
    pqe_init_1D((void **)  p_b,  sizeof(double), *p_m);
    pqe_init_1D((void **)  p_c,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_d,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_u,  sizeof(int),    *p_n);
    
    
    (*p_c)[1] = 5.0;
    (*p_c)[2] = 12.0;
    
    (*p_d)[1]= 2.0;
    (*p_d)[2] = 1.0;
    
    (*p_a)[1][1] = 2.0;
    
    (*p_a)[1][2] = 1.0;
    (*p_a)[2][1] = 1.0;
    (*p_a)[2][2] = 5.0;
    
    (*p_b)[1] = 10.0;
    (*p_b)[2] = 8.0;
    
    (*p_w0)[1] = 1000;
    (*p_w0)[2] = 1000;
    /*w0[1] = 4.0/9.0; w0[2] = 19.0/9.0;*/
    /*w0[1] = 1; w0[2] = 2;*/
    
    (*p_u)[1] = (*p_u)[2] = 5;
    
    
} /* genere_exemple_1 */


/******************************************************************/

void genere_exemple_2(int *p_n, int *p_m, double_1D *p_c, double_1D *p_d,
double_2D *p_a, double_1D *p_b, int_1D *p_u,
double_1D *p_w0)
{
    int i;
    
    *p_n = 10;
    *p_m = 10;
    
    /* initialisation des vecteurs */
    
    
    pqe_init_1D((void **)  p_w0, sizeof(double), *p_m);
    pqe_init_2D((void ***) p_a,  sizeof(double), *p_m, *p_n);
    pqe_init_1D((void **)  p_b,  sizeof(double), *p_m);
    pqe_init_1D((void **)  p_c,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_d,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_u,  sizeof(int),    *p_n);
    
    /* Données d'OR_Library pour tester ave(*p_c) pqe_main*/
    
    
    /* les (*p_c)_i, coeff des x_i*/
    (*p_c)[1] = 63.0;
    (*p_c)[2] = 15.0;
    (*p_c)[3] = 44.0;
    (*p_c)[4] = 91.0;
    (*p_c)[5] = 45.0;
    (*p_c)[6] = 50.0;
    (*p_c)[7] = 89.0;
    (*p_c)[8] = 58.0;
    (*p_c)[9] = 86.0;
    (*p_c)[10] = 82.0;
    
    /* les d_i, coeff des x_i²*/
    (*p_d)[1] = 19.0;
    (*p_d)[2] = 27.0;
    (*p_d)[3] = 23.0;
    (*p_d)[4] = 53.0;
    (*p_d)[5] = 42.0;
    (*p_d)[6] = 26.0;
    (*p_d)[7] = 33.0;
    (*p_d)[8] = 23.0;
    (*p_d)[9] = 41.0;
    (*p_d)[10] = 19.0;
    
    /* la matrice des contraintes, les coeff a[j][i]*/
    
    
    /* 1ère ligne*/
    (*p_a)[1][1] = 3.0;
    (*p_a)[1][2] = 5.0;
    (*p_a)[1][3] = 5.0;
    (*p_a)[1][4] = 6.0;
    (*p_a)[1][5] = 4.0;
    (*p_a)[1][6] = 4.0;
    (*p_a)[1][7] = 5.0;
    (*p_a)[1][8] = 6.0;
    (*p_a)[1][9] = 4.0;
    (*p_a)[1][10] = 4.0;
    /* 2ème ligne*/
    (*p_a)[2][1] = 5.0;
    (*p_a)[2][2] = 4.0;
    (*p_a)[2][3] = 5.0;
    (*p_a)[2][4] = 4.0;
    (*p_a)[2][5] = 1.0;
    (*p_a)[2][6] = 4.0;
    (*p_a)[2][7] = 4.0;
    (*p_a)[2][8] = 2.0;
    (*p_a)[2][9] = 5.0;
    (*p_a)[2][10] = 2.0;
    /* 3ème ligne*/
    (*p_a)[3][1] = 1.0;
    (*p_a)[3][2] = 5.0;
    (*p_a)[3][3] = 2.0;
    (*p_a)[3][4] = 4.0;
    (*p_a)[3][5] = 7.0;
    (*p_a)[3][6] = 3.0;
    (*p_a)[3][7] = 1.0;
    (*p_a)[3][8] = 5.0;
    (*p_a)[3][9] = 7.0;
    (*p_a)[3][10] = 6.0;
    /* 4ème ligne*/
    (*p_a)[4][1] = 3.0;
    (*p_a)[4][2] = 2.0;
    (*p_a)[4][3] = 6.0;
    (*p_a)[4][4] = 3.0;
    (*p_a)[4][5] = 2.0;
    (*p_a)[4][6] = 1.0;
    (*p_a)[4][7] = 6.0;
    (*p_a)[4][8] = 1.0;
    (*p_a)[4][9] = 7.0;
    (*p_a)[4][10] = 3.0;
    /* 5ème ligne*/
    (*p_a)[5][1] = 6.0;
    (*p_a)[5][2] = 6.0;
    (*p_a)[5][3] = 6.0;
    (*p_a)[5][4] = 4.0;
    (*p_a)[5][5] = 5.0;
    (*p_a)[5][6] = 2.0;
    (*p_a)[5][7] = 2.0;
    (*p_a)[5][8] = 4.0;
    (*p_a)[5][9] = 3.0;
    (*p_a)[5][10] = 2.0;
    /* 6ème ligne*/
    (*p_a)[6][1] = 5.0;
    (*p_a)[6][2] = 5.0;
    (*p_a)[6][3] = 2.0;
    (*p_a)[6][4] = 1.0;
    (*p_a)[6][5] = 3.0;
    (*p_a)[6][6] = 5.0;
    (*p_a)[6][7] = 5.0;
    (*p_a)[6][8] = 7.0;
    (*p_a)[6][9] = 4.0;
    (*p_a)[6][10] = 3.0;
    /* 7ème ligne*/
    (*p_a)[7][1] = 3.0;
    (*p_a)[7][2] = 6.0;
    (*p_a)[7][3] = 6.0;
    (*p_a)[7][4] = 3.0;
    (*p_a)[7][5] = 1.0;
    (*p_a)[7][6] = 6.0;
    (*p_a)[7][7] = 1.0;
    (*p_a)[7][8] = 6.0;
    (*p_a)[7][9] = 7.0;
    (*p_a)[7][10] =1.0;
    /* 8ème ligne*/
    (*p_a)[8][1] = 1.0;
    (*p_a)[8][2] = 2.0;
    (*p_a)[8][3] = 1.0;
    (*p_a)[8][4] = 7.0;
    (*p_a)[8][5] = 8.0;
    (*p_a)[8][6] = 7.0;
    (*p_a)[8][7] = 6.0;
    (*p_a)[8][8] = 5.0;
    (*p_a)[8][9] = 8.0;
    (*p_a)[8][10] = 7.0;
    /* 9ème ligne*/
    (*p_a)[9][1] = 8.0;
    (*p_a)[9][2] = 5.0;
    (*p_a)[9][3] = 2.0;
    (*p_a)[9][4] = 5.0;
    (*p_a)[9][5] = 3.0;
    (*p_a)[9][6] = 8.0;
    (*p_a)[9][7] = 1.0;
    (*p_a)[9][8] = 3.0;
    (*p_a)[9][9] = 3.0;
    (*p_a)[9][10] =5.0;
    /* 10ème ligne*/
    (*p_a)[10][1] = 1.0;
    (*p_a)[10][2] = 1.0;
    (*p_a)[10][3] = 1.0;
    (*p_a)[10][4] = 1.0;
    (*p_a)[10][5] = 1.0;
    (*p_a)[10][6] = 1.0;
    (*p_a)[10][7] = 1.0;
    (*p_a)[10][8] = 1.0;
    (*p_a)[10][9] = 1.0;
    (*p_a)[10][10] = 1.0;
    
    
    /* membre de droite des contraintes*/
    (*p_b)[1] = 380.0;
    (*p_b)[2] = 415.0;
    (*p_b)[3] = 385.0;
    (*p_b)[4] = 405.0;
    (*p_b)[5] = 470.0;
    (*p_b)[6] = 415.0;
    (*p_b)[7] = 400.0;
    (*p_b)[8] = 460.0;
    (*p_b)[9] = 400.0;
    (*p_b)[10] = 200.0;
    
    
    /* borne sup des variables x_i*/
    (*p_u)[1] = (*p_u)[2] = 1; (*p_u)[3] =1; (*p_u)[4] = 1; (*p_u)[5] = 1; (*p_u)[6] = 1; (*p_u)[7]=2; (*p_u)[8]=2;
    (*p_u)[9] = 2; (*p_u)[10] =4;
    
    
    /* un w0 de départ*/
    (*p_w0)[1] = 1; (*p_w0)[2] = 2; (*p_w0)[3] = 1; (*p_w0)[4] = 2; (*p_w0)[5] = 1; (*p_w0)[6] = 2; (*p_w0)[7] = 1; (*p_w0)[8]=2;
    (*p_w0)[9] = 1; (*p_w0)[10] = 2;
    
    
} /* genere_exemple_2 */

/******************************************************************/

void genere_exemple_3(int *p_n, int *p_m, double_1D *p_c, double_1D *p_d,
double_2D *p_a, double_1D *p_b, int_1D *p_u,
double_1D *p_w0)
{
    int i;
    
    *p_n = 10;
    *p_m = 10;
    
    /* initialisation des vecteurs */
    
    
    pqe_init_1D((void **)  p_w0, sizeof(double), *p_m);
    pqe_init_2D((void ***) p_a,  sizeof(double), *p_m, *p_n);
    pqe_init_1D((void **)  p_b,  sizeof(double), *p_m);
    pqe_init_1D((void **)  p_c,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_d,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_u,  sizeof(int),    *p_n);
    
    /* Données d'OR_Library pour tester ave(*p_c) pqe_main*/
    
    
    /* les (*p_c)_i, coeff des x_i*/
    (*p_c)[1] = 63.0;
    (*p_c)[2] = 15.0;
    (*p_c)[3] = 44.0;
    (*p_c)[4] = 91.0;
    (*p_c)[5] = 45.0;
    (*p_c)[6] = 50.0;
    (*p_c)[7] = 89.0;
    (*p_c)[8] = 58.0;
    (*p_c)[9] = 86.0;
    (*p_c)[10] = 82.0;
    
    /* les d_i, coeff des x_i²*/
    (*p_d)[1] = 19.0;
    (*p_d)[2] = 10.0;
    (*p_d)[3] = 23.0;
    (*p_d)[4] = 53.0;
    (*p_d)[5] = 42.0;
    (*p_d)[6] = 26.0;
    (*p_d)[7] = 33.0;
    (*p_d)[8] = 23.0;
    (*p_d)[9] = 41.0;
    (*p_d)[10] = 19.0;
    
    /* la matrice des contraintes, les coeff a[j][i]*/
    
    
    /* 1ère ligne*/
    (*p_a)[1][1] = 3.0;
    (*p_a)[1][2] = 5.0;
    (*p_a)[1][3] = 5.0;
    (*p_a)[1][4] = 6.0;
    (*p_a)[1][5] = 4.0;
    (*p_a)[1][6] = 4.0;
    (*p_a)[1][7] = 5.0;
    (*p_a)[1][8] = 6.0;
    (*p_a)[1][9] = 4.0;
    (*p_a)[1][10] = 4.0;
    /* 2ème ligne*/
    (*p_a)[2][1] = 5.0;
    (*p_a)[2][2] = 4.0;
    (*p_a)[2][3] = 5.0;
    (*p_a)[2][4] = 4.0;
    (*p_a)[2][5] = 1.0;
    (*p_a)[2][6] = 4.0;
    (*p_a)[2][7] = 4.0;
    (*p_a)[2][8] = 2.0;
    (*p_a)[2][9] = 5.0;
    (*p_a)[2][10] = 2.0;
    /* 3ème ligne*/
    (*p_a)[3][1] = 1.0;
    (*p_a)[3][2] = 5.0;
    (*p_a)[3][3] = 2.0;
    (*p_a)[3][4] = 4.0;
    (*p_a)[3][5] = 7.0;
    (*p_a)[3][6] = 3.0;
    (*p_a)[3][7] = 1.0;
    (*p_a)[3][8] = 5.0;
    (*p_a)[3][9] = 7.0;
    (*p_a)[3][10] = 6.0;
    /* 4ème ligne*/
    (*p_a)[4][1] = 3.0;
    (*p_a)[4][2] = 2.0;
    (*p_a)[4][3] = 6.0;
    (*p_a)[4][4] = 3.0;
    (*p_a)[4][5] = 2.0;
    (*p_a)[4][6] = 1.0;
    (*p_a)[4][7] = 6.0;
    (*p_a)[4][8] = 1.0;
    (*p_a)[4][9] = 7.0;
    (*p_a)[4][10] = 3.0;
    /* 5ème ligne*/
    (*p_a)[5][1] = 6.0;
    (*p_a)[5][2] = 6.0;
    (*p_a)[5][3] = 6.0;
    (*p_a)[5][4] = 4.0;
    (*p_a)[5][5] = 5.0;
    (*p_a)[5][6] = 2.0;
    (*p_a)[5][7] = 2.0;
    (*p_a)[5][8] = 4.0;
    (*p_a)[5][9] = 3.0;
    (*p_a)[5][10] = 2.0;
    /* 6ème ligne*/
    (*p_a)[6][1] = 5.0;
    (*p_a)[6][2] = 5.0;
    (*p_a)[6][3] = 2.0;
    (*p_a)[6][4] = 1.0;
    (*p_a)[6][5] = 3.0;
    (*p_a)[6][6] = 5.0;
    (*p_a)[6][7] = 5.0;
    (*p_a)[6][8] = 7.0;
    (*p_a)[6][9] = 4.0;
    (*p_a)[6][10] = 3.0;
    /* 7ème ligne*/
    (*p_a)[7][1] = 3.0;
    (*p_a)[7][2] = 6.0;
    (*p_a)[7][3] = 6.0;
    (*p_a)[7][4] = 3.0;
    (*p_a)[7][5] = 1.0;
    (*p_a)[7][6] = 6.0;
    (*p_a)[7][7] = 1.0;
    (*p_a)[7][8] = 6.0;
    (*p_a)[7][9] = 7.0;
    (*p_a)[7][10] =1.0;
    /* 8ème ligne*/
    (*p_a)[8][1] = 1.0;
    (*p_a)[8][2] = 2.0;
    (*p_a)[8][3] = 1.0;
    (*p_a)[8][4] = 7.0;
    (*p_a)[8][5] = 8.0;
    (*p_a)[8][6] = 7.0;
    (*p_a)[8][7] = 6.0;
    (*p_a)[8][8] = 5.0;
    (*p_a)[8][9] = 8.0;
    (*p_a)[8][10] = 7.0;
    /* 9ème ligne*/
    (*p_a)[9][1] = 8.0;
    (*p_a)[9][2] = 5.0;
    (*p_a)[9][3] = 2.0;
    (*p_a)[9][4] = 5.0;
    (*p_a)[9][5] = 3.0;
    (*p_a)[9][6] = 8.0;
    (*p_a)[9][7] = 1.0;
    (*p_a)[9][8] = 3.0;
    (*p_a)[9][9] = 3.0;
    (*p_a)[9][10] =5.0;
    /* 10ème ligne*/
    (*p_a)[10][1] = 1.0;
    (*p_a)[10][2] = 1.0;
    (*p_a)[10][3] = 1.0;
    (*p_a)[10][4] = 1.0;
    (*p_a)[10][5] = 1.0;
    (*p_a)[10][6] = 1.0;
    (*p_a)[10][7] = 1.0;
    (*p_a)[10][8] = 1.0;
    (*p_a)[10][9] = 1.0;
    (*p_a)[10][10] = 1.0;
    
    
    /* membre de droite des contraintes*/
    (*p_b)[1] = 38.0;
    (*p_b)[2] = 41.0;
    (*p_b)[3] = 35.0;
    (*p_b)[4] = 40.0;
    (*p_b)[5] = 47.0;
    (*p_b)[6] = 48.0;
    (*p_b)[7] = 46.0;
    (*p_b)[8] = 90.0;
    (*p_b)[9] = 20.0;
    (*p_b)[10] = 14.0;
    
    
    /* borne sup des variables x_i*/
    (*p_u)[1] = (*p_u)[2] = 1; (*p_u)[3] =1; (*p_u)[4] = 1; (*p_u)[5] = 1; (*p_u)[6] = 1; (*p_u)[7]=2; (*p_u)[8]=2;
    (*p_u)[9] = 2; (*p_u)[10] =4;
    
    /* un w0 de départ*/
    (*p_w0)[1] = 1; (*p_w0)[2] = 2; (*p_w0)[3] = 3; (*p_w0)[4] = 4; (*p_w0)[5] = 5; (*p_w0)[6] = 6; (*p_w0)[7] = 7; (*p_w0)[8]=8;
    (*p_w0)[9] = 9; (*p_w0)[10] = 10;
    
    
    
    /* un w0 de départ* /
    (*p_w0)[1] = 1; (*p_w0)[2] = 2; (*p_w0)[3] = 1; (*p_w0)[4] = 2; (*p_w0)[5] = 1; (*p_w0)[6] = 2; (*p_w0)[7] = 1; (*p_w0)[8]=2;
    (*p_w0)[9] = 1; (*p_w0)[10] = 2;     */
    
    
} /* genere_exemple_3 */

/******************************************************************/
/******************************************************************/

void genere_exemple_4(int *p_n, int *p_m, double_1D *p_c, double_1D *p_d,
double_2D *p_a, double_1D *p_b, int_1D *p_u,
double_1D *p_w0)
{
    int i;
    
    *p_n = 10;
    *p_m = 10;
    
    /* initialisation des vecteurs */
    
    
    pqe_init_1D((void **)  p_w0, sizeof(double), *p_m);
    pqe_init_2D((void ***) p_a,  sizeof(double), *p_m, *p_n);
    pqe_init_1D((void **)  p_b,  sizeof(double), *p_m);
    pqe_init_1D((void **)  p_c,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_d,  sizeof(double), *p_n);
    pqe_init_1D((void **)  p_u,  sizeof(int),    *p_n);
    
    /* Données d'OR_Library pour tester ave(*p_c) pqe_main*/
    
    
    /* les (*p_c)_i, coeff des x_i*/
    (*p_c)[1] = 63.0;
    (*p_c)[2] = 15.0;
    (*p_c)[3] = 44.0;
    (*p_c)[4] = 91.0;
    (*p_c)[5] = 45.0;
    (*p_c)[6] = 50.0;
    (*p_c)[7] = 89.0;
    (*p_c)[8] = 58.0;
    (*p_c)[9] = 86.0;
    (*p_c)[10] = 82.0;
    
    /* les d_i, coeff des x_i²*/
    (*p_d)[1] = 19.0;
    (*p_d)[2] = 10.0;
    (*p_d)[3] = 23.0;
    (*p_d)[4] = 53.0;
    (*p_d)[5] = 42.0;
    (*p_d)[6] = 26.0;
    (*p_d)[7] = 33.0;
    (*p_d)[8] = 23.0;
    (*p_d)[9] = 41.0;
    (*p_d)[10] = 19.0;
    
    /* la matrice des contraintes, les coeff a[j][i]*/
    
    
    /* 1ère ligne*/
    (*p_a)[1][1] = 3.0;
    (*p_a)[1][2] = 5.0;
    (*p_a)[1][3] = 5.0;
    (*p_a)[1][4] = 6.0;
    (*p_a)[1][5] = 4.0;
    (*p_a)[1][6] = 4.0;
    (*p_a)[1][7] = 5.0;
    (*p_a)[1][8] = 6.0;
    (*p_a)[1][9] = 4.0;
    (*p_a)[1][10] = 4.0;
    /* 2ème ligne*/
    (*p_a)[2][1] = 5.0;
    (*p_a)[2][2] = 4.0;
    (*p_a)[2][3] = 5.0;
    (*p_a)[2][4] = 4.0;
    (*p_a)[2][5] = 1.0;
    (*p_a)[2][6] = 4.0;
    (*p_a)[2][7] = 4.0;
    (*p_a)[2][8] = 2.0;
    (*p_a)[2][9] = 5.0;
    (*p_a)[2][10] = 2.0;
    /* 3ème ligne*/
    (*p_a)[3][1] = 1.0;
    (*p_a)[3][2] = 5.0;
    (*p_a)[3][3] = 2.0;
    (*p_a)[3][4] = 4.0;
    (*p_a)[3][5] = 7.0;
    (*p_a)[3][6] = 3.0;
    (*p_a)[3][7] = 1.0;
    (*p_a)[3][8] = 5.0;
    (*p_a)[3][9] = 7.0;
    (*p_a)[3][10] = 6.0;
    /* 4ème ligne*/
    (*p_a)[4][1] = 3.0;
    (*p_a)[4][2] = 2.0;
    (*p_a)[4][3] = 6.0;
    (*p_a)[4][4] = 3.0;
    (*p_a)[4][5] = 2.0;
    (*p_a)[4][6] = 1.0;
    (*p_a)[4][7] = 6.0;
    (*p_a)[4][8] = 1.0;
    (*p_a)[4][9] = 7.0;
    (*p_a)[4][10] = 3.0;
    /* 5ème ligne*/
    (*p_a)[5][1] = 6.0;
    (*p_a)[5][2] = 6.0;
    (*p_a)[5][3] = 6.0;
    (*p_a)[5][4] = 4.0;
    (*p_a)[5][5] = 5.0;
    (*p_a)[5][6] = 2.0;
    (*p_a)[5][7] = 2.0;
    (*p_a)[5][8] = 4.0;
    (*p_a)[5][9] = 3.0;
    (*p_a)[5][10] = 2.0;
    /* 6ème ligne*/
    (*p_a)[6][1] = 5.0;
    (*p_a)[6][2] = 5.0;
    (*p_a)[6][3] = 2.0;
    (*p_a)[6][4] = 1.0;
    (*p_a)[6][5] = 3.0;
    (*p_a)[6][6] = 5.0;
    (*p_a)[6][7] = 5.0;
    (*p_a)[6][8] = 7.0;
    (*p_a)[6][9] = 4.0;
    (*p_a)[6][10] = 3.0;
    /* 7ème ligne*/
    (*p_a)[7][1] = 3.0;
    (*p_a)[7][2] = 6.0;
    (*p_a)[7][3] = 6.0;
    (*p_a)[7][4] = 3.0;
    (*p_a)[7][5] = 1.0;
    (*p_a)[7][6] = 6.0;
    (*p_a)[7][7] = 1.0;
    (*p_a)[7][8] = 6.0;
    (*p_a)[7][9] = 7.0;
    (*p_a)[7][10] =1.0;
    /* 8ème ligne*/
    (*p_a)[8][1] = 1.0;
    (*p_a)[8][2] = 2.0;
    (*p_a)[8][3] = 1.0;
    (*p_a)[8][4] = 7.0;
    (*p_a)[8][5] = 8.0;
    (*p_a)[8][6] = 7.0;
    (*p_a)[8][7] = 6.0;
    (*p_a)[8][8] = 5.0;
    (*p_a)[8][9] = 8.0;
    (*p_a)[8][10] = 7.0;
    /* 9ème ligne*/
    (*p_a)[9][1] = 8.0;
    (*p_a)[9][2] = 5.0;
    (*p_a)[9][3] = 2.0;
    (*p_a)[9][4] = 5.0;
    (*p_a)[9][5] = 3.0;
    (*p_a)[9][6] = 8.0;
    (*p_a)[9][7] = 1.0;
    (*p_a)[9][8] = 3.0;
    (*p_a)[9][9] = 3.0;
    (*p_a)[9][10] =5.0;
    /* 10ème ligne*/
    (*p_a)[10][1] = 1.0;
    (*p_a)[10][2] = 1.0;
    (*p_a)[10][3] = 1.0;
    (*p_a)[10][4] = 1.0;
    (*p_a)[10][5] = 1.0;
    (*p_a)[10][6] = 1.0;
    (*p_a)[10][7] = 1.0;
    (*p_a)[10][8] = 1.0;
    (*p_a)[10][9] = 1.0;
    (*p_a)[10][10] = 1.0;
    
    
    /* membre de droite des contraintes*/
    (*p_b)[1] = 38.0;
    (*p_b)[2] = 41.0;
    (*p_b)[3] = 35.0;
    (*p_b)[4] = 40.0;
    (*p_b)[5] = 47.0;
    (*p_b)[6] = 48.0;
    (*p_b)[7] = 46.0;
    (*p_b)[8] = 90.0;
    (*p_b)[9] = 41.0;
    (*p_b)[10] = 14.0;
    
    
    /* borne sup des variables x_i*/
    (*p_u)[1] = (*p_u)[2] = 1; (*p_u)[3] =1; (*p_u)[4] = 1; (*p_u)[5] = 1; (*p_u)[6] = 1; (*p_u)[7]=2; (*p_u)[8]=2;
    (*p_u)[9] = 2; (*p_u)[10] =4;
    
    /* un w0 de départ*/
    (*p_w0)[1] = 1; (*p_w0)[2] = 2; (*p_w0)[3] = 3; (*p_w0)[4] = 4; (*p_w0)[5] = 5; (*p_w0)[6] = 6; (*p_w0)[7] = 7; (*p_w0)[8]=8;
    (*p_w0)[9] = 9; (*p_w0)[10] = 10;
    
    
    
    /* un w0 de départ* /
    (*p_w0)[1] = 1; (*p_w0)[2] = 2; (*p_w0)[3] = 1; (*p_w0)[4] = 2; (*p_w0)[5] = 1; (*p_w0)[6] = 2; (*p_w0)[7] = 1; (*p_w0)[8]=2;
    (*p_w0)[9] = 1; (*p_w0)[10] = 2;     */
    
    
} /* genere_exemple_4 */

/******************************************************************/


void ecrire_fichier_mps(int n, int m, double_1D c, double_1D d, double_2D a,
double_1D b, int_1D u, bool reso_ent, char_1D nom_fichier)
{
    FILE *fich;
    int   i, j, cpt;
    
    fich = fopen(nom_fichier, "w");
    
    fprintf(fich, "NAME %s\n", nom_fichier);
    fprintf(fich, "OBJSENSE\n MAX\n");
    fprintf(fich, "ROWS\n");
    fprintf(fich, " N  obj\n");
    for (j = 1; j <= m; j++) fprintf(fich, " L c%d\n", j);
    fprintf(fich, "COLUMNS");
    for (i = 1; i <= n; i++) {
        fprintf(fich, "\n x%d obj %lf\n x%d", i, c[i], i);
        cpt = 0;
        for (j = 1; j <= m; j++) {
            fprintf(fich, " c%d %lf", j, a[j][i]); cpt++;
            if ( ((cpt % 2) == 0) && (j < m))fprintf(fich, "\n x%d", i);
        }
    }
    fprintf(fich, "\nRHS\n rhs");
    cpt = 0;
    for (j = 1; j <= m; j++) {
        fprintf(fich, " c%d %lf", j, b[j]); cpt++;
        if ( ((cpt % 2) == 0) && (j < m)) fprintf(fich, "\n rhs");
    }
    fprintf(fich, "\nBOUNDS\n");
    for (i = 1; i <= n; i++) {
        if (reso_ent) fprintf(fich, " LI"); else fprintf(fich, " LO");
        fprintf(fich, " BOUND x%d 0\n", i);
        if (reso_ent) fprintf(fich, " UI"); else fprintf(fich, " UP");
        fprintf(fich, " BOUND x%d %d\n", i, u[i]);
    }
    fprintf(fich, "QMATRIX\n");
    for (i = 1; i <= n; i++) fprintf(fich, " x%d x%d %lf\n", i, i, - 2.0 * d[i]);
    fprintf(fich, "ENDATA");
    
    fclose(fich);
} /* ecrire_fichier_mps */

/******************************************************************/
/* PRE-TRAITEMENT A AJOUTER : supprimerr les contraintes trivialement vérifiées pour tout x */

bool probleme_trivial(int n, int m, double_2D a, double_1D b, int_1D u)
{
    bool trivial = VRAI;
    double somme_aji_ui;
    int    i, j;
    
    for (j = 1; ((j <= m) && (trivial)); j++) {
        somme_aji_ui = 0.0;
        for (i = 1; i <= n; i++) somme_aji_ui += a[j][i] * u[i];
        /* printf("  j = %d somme_aji_ui = %f \n", j, somme_aji_ui);  */
        if (somme_aji_ui > b[j]) trivial = FAUX;
    }
    
    return(trivial);
    
} /* probleme_trivial */

/******************************************************************/

void pqe_main_calcul_borne3(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b, 
double *p_borne3, double *p_temps, double_1D vec_dul, 
double_1D solution_y_borne3, int_1D copycstat, int_1D copyrstat, 
int_1D getcstat, int_1D getrstat)
{
    /* --------------------------------------------------------- */
    /* Resolution du pb linearise non agrege en continu : borne3 */ 
    /* --------------------------------------------------------- */
    /* On suppose l'environnement de cplex déjà ouvert           */
    /* --------------------------------------------------------- */
    
    int       nbcols, nbvar, nblignes, cpt, i, j, k, statut;
    double_1D lb, ub; /* retrait de solution_y; car mis en arguments et changement de nom pr diff de solution_y de borne2*/
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense, ctype;
    bool      bidon;
    
    nbcols = nbvar = somme_ui(n, u);
    nblignes = m;
    print("\n nbcols = %d nblignes = %d nbvar = %d\n", nbcols, nblignes, nbvar);
    pqe_common_init_temps(p_temps);
    
    pqe_init_1D((void **)  &coef_obj,   sizeof(float),  nbvar);
    pqe_init_2D((void ***) &mat_cont,   sizeof(float),  nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs,    sizeof(float),  nblignes);
    pqe_init_1D((void **)  &sense,      sizeof(char),   nblignes);
    pqe_init_1D((void **)  &lb,         sizeof(double), nbvar);
    pqe_init_1D((void **)  &ub,         sizeof(double), nbvar);
    pqe_init_1D((void **)  &ctype,      sizeof(char),   nbvar);
    /*pqe_init_1D((void **)  &solution_y, sizeof(double), nbcols);*/
    
    /* Coef obj */
    cpt = 0 ;
    for (i = 1; i<= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
        coef_obj[cpt] = (float) s[i][k];
        ctype[cpt] = 'C';
        lb[cpt] = 0.0;
        ub[cpt] = 1.0;
        cpt++;
    }
    
    /* coefficients de la matrice des contraintes */
    for (j = 1; j <= m; j++){ 
        cpt = 0;
        for (i = 1; i <= n; i++) 
        for (k = 1; k <= u[i]; k++) {
            
            mat_cont[j-1][cpt] = a[j][i];
            cpt++;
        }
    }
    
    for (j = 1; j <= m; j++) {
        vec_rhs[j-1]    = b[j];
        sense[j-1] = 'L';
    }
    
    *p_borne3 = (double) pqe_cplex_reso_lp_mixte(nblignes, nbcols, coef_obj,
                            mat_cont, vec_rhs, &statut, solution_y_borne3, 
                            vec_dul, sense, ctype, lb, ub, FAUX, VRAI, NULL,
                            copycstat, copyrstat, getcstat, getrstat) ;
    
    
    /*Essai pour voir si la résolution du linéarisé entière donne une solution        optimale alors que lorsqu'on résout le quadratique nous n'en avons pas p        our certains problèmes de grande taille*/
    
    /*
    *p_borne3 = (double) pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj, mat_cont, vec_rhs, &statut, solution_y, sense, VRAI, NULL, FAUX, NULL, &bidon);
    
    */  
    /* Decalage d'indices sur vec_dul */
    for (j = m; j >= 1; j--) {
        vec_dul[j] = vec_dul[j-1];
    }
    /* Decalage d'indices sur solution_y_borne3*/
    for (i = nbcols; i >= 1; i--) {
        solution_y_borne3[i] = solution_y_borne3[i-1];
        print("\r solution_y_borne3[%d] = %lf\n", i, solution_y_borne3[i]);
    }
    
    pqe_common_calcule_temps(p_temps);
    print("\n borne3 = %f (%.2f s) \n", *p_borne3, *p_temps);
    
    pqe_init_free_1D(coef_obj);
    pqe_init_free_2D((void *) mat_cont,  nblignes);
    pqe_init_free_1D(vec_rhs);
    pqe_init_free_1D(sense);
    pqe_init_free_1D(lb);
    pqe_init_free_1D(ub);
    pqe_init_free_1D(ctype);
    
} /* pqe_main_calcul_borne3 */

/******************************************************************/

void calcul_opt_linea(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b, 
double *p_opt_linea, double *p_temps, bool sol_adm, double_1D primaly,
bool *p_temps_limite, double *p_nb_noeuds)
{
    /* --------------------------------------------------------- */
    /* Resolution du pb linearise non agrege en continu : borne3 */ 
    /* --------------------------------------------------------- */
    /* On suppose l'environnement de cplex déjà ouvert           */
    /* --------------------------------------------------------- */
    
    int       nbcols, nbvar, nblignes, cpt, i, j, k, statut;
    double_1D lb, ub, solution_y;
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense, ctype;
    
    nbcols = nbvar = somme_ui(n, u);
    nblignes = m;
    pqe_common_init_temps(p_temps);
    
    pqe_init_1D((void **)  &coef_obj,   sizeof(float),  nbvar);
    pqe_init_2D((void ***) &mat_cont,   sizeof(float),  nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs,    sizeof(float),  nblignes);
    pqe_init_1D((void **)  &sense,      sizeof(char),   nblignes);
    pqe_init_1D((void **)  &lb,         sizeof(double), nbvar);
    pqe_init_1D((void **)  &ub,         sizeof(double), nbvar);
    pqe_init_1D((void **)  &ctype,      sizeof(char),   nbvar);
    pqe_init_1D((void **)  &solution_y, sizeof(double), nbcols);
    
    /* Coef obj */
    cpt = 0 ;
    for (i = 1; i<= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
        coef_obj[cpt] = (float) s[i][k];
        ctype[cpt] = 'C';
        lb[cpt] = 0.0;
        ub[cpt] = 1.0;
        cpt++;
    }
    
    /* coefficients de la matrice des contraintes */
    for (j = 1; j <= m; j++){ 
        cpt = 0;
        for (i = 1; i <= n; i++) 
        for (k = 1; k <= u[i]; k++) {
            
            mat_cont[j-1][cpt] = a[j][i];
            cpt++;
        }
    }
    
    for (j = 1; j <= m; j++) {
        vec_rhs[j-1]    = b[j];
        sense[j-1] = 'L';
    }
    
    /*Essai pour voir si la résolution du linéarisé entière donne une solution        optimale alors que lorsqu'on résout le quadratique nous n'en avons pas p        our certains problèmes de grande taille*/
    
    
    *p_opt_linea = (double) pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj, 
    mat_cont, vec_rhs, &statut, solution_y, sense, VRAI, NULL, 
    sol_adm, primaly, p_temps_limite, p_nb_noeuds);
    
    
    
    pqe_common_calcule_temps(p_temps);
    printf("\nopt_linea = %.2f (%.1f s) \n", *p_opt_linea, *p_temps);
    
    pqe_init_free_1D(coef_obj);
    pqe_init_free_2D((void *) mat_cont,  nblignes);
    pqe_init_free_1D(vec_rhs);
    pqe_init_free_1D(sense);
    pqe_init_free_1D(lb);
    pqe_init_free_1D(ub);
    pqe_init_free_1D(ctype);
    pqe_init_free_1D(solution_y);
    
} /* calcul_opt_linea */

/******************************************************************/

void resout_pb_aux(int n, int m, double_1D cprime, double_2D a, double_1D bprime,
double cste, double *p_borne, double_1D solution)
{
    /* ----------------------------------------------------------- */
    /* Amélioration de la borne par appel à cplex en 0-1 : borne 2 */
    /* ----------------------------------------------------------- */
    
    int       nbcols, nblignes, nbvar, cpt, i, j, k, statut;
    double    bidon2;
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense;
    bool      bidon;
    
    nbcols = nbvar = n;
    nblignes = m;
    
    
    pqe_init_1D((void **)  &coef_obj,   sizeof(float),  nbvar);
    pqe_init_2D((void ***) &mat_cont,   sizeof(float),  nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs,    sizeof(float),  nblignes);
    pqe_init_1D((void **)  &sense,      sizeof(char),   nblignes);
    
    
    
    /* Coef obj */
    cpt = 0 ;
    for (i = 1; i<= n; i++) {
        coef_obj[cpt] = (float) cprime[i];
        cpt++;
    }
    
    
    for (j = 1 ; j <= m; j++) 
    for (i=1; i <= n; i++) mat_cont[j-1][i-1] = a[j][i];
    
    for (j = 1; j <= m; j++) {
        vec_rhs[j-1] = bprime[j];
        sense[j-1]   = 'L';
    }
    
    *p_borne = cste + (double) pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj,
    mat_cont, vec_rhs, &statut, solution, 
    sense, VRAI, NULL, FAUX, NULL, &bidon, &bidon2) ;
    /* Decalage d'indices sur solution_y*/
    for(i=nbcols; i >= 1; i--){
        solution[i] = solution[i-1];
        print("\r solution[%d] = %lf\n", i, solution[i]);
    }
    
    
    pqe_init_free_1D(coef_obj);
    pqe_init_free_2D((void *) mat_cont,  nblignes);
    pqe_init_free_1D(vec_rhs);
    pqe_init_free_1D(sense);
    
} /* resout_pb_aux */

/******************************************************************/

void pqe_main_calcul_borne2(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
double_1D vec_agregation, double *p_borne2, double *p_temps, double_1D solution_y)
{
    /* ----------------------------------------------------------- */
    /* Amélioration de la borne par appel à cplex en 0-1 : borne 2 */
    /* ----------------------------------------------------------- */
    
    int       nbcols, nblignes, nbvar, cpt, i, j, k, statut;
    double    bidon2;
    double_1D A;/* ici retrait de solution_y car mise en argument*/
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense;
    bool      bidon;
    
    nbcols = nbvar = somme_ui(n, u);
    nblignes = 1;
    
    pqe_init_1D((void **)  &A,          sizeof(double), n);
    pqe_init_1D((void **)  &coef_obj,   sizeof(float),  nbvar);
    pqe_init_2D((void ***) &mat_cont,   sizeof(float),  nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs,    sizeof(float),  nblignes);
    pqe_init_1D((void **)  &sense,      sizeof(char),   nblignes);
    /*pqe_init_1D((void **)  &solution_y, sizeof(double), nbcols);*/
    
    pqe_common_init_temps(p_temps);
    
    /* Coef obj */
    cpt = 0 ;
    for (i = 1; i<= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
        coef_obj[cpt] = (float) s[i][k];
        cpt++;
    }
    
    /* coefficients de la contrainte aggregée */
    /* Aggregation des contraintes */
    for (i = 1; i <= n; i++) {
        A[i] = 0.0;
        for (j = 1; j <= m; j++) A[i] += vec_agregation[j] * a[j][i];
    }
    
    /* contrainte aggregée avec les variable yik */
    
    cpt = 0;
    for (i=1; i <= n; i++)
    for (j = 1 ; j <= u[i]; j++) 
    mat_cont[0][cpt++] = A[i];
    
    vec_rhs[0] = 0.0;
    for (j = 1; j <= m; j++) vec_rhs[0] += vec_agregation[j] * b[j];
    
    sense[0] = 'L';
    
    *p_borne2 = (double) pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj,
    mat_cont, vec_rhs, &statut, solution_y, 
    sense, VRAI, NULL, FAUX, NULL, &bidon, &bidon2) ;
    /* Decalage d'indices sur solution_y*/
    for (i = nbcols; i >= 1; i--) {
        solution_y[i] = solution_y[i-1];
        print("\r solution_y[%d] = %lf\n", i, solution_y[i]);
    }
    
    pqe_common_calcule_temps(p_temps);
    
    pqe_init_free_1D(coef_obj);
    pqe_init_free_2D((void *) mat_cont,  nblignes);
    pqe_init_free_1D(vec_rhs);
    pqe_init_free_1D(sense);
} /* pqe_main_calcul_borne2 */
    
/******************************************************************/

bool pqe_main_appartient(int i, int_1D I, int taille)
{
    int k;
    bool trouve = FAUX;

    k = 1;
    while ((trouve == FAUX) && (k <= taille)) {
        if (i == I[k]) trouve = VRAI;
	k++;
    }
    return(trouve);

} /* appartient */

/******************************************************************/

void pqe_main_calcul_borne2_I(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
                              double_1D vec_agregation, int_1D I, int taille, double *p_borne2, 
			      double *p_temps, double_1D solution_y)
{
    /* ----------------------------------------------------------- */
    /* Amélioration de la borne par appel à cplex en 0-1 : borne 2 */
    /* ----------------------------------------------------------- */
    
    int       nbcols, nblignes, nbvar, cpt, i, j, k, statut;
    double    bidon2;
    double_1D A;/* ici retrait de solution_y car mise en argument*/
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense;
    bool      bidon;
    
    nbcols = nbvar = somme_ui(n, u);
    nblignes = 2;
    
    pqe_init_1D((void **)  &A,          sizeof(double), n);
    pqe_init_1D((void **)  &coef_obj,   sizeof(float),  nbvar);
    pqe_init_2D((void ***) &mat_cont,   sizeof(float),  nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs,    sizeof(float),  nblignes);
    pqe_init_1D((void **)  &sense,      sizeof(char),   nblignes);
    /*pqe_init_1D((void **)  &solution_y, sizeof(double), nbcols);*/
    
    pqe_common_init_temps(p_temps);
    
    /* Coef obj */
    cpt = 0 ;
    for (i = 1; i<= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
        coef_obj[cpt] = (float) s[i][k];
        cpt++;
    }
    
    /* coefficients de la contrainte aggregée */
    /* Aggregation des contraintes */
    for (i = 1; i <= n; i++) {
        A[i] = 0.0;
        for (j = 1; j <= m; j++) A[i] += vec_agregation[j] * a[j][i];
    }
    
    /* contrainte aggregée avec les variable yik */
    
    cpt = 0;
    for (i=1; i <= n; i++)
       for (j = 1 ; j <= u[i]; j++) {
           mat_cont[0][cpt] = A[i];
           mat_cont[1][cpt] = (pqe_main_appartient(i, I, taille) ? 1.0 : 0.0);
	   cpt++;
       }
    
    vec_rhs[0] = 0.0;
    for (j = 1; j <= m; j++) vec_rhs[0] += vec_agregation[j] * b[j];
    vec_rhs[1] = 1.0;
    
    sense[0] = 'L';
    sense[1] = 'G';
    
    *p_borne2 = (double) pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj,
                                             mat_cont, vec_rhs, &statut, solution_y, 
                                             sense, VRAI, NULL, FAUX, NULL, &bidon, &bidon2) ;
    /* Decalage d'indices sur solution_y*/
    for (i = nbcols; i >= 1; i--) {
        solution_y[i] = solution_y[i-1];
        print("\r solution_y[%d] = %lf\n", i, solution_y[i]);
    }
    
    pqe_common_calcule_temps(p_temps);
    
    pqe_init_free_1D(coef_obj);
    pqe_init_free_2D((void *) mat_cont,  nblignes);
    pqe_init_free_1D(vec_rhs);
    pqe_init_free_1D(sense);
    
} /* pqe_main_calcul_borne2_I */

/******************************************************************/

void affichage_resu_ecran(double valeur_continue, double borne1, double borne2, double borne3, double valeur_opt,
double tps_continu, double temps1, double temps2, double temps3, double tps_opt)
{
    printf(" ------------------------------------------------------------------\n");
    printf("        | Rel. cont | borne 1   | borne 3   | borne 2   | optimum   |\n");
    printf(" ------------------------------------------------------------------\n");
    printf(" Valeur | %9.2f | %9.2f | %9.2f | %9.2f | %9.2f |\n",
    valeur_continue, borne1, borne3, borne2, valeur_opt);
    printf(" ------------------------------------------------------------------\n");
    printf(" Temps s| %9.2f | %9.2f | %9.2f | %9.2f | %9.2f |\n",
    tps_continu, temps1, temps3, temps3+temps2, tps_opt);
    printf(" ------------------------------------------------------------------\n");
    
} /* affichage_resu_ecran */

/******************************************************************/

void saut_colonne(FILE *fich, int n)
{
    int i;
    for (i = 1; i <= n; i++) fprintf(fich, "\t");
} /* saut_colonne */

/******************************************************************/

void ecrit_fichier_borne(FILE *fich, double valeur, double temps, char *nom_borne)
{
    /* On suppose le fichier fich deja ouvert en ecriture */
    
    fprintf(fich, "%s", nom_borne);
    saut_colonne(fich, 1);
    fprintf(fich, "%f", valeur);
    saut_colonne(fich, 1);
    fprintf(fich, "%.2f", temps);
    saut_colonne(fich, 3);
    
    printf("\n %s = %f (%.2f s)\n", nom_borne, valeur, temps);
} /* ecrit_fichier_borne */

/******************************************************************/

void ecrit_fichier_opt(FILE *fich, bool temps_limite, double valeur, double temps, 
                       double nb_noeuds, char *nom_opt, bool existe_pretraitement,
		       double temps_pretraitement)
{
    /* On suppose le fichier fich deja ouvert en ecriture */
    
    fprintf(fich, "%s", nom_opt);
    saut_colonne(fich, 1);
    fprintf(fich, "%d", temps_limite);
    saut_colonne(fich, 1);
    fprintf(fich, "%f", valeur);
    saut_colonne(fich, 1);
    fprintf(fich, "%.2f", temps);
    saut_colonne(fich, 1);
    fprintf(fich, "%.0f", nb_noeuds);
    if (existe_pretraitement) { 
	    saut_colonne(fich, 1);
            fprintf(fich, "%.2f", temps_pretraitement);
    }
    saut_colonne(fich, 2);
    
    printf("\n %s = %f (%.2f s) nb_noeuds = %.0f - temps_limite = %d \n", 
    nom_opt, valeur, temps, nb_noeuds, temps_limite);
} /* ecrit_fichier_opt */


/******************************************************************/

bool pqe_main_est_admissible_y(double_2D y, int n, int m, double_1D b, double_2D a, int_1D u)
{
    int i, j, k;
    bool adm;
    double ecart, somme_aji_xi, xi;
    const double epsilon_adm = 1e-5;
    
    adm = VRAI;
    
    for (j = 1; ((j <= m) && (adm)); j++) {
        somme_aji_xi = 0.0;
        for (i = 1; i <= n; i++) {
            /* Calcul de xi */
            xi = 0.0;
            for(k=1;k<=u[i];k++){
                xi += y[i][k];
            }
            somme_aji_xi += a[j][i] * xi;
        }	
        
        if (somme_aji_xi > b[j] + epsilon_adm) adm = FAUX;
        print("\n adm = %d somme_aji_xi = %.2f et b[%d] = %.2f\n", adm, somme_aji_xi, j, b[j]);
    }
    
    return adm;
    
} /* pqe_main_est_admissible_y */

/******************************************************************/

bool pqe_main_est_admissible_y1D(double_1D y, int n, int m, double_1D b, double_2D a, int_1D u)
{
    int i, j, k, cpt;
    bool adm;
    double ecart, somme_aji_xi, xi;
    double_1D  x;
    const double epsilon_adm = 1e-5;
    
    pqe_init_1D((void **) &x, sizeof(double), n);
    cpt = 1;
    for (i = 1; i <= n; i++) {
        /* Calcul de xi */
        x[i] = 0.0;
        for(k=1;k<=u[i];k++){
            x[i] += y[cpt];
            cpt++;
        }
    }
    adm = VRAI;
    for (j = 1; ((j <= m) && (adm)); j++) {
        somme_aji_xi = 0.0;
        for (i = 1; i <= n; i++) {
            somme_aji_xi += a[j][i] * x[i];
        }	
        
        print(" j = %d    somme_aji_xi = %f   b[%d] = %f \n", j, somme_aji_xi, j, b[j]);
        if (somme_aji_xi > b[j] + epsilon_adm) adm = FAUX;
        print("\n adm = %d somme_aji_xi = %.2f et b[%d] = %.2f\n", adm, somme_aji_xi, j, b[j]);
    }
    
    pqe_init_free_1D(x);
    return (adm);
    
} /* pqe_main_est_admissible_y */

/******************************************************************/

/*bool est_admissible(double_1D vec_x, int n, int m, double_1D b, double_2D a)*/
bool est_admissible(eltx_1D vec_x, int n, int m, double_1D b, double_2D a)
{
    int i, j;
    bool adm;
    double ecart, somme_aji_xi;
    
    adm = VRAI;
    for(j=1; ((j<=m) && (adm));j++){
        somme_aji_xi = 0.0;
        /*for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * vec_x[i];*/
        for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * vec_x[i].val_x;
        ecart = somme_aji_xi - b[j];
        if(ecart > 0.0) adm = FAUX;
    }
    
    return adm;
    
} /* est_admissible */

/******************************************************************/

bool est_admissible_x(double_1D vec_x, int n, int m, double_1D b, double_2D a)
{
    int i, j;
    bool adm;
    double ecart, somme_aji_xi;
    const double precision_ecart = 0.0000001;
    
    adm = VRAI;
    for(j=1; ((j<=m) && (adm));j++){
        somme_aji_xi = 0.0;
        /*for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * vec_x[i];*/
        for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * vec_x[i];
        ecart = somme_aji_xi - b[j];
        if(ecart > 0.0 + precision_ecart) {
            adm = FAUX;
            print("\n j = % d, b[%d] = %f, somme_aji_xi = %f\n", j, j, b[j], somme_aji_xi);
        }
    }
    
    return adm;
    
} /* est_admissible_x */

/******************************************************************/
/******************************************************************/

bool est_admissible_x_int(int_1D vec_x, int n, int m, double_1D b, double_2D a)
{
    int i, j;
    bool adm;
    double ecart, somme_aji_xi;
    const double precision_ecart = 0.0000001;
    
    adm = VRAI;
    for(j=1; ((j<=m) && (adm));j++){
        somme_aji_xi = 0.0;
        /*for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * vec_x[i];*/
        for(i=1; i<=n; i++) somme_aji_xi += a[j][i] * (double) vec_x[i];
        ecart = somme_aji_xi - b[j];
        if(ecart > 0.0 + precision_ecart) {
            adm = FAUX;
            printf("\n j = % d, b[%d] = %f, somme_aji_xi = %f\n", j, j, b[j], somme_aji_xi);
        }
    }
    
    return adm;
    
} /* est_admissible_x_int */



void pqe_main_primalisation(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
double_1D c, double_1D d, double constante, double_1D solution_y,
double *p_borne_inf, bool *p_optimum_trouve, double_1D primalx,
double_1D primaly, double borne2, double *p_temps)
{
    /* ------------------------------------------------ */
    /* Primalisation de la solution solution_y fournie par borne 2 */
    /* Calcul d'une solution admissible dans primalx et primaly de valeur *p_borne_inf                */
    /* si solution_y est admissible pour le pb de depart, on a trouve l'optimum */
    /* ------------------------------------------------ */
    
    int           i, j, k, cpt, indice, indice_i, indice_k, iy;
    bool          adm, trouve_x, encore, trouve_y, adm2, ilfautdetruire;
    double        somme_yik, ecart, ecartmax, somme_aji_xi, borne_inf_amelioree;
    double_1D     tab_ecart;
    double_2D     y_b2; /* solution en y fournie par le calcul de borne 2 */
    p_cellule     L, p, pmax, prec, adetruire;
    contrainte_1D tab_cont;
    coeff_1D      tab_coeff; /* les coeff des y_ik de borne_inf*/
    const double  eps_ecart = 1e-6;
    /* elt_coeff     tab_coeff[10000]; */
    
    elty_2D       vec_y;
    
    elt_y         somme_y  /* ,      vec_y[1000][1000]*/;
    eltx_1D       vec_x;
    
    *p_optimum_trouve = FAUX;
    
    pqe_init_1D((void **)  &vec_x,    sizeof(elt_x), n);
    print("\n vec_x = %p\n", vec_x);
    
    pqe_common_init_temps(p_temps);
    pqe_init_1D((void **)  &tab_ecart, sizeof(double), m);
    print("\n tab_ecart = %p\n", tab_ecart);
    /* --------- Ecriture de la solution en y_ik --------- */
    
    pqe_init_2D((void ***) &y_b2, sizeof(double), n, pqe_main_max_ui(n, u));
    print("\n y_b2 = %p pqe_main_max_ui(n,u) = %d \n", y_b2, pqe_main_max_ui(n, u));
    /*
    y_b2 = (double **) malloc(sizeof(double *)*(n + 1));
    for(i = 1; i < 1 + n; i++) y_b2[i] = (double *) malloc(sizeof(double)*(1 + u[i]));
    */
    
    pqe_init_1D((void **)  &tab_coeff, sizeof(elt_coeff), somme_ui(n,u)); 
    print("\n tab_coeff = %p\n", tab_coeff);
    
    /*
    print("\n somme_ui = %d \n", somme_ui(n,u));
    for(i=1; i<=n; i++) printf("u[%d]=%d ", i, u[i]);
    print("\n");
    */
    
    
    
    cpt = 1;
    for(i = 1; i <= n; i++)
    for(k = 1; k <= u[i]; k++) {
        y_b2[i][k] = solution_y[cpt];
        print("\r y_b2[%d][%d]= %lf\n ", i, k, y_b2[i][k]);
        cpt++;
    }
    
    /*--------- Ecriture de la solution en x_i ---------- */
    
    for(i = 1; i <= n; i++) {
        somme_yik = 0.0;	
        for(k = 1; k <= u[i]; k++){
            somme_yik = somme_yik + y_b2[i][k];
        }
        /* vec_x[i] = somme_yik;
        print("\r vec_x[%d]= %lf\n", i, vec_x[i]);*/
        vec_x[i].val_x = somme_yik;
        vec_x[i].modifier = FAUX;
        print("\r vec_x[%d].val_x= %lf\n", i, vec_x[i].val_x);
    }
    
    /* admissibilite en x et creation de la liste des ecarts */
    
    adm = VRAI;
    L   = NULL;
    
    for(j = 1; j <= m; j++) {      
        somme_aji_xi = 0.0;  
        /*for(i = 1; i <= n; i++) somme_aji_xi += a[j][i] * vec_x[i];*/
        for(i = 1; i <= n; i++) somme_aji_xi += a[j][i] * vec_x[i].val_x;
        ecart = somme_aji_xi - b[j];
        if(ecart > 0.0) {
            print(" somme_aji_xi = %lf bj = %lf\n", somme_aji_xi, b[j]); 
            adm         = FAUX; 
            p           = (p_cellule) malloc (sizeof(cellule));
            p->ecart    = ecart;
            p->ind_cont = j;
            p->suivant  = L;
            L           = p;
        }
    }
    
    if (!adm) {
        print("\r Solution non admissible\n");
        p  = L;
        while (p != NULL) {
            print(" Contrainte %d : ecart = %lf\n", p->ind_cont, p->ecart);
            p = p->suivant;
        }
        
        /* initialisation du tableau des contraintes */
        pqe_init_1D((void **)  &tab_cont, sizeof(contrainte), m);
        
        for (j = 1; j <= m; j++) {
            tab_cont[j].triee = FAUX;
            tab_cont[j].tab   = NULL;
        }
        
        while (!adm) {
            print("\n entree ds la boucle while !adm\n");
            /* Recherche de la contrainte j la plus violée */
            
            p = L; 
            pmax = L;
            ecartmax = L->ecart;
            
            while (p != NULL) {
                if (p->ecart > ecartmax) {
                    pmax = p;
                    ecartmax = p->ecart;
                }
                p = p->suivant;
            }
            
            j = pmax->ind_cont;
            
            if (!tab_cont[j].triee) {
                /* on trie alors les coeff de la contrainte j */
                pqe_init_1D((void **)  &(tab_cont[j].tab), sizeof(elt), n);
                
                for (i = 1; i <= n; i++) {
                    tab_cont[j].tab[i].cle = a[j][i];
                    tab_cont[j].tab[i].ind = i;
                }
                
                pqe_main_tri_rapide(1, n, tab_cont[j].tab);
                
                tab_cont[j].triee = VRAI;
            }
            
            trouve_x = FAUX;
            i = n;
            while (!trouve_x) {
                indice = tab_cont[j].tab[i].ind;
                
                if (vec_x[indice].val_x > 0) {
                    print("l'indice de la var non nulle de plus fort coeff est %d ", indice);
                    print("dans la contrainte %d\n", j);
                    (vec_x[indice].val_x)--;
                    (vec_x[indice].modifier) = VRAI;
                    trouve_x = VRAI;
                }
                else i--;
            }
            
            /* printf("\n nouvelle valeur de vec_x[%d] : %.1f \n", indice,vec_x[indice]); */
            print("\n nouvelle valeur de vec_x[%d].val_x : %.1f \n", indice,vec_x[indice].val_x); 
            print("\n et modifier vaut %d \n", vec_x[indice].modifier);
            
            /* Mise à jour de la liste des ecarts */
            
            p = L;
            
            while (p != NULL) {
                print("entree dans la liste tt que non nulle\n");
                (p->ecart) -= a[p->ind_cont][indice];
                p = p->suivant;
            }
            
            /* Retrait des maillons ayant un ecart <= 0.0 */
            /* on commence par retirer les elts nul en debut de liste */
            encore = VRAI;
            while (encore) {
                print("while encore\n");
                if (L == NULL) encore = FAUX;
                else {
                    if (L->ecart <= 0.0) { print("suppression en tete \n "); L = L->suivant; }
                    else encore = FAUX;
                }
            }
            print("sortie de while encore\n");
            if (L != NULL) {
                /* le premier elt est >0 et on retire les elts <= 0 du reste de la liste */
                print("entrer dans liste non vide\n");
                prec = L;
                p = L->suivant;
                while (p != NULL) {
                    print("while p!= null\n");
                    if (p->ecart <= 0.0) {
                        print("entrer ds le si\n");
                        prec->suivant = p->suivant;
                        adetruire = p; 
                        p = p->suivant;
                        free(adetruire);
                    }
                    else {
                        prec = p;
                        p = p->suivant;
                    }
                    print("top avant free\n");
                    print("\n adetruire vaut %p\n", adetruire);
                    print("top après free\n");
                }
                
            }
            print("sortie de la liste\n");
            adm = (L == NULL);
            
            print("\n adm vaut %d et l'admissibilite : %d\n", adm, est_admissible(vec_x, n, m, b, a));
            
        } /* while (!adm) */
        
        *p_borne_inf = constante;
        for (i = 1; i <= n; i++) {
            (*p_borne_inf) += c[i] * vec_x[i].val_x - d[i] * vec_x[i].val_x * vec_x[i].val_x;
        }
        
        print("\n borne_inf = %f\n", *p_borne_inf);
        
        /* libereration de tab_cont */
        print("\n avant lib de tab_cont\n");
        for (j = 1; j <= m; j++) 
        if (tab_cont[j].triee) pqe_init_free_1D((void *) tab_cont[j].tab);
        pqe_init_free_1D(tab_cont);
        print("\n apres lib de tab_cont\n");
        
        
        /*------------------------------------------------------*/
        /* Phase de remplissage de la sol relative à borne_inf  */
        /*------------------------------------------------------*/
        
        /* ecriture de la sol issue de borne_inf en y_ik*/
        pqe_init_2D((void ***) &vec_y, sizeof(elt_y), n, pqe_main_max_ui(n, u));
        
        
        print("\n vec_y = %p\n", vec_y);
        print("\n vec_y[%d] = %p\n", n, vec_y[n]);
        /*
        vec_y = (elt_y**)malloc(sizeof(elt_y*)*(n + 1));
        for(i = 0; i < 1 + n; i++) vec_y[i] = (elt_y*)malloc(sizeof(elt_y)*(1 + u[i]));
        */
        
        for(i=1; i<=n; i++){
            if (vec_x[i].val_x ==u[i]){
                for(k=1; k <= u[i];k++) {
                    vec_y[i][k].val_y = 1;
                    vec_y[i][k].modifier = vec_x[i].modifier;
                }/*fin for k*/
                
            }/*fin alors si x_i = u_i*/
            else { 
                if (vec_x[i].val_x != 0) {
                    for(k=1; k<= vec_x[i].val_x;k++){
                        vec_y[i][k].val_y = 1;
                        vec_y[i][k].modifier = vec_x[i].modifier;
                    }/*fin k jusqu'à val_x*/
                    for(k= (vec_x[i].val_x)+1; k<= u[i];k++){
                        vec_y[i][k].val_y = 0;
                        vec_y[i][k].modifier = vec_x[i].modifier;
                    }/*fin k jusqu'à u_i*/
                }/* fin si x_i != 0*/
                else {
                    for(k=1; k<= u[i]; k++){
                        vec_y[i][k].val_y = 0;
                        vec_y[i][k].modifier = vec_x[i].modifier;
                    }/*fin k dans x_i =0*/
                }/*fin else dans si x_i = 0*/
            }/*fin else dans si x_i=u_i*/
        }/* fin for i ecriture borne_in en y*/
        
        /* affichage de borne_inf en y_ik pour vérifier que nous avons la meme chose qu'en x*/
        for(i=1; i<=n;i++){
            for(k=1; k<=u[i];k++){
                print("\nvec_y[%d][%d] = %d et a ete modifier = %d\n", i, k, vec_y[i][k].val_y, vec_y[i][k].modifier);
            }
        }/*fin for i affichage y*/
        
        /*Ecriture du tableau des coefficients des y_ik pour les trier: ordre décroissant*/
        /* et nous prenons la var y_ik non modfiee qui a le plus fort coeff parmi les var à 0*/
        
        
        
        /*remplissons le tableau des coefficients des y_ik*/
        
        cpt = 1;
        for(i=1; i<= n;i++) {	      
            for(k=1; k<=u[i];k++){
                tab_coeff[cpt].cle = s[i][k];
                tab_coeff[cpt].ind_i = i; 
                tab_coeff[cpt].ind_k = k;
                print("\ncoeff[%d] =  %2.lf est le coeff de la var y[%d][%d]\n",cpt, tab_coeff[cpt].cle,tab_coeff[cpt].ind_i,tab_coeff[cpt].ind_k);
                cpt++;
            }
            
        }
        
        tri_rapide_coeff (1, somme_ui(n,u),tab_coeff);
        /*i = somme_ui(n,u);*/
        /*printf("\nle plus fort coeff est %lf et correspond à la variable y[%d][%d]\n", tab_coeff[i].cle, tab_coeff[i].ind_i, tab_coeff[i].ind_k);*/
        
        
        /*--------------------------------------------------------------------------------*/
        /* recherche de la var y_ik non modifee qui vaut 0 et qui a le plus fort coeff obj*/
        /*--------------------------------------------------------------------------------*/
        
        
        for (j = 1; j <= m; j++) {
            somme_aji_xi = 0.0;  
            /*for(i = 1; i <= n; i++) somme_aji_xi += a[j][i] * vec_x[i];*/
            for(i = 1; i <= n; i++) somme_aji_xi += a[j][i] * vec_x[i].val_x;
            tab_ecart[j] = somme_aji_xi - b[j];
            /* printf("tab_ecart[%d] = %f\n", j, tab_ecart[j]); */
            if (tab_ecart[j] > 0.0 + eps_ecart) printf("\n ERREUR !!! ecart > 0 : %f\n", tab_ecart[j]);
        }
        
        borne_inf_amelioree = (*p_borne_inf);
        
        for (iy = somme_ui(n,u); iy >= 1; iy--) {
            indice_i = tab_coeff[iy].ind_i;
            indice_k = tab_coeff[iy].ind_k;
            
            /* printf("\n top 4 i = %d k = %d valy = %d modifier = %d \n", indice_i, indice_k, 
            vec_y[indice_i][indice_k].val_y, vec_y[indice_i][indice_k].modifier);
            */
            if ((vec_y[indice_i][indice_k].val_y == 0) &&(!(vec_y[indice_i][indice_k].modifier))) {
                /* Examen de iy */
                adm2 = VRAI;
                /* 	    printf("adm2 vaut ici %d\n", adm2); */
                for (j = 1; ((j <= m) && (adm2)); j++) {
                    
                    /*
                    printf("a[%d][%d] = %f \n", j, indice_i, a[j][indice_i]);
                    printf("tab_ecart[%d] = %f \n", j, tab_ecart[j]);
                    */
                    if ((tab_ecart[j] + a[j][indice_i]) > 0.0) adm2 = FAUX; 
                }
                /* printf("adm2 vaut ici %d\n", adm2); */
                if (adm2) {
                    
                    print("la var nulle non modifiee est y[%d][%d] ", indice_i, indice_k);
                    print("de coefficient %lf\n", tab_coeff[iy].cle);
                    for (j = 1; j <= m; j++) tab_ecart[j] += a[j][indice_i];
                    (vec_y[indice_i][indice_k].val_y)++;
                    (vec_y[indice_i][indice_k].modifier) = VRAI;
                    borne_inf_amelioree += c[indice_i] - 2.0 * d[indice_i] * vec_x[indice_i].val_x 
                    - d[indice_i];
                    (vec_x[indice_i].val_x) ++;
                }
            }
        } /*fin du pour iy */
        print("\n borne_inf = %f et borne_inf_amelioree = %f\n", *p_borne_inf, borne_inf_amelioree);
        *p_borne_inf = borne_inf_amelioree;
        
        /* M-a-j de primaly et primalx */
        
        cpt = 0;
        for(i = 1; i <= n; i++) {
            somme_yik = 0.0;
            for(k = 1; k <= u[i]; k++) {
                primaly[cpt] = vec_y[i][k].val_y;
                cpt++;
                somme_yik += vec_y[i][k].val_y;
            }
            primalx[i-1] = somme_yik;
        }
        
        
        
        
        pqe_init_free_2D((void *) vec_y, n);
        
        
        
        
        
        
    }
    
    else {/* remarque pourquoi ne pas dire à ce niveau que Borne_inf = Borne2 = Opt?*/
        print("on est a l'optimum avec borne 2\n");
        /**p_optimum_trouve = VRAI;*/
        *p_optimum_trouve = VRAI;
        *p_borne_inf      = borne2;
    }
    
    if (*p_borne_inf == borne2) *p_optimum_trouve = VRAI;
    
    pqe_init_free_1D((void *) tab_coeff); 
    
    /* for(i = 1; i < 1 + n; i++) free(y_b2[i]); */
    pqe_init_free_2D((void *) y_b2, n);
    /* free(y_b2); / * faire ça proprement */
    pqe_init_free_1D((void *) tab_ecart);
    pqe_init_free_1D(vec_x);
    
    /*printf("\n primalx");
    for(i=0; i<n; i++) printf(" [%d] = %.2f ", i, primalx[i]);*/
    
    pqe_common_calcule_temps(p_temps);
    
} /* pqe_main_primalisation */

/******************************************************************/

void heuri_D(int n, int m, double_2D a, double_1D b, double_1D c, double_1D d, int_1D u, 
double_1D w, double_1D xD, double *infD, double *p_temps) 
{
    int       i, j, k, choix, card_I;
    double_1D R, E, b_copie;
    int_1D    xinf;
    double    ratio, valinf, Rmax, theta, minbsura;
    bool_1D   I;
    bool      fin;
    
    print("avant initialisation variables\n");
    pqe_init_1D((void **) &R, sizeof(double), n);
    pqe_init_1D((void **) &E, sizeof(double), n);
    pqe_init_1D((void **) &b_copie, sizeof(double), m);
    pqe_init_1D((void **) &xinf, sizeof(int), n);
    pqe_init_1D((void **) &I, sizeof(int), n);
    print("apres initialisation variables\n");
    pqe_common_init_temps(p_temps);
    
    *infD = 0.0;
    for (choix = 1; choix <= 3; choix++) {
        print("choix = %d\n", choix);
        /* Step 1 */
        for (i = 1; i <= n; i++) {
            print("i = %d\n", i);		
            xinf[i] = 0;
            print("apres calcul xinf\n");
            print("xinf[%d] = %d", i, xinf[i]);
            I[i]    = VRAI;
            print("apres  I[%d] \n", i);
        }
        
        card_I = n;
        print("avant recopie de b\n");
        for (j = 1; j <= m; j++) b_copie[j] = b[j];
        print("apres recopie de b\n");
        
        valinf = 0.0;
        
        fin = FAUX;
        while (!(fin)) {
            switch (choix) {
                case 1 : for (i = 1; i <= n; i++) {		 
                    if (I[i]) {
                        E[i] = 0.0;
                        for (j = 1; j <= m; j++) E[i] += a[j][i];
                    }
                }
                break; 
                case 2 : for (i = 1; i <= n; i++) {
                    if (I[i]) {
                        E[i] = 0.0;
                        for (j = 1; j <= m; j++) E[i] += a[j][i] * w[j];
                    }
                }
                break; 
                case 3 : for (i = 1; i <= n; i++) {
                    if (I[i]) {
                        E[i] = 0.0;
                        for (j = 1; j <= m; j++) E[i] += (a[j][i] / b_copie[j]);
                    }
                }
            }		
            
            k = 0;
            Rmax = - INFINI;
            print("on a fini les choix\n");
            for (i = 1; i <= n; i++) {
                if (I[i]) {
                    R[i] = (c[i] - d[i] * (1.0 + 2.0 * (double) xinf[i])) / E[i];
                    if (R[i] > Rmax) {
                        Rmax = R[i];
                        k = i;
                    }
                }
            }
            print("avant minbsura\n");
            minbsura = b_copie[1] / a[1][k];
            for (j = 2; j <= m; j++) {
                ratio = b_copie[j] / a[j][k];
                if (ratio < minbsura) minbsura = ratio;
            }
            
            theta = (minbsura < 1.0 ? minbsura : 1.0);
            print("\n choix = %d, k = %d minbsura = %f, theta = %f, card_I = %d", choix, k, minbsura, theta, card_I);
            
            if (theta < 1.0) { 
                I[k] = FAUX; 
                card_I--;
                if (card_I <= 0) fin = VRAI;
            }
            else {
                /* theta == 1 */
                print("\n xinf[k] = %d u[k] = %d, a[1][3] = %f b_copie[1] = %f", xinf[k], u[k], a[1][3], b_copie[1]);
                xinf[k]++;
                if (xinf[k] == u[k]) {
                    I[k] = FAUX; 
                    card_I--;
                    if (card_I <= 0) fin = VRAI;
                }
                /* else { */
                print(" cas theta == 1 et xinf_k different de u_k\n");
                for (j = 1; ((j <= m) && (fin == FAUX)); j++) {
                    b_copie[j] -= a[j][k];	
                    print("b_copie[%d]=%f ", j, b_copie[j]);
                    if (b_copie[j] <= epsilon) fin = VRAI;
                }
                /*  } */
            }
            
            if (est_admissible_x_int(xinf, n, m, b, a)) {print("\n xinf admissible");}
            else printf("\n xinf non admissible");
        } /* while non fin */
        
        if (est_admissible_x_int(xinf, n, m, b, a)) {print("\n xinf admissible");}
        else printf("\n xinf non admissible");
        print("\n choix = %d \n", choix);
        valinf = 0.0;
        for (i = 1; i <= n; i++) {
            valinf += (c[i] * ((double) xinf[i])) - d[i] * ((double) xinf[i] * xinf[i]);
            print("xinf[%d] = %d   ", i, xinf[i]); 
        }
        
        printf("\n choix = %d valinf = %lf\n", choix, valinf); 
        
        if (valinf > *infD) {
            *infD = valinf;
            for (i = 1; i <= n; i++) xD[i] = (double) xinf[i];
        } 
    } /* pour choix = 1, 2, 3 */
    
    pqe_init_free_1D(R);
    pqe_init_free_1D(E);
    pqe_init_free_1D(xinf);
    pqe_init_free_1D(I);
    pqe_init_free_1D(b_copie);
    
    pqe_common_calcule_temps(p_temps);
} /* heuri_D */


/******************************************************************/

void pqe_main_heur_cube(double_1D solution_y_borne3, int n, int m, double_2D a,
                        double_1D b, double_1D c, double_1D d, int_1D u, double *p_borne_inf2,
                        double_1D sol_inf2, double_1D sol_x_b3, double *p_temps_inf2)
{
    int       i, j, k, cpt;
    int_1D    nl;
    double    somme_yik, cste;
    double_1D cprime, bprime;
    double_2D y_b3;

    /*------------------------------------------------------------*/
    /*   Le cube UNITE                                            */
    /*------------------------------------------------------------*/
        
    pqe_init_1D((void **)  &nl,       sizeof(int),    n);
    pqe_init_1D((void **)  &bprime,   sizeof(double), m);
    pqe_init_1D((void **)  &cprime,   sizeof(double), n);

    pqe_common_init_temps(p_temps_inf2);
    y_b3 = (double **) malloc(sizeof(double *)*(n + 1));
    for (i = 1; i < 1 + n; i++) y_b3[i] = (double *) malloc(sizeof(double)*(1 + u[i]));
        
    /*-------- Ecriture de la solution en y_ik ----------*/
    cpt = 1;
    for (i = 1; i <= n; i++)
    for (k = 1; k <= u[i]; k++) {
        y_b3[i][k] = solution_y_borne3[cpt];
            /* print("\r y_b3[%d][%d]= %f    cpt = %d solution_y_borne3[%d] = %f \n ", 
            i, k, y_b3[i][k], cpt, cpt, solution_y_borne3[cpt]); */
        cpt++;
    }
        
    /*--------- Ecriture de la solution en x_i ---------- */
        
    for (i = 1; i <= n; i++) {
        somme_yik = 0.0;	
        for (k = 1; k <= u[i]; k++){
            print("\n y_b3[%d][%d] = %f ", i, k, y_b3[i][k]);
            somme_yik = somme_yik + y_b3[i][k];
        }
        print("\n");
            
        sol_x_b3[i] = somme_yik;
        print("\r sol_x_b3[%d]= %lf\n", i, sol_x_b3[i]);
        nl[i] = floor(sol_x_b3[i]);
        if (nl[i] == u[i]) (nl[i])--;
    }

    /* Calcul des coefficients du problème auxiliaire */
        
    cste = 0.0;
    for (i = 1; i <= n; i++) {
        cste += (c[i] * (double) nl[i]) - (d[i] * (double) nl[i] * (double) nl[i]);
        cprime[i] = c[i] - 2.0 * d[i] * ((double) nl[i]) - d[i];
    }
        
    for (j = 1; j <= m; j++) {
        bprime[j] = b[j];
        for (i = 1; i <= n; i++) {
            bprime[j] -= a[j][i] * ((double) nl[i]);
        }
    }
        
    print("\n cste = %f \n", cste);
    resout_pb_aux(n, m, cprime, a, bprime, cste, p_borne_inf2, sol_inf2);
        
    /* sol_inf2 est une solution en 01, on la ré-écrit en entier */
    for (i = 1; i <= n; i++) sol_inf2[i] += (double) nl[i];
        
    /* print("\n Admissibilite de la solution du cube : %d \n", est_admissible_x(sol_inf2, n, m, b, a)); */
    if (! est_admissible_x(sol_inf2, n, m, b, a)) {
        for (i = 1; i <= n; i++) sol_inf2[i] = 0.0;
	*p_borne_inf2 = 0.0;
    }

    pqe_common_calcule_temps(p_temps_inf2);
        
    for(i = 1; i < 1 + n; i++) free(y_b3[i]);
    free(y_b3);
    pqe_init_free_1D(nl);
    pqe_init_free_1D(bprime);
    pqe_init_free_1D(cprime);

} /* pqe_main_heur_cube */

/******************************************************************/

void pqe_main_ecrit_carac(FILE *fich, char *entete, int n, int m, int_1D u)
{
    /* on suppose que le fichier fichresu est déja ouvert */

    int nb_var_01 = 0;
    int i;
    double u_moy, u_max, etu, tmp; /* etu = ecart type des ui */

    u_moy = 0.0;
    u_max = u[1];
    for (i = 1; i <= n; i++) {
	u_moy += u[i];
        if (u[i] == 1) nb_var_01++;
	if (u[i] > u_max) u_max = u[i];
    }
    u_moy /= (double) n;

    tmp = 0.0;
    for (i = 1; i <= n; i++) {
        tmp += (u[i] - u_moy) * (u[i] - u_moy);	
    }
    tmp /= (double) n;
    etu = sqrt(tmp);

    fprintf(fich, "%s", entete);
    saut_colonne(fich, 1);
    fprintf(fich, "%d", n);
    saut_colonne(fich, 1);
    fprintf(fich, "%d", m);
    saut_colonne(fich, 1);
    fprintf(fich, "%d", nb_var_01);
    saut_colonne(fich, 1);
    fprintf(fich, "%.1f", u_moy);
    saut_colonne(fich, 1);
    fprintf(fich, "%.0f", u_max);
    saut_colonne(fich, 1);
    fprintf(fich, "%.1f", etu);
    saut_colonne(fich, 1);

} /* pqe_main_ecrit_carac */

/******************************************************************/

void generation_nonsep(int n, int m, double_2D q, double_1D c, double_1D b,
                       double_2D a, int_1D u)
/* on suppose n et m connues, et les differents vecteurs déjà alloués */
{
    int i, j, de;
    double somme_qij, somme_aji_ui;
    FILE *pourvoir;

                 for (i = 1; i <= n; i++) {
                    c[i] = 4.0 + (double) pqe_common_irand((int) (2*n*100)-3);
                    for (j = i+1; j <= n; j++) {
                        
                         q[i][j] = (double) pqe_common_irand(101);
                         de = pqe_common_irand(3);
                         if (de==0) q[i][j] = -q[i][j]; 
                         q[j][i] = q[i][j];
                    } 
                    u[i] = 1 + pqe_common_irand(10);
                    
                    /*ecriture d'instances avec u_i = c_i/2d_i exactement*/
                    /* u[i] = 1.0 +  pqe_common_irand(c[i]/(2.0 * d[i]));*/
                    if (u[i] == 0) printf("\n U%d nul !!!!! \n",i);
                }
                 for (i = 1; i <= n; i++) {
                     somme_qij = 0.0;
                     for (j = 1; j <= n; j++)
                         if (j != i) somme_qij += fabs(q[i][j]);
                     q[i][i] = 1.0 + somme_qij + (double) pqe_common_irand((int) somme_qij);
                 }
                u[0] = 0;
                /*suite de la generation des problemes n=m et m=5%n*/			
                
                for (j = 1; j <= m; j++) {
                    somme_aji_ui = 0.0;
                    for (i = 1; i <= n; i++) {
                        a[j][i] = 1.0 + pqe_common_irand(100);
                        somme_aji_ui += a[j][i] * u[i];
                    }
                    b[j] = 50.0 + pqe_common_irand(somme_aji_ui - 50);
                }
    pourvoir = fopen("pourvoir.lp", "w");
    fprintf(pourvoir, "maximize\n obj: ");
    for (i = 1; i <= n; i++) {if (c[i] >= 0.0) fprintf(pourvoir, "+ "); fprintf(pourvoir, " %.2f x%d\n", c[i], i);}
    fprintf(pourvoir, " + [ ");

    for (i = 1; i <= n; i++) {
        fprintf(pourvoir, "\n");
       fprintf(pourvoir, "- %.2f x%d ^2 ", 2.0* q[i][i], i);
       for (j = i+1; j<= n; j++) {
           if (q[i][j] <= 0) fprintf(pourvoir, " + "); else fprintf(pourvoir, " - ");
           fprintf(pourvoir, " %.2f x%d * x%d ", fabs(4.0 * q[i][j]), i, j);
       }
     }


    fprintf(pourvoir, "] / 2\n");
    fprintf(pourvoir, "\nSubject to\n");
    for (j = 1; j <= m; j++){
        for (i = 1; i<=n; i++) fprintf(pourvoir, " + %.2f x%d ", a[j][i], i);
        fprintf(pourvoir, " <= %.2f \n", b[j]);
    }
    fprintf(pourvoir, "\nbounds\n");
    for (i=1; i <=n; i++)
        fprintf(pourvoir, "0 <= x%d <= %d\n", i, u[i]);
    fprintf(pourvoir, "\ngenerals\n");
    for (i=1; i <=n; i++) fprintf(pourvoir, " x%d\n", i);
    fprintf(pourvoir, "\nend\n");
    fclose(pourvoir);
                
} /* generation_nonsep */

/******************************************************************/

void jeux_essai(int n, int m, int germe_inf, int germe_sup)
{
    /* on ne peut lancer dans le meme appel de jeux_essai notre BB (variable BB_borne2 à vrai)
     * et l'optimum par linea MKP */

    int       nbvar, nblignes, nbcols, cpt, statut, germe, nprime, mprime, nsauv;
    int       i, j, k, l, nbiter, plus_grand, var_frac_racine, umax;
    int       nb_suppr;
    int_1D    u, uprime, correspondance, solution, pr_indices, copycstat, copyrstat, getcstat, getrstat;
    int_1D    ind_cont_suppr;
    double    nb_noeuds, nb_noeuds_sp, nb_noeuds_linea, nb_noeuds_linea_sp, somme_yik, bmin;
    /*somme_yik pour ecrire en x_i la solution en y fournie par borne3*/ 
    double    cmin,somme_ponderee, somme_aji_ui, borne_inf;
    double    somme_sur_j_aji, borne1, borne2, temps1, temps2;
    double    valeur_continue, tps_continu, valeur_opt, tps_opt, tps_opt_sp;
    double    borne3, temps3, opt_linea, tps_opt_linea, tps_opt_linea_sp, tps_prim;
    double    borne_inf2, tps_inf2, constante, val_frac_racine, temps_notre_bb;
    double    temps_bb_b1, nb_noeuds_b1, infD, tempsD, temps_fix01, temps_paq0;
    double    temps_procedure_P, demi_somme;
    double_1D b, c, d, w0, w_etoile, vec_dul, A, solution_y, lb, ub; /* vec_x;*/
    double_1D bprime, cprime, dprime;
    double_1D primalx, primaly, solution_y_borne3, sol_x_b3; 
              /* sol_x_b3 : sol en x que l on va calculer a partir de y_b3*/
    double_1D sol_inf2, x_opt, sol_x_b1, xD;
    /*
    double_1D u_cube; 
               * nouvelle borne pour les x_i declares en double car floor fournit 
	       * des doubles mais ce sont des int */
    double_2D a, q, aprime, s, sprime, y_cont; 
    float     borne_01;
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont, new_mat;
    char_1D   sense, ctype;
    char      rep;
    bool      reso_ent, temps_limite, temps_limite_sp, temps_limite_c, trouve; 
    bool      temps_limite_linea, temps_limite_linea_sp, sol_adm;
    bool      notre_temps_limite_atteint, temps_limite_atteint_b1;
    bool      optimum_trouve     = FAUX;
    bool      lecture_ds_fichier = FAUX;
    bool      correlation        = FAUX; /*MODIF pour cube unite sur instances correlees*/
    bool      b_procedure_P      = VRAI;
    bool      b_fixations_01     = FAUX;
    bool      b_fixations_paquets = VRAI;
    bool      b_borne1           = FAUX;
    bool      b_borne_cont       = VRAI;
    bool      BB_borne1          = FAUX;
    bool      b_borne2           = VRAI;
    bool      BB_borne2          = VRAI;
    bool      BB_linea_ss_prima  = FAUX;
    bool      BB_linea_avec_prima = FAUX;
    bool      BB_quad_avec_prima = FAUX;
    bool      BB_quad_ss_prima   = FAUX;
    char      nom_fichier[100], nom_fich_lu[30];
    char      nom_fichier_excel[100];
    char      ma_commande[1000];
    st_minilib2 *minilib2;
    st_fix    *varfix;
    FILE      *fichresu, *fich_lu;
    
    char rep2[4];


    /*    varie_n[0] = 5;
    varie_n[1] = 10;
    varie_n[2] = 20;
    varie_n[3] = 50;
    varie_n[4] = 70;
    varie_n[5] = 100;
    varie_n[6] = 200;
    varie_n[7] = 500;
    varie_n[8] = 1000;
    varie_n[9] = 1500;
    varie_n[10] = 2000;
    */
    /* Ecriture des resultats dans fichiers Pour instances n=m et m=5%n*/

    sprintf(nom_fichier_excel, "resu_%d_%d_%d_%d.xls", n, m, germe_inf, germe_sup);
    fichresu = fopen(nom_fichier_excel,"w");
    
    
    /* Ecriture des fichiers resu pour instances corrélées: resu_core_n_m_germeinf_germesup.xls */
    /*
    sprintf(nom_fichier_excel, "resu_core_%d_%d_%d_%d.xls", n, m, germe_inf, germe_sup);
    fichresu = fopen(nom_fichier_excel, "w");
    */  
    /*    for (l = 6; l <= 6; l++) {
    n = varie_n[l];
    m = 5;
    */

    pqe_init_2D((void ***) &a,        sizeof(double), m, n);
    pqe_init_2D((void ***) &q,        sizeof(double), n, n);
    pqe_init_1D((void **)  &b,        sizeof(double), m);
    pqe_init_1D((void **)  &c,        sizeof(double), n);
    pqe_init_1D((void **)  &d,        sizeof(double), n);
    pqe_init_1D((void **)  &u,        sizeof(int),    n);
    pqe_init_1D((void **)  &vec_dul,  sizeof(double), m);
    pqe_init_1D((void **)  &w0,       sizeof(double), m);
    pqe_init_1D((void **)  &w_etoile, sizeof(double), m);
    pqe_init_1D((void **)  &sol_x_b3, sizeof(double), n);
    /* pqe_init_1D((void **)  &u_cube,   sizeof(double), n); */
    pqe_init_1D((void **)  &sol_inf2, sizeof(double), n);
    pqe_init_1D((void **)  &x_opt,    sizeof(double), n);
    pqe_init_1D((void **)  &xD,       sizeof(double), n);
    pqe_init_1D((void **)  &sol_x_b1, sizeof(double), n);
    pqe_init_1D ((void **)  &minilib2,     sizeof(st_minilib2), n);
    pqe_init_1D ((void **)  &cprime,       sizeof(double), n);
    pqe_init_1D ((void **)  &dprime,       sizeof(double), n);
    pqe_init_1D ((void **)  &uprime,       sizeof(int),    n);
    pqe_init_1D ((void **)  &bprime,       sizeof(double), m);
    pqe_init_2D ((void ***) &aprime,       sizeof(double), m, n);
    pqe_init_1D ((void **)  &correspondance,    sizeof(int),    n);
    pqe_init_1D ((void **)  &varfix,            sizeof(st_fix), n);
    pqe_init_1D((void **)  &ind_cont_suppr, sizeof(int), m);
				
    //generation_nonsep(n, m, q, c, b, a, u);

   // exit(1);

    for (germe = germe_inf; germe<= germe_sup; germe += 10) {
        /*optimum_trouve = FAUX Ne faut-il pas le placer ici ????????*/
        
        fprintf(fichresu, "%d\t%d\t%d\t", n, m, germe);
        
        
      /* ancien lecture fichier  if (lecture_ds_fichier) {
            printf(" Nom du fichier : ");
            scanf("%s", nom_fich_lu);
            fich_lu = fopen(nom_fich_lu, "r");
            printf("\n fichier ouvert\n");
            fscanf(fich_lu, "%d", &n);
            printf("n = %d \n", n);
            fscanf(fich_lu, "%d", &m);
            printf("m = %d \n", m);
            fscanf(fich_lu, "%lf", &constante);
            printf("n = %d, m = %d, constante = %lf", n, m, constante);
            for (i = 1; i <= n; i++) fscanf(fich_lu, "%lf", &c[i]);
            for (i = 1; i <= n; i++) fscanf(fich_lu, "%lf", &d[i]);
            for (i = 1; i <= n; i++) fscanf(fich_lu, "%d", &u[i]);
            for (j = 1; j <= m; j++) fscanf(fich_lu, "%lf", &b[j]);
            for (j = 1; j <= m; j++) 
            for (i = 1; i <= n; i++) 
            fscanf(fich_lu, "%lf", &a[j][i]);
            for (j = 1; j <= m; j++) w0[j] = 1.0;
            fclose(fich_lu);
        }*/
/* nouveau lecture fichier du 17/11/2008*/
	if (lecture_ds_fichier){
            printf(" Nom du fichier : ");
            scanf("%s", nom_fich_lu);
            fich_lu = fopen(nom_fich_lu, "r");
            printf("\n fichier ouvert\n");
            fscanf(fich_lu, "%d", &n);
            printf("n = %d \n", n);
            fscanf(fich_lu, "%d", &m);
            printf("m = %d \n", m);
            for (i = 1; i <= n; i++) fscanf(fich_lu, "%lf", &c[i]);
            for (j = 1; j <= m; j++) 
            for (i = 1; i <= n; i++) 
            fscanf(fich_lu, "%lf", &a[j][i]);
	    for (j = 1; j <= m; j++) fscanf(fich_lu, "%lf", &b[j]); 
	    for (j = 1; j <= m; j++) w0[j] = 1.0;
	    for (i = 1; i <= n; i++) {
		d[i] = 1.0 + (double) pqe_common_irand(( floor(c[i]/2.0) )-1) ;
                u[i] = 1 + pqe_common_irand(floor(c[i]/(2.0 * d[i])));
                if (u[i] == 0) printf("\n U%d nul !!!!! \n",i);
            }	
            fclose(fich_lu);


	}/*fin lecture fichier*/
        else {
            /* ---------------------- */
            /* Génération du problème */
            /* ---------------------- */
            
            pqe_common_init_rand(germe);
            
            if (correlation) {
                
                /* --------------------------------------------------------------------------*/
                /* CORRELATION DES COEFF LINEAIRES DE LA FONTION OBJECTIF ET DES CONTRAINTES*/
                /*		         	i e Somme_sur_j_aji = c_i                   */
                /*--------------------------------------------------------------------------*/
                /* 1. Génération des var a_ji, b_j et c_i w_0 choisi arbitrairement*/
                
                
                for (j = 1; j <= m; j++) {
                    for (i = 1; i <= n; i++) {
                        a[j][i] = 1.0 + pqe_common_irand(100);
                    }
                }
		printf("\n Correlation : d = ");
                for (i=1; i<= n; i++) {
                    
                    somme_sur_j_aji= demi_somme = 0.0;
                    for (j=1; j<= m; j++) {
                        somme_sur_j_aji += a[j][i];
			if (j%((int) (n/2)) == 0) demi_somme += a[j][i];
                    }
                    c[i] = somme_sur_j_aji;  
                    d[i] = 20 * floor(demi_somme);
		    printf(" %.0f ", d[i]);
                }
                
                /* 2. Génération des variables d_i et u_i qui dépendent de c_i, on pourra facilement rajouter d_i = cste ici*/
                
                /*
                for (i = 1; i <= n; i++) {
                d[i] = 1.0 + (double) pqe_common_irand(( floor(c[i]/2.0) )-1) ; 
                u[i] = 1 + pqe_common_irand(floor(c[i]/(2.0 * d[i])));
                if (u[i] == 0) printf("\n U%d nul !!!!! \n",i);
            }			
                */
                /* Cas di = cste = cmin/2; */
                cmin = c[1];
                for (i = 2; i <= n; i++) {
                    if (c[i] < cmin) cmin = c[i];
                }
		printf("\nCorrelation : u = ");
                for (i= 1; i <= n; i++) { 
                    /* d[i] = floor(cmin/2.0); */
                    /* d[i] = floor(cmin/10.0); */
                    /* d[i] = floor(c[i]/20.0); */
                    /* d[i] = floor(.25*c[i]); */
		    if (floor(c[i]/(2.0 * d[i])) < 1.0) u[i] = 1;
		    else u[i] = 1 + pqe_common_irand(floor(c[i]/(2.0 * d[i])));
                    if (u[i] == 0) printf("\n u%d nul !!!!! \n",i);
		    printf(" %.0f ", u[i]);
                }
		/* scanf("%s", rep2);*/
		printf("\n");
		printf("\n Correlation : c = ");
                for (i= 1; i <= n; i++) { 
		    printf(" %.0f ", c[i]);
                }
		printf("\n");
		printf("\n Correlation : u = ");
                for (i= 1; i <= n; i++) { 
		    printf(" %d ", u[i]);
                }
		printf("\n");
                u[0] = 0;
                for (j = 1; j <= m; j++) {
                    somme_aji_ui = 0.0;
                    for (i = 1; i <= n; i++) {
                        somme_aji_ui += a[j][i] * u[i];
                    }
                    b[j] = 50.0 + pqe_common_irand(somme_aji_ui - 50);
                    w0[j] = 1.0;
                }
                constante = 0.0;
		/*****************************************************************/
		/* comptage des variables 0-1                                    */
		/*****************************************************************/
		/* voir avec eric
                for(i=1; i <=n; i++){
		   }
		*/
            }
            else {    
                /* Problemes n=m et m=5%n */  
                
                 for (i = 1; i <= n; i++) {
                    c[i] = 4.0 + (double) pqe_common_irand(97);
                    d[i] = 1.0 + (double) pqe_common_irand(( floor(c[i]/2.0) )-1) ;
                 
                    u[i] = 1 + pqe_common_irand(floor(c[i]/(2.0 * d[i])));
                    
                    /*ecriture d'instances avec u_i = c_i/2d_i exactement*/
                    /* u[i] = 1.0 +  pqe_common_irand(c[i]/(2.0 * d[i]));*/
                    if (u[i] == 0) printf("\n U%d nul !!!!! \n",i);
                }
                u[0] = 0;
                /*suite de la generation des problemes n=m et m=5%n*/			
                
                for (j = 1; j <= m; j++) {
                    somme_aji_ui = 0.0;
                    for (i = 1; i <= n; i++) {
                        a[j][i] = 1.0 + pqe_common_irand(100);
                        somme_aji_ui += a[j][i] * u[i];
                    }
                    b[j] = 50.0 + pqe_common_irand(somme_aji_ui - 50);
                    w0[j] = 1.0;
                }
                
                constante = 0.0;
            }


        }
        if (probleme_trivial (n, m, a, b, u)) {
            printf("\n Problème trivial : contraintes inutiles : b trop grand");
            printf("\n Continuer (o/n) ? : ");
            scanf("%c", &rep);
            if (rep != 'o') exit(1);
        }
        
	pqe_main_ecrit_carac(fichresu, "Avant_pretraitement", n, m, u);

        pqe_cplex_init();

        if (BB_borne2) pqe_bb_eliminer_contraintes(n, m, &mprime, a, b, &nb_suppr, ind_cont_suppr);
	else mprime = m;
	/* dorenavant, le nombre de contraintes est mpriime */

	if (BB_borne2) {
            pqe_bb_init_minilib2(n, u, minilib2);
	    if (b_procedure_P) pqe_bb_procedure_P(n, mprime, a, b, u, minilib2, &temps_procedure_P);
	}

        
        /* -------------------------------------------------------- */
        /* Calcul de la relaxation continue du pb initial par cplex */
        /* -------------------------------------------------------- */
        
        if (b_borne_cont) {
            reso_ent = FAUX;
            sprintf(nom_fichier, "pqc_%d_%d_%d.mps", n, m,&germe);
            ecrire_fichier_mps(n, m, c, d, a, b, u, reso_ent, nom_fichier);
            pqe_cplex_reso_pq_fichier(&valeur_continue, &tps_continu, nom_fichier, 
            reso_ent, NULL, &temps_limite_c, FAUX, n, NULL, NULL);
            sprintf(ma_commande, "rm -f %s", nom_fichier);
            system(ma_commande); 
            ecrit_fichier_borne(fichresu, valeur_continue, tps_continu, "Rel_cont");
        }
        
        /* ------------- */
        /* Linearisation */
        /* ------------- */
        
        s = (double **) malloc(sizeof(double *) * (1+n));
        for (i = 0; i < 1 + n; i++) s[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
        
	umax = u[1];
        for (i = 2; i <=  n; i++) 
		if (umax < u[i]) umax = u[i];
        pqe_init_2D ((void ***) &sprime,       sizeof(double), n, umax);
        /* sprime = (double **) malloc(sizeof(double *) * (1+n));
        for (i = 0; i < 1 + n; i++) sprime[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
	nsauv = n; print("u[0] = %d\n", u[0]);
	*/
        
        y_cont = (double **) malloc(sizeof(double *) * (1+n));
        for (i = 0; i < 1 + n; i++) y_cont[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
        
        pqe_main_linea_fction_quadra(n, u, c, d, s);
        
        
        /* --------------------------------------------------------- */
        /* Resolution du pb linearise non agrege en continu : borne3 */ 
        /* --------------------------------------------------------- */
        
        pqe_init_1D((void **)   &solution_y_borne3, sizeof(double), somme_ui(n,u));
        /* 
        pqe_init_1D((void **)   &copycstat, sizeof(int), somme_ui(n,u));
        pqe_init_1D((void **)   &copyrstat, sizeof(int), m);
        pqe_init_1D((void **)   &getcstat, sizeof(int), somme_ui(n,u));
        pqe_init_1D((void **)   &getrstat, sizeof(int), m);
        */
        copycstat = (int *) malloc(sizeof(int) * (1 +  somme_ui(n,u)));
        copyrstat = (int *) malloc(sizeof(int) * (1 +  mprime));
        getcstat =  (int *) malloc(sizeof(int) * (1 +  somme_ui(n,u)));
        getrstat =  (int *) malloc(sizeof(int) * (1 +  mprime));
        
        pqe_main_calcul_borne3(n, mprime, u, s, a, b, &borne3, &temps3, vec_dul, solution_y_borne3,
        NULL, NULL, getcstat, getrstat);
        for (i = somme_ui(n, u); i >= 1; i--) {
            print("\n solution_y_borne3[%d] = %lf\n", i, solution_y_borne3[i]);
        }
        
        for (i = 0; i < somme_ui(n, u); i++) copycstat[i] = getcstat[i];
        for (j = 0; j < mprime; j++) copyrstat[j] = getrstat[j];    
        
        ecrit_fichier_borne(fichresu, borne3, temps3, "Borne3");
        
        /*------------------------------------------------------------*/
        /*   Le cube UNITE                                            */
        /*------------------------------------------------------------*/

        pqe_main_heur_cube(solution_y_borne3, n, mprime, a, b, c, d, u, &borne_inf2,
                           sol_inf2, sol_x_b3, &tps_inf2);
	
        ecrit_fichier_borne(fichresu, borne_inf2, tps_inf2, "Borne_inf2");

        /* ------------------------------------------- */
        /* Resolution du pb linearise agrege : borne 1 */ 
        /* ------------------------------------------- */
        
        if (b_borne1) {
            pqe_main_calcule_w_etoile(n, mprime, w0, s, a, b, u, w_etoile, &borne1, &nbiter, &temps1, y_cont);
        
            ecrit_fichier_borne(fichresu, borne1, temps1, "Borne1");
        
            heuri_D(n, mprime, a, b, c, d, u, w_etoile, xD, &infD, &tempsD);
        
            ecrit_fichier_borne(fichresu, infD, tempsD, "Heur_D");
	}
        /* 
        *
        print("\n TEST D'ADMISSIBILITE de y : %d\n", pqe_main_est_admissible_y(y_cont, n, mprime, b, a, u));
        *
        * for (j = 1; j <= mprime; j++) 
        printf("w_etoile[%d] (%.6f) / vec_dul[%d] (%.6f) = %.6f \n", 
        j, w_etoile[j], j, vec_dul[j-1], w_etoile[j]/vec_dul[j-1]); 
        */
        
        
        /* ----------------------------------------------------------- */
        /* Amélioration de la borne par appel à cplex en 0-1 : borne 2 */
        /* ----------------------------------------------------------- */
        
        pqe_init_1D((void **)   &solution_y, sizeof(double), somme_ui(n,u));
        if (b_borne2) {
            pqe_main_calcul_borne2(n, mprime, u, s, a, b, vec_dul, &borne2, &temps2, solution_y); 
            ecrit_fichier_borne(fichresu, borne2, temps3 + temps2, "Borne2");
        }
        
        
        
        /* pqe_main_calcul_borne2(n, mprime, u, s, a, b, w_etoile, &borne2, &temps2, solution_y);
        for (i=1; i<= somme_ui(n,u);i++) print("\r solution_y[%d] = %lf\n", i, solution_y[i]);
        fprintf(fichresu, "%d\t%d\t%d\t%lf\t%lf\t%d\t\t%lf\t%lf\t\t\%lf\t%lf\t\t%lf\t%lf\t%.0f\t%lf\t%lf",
        n, mprime, germe, borne1, temps1, nbiter, borne2, temps1+temps2,
        valeur_continue, tps_continu, valeur_opt, tps_opt, nb_noeuds, borne3, temps3);
        if (temps_limite) fprintf(fichresu, "\tTEMPS LIMITE\n");
        else fprintf(fichresu, "\n");*/
        
        
        /* ------------------------------------------------ */
        /* Primalisation de la solution fournie par borne 2 */
        /* Calcul d'une solution admissible                 */
        /* ------------------------------------------------ */
        
        pqe_init_1D((void **) &primalx, sizeof(double), n);
        pqe_init_1D((void **) &primaly, sizeof(double), somme_ui(n,u));
        
        
        if (b_borne2) {
            pqe_main_primalisation (n, mprime, u, s, a, b, c, d, constante, solution_y, &borne_inf, 
            &optimum_trouve, primalx, primaly, borne2, &tps_prim);
        }
        
        /* à mon avis c'est là qu'il y a un problème car ce n'est pas borne2 qui est égale à Borne_inf mais B_inf = val_opt*/
        /*  if ((optimum_trouve) || (borne2 == borne_inf)) {
        if (optimum_trouve) printf("\n Le calcul de borne 2 a permis de trouver directement l'optimum\n");
        else printf("\n Le calcul de borne2 puis la pqe_main_primalisation ont permis de trouver l'optimum\n");*/
        if (optimum_trouve) {
            printf("\n Le calcul de borne 2 a permis de trouver directement l'optimum\n");
            valeur_opt = borne2;
            tps_opt = tps_opt_sp = 0.0;
            nb_noeuds = nb_noeuds_sp = 0.0;
            
	    if (BB_borne2)
	        ecrit_fichier_opt(fichresu, FAUX, valeur_opt, tps_opt, nb_noeuds, 
			       "opt_notre_bb=borne2", FAUX, 0.0);
            if (BB_quad_avec_prima) 
	        ecrit_fichier_opt(fichresu, FAUX, valeur_opt, tps_opt, nb_noeuds, 
			       "opt_q_avec_prim=borne2", FAUX, 0.0);
            if (BB_quad_ss_prima)
                ecrit_fichier_opt(fichresu, FAUX, valeur_opt, tps_opt, nb_noeuds, 
			       "opt_q_sans_prim=borne2", FAUX, 0.0);
            if (BB_linea_avec_prima) 
	        ecrit_fichier_opt(fichresu, FAUX, valeur_opt, tps_opt, nb_noeuds, 
			       "opt_lin_avec_prim=borne2", FAUX, 0.0);
            if (BB_linea_ss_prima) 
		ecrit_fichier_opt(fichresu, FAUX, valeur_opt, tps_opt, 
			    nb_noeuds, "opt_lin_sans_prim=borne2", FAUX, 0.0);
            /*optimum_trouve = VRAI;*/
        }
        else { 
            /*---------------------------------------------------------------------------------*/
            /*---------------------- Notre beau B&B -------------------------------------------*/
            /*---------------------------------------------------------------------------------*/
            
            /* for(i=1;i<=n;i++) print("\r sol_x_b3[%d]= %lf\n", i, sol_x_b3[i]);*/
            /*mise ne commentaire de notre bb pour faire tourner djerdjour*/
            
            if (BB_borne2) {
		
		if (b_fixations_paquets){ 
                    pqe_bb_fixations_paquets_0 (n, mprime, germe, a, b, c, d,
				            u, s, sol_inf2, borne_inf2, vec_dul, &temps_paq0, minilib2);
		    printf (" temps fixations par paquets à 0 : %.1f s\n", temps_paq0); 
		}
		
		if (b_fixations_01) {
		    pqe_bb_fixations_01(n, mprime, germe, a, b, c, d, u, sol_inf2, borne_inf2, &temps_fix01, minilib2);

		    printf (" temps fixations variable par variable : %.1f s\n", temps_fix01); 
		}
		
		
		pqe_bb_reecriture(n, mprime, germe, a, b, c, d, u, minilib2, &nprime, aprime, bprime,  cprime, dprime, 
				  uprime, &constante, correspondance, varfix);
		printf("\n valeur de la constante apres fixations : %f", constante);
        pqe_main_linea_fction_quadra(nprime, uprime, cprime, dprime, sprime);
        pqe_main_calcul_borne3(nprime, mprime, uprime, sprime, aprime, bprime, &borne3, &temps3, vec_dul, solution_y_borne3,
        NULL, NULL, getcstat, getrstat);
        pqe_main_heur_cube(solution_y_borne3, nprime, mprime, aprime, bprime, cprime, dprime, uprime, &borne_inf2,
                           sol_inf2, sol_x_b3, &tps_inf2);
            pqe_main_calcul_borne2(nprime, mprime, uprime, sprime, aprime, bprime, vec_dul, &borne2, &temps2, solution_y); 
			 

                pqe_bb_trouver_var_frac(sol_x_b3, nprime, &var_frac_racine, &val_frac_racine); 
	        pqe_main_ecrit_carac(fichresu, "Apres_pretraitement", nprime, mprime, uprime);
                printf("\n Lancement du Branch-and-Bound : n = %d, m = %d\n", nprime, mprime);
                pqe_bb_branch_and_bound(nprime, mprime, germe, cprime, dprime, bprime, aprime, uprime, borne_inf2, sol_inf2, 
                              borne2, var_frac_racine, val_frac_racine, &optimum_trouve, 
                              &valeur_opt, x_opt, &temps_notre_bb, copycstat, copyrstat, 
                              vec_dul, &nb_noeuds, T_LIMITE, &notre_temps_limite_atteint);
		/*
		printf("\n xopt");
		for (i = 1; i <= nprime; i++) printf(" [%d]=%.0f", i, x_opt[i]);
		
		bmin = b[1];
		printf("\n b");
		for (j = 1; j <= mprime; j++) {printf(" [%d]=%.0f", j, b[j]); if (b[j] < bmin) bmin = b[j];}
		printf("\n bmin = %f", bmin);
		*/
			  /*suite au prob 1000 1000 10 10 pour lequel opt = 233 (car BI passe a 233) 
			  alors que opt est egale a 229 cf Clpex et Djerdjour*/	
              print("est admissible = %d\n", est_admissible_x(x_opt, n, mprime, b, a));
            
                ecrit_fichier_opt(fichresu, notre_temps_limite_atteint, valeur_opt+constante, 
                      temps_notre_bb+temps_paq0+temps_fix01+temps_procedure_P, 
                      nb_noeuds, "opt_notre_bb", VRAI, temps_procedure_P + temps_paq0+temps_fix01);
            
                /* attention valeur_opt est utilisee dans notre B&B et dans le pqe_reso_fichier de cplex !!!!!*/
                print("\n OUF ! on est sorti du b&b\n");
                print("optimum = %lf et temps_opt de notre B&B = %lf s\n", valeur_opt, temps_notre_bb);
            } 
            
            /*---------------------------------------------------------------------------------*/
            /*---------------------- B&B avec BORNE1-------------------------------------------*/
            /*---------------------------------------------------------------------------------*/
           
	    if (BB_borne1) { 
                /* remarque: pas obligatoire ici de recalculer var frac a partir de bonre1*/
	    
                for(i = 1; i <= n; i++) {
                    somme_yik = 0.0;	
                    for(k = 1; k <= u[i]; k++){
                        somme_yik = somme_yik + y_cont[i][k];
                    }
                    sol_x_b1[i] = somme_yik;
                }
                pqe_bb_trouver_var_frac(sol_x_b1, n, &var_frac_racine, &val_frac_racine); 
            
                /*remarque: il faut retirer copystat.. mais attention il faut changer la def des na ou mettre null a la place*/
                /*retrait de multiplicateur egalement*/
            
            
            
                pqe_bb_branch_and_bound_borne1(n, mprime, germe, c, d, b, a, u, borne_inf2, sol_inf2, 
                     borne1, var_frac_racine, val_frac_racine,
                     &optimum_trouve, &valeur_opt,  x_opt, &temps_bb_b1,
                     /* copycstat,  copyrstat, w_etoile, &nb_noeuds_b1, T_LIMITE, */
                     copycstat,  copyrstat, w0, &nb_noeuds_b1, T_LIMITE,
                     &temps_limite_atteint_b1);
            
                ecrit_fichier_opt(fichresu, temps_limite_atteint_b1, valeur_opt, temps_bb_b1, 
                                  nb_noeuds_b1, "opt_bb_b1", FAUX, 0.0);
            
            } 
	    
            
            pqe_init_1D((void **) &pr_indices, sizeof(int), n);
            for (i = 1; i <= n; i++)  pr_indices[i - 1] = i - 1; 
            
            /* ---------------------------------------- */
            /* Calcul de la solution optimale par cplex */
            /* ---------------------------------------- */
            
            reso_ent = VRAI;
	    if (BB_quad_ss_prima || BB_quad_avec_prima) {
                sprintf(nom_fichier, "pqe_%d_%d_%d.mps", n, m, germe);
                /*MODIF du 01/10/05 pour borne_in_cube_unite*/
                ecrire_fichier_mps(n, m, c, d, a, b, u, reso_ent, nom_fichier);
	    }
            
            
            temps_limite_sp = FAUX;
            sol_adm = FAUX;
            
	    if (BB_quad_ss_prima) {
                pqe_cplex_reso_pq_fichier(&valeur_opt, &tps_opt_sp, nom_fichier, reso_ent, 
                         &nb_noeuds_sp, &temps_limite_sp, sol_adm, n, primalx, pr_indices);
            
                ecrit_fichier_opt(fichresu, temps_limite_sp, valeur_opt, tps_opt_sp, 
                    nb_noeuds_sp, "opt_q_sans_prim", FAUX, 0.0);
	    }
            temps_limite = FAUX;
            sol_adm = VRAI;
            
	    if (BB_quad_avec_prima) {
                pqe_cplex_reso_pq_fichier(&valeur_opt, &tps_opt, nom_fichier, reso_ent, 
                    &nb_noeuds, &temps_limite, sol_adm, n, primalx, pr_indices);
                ecrit_fichier_opt(fichresu, temps_limite, valeur_opt, tps_opt, 
                    nb_noeuds, "opt_q_avec_prim", FAUX, 0.0);
	    }
            
	    if (BB_quad_ss_prima || BB_quad_avec_prima) {
                sprintf(ma_commande, "rm -f %s", nom_fichier);
                system(ma_commande);
	    }
            
            temps_limite_linea_sp = FAUX;
            sol_adm = FAUX;
            if (BB_linea_ss_prima){
                calcul_opt_linea(n, mprime, u, s, a, b, &opt_linea, &tps_opt_linea_sp, sol_adm, 
                    primaly, &temps_limite_linea_sp, &nb_noeuds_linea_sp);
                ecrit_fichier_opt(fichresu, temps_limite_linea_sp, opt_linea, tps_opt_linea_sp, 
                    nb_noeuds_linea_sp, "opt_lin_sans_prim", FAUX, 0.0);
	    }
            
            temps_limite_linea = FAUX;
            sol_adm = VRAI;
            if(BB_linea_avec_prima){
                calcul_opt_linea(n, mprime, u, s, a, b, &opt_linea, &tps_opt_linea, sol_adm, 
                    primaly, &temps_limite_linea, &nb_noeuds_linea);
                ecrit_fichier_opt(fichresu, temps_limite_linea, opt_linea, tps_opt_linea, 
                    nb_noeuds_linea, "opt_lin_avec_prim", FAUX, 0.0);
	    }
            
            pqe_init_free_1D(pr_indices);
        }
        
        /* ---------------------------------------- */
        /* Ecriture des résultats dans fichier xls  */
        /* ---------------------------------------- */
        
        print("\nProbleme n = %d m = %d germe = %d, fichier %s\n continu = %.2f, entier = %2.f\n",
            n, mprime, germe, nom_fichier,  valeur_continue, valeur_opt);
        print("\ntemps continu = %2.f s - temps entier = %2.f s, nb_noeuds = %.0f\n\n",
            tps_continu, tps_opt, nb_noeuds);
        /* fprintf(fichresu, "%d\t%d\t%d\t%lf\t%lf\t%d\t\t%lf\t%lf\t\t\%lf\t%lf\t\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\t%d",
        n, mprime, germe, borne1, temps1, nbiter, borne2, temps3+temps2,
        valeur_continue, tps_continu, valeur_opt, tps_opt, nb_noeuds, tps_opt_linea, 
        borne3, temps3, borne_inf, tps_opt_sp, nb_noeuds_sp, tps_opt_linea_sp, optimum_trouve);
        fprintf(fichresu, "%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%d\t\t%lf\t%lf\t\t%d\t%lf\t%lf\t\t%lf\t%lf\t%.0f\t%lf\t\t%lf\t%lf\t%.0f\t%lf",
        n, mprime, germe, borne3, temps3, borne1, temps1, nbiter, borne2, temps3+temps2,optimum_trouve,
        valeur_continue, tps_continu, valeur_opt, tps_opt, nb_noeuds, tps_opt_linea, 
        borne_inf, tps_opt_sp, nb_noeuds_sp, tps_opt_linea_sp);
        
        
        
        if (temps_limite) fprintf(fichresu, "\tTEMPS LIMITE\n");
        if (temps_limite_sp) fprintf(fichresu, "\tTEMPS LIMITE SP\n");
        if (temps_limite_linea) fprintf(fichresu, "\tTEMPS LIMITE LINEA\n");
        if (temps_limite_linea_sp) fprintf(fichresu, "\tTEMPS LIMITE LINEA SP\n");
        */
        
        fprintf(fichresu, "\n");
        
        /* ------------------------------------------- */
        /* Affichage des résulats numériques à l'écran */
        /* ------------------------------------------- */
        
        /* affichage_resu_ecran(valeur_continue, borne1, borne2, borne3, valeur_opt,
        tps_continu, temps1, temps2, temps3, tps_opt);
	*/
        
        
        /* Liberation des tableaux dynamiques */
        /*	printf("\n la1\n");A */
    
        pqe_init_free_2D((void *) sprime,  n); 

    
        /* for (j = 1; j <= nsauv; j++) { free(sprime[j]); printf("top 2.%d  nsauv = %d \n", j, nsauv); }
	free(sprime); */
        pqe_init_free_2D((void *) s,  n); 
        pqe_init_free_1D(solution_y);
        pqe_init_free_1D(primalx);
        pqe_init_free_1D(primaly);
        pqe_init_free_1D(solution_y_borne3);
        /*
        pqe_init_free_2D((void *) y_b3,  n);
        y_b3 = (double **) malloc(sizeof(double *)*(n + 1));
        for(i = 1; i < 1 + n; i++) y_b3[i] = (double *) malloc(sizeof(double)*(1 + u[i]));
        */      
        pqe_cplex_fin();
        /* scanf("%c", &rep); */
    }  /* boucle sur le germe */
    
    pqe_init_free_1D(w0);
    pqe_init_free_1D(w_etoile);
    pqe_init_free_2D((void *) a,  m);
    pqe_init_free_2D((void *) q,  n);
    pqe_init_free_2D((void *) aprime,  m);
    pqe_init_free_1D(b);
    pqe_init_free_1D(bprime);
    pqe_init_free_1D(c);
    pqe_init_free_1D(cprime);
    pqe_init_free_1D(sol_inf2);
    pqe_init_free_1D(x_opt);
    pqe_init_free_1D(xD);
    pqe_init_free_1D(d);
    pqe_init_free_1D(dprime);
    pqe_init_free_1D(u);
    pqe_init_free_1D(uprime);
    pqe_init_free_1D(vec_dul);
    pqe_init_free_1D(sol_x_b3);
    /* pqe_init_free_1D(u_cube); */
    pqe_init_free_1D(sol_x_b1);
    pqe_init_free_1D(minilib2);
    pqe_init_free_1D(correspondance);
    pqe_init_free_1D(varfix);
    pqe_init_free_1D(ind_cont_suppr);
    
    
    /* } boucle sur varie_n */
    fclose(fichresu);
    
} /* jeux_essai */

/******************************************************************/

int main(int argc, char **argv)
{
    printf("\n usage : %s n m germe_inf germe_sup", argv[0]);
    
    jeux_essai(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
      //jeux_essai(250, 10, atoi(argv[3]), atoi(argv[4]));
}    
    /*
    int       n, m, nbvar, nblignes, nbcols, cpt, statut;
    int       i, j, k;
    int_1D    u;
    double    somme_ponderee;
    double_1D b, c, d, w0, w_etoile, A, solution_y;
    double_2D a, s;
    float     borne_01;
    float_1D  coef_obj, vec_rhs;
    float_2D  mat_cont;
    char_1D   sense, nom_fichier;
    char      rep;
    bool      reso_ent;
    
    
    / * affichage = VRAI;  * /
    / * Generation des donnees * /
    
    if (argc>1) {
    print("atoi(argv[1] = %d",atoi(argv[1]));
    switch (atoi(argv[1])) {
    case 1 : genere_exemple_1(&n, &m, &c, &d, &a, &b, &u, &w0); break;
    case 2 : genere_exemple_2(&n, &m, &c, &d, &a, &b, &u, &w0); break;
    case 3 : genere_exemple_3(&n, &m, &c, &d, &a, &b, &u, &w0); break;
    case 4 : genere_exemple_4(&n, &m, &c, &d, &a, &b, &u, &w0); break;
    default : printf("erreur de syntaxe\n");
}
}
    else genere_exemple_1(&n, &m, &c, &d, &a, &b, &u, &w0);
    
    if (probleme_trivial (n, m, a, b, u)) {
    if (probleme_trivial (n, m, a, b, u)) {
    printf("\n Problème trivial : contraintes inutiles : b trop grand");
    printf("\n Continuer (o/n) ? : ");
    scanf("%c", &rep);
    if (rep != 'o') exit(1);
}
    
    reso_ent = VRAI;
    nom_fichier = "pqe.mps";
    ecrire_fichier_mps(n, m, c, d, a, b, u, reso_ent, nom_fichier);
    
    nbcols = nbvar = somme_ui(n, u);
    nblignes = 1;
    
    
    
    pqe_init_1D((void **)  &A, sizeof(double), n);
    pqe_init_1D((void **)  &coef_obj, sizeof(float), nbvar);
    pqe_init_2D((void ***) &mat_cont, sizeof(float), nblignes, nbcols);
    pqe_init_1D((void **)  &vec_rhs, sizeof(float), nblignes);
    pqe_init_1D((void **)  &sense, sizeof(char), nblignes);
    
    pqe_init_1D((void **) &solution_y, sizeof(double), nbcols);
    pqe_init_1D((void **) &w_etoile, sizeof(double), m);
    s = (double **) malloc(sizeof(double *) * (1+n));
    for (i = 0; i < 1 + n; i++) s[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
    
    
    
    pqe_cplex_init();
    
    pqe_main_linea_fction_quadra(n, u, c, d, s );
    
    pqe_main_calcule_w_etoile(n, m, w0, s, a, b, u, w_etoile);
    
    / * ------------------------------------------------- * /
    / * Amélioration de la borne par appel à cplex en 0-1 * /
    / * ------------------------------------------------- * /
    
    / * Coef obj * /
    
    cpt = 0 ;
    for (i = 1; i<= n; i++)
    for (k = 1 ; k <= u[i]; k++) {
    coef_obj[cpt] = (float) s[i][k];
    cpt++;
}
    
    / * coefficients de la contrainte aggregée * /
    
    / * Aggregation des contraintes * /
    for (i = 1; i <= n; i++) {
    A[i] = 0.0;
    for (j = 1; j <= m; j++) A[i] += w_etoile[j] * a[j][i];
}
    
    / * contrainte aggregée avec les variable yik * /
    
    cpt = 0;
    for (i=1; i <= n; i++) 
    for (j = 1 ; j <= u[i]; j++) 
    mat_cont[0][cpt++] = A[i];
    
    vec_rhs[0]    = 0.0;
    for (j = 1; j <= m; j++) vec_rhs[0] += w_etoile[j] * b[j];
    
    sense[0] = 'L';
    
    borne_01 = pqe_cplex_reso_lp01(nblignes, nbcols, coef_obj,
    mat_cont, vec_rhs,
    &statut, solution_y, sense,
    VRAI, NULL, FAUX, NULL, &bidon) ;
    
    printf("\n Valeur de la borne en 01 : %f\n", borne_01);
    scanf("%c", &rep);
    
    
    
    / *
    A[1] = 3;
    A[2] = 11;
    
    B = 26
    p = (int *) malloc(sizeof(int) * (1+n));
    y = (double **) malloc(sizeof(double *) * (1+n));
    for (i = 0; i < 1 + n; i++) y[i] = (double *) malloc(sizeof(double) * (1 + u[i]));
    print("\n top\n");
    pqe_main_reso_kp_y(n, u, A, B, s, y, &ietoile, &pietoile, p);
    
    for (i = 1; i <= n; i++) {
    for (k = 1; k <= u[i]; k++) print("y[%d][%d] = %f\n",i, k, y[i][k]);
    print("p[%d] = %d \n", i, p[i]);
}
    
    print("ietoile = %d, pietoile = %d\n", ietoile, pietoile);
    
    scanf("%c", &c);
    
    
    ancien test :
    
    st_tri_kp *tri;
    int nbvar = 10, i , ind_critique;
    double p[11], w[11], solution[11], cap = 26.0, valopt;
    FILE   *fich;
    
    
    p[1] = 3.0; p[2] = -1.0; p[3] = -5.0; p[4] = -9.0; p[5] = -13.0; 
    p[6] = 11.0; p[7] = 9.0; p[8] = 7.0; p[9] = 5.0; p[10] = 3.0;
    
    w[1] = 3.0; w[2] = 3.0; w[3] = 3.0; w[4] = 3.0; w[5] = 3.0;
    w[6] = 11.0; w[7] = 11.0; w[8] = 11.0; w[9] = 11.0; w[10] = 11.0;
    
    tri = (st_tri_kp *) malloc ((nbvar+1) * sizeof(st_tri_kp));
    
    for (i = 1; i <= nbvar; i++) 
    { 
    tri[i].cle = p[i]/w[i];  
    tri[i].ind = i; 
    solution[i] = 0;
}  
    
    valopt = 0.0;
    reso_kp (nbvar, p, w, cap, 1, nbvar, &valopt, solution, &ind_critique, tri);
    
    fich = fopen("w","resukp.txt");
    print("\n valopt = %lf \n", valopt);
    print("\n ind critique : %d\n", ind_critique);
    
    for (i=1 ; i <= nbvar; i++) print("x[%d] = %f \n", i, (float) solution[i]);
    
    fclose(fich);
    free(tri);
} 
    */



