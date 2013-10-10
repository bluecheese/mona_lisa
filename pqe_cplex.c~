/***************************************************************************/
/*                                                                         */
/*                    PROGRAMMATION         QUADRATIQUE                    */
/*                   Methode de Decomposition Lagrangienne                 */
/*                                                                         */
/*  Fichier : knap_cplex.c                                                 */
/*                                                                         */
/*  Eric SOUTIF (e-mail : soutif@cnam.fr)                                  */
/*  Dominique QUADRI (e-mail : quadri@lri.fr)                              */
/*                                                                         */
/***************************************************************************/



#include <stdio.h>
#include <string.h>
#include <math.h>
/* #include "C:\ILOG\Cplex90\include\ilcplex\cplex.h"
#if defined(SYSTEME_ESRA)
#include "/usr/local/ilog/cplex90/include/ilcplex/cplex.h"
#endif
*/
#include "pqe_entete_cplex.h"



/***********************/
/* Types et constantes */
/***********************/

#include "pqe_const.h"
#include "pqe_types.h"

#include "pqe_varext.h"

/* --------------------------------- */
/* Prototypes des fonctions externes */
/* --------------------------------- */

#include "pqe_prototypes.h"

/* --------------------------------- */
/* Prototypes des fonctions internes */
/* --------------------------------- */

void pqe_cplex_masquer_affichage (CPXENVptr);
void pqe_cplex_affichage_par_defaut (CPXENVptr, bool);
void pqe_cplex_arret_sur_erreur (int, char *);
void pqe_cplex_verif_statut (int, char *);
void pqe_cplex_parametrage(CPXENVptr, CPXLPptr, int, char_1D, bool, int, int, char *, int_1D, double_1D);
void pqe_cplex_set_pb_data_lin2 (int, int, int, float_1D, float_1D, st_tri_1D, float_1D, float_2D, int,
            char_1D, int *, int *, int *, double_1D, double_1D, char_1D, int_1D, int_1D, int_1D, double_1D, 
            double_1D, double_1D, char_1D, double_1D, double_1D, bival_1D, int_1D, int_1D, int, int_1D, double_1D);
void pqe_cplex_set_pb_data_lkp2(float coef_obj[], float coef_contrainte1[], float coef_contrainte2[],
                             int n, int capa1, int capa2, char probname[], 
                             int *p_mac, int *p_mar, int *p_objsen, double objx[], double rhsx[],
                             char senx[], int matbeg[], int matcnt[], int matind[], double matval[],
                             double bdl[], double bdu[], char **p_dataname, char **p_objname, char **p_rhsname,
                             char **p_rngname, char **p_bndname, char *cname[], char cstore[], char *rname[],
                             char rstore[], int *p_macsz, int *p_marsz, int *p_matsz, unsigned *p_cstorsz,
                             unsigned *p_rstorsz);
			     

void pqe_cplex_set_pb_data_LPi_stocha (bool prime, int i, int nbvar, int COLSPACE, int ROWSPACE, float_1D sur_phi, st_tri_1D ordre_poids,
                                    float_1D q_lin, float_2D q_qua, float_1D cout, float CAP, float_1D moyenne,
                                    float k, int *p_objsen, 
                                    double_1D coef_obj, double_1D rhs, char_1D sense, int_1D matbeg, 
                                    int_1D matcnt, int_1D matind, double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal);

void pqe_cplex_set_pb_data_lin2_stocha (int nbvar, int COLSPACE, int ROWSPACE, float_1D coef_c1z,
                                     float_1D sous_phi, st_tri_1D ordre_poids,
                                     float_1D q_lin, float_2D q_qua, float_1D cout, float CAP,
                                     float_1D moyenne, float k,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement, int_1D pr_indices,
                                     double_1D pr_value);
int pqe_cplex_stocha_callback (CPXENVptr env, void *cbdata, int wherefrom, void *cbhandle);

/******************************************************************************/

/* Part I definitions */

static CPXENVptr     env;
static CPXLPptr      lp = NULL;
static int           status;

#define maxn    500

#define MACSZ   maxn  /* anciennement 4 puis maxn */
#define MARSZ   maxn  /* anciennement 4 puis maxn */
#define MATSZ   maxn  /* anciennement 12 puis maxn - A revoir */
#define CSTORSZ 500
#define RSTORSZ 1000

/* Part I arguments */

static char     probname[16];
static int      mac;
static int      mar;
static int      objsen;
static double   objx[MACSZ];
static double   rhsx[MARSZ];
static char     senx[MARSZ];
static int      matbeg[MACSZ];
static int      matcnt[MACSZ];
static int      matind[MATSZ];
static double   matval[MATSZ];
static double   bdl[MACSZ];
static double   bdu[MACSZ];
static char     *dataname;
static char     *objname;
static char     *rhsname;
static char     *rngname;
static char     *bndname;
static char     cstore[CSTORSZ];
static char     *cname[MACSZ];
static char     rstore[RSTORSZ];
static char     *rname[MARSZ];
static int      macsz;
static int      marsz;
static int      matsz;
static unsigned cstorsz;
static unsigned rstorsz;
static double   obj;

/* Autres, notamment pour sauvegarder la base */
int      cstat[maxn];
int      rstat[maxn];



/*****************************************************************************/

void pqe_cplex_masquer_affichage (CPXENVptr env)
{
	CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
}

/*****************************************************************************/

void pqe_cplex_affichage_par_defaut (CPXENVptr env, bool affichage)
{
	if (affichage) CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
	else           CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
} /* pqe_cplex_affichage_on */

/*****************************************************************************/

void pqe_cplex_arret_sur_erreur (int statut, char *message)
{
	char ligne[1000];

	sprintf(ligne, "%s - statut = %d", message, statut);
	pqe_common_problemo(ligne);
	pqe_cplex_fin();

} /* pqe_cplex_arret_sur_erreur */

/*****************************************************************************/

void pqe_cplex_verif_statut (int statut, char *message)
{
	if (statut) pqe_cplex_arret_sur_erreur(statut, message);

} /* pqe_cplex_verif_statut */

/*****************************************************************************/

void pqe_cplex_init(void)
{
	env = CPXopenCPLEX (&status);
	pqe_cplex_verif_statut (status, "Erreur initialisation Cplex");
	pqe_cplex_affichage_par_defaut (env, affichage);
	CPXsetintparam (env, CPX_PARAM_NODELIM, N_LIMITE);
/*
	CPXsetintparam (env, CPX_PARAM_ROWREADLIM, 100000);
	CPXsetintparam (env, CPX_PARAM_COLREADLIM, 100000);
*/
} /* pqe_cplex_init */

/*****************************************************************************/

void pqe_cplex_init_tab(int COLSPACE, int ROWSPACE, int NZSPACE, int AR_COLSPACE,
                          int AR_ROWSPACE, int AR_NZSPACE, double_1D *p_coef_obj, 
                          double_1D *p_lb, double_1D *p_ub, double_1D *p_rhs, 
                          char_1D *p_ctype, char_1D *p_sense, int_1D *p_matbeg, int_1D *p_matcnt, 
                          int_1D *p_colindex, int_1D *p_priority, int_1D *p_matind, 
                          double_1D *p_matval, double_1D *p_cprim, double_1D *p_rprim, 
                          int_1D *p_pr_indices, double_1D *p_pr_value, double_1D *p_ar_rhs,  
                          char_1D *p_ar_sense, int_1D *p_ar_matbeg, int_1D *p_ar_matind, 
                          double_1D *p_ar_matval, char_2D *p_newrowname, char_2D *p_colname,  
                          char_2D *p_rowname,  char_1D *p_colnamestore, char_1D *p_rownamestore, 
                          char_1D *p_probname, double_1D *p_x, double_1D *p_pi,  
                          double_1D *p_slack, double_1D *p_dj)  
{
	pqe_init_1D((void **) p_coef_obj,   sizeof(double), COLSPACE);
	pqe_init_1D((void **) p_lb,         sizeof(double), COLSPACE);
	pqe_init_1D((void **) p_ub,         sizeof(double), COLSPACE);
	pqe_init_1D((void **) p_rhs,        sizeof(double), ROWSPACE);
	pqe_init_1D((void **) p_ctype,      sizeof(char),   COLSPACE);
	pqe_init_1D((void **) p_sense,      sizeof(char),   ROWSPACE);
	pqe_init_1D((void **) p_matbeg,     sizeof(int),    COLSPACE);
	pqe_init_1D((void **) p_matcnt,     sizeof(int),    COLSPACE);
	pqe_init_1D((void **) p_colindex,   sizeof(int),    COLSPACE);
	pqe_init_1D((void **) p_priority,   sizeof(int),    COLSPACE);
	pqe_init_1D((void **) p_matind,     sizeof(int),    NZSPACE);
	pqe_init_1D((void **) p_matval,     sizeof(double), NZSPACE);
	pqe_init_1D((void **) p_cprim,      sizeof(double), COLSPACE);
	pqe_init_1D((void **) p_rprim,      sizeof(double), COLSPACE);
	pqe_init_1D((void **) p_pr_indices, sizeof(int),    COLSPACE);
	pqe_init_1D((void **) p_pr_value,   sizeof(double), COLSPACE);

	pqe_init_1D((void **) p_ar_rhs,     sizeof(double), AR_ROWSPACE);
	pqe_init_1D((void **) p_ar_sense,   sizeof(char),   AR_ROWSPACE);
	pqe_init_1D((void **) p_ar_matbeg,  sizeof(int),    AR_COLSPACE);
	pqe_init_1D((void **) p_ar_matind,  sizeof(int),    AR_NZSPACE);
	pqe_init_1D((void **) p_ar_matval,  sizeof(double), AR_NZSPACE);
	pqe_init_1D((void **) p_newrowname, sizeof(char *), 1);

	pqe_init_1D((void **) p_colname,      sizeof(char *), COLSPACE);
	pqe_init_1D((void **) p_rowname,      sizeof(char *), ROWSPACE);
	pqe_init_1D((void **) p_colnamestore, sizeof(char),   COLSPACE * 6);
	pqe_init_1D((void **) p_rownamestore, sizeof(char),   ROWSPACE * 10);
	pqe_init_1D((void **) p_probname,     sizeof(char),   20);

	pqe_init_1D((void **) p_x,     sizeof(double), COLSPACE + ROWSPACE);
	pqe_init_1D((void **) p_pi,    sizeof(double), COLSPACE + ROWSPACE);
	pqe_init_1D((void **) p_slack, sizeof(double), COLSPACE + ROWSPACE);
	pqe_init_1D((void **) p_dj,    sizeof(double), COLSPACE + ROWSPACE);

} /* pqe_cplex_init_tab */

/*****************************************************************************/

void pqe_cplex_free_tab(int COLSPACE, int ROWSPACE, double_1D coef_obj, double_1D lb,
                          double_1D ub, double_1D rhs, char_1D ctype, char_1D sense, 
                          int_1D matbeg, int_1D matcnt, int_1D colindex, int_1D priority, 
                          int_1D matind, double_1D matval, double_1D cprim, double_1D rprim, 
                          int_1D pr_indices, double_1D pr_value, double_1D ar_rhs,  
                          char_1D ar_sense, int_1D ar_matbeg, int_1D ar_matind, 
                          double_1D ar_matval, char_2D newrowname, char_2D colname, 
                          char_2D rowname,  char_1D colnamestore, char_1D rownamestore, 
                          char_1D probname, double_1D x, double_1D pi,  double_1D slack, 
                          double_1D dj)  
{
	pqe_init_free_1D(coef_obj);
	pqe_init_free_1D(lb);
	pqe_init_free_1D(ub);
	pqe_init_free_1D(rhs);
	pqe_init_free_1D(ctype);
	pqe_init_free_1D(sense);
	pqe_init_free_1D(matbeg);
	pqe_init_free_1D(matcnt);
	pqe_init_free_1D(colindex);
	pqe_init_free_1D(priority);
	pqe_init_free_1D(matind);
	pqe_init_free_1D(matval);
	pqe_init_free_1D(cprim);
	pqe_init_free_1D(rprim);
	pqe_init_free_1D(pr_indices);
	pqe_init_free_1D(pr_value);
	pqe_init_free_1D(ar_rhs);
	pqe_init_free_1D(ar_sense);
	pqe_init_free_1D(ar_matbeg);
	pqe_init_free_1D(ar_matind);
	pqe_init_free_1D(ar_matval);

/*
	pqe_init_free_2D((void **) newrowname, 1);
	pqe_init_free_2D((void **) colname, COLSPACE);
	pqe_init_free_2D((void **) rowname, ROWSPACE);
*/
	pqe_init_free_1D(newrowname);
	pqe_init_free_1D(colname);
	pqe_init_free_1D(rowname);

	pqe_init_free_1D(colnamestore);
	pqe_init_free_1D(rownamestore);
	pqe_init_free_1D(probname);
	pqe_init_free_1D(x);
	pqe_init_free_1D(pi);
	pqe_init_free_1D(slack);
	pqe_init_free_1D(dj);

} /* pqe_cplex_free_tab */

/*****************************************************************************/

void pqe_cplex_relaxation_continue(CPXENVptr env, CPXLPptr lp, double *p_temps_UB, float *p_meilleur_dual,
                                     double_1D x, double_1D pi, double_1D slack, double_1D dj, bool *p_faisable)
{
	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	int    solstat, status; 
	double val_UB;

	pqe_cplex_masquer_affichage(env);
	print("\r Valeur de la relaxation continue :");
	pqe_common_init_temps (p_temps_UB);
	status = CPXprimopt (env, lp);
	if (status) pqe_common_problemo("relaxation continue");
	*p_faisable = (CPXgetstat(env, lp) != CPX_STAT_INFEASIBLE);
	if (! (*p_faisable)) print("Infaisabilite prouv\8Ee\n");
	pqe_common_calcule_temps(p_temps_UB);
	if ((status = CPXsolution (env, lp, &solstat, &val_UB, x, pi, slack, dj)))
		pqe_common_problemo("Cpxsolution");
	if (*p_faisable) print (" %lf       ", val_UB);
	if (affichage) fflush(stdout);
	*p_meilleur_dual = (float) val_UB;
	pqe_cplex_affichage_par_defaut(env, affichage);

} /* pqe_cplex_relaxation_continue */

/*****************************************************************************/

void pqe_cplex_parametrage(CPXENVptr env, CPXLPptr lp, int COLSPACE, char_1D ctype, bool sol_adm, int type_BB,
                             int var_select, char *nom_prob, int_1D pr_indices, double_1D pr_value)
{
	int status;

	print("\n Resolution en 0-1 :\n");
	status = CPXcopyctype (env, lp, ctype);
	if (status) pqe_common_problemo("Pb chargement ctype");

	/* ----------------------------------------------------------------------- */
	/* Chargement d'une bonne solution admissible et d'un ordre de branchement */
	/* ----------------------------------------------------------------------- */

	if (sol_adm)
	{
		//CPXsetintparam (env, CPX_PARAM_MIPSTART, CPX_ON); valable pour Cplex90 pas pour Cplex101
		CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON); 
		status = CPXcopymipstart (env, lp, COLSPACE, pr_indices, pr_value);
		if (status != 0) pqe_common_problemo("CPXcopymipstart");
	}

	switch (type_BB)
	{
		case BB_PROF :       CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_DFS); break;
		case BB_BEST_BOUND : CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTBOUND); break;
		case BB_BEST_EST :   CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTEST); break;
	}

	switch (var_select)
	{
		case VAR_SEL_POUSSEE : CPXsetintparam (env, CPX_PARAM_VARSEL, CPX_VARSEL_STRONG);
		case VAR_SEL_DEFAUT :  break;
	}

    if (mode_debug) CPXlpwrite (env, lp, nom_prob); 

	/* ----------------------------------------------------- */
	/* le temps de resolution est limite a T_LIMITE secondes */
	/* ----------------------------------------------------- */

	/* 
		CPXsetintparam (env, CPX_PARAM_STARTALG, CPXALG_BARRIER); 
		CPXsetintparam (env, CPX_PARAM_SUBALG, CPXALG_DUAL); 
	*/

	CPXsetdblparam (env, CPX_PARAM_TILIM, T_LIMITE);

} /* pqe_cplex_parametrage */

/*****************************************************************************/

void pqe_cplex_recuperation_solution(CPXENVptr env, CPXLPptr lp, bool *p_temps_limite, double *p_obj,
                                       double_1D v_sol, double *p_nb_noeuds)
{
	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */
	
	int  mipstat;

	mipstat = CPXgetstat(env,lp);

	if (   ( mipstat != CPXMIP_OPTIMAL ) 
	    && ( mipstat != CPXMIP_OPTIMAL_TOL )
		&& ( mipstat != CPXMIP_NODE_LIM_FEAS) 
		&& ( mipstat != CPXMIP_TIME_LIM_FEAS) 
		&& ( mipstat != CPXMIP_TIME_LIM_INFEAS) ) pqe_cplex_arret_sur_erreur(mipstat, "recuperation_sol");

	if (mipstat == CPXMIP_NODE_LIM_FEAS) 
	{
		pqe_common_problemo("Limite du nb de noeuds atteinte, il existe une solution entiere...");
	}
	
	*p_temps_limite = ((mipstat == CPXMIP_TIME_LIM_FEAS) || (mipstat == CPXMIP_TIME_LIM_INFEAS));

	CPXgetmipobjval(env, lp, p_obj);
	CPXgetmipx (env, lp, v_sol, 0, CPXgetnumcols(env,lp) - 1);
	(*p_nb_noeuds) = CPXgetnodecnt(env, lp);

} /* pqe_cplex_recuperation_solution */

/*****************************************************************************/

void pqe_cplex_set_souskpi_data(double coef_obj[maxn], double coef_contrainte[maxn], int n, int capa,
                             char probname[], 
                             int *p_mac, int *p_mar, int *p_objsen, double objx[], double rhsx[],
                             char senx[], int matbeg[], int matcnt[], int matind[], double matval[],
                             double_1D lb, double_1D ub, char **p_dataname, char **p_objname, char **p_rhsname,
                             char **p_rngname, char **p_bndname, char *cname[], char cstore[], char *rname[],
                             char rstore[], int *p_macsz, int *p_marsz, int *p_matsz, unsigned *p_cstorsz,
                             unsigned *p_rstorsz, int ii)
{
    int        i, j, k, cpt, amax;
    static int num_appel = 0;
    char chaine[5];

    strcpy(probname, "souskpi");

    *p_mac = n;          /* Nb effectif de variables */
    *p_mar = n + 2;      /* Nb de contraintes */
    *p_objsen = CPX_MIN; /* Minimisation */

    /* -------------------------------- */
    /* Coefficients de la fonction obj. */
    /* -------------------------------- */

    for (i=0; i<n; i++) objx[i] = coef_obj[i];
    num_appel ++;
		

    /* ---------------------- */
    /* Vecteur second-membre) */
    /* ---------------------- */
    rhsx[0] = capa; senx[0] = 'L';
  
    for (k = 1; k <= n; k++) {
        rhsx[k] = capa - coef_contrainte[k - 1];
        senx[k] = 'G';
    }
    
    rhsx[n + 1] = 0.0; senx[n + 1] = 'E';

    for (i=0; i < n; i++) {
        lb[i] = 0.0;
        ub[i] = 1.0;
    }

    /* ----------------------------------------------------------------------------------------------------- */
    /* Matrice des contraintes         :                                                                     */
    /* Elle est lue colonne par colonne :                                                                    */
    /* matbeg[j] est l'indice dans matind et matval du premier element non nul (eventuel) de la colonne j    */
    /* matcnt[j] est le nombres de coefficients non nuls de la colonne j                                     */
    /* matind[k] est le numero de la ligne a laquelle appartient le kieme element non nul de la matrice, les */
    /*           elements etant lus colonne apres colonne                                                    */ 
    /* matval[k] est la valeur du kieme element non nuls de la matrice                                       */
    /* ----------------------------------------------------------------------------------------------------- */
        
    cpt = 0;
    amax = coef_contrainte[0];
    for (k = 1; k < n; k++) amax = (coef_contrainte[k] > amax ? coef_contrainte[k] : amax);

    /* print(" \n tip 1 \n"); */

    for (j = 0; j <= n - 1; j++)
    {
        /* print(" j=%d, n = %d, ii = %d \n",j,n,ii); */
        matbeg[j] = cpt;
        matcnt[j] = (j == ii ? n + 2 : n + 1);
      
        for (k = 0; k <= n; k++) {
            /* print(" tip%d ",k); */
            matind[cpt]   = k;
            if (k != j + 1) matval[cpt] = coef_contrainte[j];
            else matval[cpt] = amax - coef_contrainte[j];
            cpt++;
        }
        
        if (j == ii) {
            matind[cpt] = n + 1;
            matval[cpt++] = 1.0;
        }
    }   
  
    /*  print(" \n tip 3 \n"); */
  
    *p_dataname = NULL;
    *p_objname = NULL;
    *p_rhsname = NULL;
    *p_rngname = NULL;
    *p_bndname = NULL;

    cpt = 0;
    for (i=0; i<n; i++)
    {
	cname[i] = &cstore[cpt];
	cstore[cpt++] = 'x';
        pqe_common_intochar(chaine,i+1);
        j = 0;
	while (chaine[j] != '\0') {cstore[cpt++] = chaine[j++];};
	cstore[cpt++] = '\0';
    }

    rstore[0] = 'c'; rstore[1] = '1', rstore[2] = '\0';
    rname[0] = &rstore[0];

    *p_macsz = MACSZ;
    *p_marsz = MARSZ;
    *p_matsz = MATSZ;
    *p_cstorsz = CSTORSZ;
    *p_rstorsz = RSTORSZ;
    
    /*    print(" \n sortie de souskpidata \n"); */

} /* pqe_cplex_set_souskpi_data */

/*****************************************************************************/

void pqe_cplex_set_pb_data(float coef_obj[maxn], float coef_contrainte[maxn], int n, int capa, char probname[],
                             int *p_mac, int *p_mar, int *p_objsen, double objx[], double rhsx[],
                             char senx[], int matbeg[], int matcnt[], int matind[], double matval[],
                             double bdl[], double bdu[], char **p_dataname, char **p_objname, char **p_rhsname,
                             char **p_rngname, char **p_bndname, char *cname[], char cstore[], char *rname[],
                             char rstore[], int *p_macsz, int *p_marsz, int *p_matsz, unsigned *p_cstorsz,
                             unsigned *p_rstorsz)
{
	int i,j,cpt;
	static int num_appel = 0;
	char chaine[5];

	strcpy(probname, "lkp");

	*p_mac = n;          /* Nb effectif de variables */
	*p_mar = 1;          /* Nb de contraintes */
	*p_objsen = CPX_MAX; /* Maximisation */

	/* -------------------------------- */
	/* Coefficients de la fonction obj. */
	/* -------------------------------- */

	for (i=0; i<n; i++) objx[i] = coef_obj[i];
	num_appel ++;
		

	/* ------------------------- */
	/* Vecteur b (second-membre) */
	/* ------------------------- */
	rhsx[0] = capa; senx[0] = 'L';

	/* ------------------------- */
	/* Matrice A des contraintes */
	/* ------------------------- */
	cpt = 0;
	for (i=0; i<n; i++)
	{
		matcnt[i] = 0;
		matbeg[i] = cpt;
		if (coef_contrainte[i] != 0)
		{
			matcnt[i] ++;
			matind[cpt] = 0;
			matval[cpt] = coef_contrainte[i];
			cpt ++;
		}

		bdl[i] = 0.0;
		bdu[i] = 1.0;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	*p_dataname = NULL;
	*p_objname = NULL;
	*p_rhsname = NULL;
	*p_rngname = NULL;
	*p_bndname = NULL;

	cpt = 0;
	for (i=0; i<n; i++)
	{
		cname[i] = &cstore[cpt];
		cstore[cpt++] = 'y';
		pqe_common_intochar(chaine,i+1);
		j = 0;
		while (chaine[j] != '\0') {cstore[cpt++] = chaine[j++];};
		cstore[cpt++] = '\0';
	}

	rstore[0] = 'c'; rstore[1] = '1', rstore[2] = '\0';
	rname[0] = &rstore[0];

	*p_macsz = MACSZ;
	*p_marsz = MARSZ;
	*p_matsz = MATSZ;
	*p_cstorsz = CSTORSZ;
	*p_rstorsz = RSTORSZ;

} /* pqe_cplex_set_pb_data */

/*****************************************************************************/

void pqe_cplex_set_pb_data_lkp2(float coef_obj[maxn], float coef_contrainte1[maxn], float coef_contrainte2[maxn],
                             int n, int capa1, int capa2, char probname[], 
                             int *p_mac, int *p_mar, int *p_objsen, double objx[], double rhsx[],
                             char senx[], int matbeg[], int matcnt[], int matind[], double matval[],
                             double bdl[], double bdu[], char **p_dataname, char **p_objname, char **p_rhsname,
                             char **p_rngname, char **p_bndname, char *cname[], char cstore[], char *rname[],
                             char rstore[], int *p_macsz, int *p_marsz, int *p_matsz, unsigned *p_cstorsz,
                             unsigned *p_rstorsz)
{
	int i,j,cpt;
	static int num_appel = 0;
	char chaine[5];

	strcpy(probname, "2lkp");

	*p_mac = n;          /* Nb effectif de variables */
	*p_mar = 2;          /* Nb de contraintes */
	*p_objsen = CPX_MAX; /* Maximisation */

	/* -------------------------------- */
	/* Coefficients de la fonction obj. */
	/* -------------------------------- */

	for (i=0; i<n; i++) objx[i] = coef_obj[i];
	num_appel ++;
		

	/* ------------------------- */
	/* Vecteur b (second-membre) */
	/* ------------------------- */
	rhsx[0] = capa1; senx[0] = 'L';
	rhsx[1] = capa2; senx[1] = 'G';

	/* ------------------------- */
	/* Matrice A des contraintes */
	/* ------------------------- */
	cpt = 0;
	for (i=0; i<n; i++)
	{
		matcnt[i] = 0;
		matbeg[i] = cpt;
		if (coef_contrainte1[i] != 0)
		{
			matcnt[i] ++;
			matind[cpt] = 0;
			matval[cpt] = coef_contrainte1[i];
			cpt ++;
		}
		if (coef_contrainte2[i] != 0)
        {
            matcnt[i] ++;
            matind[cpt] = 1;
            matval[cpt] = coef_contrainte2[i];
            cpt ++;
        }
                

		bdl[i] = 0.0;
		bdu[i] = 1.0;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	*p_dataname = NULL;
	*p_objname = NULL;
	*p_rhsname = NULL;
	*p_rngname = NULL;
	*p_bndname = NULL;

	cpt = 0;
	for (i=0; i<n; i++)
	{
		cname[i] = &cstore[cpt];
		cstore[cpt++] = 'y';
		pqe_common_intochar(chaine,i+1);
		j = 0;
		while (chaine[j] != '\0') {cstore[cpt++] = chaine[j++];};
		cstore[cpt++] = '\0';
	}

	rstore[0] = 'c'; rstore[1] = '1', rstore[2] = '\0';
	rname[0] = &rstore[0];

	*p_macsz = MACSZ;
	*p_marsz = MARSZ;
	*p_matsz = MATSZ;
	*p_cstorsz = CSTORSZ;
	*p_rstorsz = RSTORSZ;

} /* pqe_cplex_set_pb_data_lkp2 */

/*****************************************************************************/

void pqe_cplex_set_pb_data_demi_lin (int nbvar, float_1D coef_c1z, st_tri_1D ordre_poids,
                                     float_1D q_lin, float_2D q_qua, int capacite,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement)
{
	char part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE];
	int  i, j, cpt2, indice, minij, maxij, nb;

	pqe_common_intochar(part1,nbvar);
    pqe_common_intochar(part2,(int) (densite*100));
    pqe_common_intochar(part3,seed);
	
	strcpy(probname,"newlin");
	strcat(probname,part1);
	strcat(probname,"D");
	strcat(probname,part2);
	strcat(probname,"G");
	strcat(probname,part3);

	*numcols_p = 2 * nbvar;      /* zn n'est pas inutile          */
	*numrows_p = 2 * nbvar + 1;  /* c1zn et c2zn sont necessaires */
	*objsen_p  = -1;             /* Maximisation                  */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	for (i = 0; i < nbvar; i++) obj[i] = q_lin[ordre_poids[i + 1].indice];
	for (i = nbvar; i <= 2 * nbvar - 1; i++) obj[i] = 0.5;

	/* ----------------------------------------------------------------------------------------------------- */
	/* Matrice des contraintes         :                                                                     */
	/* Elle est lue colonne par colonne :                                                                    */
	/* matbeg[j] est l'indice dans matind et matval du premier element non nul (eventuel) de la colonne j    */
	/* matcnt[j] est le nombres de coefficients non nuls de la colonne j                                     */
	/* matind[k] est le numero de la ligne a laquelle appartient le kieme element non nul de la matrice, les */
	/*           elements etant lus colonne apres colonne                                                    */ 
	/* matval[k] est la valeur du kieme element non nuls de la matrice                                       */
	/* ----------------------------------------------------------------------------------------------------- */

	/* -------------------- */
	/* n premieres colonnes */
	/* -------------------- */

	for (j = 0; j <= nbvar - 1; j++)
	{
		indice = ordre_poids[j + 1].indice;

		matbeg[j] = j * (nbvar + 1);
		matcnt[j] = nbvar + 1;

		matind[matbeg[j]] = 0;
		matval[matbeg[j]] = a[indice];                   /* Contrainte de capacite */

		matind[matbeg[j] + 1] = j + 1;
		matval[matbeg[j] + 1] = - coef_c1z[indice];        /* Contrainte c1z */

		cpt2 = matbeg[j] + 2;

		for (i = 1; i < j + 1; i++)
		{
			/* Contrainte c2z */

			minij = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
			maxij = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

			matind[cpt2] = nbvar + i;
			matval[cpt2] = - q_qua[minij][maxij];
			cpt2++;
		}

		for (i = j + 2; i <= nbvar; i++)
		{
			/* Contrainte c2z */

			minij = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
			maxij = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

			matind[cpt2] = nbvar + i;
			matval[cpt2] = - q_qua[minij][maxij];
			cpt2++;
		}
	}

	/* -------------------- */
	/* Colonnes de n a 2n-1 */
	/* -------------------- */

	nb = nbvar * (nbvar + 1);

	for (j = nbvar; j <= 2 * nbvar - 1; j++)
	{
		matbeg[j] = nb + 2 * (j - nbvar);
		matcnt[j] = 2;

		matind[matbeg[j]] = j - nbvar + 1;
		matval[matbeg[j]] = 1.0;

		matind[matbeg[j] + 1] = j + 1;
		matval[matbeg[j] + 1] = 1.0;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'I';
	}

	for (j = nbvar; j <= 2 * nbvar - 1; j++)
	{
		lb[j] = 0.0;
		/* ub[j] = INFBOUND; */
		ub[j] = coef_c1z[ordre_poids[j - nbvar + 1].indice];
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = (float) capacite;

	for (i = 1; i <= 2 * nbvar; i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}

	/* ----------------------------------------------------------------------------- */
	/* Preparation du vecteur contenant la meilleure solution admissible deja connue */
	/* Preparation des vecteurs colindex et priority pour l'ordre de branchement     */
	/* ----------------------------------------------------------------------------- */

	for (i = 0; i <= nbvar - 1; i++)
	{
		indice = ordre_poids[i + 1].indice;
		/* ----------------- */
		/* Variables x1 a xn */
		/* ----------------- */
		colindex[i] = i;

		switch (type_branchement)
		{
			case BRAN_NAT :   priority[i] = nbvar - i + 60;
			case BRAN_PRI_1 : priority[i] = (b_primal[indice] == 1 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_PRI_0 : priority[i] = (b_primal[indice] == 0 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_RAT   : priority[i] = nbvar - i + 60; /* Suppose ordre_var == VAR_RAT_DEC */
		}

		indice = ordre_poids[i + 1].indice;
		cprim[i] = (double) b_primal[indice];
		rprim[i] = 0.0;
	}

	for (i = nbvar; i <= 2 * nbvar - 1; i++)
	{
		/* ----------------- */
		/* Variables z1 a zn */
		/* ----------------- */
		indice = ordre_poids[i - nbvar + 1].indice;
		rprim[i] = 0.0;
		cprim[i] = 0.0;
	}

} /* pqe_cplex_set_pb_data_demi_lin */

/*****************************************************************************/

void pqe_cplex_set_pb_data_lin_croisee (int nbvar, float_1D Lmax, float_1D sur_phi, float_1D Dmax, int COLSPACE, int ROWSPACE,
                                     st_tri_1D ordre_poids, float_1D q_lin, float_2D q_qua, int capacite,
                                     char_1D probname, int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement)
{
/*
	char part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE];
	int  i, j, compteur, cpt, indice, gauche, droite, nb;

	pqe_common_intochar(part1,nbvar);
    pqe_common_intochar(part2,(int) (densite*100));
    pqe_common_intochar(part3,seed);
	
	strcpy(probname,"lin_cr");
	strcat(probname,part1);
	strcat(probname,"D");
	strcat(probname,part2);
	strcat(probname,"G");
	strcat(probname,part3);

	*objsen_p = CPX_MAX;               Maximisation                    

	 Fonction objectif 

	obj[0] = 1.0;
	for (i = 1; i <= nbvar; i++) obj[i] = q_lin[ordre_poids[i].indice];
	for (i = nbvar + 1; i < COLSPACE; i++) obj[i] = 0.0;

	 Matrice des contraintes                                                                  
	 Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,      
	 la matrice des contraintes etant lus colonne par colonne. matcnt[j] est le nombre de     
	 coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      
	 non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] 
	 indique le numero de la ligne correspondant au coefficient.                             

	compteur = 0;  Indice dans matval et matind 

	 --------- 
	 Colonne 0 
	 --------- 

	matbeg[0] = 0;
	matcnt[0] = 3;
	matval[0] = matval[1] = matval[2] = 1.0;
	matind[compteur++] = ROWSPACE - 3;
	matind[compteur++] = ROWSPACE - 2;
	matind[compteur++] = ROWSPACE - 1;

	 --------- 
	 Colonne 1 
	 --------- 

	indice = ordre_poids[1].indice;
	matbeg[1] = compteur;
	matcnt[1] = 2 * nbvar;

	matval[compteur] = a[indice];
	matind[compteur++] = 0;

	matval[compteur] = -Lmax[indice];
	matind[compteur++] = 1;

	for (i = 2; i <= nbvar; i++)
	{
		gauche = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
		droite = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

		matval[compteur] = -q_qua[gauche][droite];
		matind[compteur++] = 3 * nbvar + i - 4;
	}


	 ---------------- 
	 Colonnes 2 a n-1 
	 ---------------- 

	for (j = 1; j <= nbvar; j++)
	{
		matcnt[j] = nbvar + nb
	}
	
	for (j = 0; j <= nbvar - 2; j++)
	{
		indice = ordre_poids[j + 1].indice;

		matbeg[j] = j * (j + 3) / 2;
		matcnt[j] = 2 + j;

		matind[matbeg[j]] = 0;
		matval[matbeg[j]] = a[indice];                    Contrainte de capacite 

		matind[matbeg[j] + 1] = j + 1;
		matval[matbeg[j] + 1] = - coef_c1z[indice];         Contrainte c1z 

		for (i = 2; i < matcnt[j]; i++)
		{
			 Contrainte c2z 

			minij = (ordre_poids[i-1].indice < indice ? ordre_poids[i-1].indice : indice);
			maxij = (ordre_poids[i-1].indice > indice ? ordre_poids[i-1].indice : indice); 

			matind[matbeg[j] + i] = nbvar + i - 2;
			matval[matbeg[j] + i] = - q_qua[minij][maxij];
		}
	}

	 n-1 eme colonne 

	j = nbvar - 1;
	indice = ordre_poids[nbvar].indice;
	matbeg[nbvar - 1] = (nbvar - 1) * (nbvar + 2) / 2;
	matcnt[nbvar - 1] = nbvar;
	matind[matbeg[nbvar - 1]] = 0;
	matval[matbeg[nbvar - 1]] = a[ordre_poids[nbvar].indice];

	for (i = 1; i < matcnt[j]; i++)
	{
		 Contrainte c2z 

		minij = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
		maxij = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

		matind[matbeg[j] + i] = nbvar + i - 1;
		matval[matbeg[j] + i] = - q_qua[minij][maxij];
	}

	 Colonnes de n a 2n-2 

	nb = (nbvar * nbvar + 3 * nbvar - 2) / 2;

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		matbeg[j] = nb + 2 * (j - nbvar);
		matcnt[j] = 2;

		matind[matbeg[j]] = j - nbvar + 1;
		matval[matbeg[j]] = 1.0;

		matind[matbeg[j] + 1] = j;
		matval[matbeg[j] + 1] = 1.0;
	}

	 ------ 
	 Bornes 
	 ------ 

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'I';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		lb[j] = 0.0;
		ub[j] = INFBOUND;
		ctype[j] = 'C';
	}

	 ---------------- 
	 Membre de droite 
	 ---------------- 

	sense[0] = 'L';
	rhs[0]   = (float) capacite;

	for (i = 1; i <= 2 * nbvar - 2; i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}

	 ----------------------------------------------------------------------------- 
	 Preparation du vecteur contenant la meilleure solution admissible deja connue 
	 Preparation des vecteurs colindex et priority pour l'ordre de branchement     
	 ----------------------------------------------------------------------------- 

	for (i = 0; i <= nbvar - 1; i++)
	{
		indice = ordre_poids[i + 1].indice;
		 ----------------- 
		 Variables x1 a xn 
		 ----------------- 
		colindex[i] = i;

		switch (type_branchement)
		{
			case BRAN_NAT :   priority[i] = nbvar - i + 60;
			case BRAN_PRI_1 : priority[i] = (b_primal[indice] == 1 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_PRI_0 : priority[i] = (b_primal[indice] == 0 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_RAT   : priority[i] = nbvar - i + 60;  Suppose ordre_var == VAR_RAT_DEC 
		}

		cprim[i] = (double) b_primal[indice];
		rprim[i] = 0.0;
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		 ------------------- 
		 Variables z1 a zn-1 
		 ------------------- 
		indice = ordre_poids[i - nbvar + 1].indice;
		rprim[i] = 0.0;
		cprim[i] = 0.0;


		if (b_primal[indice] > 0)
		{
			for (j = i - nbvar + 2; j <= nbvar; j++)
			{
				minij = (ordre_poids[j].indice < indice ? ordre_poids[j].indice : indice);
				maxij = (ordre_poids[j].indice > indice ? ordre_poids[j].indice : indice); 
				cprim[i] += q_qua[minij][maxij] * b_primal[minij] * b_primal[maxij];
			}
		}

	}
*/
} /* pqe_cplex_set_pb_data_lin_croisee */
	
/*****************************************************************************/

void pqe_cplex_set_pb_data_lin_classique (int nbvar, float_1D q_lin, float_2D q_qua, int capacite,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement, int_1D pr_indices,
                                     double_1D pr_value)
{
	char part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE];
	int  i, j, cpt, cpt2, col, saut;

	pqe_common_intochar(part1,nbvar);
    pqe_common_intochar(part2,(int) (densite*100));
    pqe_common_intochar(part3,seed);
	
	strcpy(probname,"linclass");
	strcat(probname,part1);
	strcat(probname,"D");
	strcat(probname,part2);
	strcat(probname,"G");
	strcat(probname,part3);

	*numcols_p = nbvar * (nbvar + 1) / 2;   /* Nb variables   */
	*numrows_p = nbvar * nbvar - nbvar + 1; /* Nb contraintes */
	*objsen_p = CPX_MAX;                    /* Maximisation   */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;
	for (i = 1; i <= nbvar; i++) obj[cpt++] = q_lin[i];
	for (i = 1; i <= nbvar; i++) 
		for (j = i + 1; j <= nbvar; j++)
			obj[cpt++] = q_qua[i][j];

	/* ------------------------------------------- */
	/* Matrice des contraintes                     */
	/* Sigma ai <= b                               */
	/* wij <= xi pour i de 1 a n-1 et j de i+1 a n */
	/* wij <= xj pour i de 1 a n-1 et j de i+1 a n */
	/* ------------------------------------------- */

	/* ---------------- */
	/* Colonnes 0 a n-1 */
	/* ---------------- */
	cpt = 0;
	col = 0;

	for (i = 1; i <= nbvar; i++)
	{
		matbeg[col]   = cpt;
		matcnt[col++] = nbvar;  /* 1 + (nbvar - i) + (i - 1) */		

		matind[cpt]   = 0;
		matval[cpt++] = a[i];                   /* Contrainte de capacite */
		
		for (j = 1; j <= nbvar - i; j++)
		{
			matind[cpt]   = j + (i-1)*nbvar - i*(i-1)/2;
			matval[cpt++] = -1.0;
		}
	
		saut =  nbvar * (nbvar - 1) / 2 + i -1;
		for (j = 1; j <= i - 1; j++)
		{
			matind[cpt]   = saut;
			matval[cpt++] = -1.0;
			saut += nbvar - 1 - j;
		}
	}

	/* ------------------ */
	/* Colonnes suivantes */
	/* ------------------ */

	saut = nbvar * (nbvar - 1) / 2;
	cpt2 = 1;
	for (i = 1; i <= nbvar; i++) 
		for (j = i + 1; j <= nbvar; j++)
		{
			matbeg[col]    = cpt;
			matcnt[col++] = 2;

			matind[cpt] =   cpt2;
			matval[cpt++] = 1.0;
			
			matind[cpt] =   cpt2 + saut;
			matval[cpt++] = 1.0;
	
			cpt2++;
		}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (i = 1; i <= nbvar; i++)
	{
		lb[i - 1] = 0.0;
		ub[i - 1] = 1.0;
		ctype[i - 1] = 'B';
	}

	cpt2 = nbvar;
	for (i = 1; i <= nbvar; i++) 
		for (j = i + 1; j <= nbvar; j++)
		{
			lb[cpt2] = 0.0;
			ub[cpt2] = CPX_INFBOUND;
			ctype[cpt2] = 'C';
			cpt2++;
		}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = (float) capacite;

	for (i = 1; i <= nbvar * (nbvar - 1); i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}

} /* pqe_cplex_set_pb_data_lin_classique */

/*****************************************************************************/

void pqe_cplex_set_pb_data_sous_phii (int i, int nbvar, int COLSPACE, int ROWSPACE, float_1D sur_phi, st_tri_1D ordre_poids,
                                    float_1D q_lin, float_2D q_qua, int_1D a, int capacite, int *p_objsen, 
                                    double_1D coef_obj, double_1D rhs, char_1D sense, int_1D matbeg, 
                                    int_1D matcnt, int_1D matind, double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal)
{
	int    k, j, l, cpt, mini, maxi, indice, ind_i, compteur;

	*p_objsen = CPX_MIN;  /* Minimisation */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;
	ind_i = ordre_poids[i].indice;
	for (j = 1; j <= i; j++) coef_obj[cpt++] = 0.0;
	for (j = i + 1; j <= nbvar; j++) 
	{
		mini = (ordre_poids[j].indice < ind_i ? ordre_poids[j].indice : ind_i);
		maxi = (ordre_poids[j].indice > ind_i ? ordre_poids[j].indice : ind_i); 
		coef_obj[cpt++] = q_qua[mini][maxi];
	}
	while (cpt < COLSPACE) coef_obj[cpt++] = 0.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* ---------------- */
	/* Colonnes 0 a i-2 */
	/* ---------------- */

	for (k = 0; k <= i - 2; k++)
	{
		indice = ordre_poids[k + 1].indice;
		matbeg[k] = compteur;
		matcnt[k] = 3 + k;
		matval[compteur]   = a[indice];
		matind[compteur++] = 0;
		matval[compteur]   = -sur_phi[indice];
		matind[compteur++] = k + 2; 
		for (l = 1; l <= k; l++)
		{
			mini = (ordre_poids[l].indice < indice ? ordre_poids[l].indice : indice);
			maxi = (ordre_poids[l].indice > indice ? ordre_poids[l].indice : indice); 
			matval[compteur]   = -q_qua[mini][maxi];
			matind[compteur++] = l + nbvar - 1;
		}
		matval[compteur]   = q_lin[indice];
		matind[compteur++] = ROWSPACE-1;
	}

	/* ----------- */
	/* Colonne i-1 */
	/* ----------- */

	matbeg[i - 1] = compteur;
	matcnt[i - 1] = i + 2;
	matval[compteur]   = a[ind_i];
	matind[compteur++] = 0; 
	matval[compteur]   = 0.0; /* 1.0 si on veut xi = 1, 0.0 sinon */
	matind[compteur++] = 1; 
	for (l = 1; l <= i-1; l++)
	{
		mini = (ordre_poids[l].indice < ind_i ? ordre_poids[l].indice : ind_i);
		maxi = (ordre_poids[l].indice > ind_i ? ordre_poids[l].indice : ind_i); 
		matval[compteur]   = -q_qua[mini][maxi];
		matind[compteur++] = l + nbvar - 1;
	}
	matval[compteur]   = q_lin[ind_i];
	matind[compteur++] = ROWSPACE-1;

	/* ---------------- */
	/* Colonnes i a n-1 */
	/* ---------------- */

	for (k = i; k <= nbvar - 1; k++)
	{
		indice = ordre_poids[k + 1].indice;
		matbeg[k] = compteur;
		matcnt[k] = (k == (nbvar-1) ? k + 1 : k + 2);
		matval[compteur]   = a[indice];
		matind[compteur++] = 0;
		if (k < nbvar - 1)
		{
			matval[compteur]   = -sur_phi[indice];
			matind[compteur++] = k + 1; 
		}
		for (l = 1; l <= i-1; l++)
		{
			mini = (ordre_poids[l].indice < indice ? ordre_poids[l].indice : indice);
			maxi = (ordre_poids[l].indice > indice ? ordre_poids[l].indice : indice); 
			matval[compteur]   = -q_qua[mini][maxi];
			matind[compteur++] = l + nbvar - 1;
		}
		for (l = i + 1; l <= k; l++)
		{
			mini = (ordre_poids[l].indice < indice ? ordre_poids[l].indice : indice);
			maxi = (ordre_poids[l].indice > indice ? ordre_poids[l].indice : indice); 
			matval[compteur]   = -q_qua[mini][maxi];
			matind[compteur++] = l + nbvar - 2;
		}
		mini = (indice < ind_i ? indice : ind_i);
		maxi = (indice > ind_i ? indice : ind_i); 
		matval[compteur]   = q_lin[indice] + q_qua[mini][maxi];
		matind[compteur++] = ROWSPACE-1;
	}

	/* ----------------- */
	/* Colonnes n a 2n-3 */
	/* ----------------- */

	for (k = nbvar; k <= 2 * nbvar - 3; k++)
	{
		matbeg[k] = compteur;
		matcnt[k] = 3;
		matval[compteur]   = 1.0;
		matind[compteur++] = k - nbvar + 2;
		matval[compteur]   = 1.0;
		matind[compteur++] = k;
		matval[compteur]   = 1.0;
		matind[compteur++] = ROWSPACE-1;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'B';
	}

	for (j = nbvar; j <= 2 * nbvar - 3; j++)
	{
		lb[j] = -CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = (float) capacite;
	sense[1] = 'E';
	rhs[1]   = 0.0;   /* 1.0 si on veut xi = 1 , 0.0 sinon */
	for (k = 2; k <= 2 * nbvar - 3; k++)
	{
		sense[k] = 'L';
		rhs[k]   = 0.0;
	}
	sense[ROWSPACE-1] = 'G';
	rhs[ROWSPACE-1]   = meilleur_primal;

} /* pqe_cplex_set_pb_data_sous_phii */

/*****************************************************************************/

void pqe_cplex_set_pb_data_LPi_stocha (bool prime, int i, int nbvar, int COLSPACE, int ROWSPACE, float_1D sur_phi, st_tri_1D ordre_poids,
                                    float_1D q_lin, float_2D q_qua, float_1D cout, float CAP, float_1D moyenne,
                                    float k, int *p_objsen, 
                                    double_1D coef_obj, double_1D rhs, char_1D sense, int_1D matbeg, 
                                    int_1D matcnt, int_1D matind, double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal)
{
	/* ------------------------------------------- */
	/* prime indique s'il s'agit de LP'i ou de LPi */
	/* ------------------------------------------- */

	int    p, j, l, cpt, mini, maxi, indice, ind_i, compteur;
	int    extra_ligne, controle1, controle2;

	*p_objsen = CPX_MIN;  /* Minimisation */
	extra_ligne = (prime ? 1 : 0);


	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;
	ind_i = ordre_poids[i].indice;
	for (j = 1; j <= i; j++) coef_obj[cpt++] = 0.0;
	for (j = i + 1; j <= nbvar; j++) 
	{
		mini = (ordre_poids[j].indice < ind_i ? ordre_poids[j].indice : ind_i);
		maxi = (ordre_poids[j].indice > ind_i ? ordre_poids[j].indice : ind_i); 
		coef_obj[cpt++] = q_qua[mini][maxi];
	}
	for (j = 1; j <= nbvar - 1; j++) coef_obj[cpt++] = 0.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* -------------------- */
	/* Colonnes 0 a nbvar-1 */
	/* -------------------- */

	for (p = 0; p <= nbvar - 1; p++)
	{
		controle1 = compteur;
		indice = ordre_poids[p + 1].indice;
		matbeg[p] = compteur;
		matcnt[p] = (p == (nbvar-1) ? p + 3 : p + 4);
		if ((prime) && (p == i - 1)) 
		{
			matcnt[p] = matcnt[p] + 1;
			matval[compteur] = 1.0;
			matind[compteur++] = 0;
		}
		matval[compteur]   = cout[indice];
		matind[compteur++] = extra_ligne;

        matval[compteur]   = moyenne[indice];
        matind[compteur++] = 1 + extra_ligne;
		if (p < nbvar - 1)
		{
			matval[compteur]   = -sur_phi[indice];
			matind[compteur++] = p + 2 + extra_ligne; 
		}
		for (l = 1; l <= p; l++)
		{
			mini = (ordre_poids[l].indice < indice ? ordre_poids[l].indice : indice);
			maxi = (ordre_poids[l].indice > indice ? ordre_poids[l].indice : indice); 
			matval[compteur]   = -q_qua[mini][maxi];
			matind[compteur++] = l + nbvar + extra_ligne;
		}
		matval[compteur]   = (double) q_lin[indice];
        print("\n matval[compteur]   = %lf", matval[compteur]);
		matind[compteur++] = 2 * nbvar + extra_ligne;
		controle2 = compteur;
		if (controle2 - controle1 != matcnt[p]) pqe_common_problemo("cplex_set_pb_data_LPi_stocha");
	}

	/* ----------------- */
	/* Colonnes n a 2n-2 */
	/* ----------------- */

	for (p = nbvar; p <= 2 * nbvar - 2; p++)
	{
		matbeg[p] = compteur;
		matcnt[p] = 3;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + p - nbvar + 2;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + p + 1;
		matval[compteur]   = 1.0;
		matind[compteur++] = 2 * nbvar + extra_ligne;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'B';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		lb[j] = -CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	if (prime)
	{
		sense[0] = 'E';
		rhs[0]   = 0.0;
	}
	sense[extra_ligne] = 'L';
	rhs[extra_ligne]   = CAP;
    sense[extra_ligne + 1] = 'G';
    rhs[extra_ligne + 1] = k;
	for (p = extra_ligne + 2; p <= 2 * nbvar + extra_ligne - 1; p++)
	{
		sense[p] = 'L';
		rhs[p]   = 0.0;
	}
	sense[2 * nbvar + extra_ligne] = 'G';
	rhs[2 * nbvar + extra_ligne]   = meilleur_primal;

} /* pqe_cplex_set_pb_data_LPi_stocha */

/*****************************************************************************/

void pqe_cplex_set_pb_data_LPi (bool prime, int i, int nbvar, int COLSPACE, int ROWSPACE, float_1D sur_phi,
                                    st_tri_1D ordre_poids, float_1D q_lin, float_2D q_qua, int_1D a, 
                                    int capacite, int *p_objsen, double_1D coef_obj, double_1D rhs, 
                                    char_1D sense, int_1D matbeg, int_1D matcnt, int_1D matind, 
                                    double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal, int sens)
{
	/* ------------------------------------------------------------------------ */
	/* prime indique s'il s'agit de LP'i ou de LPi : xi = 1 ou 0 ajout\8Ee ou non */
	/* ------------------------------------------------------------------------ */

	int    k, j, l, cpt, mini, maxi, indice, ind_i, compteur;
	int    extra_ligne, controle1, controle2;

	*p_objsen = (sens == 1 ? CPX_MAX : CPX_MIN);  /* Minimisation ou maximisation */
	extra_ligne = (prime ? 1 : 0);


	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;
	ind_i = ordre_poids[i].indice;
	for (j = 1; j <= i; j++) coef_obj[cpt++] = 0.0;
	for (j = i + 1; j <= nbvar; j++) 
	{
		mini = (ordre_poids[j].indice < ind_i ? ordre_poids[j].indice : ind_i);
		maxi = (ordre_poids[j].indice > ind_i ? ordre_poids[j].indice : ind_i); 
		coef_obj[cpt++] = q_qua[mini][maxi];
	}
	for (j = 1; j <= nbvar - 1; j++) coef_obj[cpt++] = 0.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* -------------------- */
	/* Colonnes 0 a nbvar-1 */
	/* -------------------- */

	for (k = 0; k <= nbvar - 1; k++)
	{
		controle1 = compteur;
		indice = ordre_poids[k + 1].indice;
		matbeg[k] = compteur;
		matcnt[k] = (k == (nbvar-1) ? k + 2 : k + 3);
		if ((prime) && (k == i - 1)) 
		{
			matcnt[k] = matcnt[k] + 1;
			matval[compteur] = 1.0;
			matind[compteur++] = 0;
		}
		matval[compteur]   = a[indice];
		matind[compteur++] = extra_ligne;
		if (k < nbvar - 1)
		{
			matval[compteur]   = -sur_phi[indice];
			matind[compteur++] = k + 1 + extra_ligne; 
		}
		for (l = 1; l <= k; l++)
		{
			mini = (ordre_poids[l].indice < indice ? ordre_poids[l].indice : indice);
			maxi = (ordre_poids[l].indice > indice ? ordre_poids[l].indice : indice); 
			matval[compteur]   = -q_qua[mini][maxi];
			matind[compteur++] = l + nbvar + extra_ligne - 1;
		}
		matval[compteur]   = q_lin[indice];
		matind[compteur++] = 2 * nbvar + extra_ligne - 1;
		controle2 = compteur;
		if (controle2 - controle1 != matcnt[k]) pqe_common_problemo("cplex_set_pb_data_LPi");
	}

	/* ----------------- */
	/* Colonnes n a 2n-2 */
	/* ----------------- */

	for (k = nbvar; k <= 2 * nbvar - 2; k++)
	{
		matbeg[k] = compteur;
		matcnt[k] = 3;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k - nbvar + 1;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k;
		matval[compteur]   = 1.0;
		matind[compteur++] = 2 * nbvar + extra_ligne - 1;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
        lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'B';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		lb[j] = -CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	if (prime)
	{
		sense[0] = 'E';
		rhs[0]   = (sens == 1 ? 1.0 : 0.0);
	}
	sense[extra_ligne] = 'L';
	rhs[extra_ligne]   = (float) capacite;
	for (k = extra_ligne + 1; k <= 2 * nbvar + extra_ligne - 2; k++)
	{
		sense[k] = 'L';
		rhs[k]   = 0.0;
	}
	sense[2 * nbvar + extra_ligne - 1] = 'G';
	rhs[2 * nbvar + extra_ligne - 1]   = meilleur_primal;

} /* pqe_cplex_set_pb_data_LPi */

/*****************************************************************************/

void pqe_cplex_set_pb_data_max_LPi (bool prime, int i, int N, int COLSPACE, int ROWSPACE, double_1D sur_phi,
                                    double_1D sous_phi, double_1D rlin, double_2D rqua, 
                                    double_1D weight, 
                                    int CAPA, int *p_objsen, double_1D coef_obj, double_1D rhs, 
                                    char_1D sense, int_1D matbeg, int_1D matcnt, int_1D matind, 
                                    double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal, int sens)
{
	/* ------------------------------------------------------------------------ */
	/* prime indique s'il s'agit de LP'i ou de LPi : xi = 1 ou 0 ajout\8Ee ou non */
	/* ------------------------------------------------------------------------ */

	int    k, j, l, cpt, compteur;
	int    extra_ligne, controle1, controle2;

	*p_objsen = (sens == 1 ? CPX_MAX : CPX_MIN);  /* Minimisation ou maximisation */
	extra_ligne = (prime ? 1 : 0);


	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;

	for (j = 1; j <= i; j++) coef_obj[cpt++] = 0.0;
	for (j = i + 1; j <= N; j++) {
		coef_obj[cpt++] = rqua[i][j];
	}
	for (j = 1; j <= N - 1; j++) coef_obj[cpt++] = 0.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* -------------------- */
	/* Colonnes 0 a nbvar-1 */
	/* -------------------- */

	for (k = 0; k <= N - 1; k++)
	{
		controle1 = compteur;
		matbeg[k] = compteur;
		matcnt[k] = (k == (N-1) ? k + 2 : k + 3);
		if ((prime) && (k == i - 1)) {
			matcnt[k]++;
			matval[compteur] = 1.0;
			matind[compteur++] = 0;
		}
		matval[compteur]   = weight[k+1];
		matind[compteur++] = extra_ligne;
		if (k < N - 1)
		{
			matval[compteur]   = sous_phi[k+1] - sur_phi[k+1];
			matind[compteur++] = k + 1 + extra_ligne; 
		}
                
		for (l = 1; l <= k; l++) {
			matval[compteur]   = -rqua[l][k+1];
			matind[compteur++] = l + N + extra_ligne - 1;
		}
		matval[compteur]   = rlin[k+1] + (k < N - 1 ? sous_phi[k+1] : 0.0);
		matind[compteur++] = 2 * N + extra_ligne - 1;
		controle2 = compteur;
		if (controle2 - controle1 != matcnt[k]) pqe_common_problemo("cplex_set_pb_data_maxLPi");
	}

	/* ----------------- */
	/* Colonnes n a 2n-2 */
	/* ----------------- */

	for (k = N; k <= 2 * N - 2; k++)
	{
		matbeg[k] = compteur;
		matcnt[k] = 3;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k - N + 1;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k;
		matval[compteur]   = 1.0;
		matind[compteur++] = 2 * N + extra_ligne - 1;
	}

    /* ------ */
    /* Bornes */
    /* ------ */

    for (j = 0; j <= N - 1; j++)
    {
        lb[j] = 0.0;
        ub[j] = 1.0;
        ctype[j] = 'C';
    }

    for (j = N; j <= 2 * N - 2; j++)
    {
        lb[j] = -CPX_INFBOUND;
        ub[j] = CPX_INFBOUND;
        ctype[j] = 'C';
    }

    /* ---------------- */
    /* Membre de droite */
    /* ---------------- */

    if (prime) {
        sense[0] = 'E';
        rhs[0]   = (sens == 1 ? 1.0 : 0.0);
    }
    
    sense[extra_ligne] = 'L';
    rhs[extra_ligne]   = (double) CAPA;
    
    for (k = extra_ligne + 1; k <= N + extra_ligne - 1; k++) {
        sense[k] = 'L';
        rhs[k]   = 0.0;
    }
    
    for (k = N + extra_ligne; k <= 2 * N + extra_ligne - 2; k++) {
        sense[k] = 'L';
        rhs[k]   = - sous_phi[k + 1 - N - extra_ligne];
    }
        
    sense[2 * N + extra_ligne - 1] = 'G';
    rhs[2 * N + extra_ligne - 1]   = meilleur_primal;

} /* pqe_cplex_set_pb_data_max_LPi */

/*****************************************************************************/

void pqe_cplex_set_pb_data_lin3 (int N, int COLSPACE, int ROWSPACE, double_1D sur_phi,
                                    double_1D sous_phi, double_1D rlin, double_2D rqua, 
                                    double_1D weight, 
                                    int CAPA, int *p_objsen, double_1D coef_obj, double_1D rhs, 
                                    char_1D sense, int_1D matbeg, int_1D matcnt, int_1D matind, 
                                    double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal, int sens)
{
	int    k, j, l, cpt, compteur;
	int    controle1, controle2;

	*p_objsen = (sens == 1 ? CPX_MAX : CPX_MIN);  /* Minimisation ou maximisation */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	cpt = 0;

	for (j = 1; j <= N; j++) coef_obj[cpt++] = rlin[j] + (j < N ? sous_phi[j] : 0.0);
	for (j = 1; j <= N - 1; j++) coef_obj[cpt++] = 1.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* -------------------- */
	/* Colonnes 0 a nbvar-1 */
	/* -------------------- */

	for (k = 0; k <= N - 1; k++)
	{
		controle1 = compteur;
		matbeg[k] = compteur;
		matcnt[k] = (k == (N-1) ? k + 1 : k + 2);
		matval[compteur]   = weight[k+1];
		matind[compteur++] = 0;
		if (k < N - 1)
		{
			matval[compteur]   = sous_phi[k+1] - sur_phi[k+1];
			matind[compteur++] = k + 1; 
		}
                
		for (l = 1; l <= k; l++) {
			matval[compteur]   = -rqua[l][k+1];
			matind[compteur++] = l + N - 1;
		}
		controle2 = compteur;
		if (controle2 - controle1 != matcnt[k]) pqe_common_problemo("cplex_set_pb_data_maxLPi");
	}

	/* ----------------- */
	/* Colonnes n a 2n-2 */
	/* ----------------- */

	for (k = N; k <= 2 * N - 2; k++)
	{
		matbeg[k] = compteur;
		matcnt[k] = 2;
		matval[compteur]   = 1.0;
		matind[compteur++] = k - N + 1;
		matval[compteur]   = 1.0;
		matind[compteur++] = k;
	}

    /* ------ */
    /* Bornes */
    /* ------ */

    for (j = 0; j <= N - 1; j++)
    {
        lb[j] = 0.0;
        ub[j] = 1.0;
        ctype[j] = 'B';
    }

    for (j = N; j <= 2 * N - 2; j++)
    {
        lb[j] = -CPX_INFBOUND;
        ub[j] = CPX_INFBOUND;
        ctype[j] = 'C';
    }

    /* ---------------- */
    /* Membre de droite */
    /* ---------------- */

    sense[0] = 'L';
    rhs[0]   = (double) CAPA;
    
    for (k = 1; k <= N - 1; k++) {
        sense[k] = 'L';
        rhs[k]   = 0.0;
    }
    
    for (k = N ; k <= 2 * N - 2; k++) {
        sense[k] = 'L';
        rhs[k]   = - sous_phi[k + 1 - N ];
    }
      

} /* pqe_cplex_set_pb_data_lin3 */


/*****************************************************************************/

void pqe_cplex_set_pb_data_min_LPi (bool prime, int i, int N, int COLSPACE, int ROWSPACE, double_1D sur_phi,
                                    double_1D sous_phi, double_1D rlin, double_2D rqua, 
                                    double_1D weight, 
                                    int CAPA, int *p_objsen, double_1D coef_obj, double_1D rhs, 
                                    char_1D sense, int_1D matbeg, int_1D matcnt, int_1D matind, 
                                    double_1D matval, double_1D lb, double_1D ub, 
                                    char_1D ctype, float meilleur_primal, int sens)
{
    /* ------------------------------------------------------------------------ */
    /* prime indique s'il s'agit de LP'i ou de LPi : xi = 1 ou 0 ajout\8Ee ou non */
    /* ------------------------------------------------------------------------ */

    int    k, j, l, cpt, compteur;
    int    extra_ligne, controle1, controle2;
    double amax; 
    
    for (amax = weight[1], k = 2; k <= N; k++) amax = (weight[k] > amax ? weight[k] : amax);
    *p_objsen = (sens == 1 ? CPX_MAX : CPX_MIN);  /* Minimisation ou maximisation */
    extra_ligne = (prime ? 1 : 0);

    /* ----------------- */
    /* Fonction objectif */
    /* ----------------- */

    cpt = 0;

    for (j = 1; j <= i; j++) coef_obj[cpt++] = 0.0;
    for (j = i + 1; j <= N; j++) {
        coef_obj[cpt++] = rqua[i][j];
    }
    for (j = 1; j <= N - 1; j++) coef_obj[cpt++] = 0.0;

	/* ---------------------------------------------------------------------------------------- */
	/* Matrice des contraintes                                                                  */
	/* Seuls les coeff non nuls sont stockes, l'un apres l'autre, dans le tableau matval,       */
	/* la matrice des contraintes etant lue colonne par colonne. matcnt[j] est le nombre de     */
	/* coeff non nuls de la colonne j et matbeg[j] l'indice dans matval du premier element      */
	/* non nul de la colonne j. Si on considere un coefficient stocke dans matval[k], matind[k] */
	/* indique le numero de la ligne correspondant au coefficient.                              */
	/* ---------------------------------------------------------------------------------------- */

	compteur = 0;  /* Indice dans matval et matind */

	/* ---------------- */
	/* Colonnes 0 a N-1 */
	/* ---------------- */

	for (k = 0; k <= N - 1; k++)
	{
		controle1 = compteur;
		matbeg[k] = compteur;
		matcnt[k] = (k == (N-1) ? k + 2 + N: k + 3 + N);
		if ((prime) && (k == i - 1)) {
			matcnt[k]++;
			matval[compteur] = 1.0;
			matind[compteur++] = 0;
		}
		matval[compteur]   = weight[k+1];
		matind[compteur++] = extra_ligne;
		if (k < N - 1)
		{
			matval[compteur]   = sous_phi[k+1] - sur_phi[k+1];
			matind[compteur++] = k + 1 + extra_ligne; 
		}
                
		for (l = 1; l <= k; l++) {
			matval[compteur]   = -rqua[l][k+1];
			matind[compteur++] = l + N + extra_ligne - 1;
		}
		matval[compteur]   = rlin[k+1] + (k < N - 1 ? sous_phi[k+1] : 0.0);
		matind[compteur++] = 2 * N + extra_ligne - 1;
        
        for (l = 1; l <= N; l++) {
            matval[compteur]   = (l == (k+1) ? amax - weight[k+1] : weight[k+1]);
            matind[compteur++] = extra_ligne + 2 * N + l - 1;
        }
        
		controle2 = compteur;
		if (controle2 - controle1 != matcnt[k]) pqe_common_problemo("cplex_set_pb_data_minLPi");
	}

	/* ----------------- */
	/* Colonnes n a 2n-2 */
	/* ----------------- */

	for (k = N; k <= 2 * N - 2; k++)
	{
		matbeg[k] = compteur;
		matcnt[k] = 3;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k - N + 1;
		matval[compteur]   = 1.0;
		matind[compteur++] = extra_ligne + k;
		matval[compteur]   = 1.0;
		matind[compteur++] = 2 * N + extra_ligne - 1;
	}

    /* ------ */
    /* Bornes */
    /* ------ */

    for (j = 0; j <= N - 1; j++)
    {
        lb[j] = 0.0;
        ub[j] = 1.0;
        ctype[j] = 'C';
    }

    for (j = N; j <= 2 * N - 2; j++)
    {
        lb[j] = -CPX_INFBOUND;
        ub[j] = CPX_INFBOUND;
        ctype[j] = 'C';
    }

    /* ---------------- */
    /* Membre de droite */
    /* ---------------- */

    if (prime) {
        sense[0] = 'E';
        rhs[0]   = (sens == 1 ? 1.0 : 0.0);
    }
    
    sense[extra_ligne] = 'L';
    rhs[extra_ligne]   = (double) CAPA;
    
    for (k = extra_ligne + 1; k <= N + extra_ligne - 1; k++) {
        sense[k] = 'L';
        rhs[k]   = 0.0;
    }
    
    for (k = N + extra_ligne; k <= 2 * N + extra_ligne - 2; k++) {
        sense[k] = 'L';
        rhs[k]   = - sous_phi[k + 1 - N - extra_ligne];
    }
        
    sense[2 * N + extra_ligne - 1] = 'G';
    rhs[2 * N + extra_ligne - 1]   = meilleur_primal;

    for (k = 2*N + extra_ligne; k <= 3*N + extra_ligne - 1; k++) {
        sense[k] = 'G';
        rhs[k]   = CAPA - weight[k - 2*N - extra_ligne + 1];
    }

} /* pqe_cplex_set_pb_data_min_LPi */

/*****************************************************************************/

void pqe_cplex_set_pb_data_lin1 (int nbvar, float_1D coef_c1z, st_tri_1D ordre_poids,
                                     float_1D q_lin, float_2D q_qua, int capacite,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement, int_1D pr_indices,
                                     double_1D pr_value)
{
	char part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE];
	int  i, j, indice, minij, maxij, nb;
	double val;

	pqe_common_intochar(part1,nbvar);
    pqe_common_intochar(part2,(int) (densite*100));
    pqe_common_intochar(part3,seed);
	
	strcpy(probname,"newlin");
	strcat(probname,part1);
	strcat(probname,"D");
	strcat(probname,part2);
	strcat(probname,"G");
	strcat(probname,part3);

	*numcols_p = 2 * nbvar - 1;  /* zn est inutile, toujours a zero */
	*numrows_p = 2 * nbvar - 1;  /* c1zn et c2zn sont inutiles      */
	*objsen_p = CPX_MAX;              /* Maximisation                    */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	for (i = 0; i < nbvar; i++) obj[i] = q_lin[ordre_poids[i + 1].indice];
	for (i = nbvar; i <= 2 * nbvar - 2; i++) obj[i] = 1.0;

	/* ----------------------- */
	/* Matrice des contraintes */
	/* ----------------------- */

	for (j = 0; j <= nbvar - 2; j++)
	{
		indice = ordre_poids[j + 1].indice;

		matbeg[j] = j * (j + 3) / 2;
		matcnt[j] = 2 + j;

		matind[matbeg[j]] = 0;
		matval[matbeg[j]] = a[indice];                   /* Contrainte de capacite */

		matind[matbeg[j] + 1] = j + 1;
		matval[matbeg[j] + 1] = - coef_c1z[indice];        /* Contrainte c1z */

		for (i = 2; i < matcnt[j]; i++)
		{
			/* Contrainte c2z */

			minij = (ordre_poids[i-1].indice < indice ? ordre_poids[i-1].indice : indice);
			maxij = (ordre_poids[i-1].indice > indice ? ordre_poids[i-1].indice : indice); 

			matind[matbeg[j] + i] = nbvar + i - 2;
			matval[matbeg[j] + i] = - q_qua[minij][maxij];
		}
	}

	/* n-1 eme colonne */

	j = nbvar - 1;
	indice = ordre_poids[nbvar].indice;
	matbeg[nbvar - 1] = (nbvar - 1) * (nbvar + 2) / 2;
	matcnt[nbvar - 1] = nbvar;
	matind[matbeg[nbvar - 1]] = 0;
	matval[matbeg[nbvar - 1]] = a[ordre_poids[nbvar].indice];

	for (i = 1; i < matcnt[j]; i++)
	{
		/* Contrainte c2z */

		minij = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
		maxij = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

		matind[matbeg[j] + i] = nbvar + i - 1;
		matval[matbeg[j] + i] = - q_qua[minij][maxij];
	}

	/* Colonnes de n a 2n-2 */

	nb = (nbvar * nbvar + 3 * nbvar - 2) / 2;

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		matbeg[j] = nb + 2 * (j - nbvar);
		matcnt[j] = 2;

		matind[matbeg[j]] = j - nbvar + 1;
		matval[matbeg[j]] = 1.0;

		matind[matbeg[j] + 1] = j;
		matval[matbeg[j] + 1] = 1.0;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'I';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		/* lb[j] = 0.0; */
		lb[j] = -CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = (float) capacite;

	for (i = 1; i <= 2 * nbvar - 2; i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}

	/* ----------------------------------------------------------------------------- */
	/* Preparation du vecteur contenant la meilleure solution admissible deja connue */
	/* Preparation des vecteurs colindex et priority pour l'ordre de branchement     */
	/* ----------------------------------------------------------------------------- */

	for (i = 1; i <= nbvar; i++)
	{
		indice = ordre_poids[i].indice;
		pr_indices[i - 1] = i - 1;
		pr_value[i - 1] = (double) b_primal[indice];
		pr_indices[i + nbvar - 1] = i + nbvar - 1;
		val = 0;
		for (j = i + 1; j <= nbvar; j++)
		{
			minij = (ordre_poids[j].indice < indice ? ordre_poids[j].indice : indice);
			maxij = (ordre_poids[j].indice > indice ? ordre_poids[j].indice : indice); 
			val += q_qua[minij][maxij] * b_primal[minij] * b_primal[maxij];	
		}
		/* pr_value[i + nbvar - 1] = val; */
		pr_value[i + nbvar - 1] = CPX_INFBOUND;
	}

	for (i = 0; i <= nbvar - 1; i++)
	{
		indice = ordre_poids[i + 1].indice;
		/* ----------------- */
		/* Variables x1 a xn */
		/* ----------------- */
		colindex[i] = i;

		switch (type_branchement)
		{
			case BRAN_NAT :   priority[i] = nbvar - i + 60;
			case BRAN_PRI_1 : priority[i] = (b_primal[indice] == 1 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_PRI_0 : priority[i] = (b_primal[indice] == 0 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_RAT   : priority[i] = nbvar - i + 60; /* Suppose ordre_var == VAR_RAT_DEC */
		}

		cprim[i] = (double) b_primal[indice];
		rprim[i] = 0.0;
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		/* ------------------- */
		/* Variables z1 a zn-1 */
		/* ------------------- */
		indice = ordre_poids[i - nbvar + 1].indice;
		rprim[i] = 0.0;
		cprim[i] = 0.0;
	}

} /* pqe_cplex_set_pb_data_lin1 */

/*****************************************************************************/

void pqe_cplex_set_pb_data_lin2 (int nbvar, int COLSPACE, int ROWSPACE, float_1D coef_c1z,
                                     float_1D sous_phi, st_tri_1D ordre_poids,
                                     float_1D q_lin, float_2D q_qua, int capacite,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement, int_1D pr_indices,
                                     double_1D pr_value)
{
	char part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE];
	int  i, decalage, j, indice, minij, maxij, nb;
	double val;
	bool   formule_brute = VRAI; /* vrai si zi <= phi_i(x) - sous-phi_i       */
	                             /* faux si zi <= phi_i(x) - sous-phi_i * x_i */

	pqe_common_intochar(part1,nbvar);
    pqe_common_intochar(part2,(int) (densite*100));
    pqe_common_intochar(part3,seed);
	strcpy(probname,"linMM");
	strcat(probname,part1);
	strcat(probname,"D");
	strcat(probname,part2);
	strcat(probname,"G");
	strcat(probname,part3);

	*numcols_p = COLSPACE;  /* zn est inutile, toujours a zero */
	*numrows_p = ROWSPACE;  /* c1zn et c2zn sont inutiles      */
	*objsen_p = CPX_MAX;              /* Maximisation                    */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	for (i = 0; i < nbvar; i++)
	{
		indice = ordre_poids[i + 1].indice;
		obj[i] = q_lin[indice];
		if (i != nbvar - 1) obj[i] += sous_phi[indice];
	}
	for (i = nbvar; i <= 2 * nbvar - 2; i++) obj[i] = 1.0;

	/* ----------------------- */
	/* Matrice des contraintes */
	/* ----------------------- */

	for (j = 0; j <= nbvar - 2; j++)
	{
		indice = ordre_poids[j + 1].indice;

		matbeg[j] = (formule_brute ? j * (j + 3) / 2 : j * (j + 5) / 2);
		matcnt[j] = (formule_brute ? 2 + j           : 3 + j);

		matind[matbeg[j]] = 0;
		matval[matbeg[j]] = a[indice];                   /* Contrainte de capacite */

		matind[matbeg[j] + 1] = j + 1;
		matval[matbeg[j] + 1] = - coef_c1z[indice];        /* Contrainte c1z */

		if (!(formule_brute))
		{
			matind[matbeg[j] + 2] = j + nbvar;
			matval[matbeg[j] + 2] = sous_phi[indice];
		}

		decalage = (formule_brute ? 0 : 1);

		for (i = 2 + decalage; i < matcnt[j]; i++)
		{
			/* Contrainte c2z */

			minij = (ordre_poids[i-decalage-1].indice < indice ? ordre_poids[i-decalage-1].indice : indice);
			maxij = (ordre_poids[i-decalage-1].indice > indice ? ordre_poids[i-decalage-1].indice : indice); 

			matind[matbeg[j] + i] = nbvar + i - decalage - 2;
			matval[matbeg[j] + i] = - q_qua[minij][maxij];
		}
	}

	/* n-1 eme colonne */

	j = nbvar - 1;
	indice = ordre_poids[nbvar].indice;
	matbeg[nbvar - 1] = (nbvar - 1) * (nbvar + 2) / 2;
	if (!(formule_brute)) matbeg[nbvar - 1] += nbvar - 1;
	matcnt[nbvar - 1] = nbvar;
	matind[matbeg[nbvar - 1]] = 0;
	matval[matbeg[nbvar - 1]] = a[ordre_poids[nbvar].indice];

	for (i = 1; i < matcnt[j]; i++)
	{
		/* Contrainte c2z */

		minij = (ordre_poids[i].indice < indice ? ordre_poids[i].indice : indice);
		maxij = (ordre_poids[i].indice > indice ? ordre_poids[i].indice : indice); 

		matind[matbeg[j] + i] = nbvar + i - 1;
		matval[matbeg[j] + i] = - q_qua[minij][maxij];
	}

	/* Colonnes de n a 2n-2 */

	nb = (nbvar * nbvar + 3 * nbvar - 2) / 2;
	if (!(formule_brute)) nb += (nbvar - 1);
  

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		matbeg[j] = nb + 2 * (j - nbvar);
		matcnt[j] = 2;

		matind[matbeg[j]] = j - nbvar + 1;
		matval[matbeg[j]] = 1.0;

		matind[matbeg[j] + 1] = j;
		matval[matbeg[j] + 1] = 1.0;
	}

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'B';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		/* lb[j] = 0.0; */
		lb[j] = -CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = (float) capacite;

	for (i = 1; i <= nbvar - 1; i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}
	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		sense[i] = 'L';
		rhs[i]   = (formule_brute ? - sous_phi[ordre_poids[i - nbvar + 1].indice] : 0.0);
	}

	/* ----------------------------------------------------------------------------- */
	/* Preparation du vecteur contenant la meilleure solution admissible deja connue */
	/* Preparation des vecteurs colindex et priority pour l'ordre de branchement     */
	/* ----------------------------------------------------------------------------- */

	for (i = 1; i <= nbvar; i++)
	{
		indice = ordre_poids[i].indice;
		pr_indices[i - 1] = i - 1;
		pr_value[i - 1] = (double) b_primal[indice];
		if (i < nbvar) pr_indices[i + nbvar - 1] = i + nbvar - 1;
		val = 0.0;
		for (j = i + 1; j <= nbvar; j++)
		{
			minij = (ordre_poids[j].indice < indice ? ordre_poids[j].indice : indice);
			maxij = (ordre_poids[j].indice > indice ? ordre_poids[j].indice : indice); 
			val  += q_qua[minij][maxij] * b_primal[minij] * b_primal[maxij];	
		}
		/* if (i < nbvar) pr_value[i + nbvar - 1] = (val - sous_phi[indice]) * b_primal[indice]; */
		if (i < nbvar) pr_value[i + nbvar - 1] = CPX_INFBOUND; 
	}
	pr_indices[2*nbvar-1]=2*nbvar-1;
	pr_value[2*nbvar-1]=0.0;

	for (i = 0; i <= nbvar - 1; i++)
	{
		indice = ordre_poids[i + 1].indice;
		/* ----------------- */
		/* Variables x1 a xn */
		/* ----------------- */
		colindex[i] = i;

		switch (type_branchement)
		{
			case BRAN_NAT :   priority[i] = nbvar - i + 60;
			case BRAN_PRI_1 : priority[i] = (b_primal[indice] == 1 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_PRI_0 : priority[i] = (b_primal[indice] == 0 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_RAT   : priority[i] = nbvar - i + 60; /* Suppose ordre_var == VAR_RAT_DEC */
		}

		cprim[i] = (double) b_primal[indice];
		rprim[i] = 0.0;
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		/* ------------------- */
		/* Variables z1 a zn-1 */
		/* ------------------- */
		indice = ordre_poids[i - nbvar + 1].indice;
		rprim[i] = 0.0;
		cprim[i] = 0.0;
	}

} /* pqe_cplex_set_pb_data_lin2 */
	
/*****************************************************************************/

void pqe_cplex_set_pb_data_lin2_stocha (int nbvar, int COLSPACE, int ROWSPACE, float_1D coef_c1z,
                                     float_1D sous_phi, st_tri_1D ordre_poids,
                                     float_1D q_lin, float_2D q_qua, float_1D cout, float CAP,
                                     float_1D moyenne, float k,
                                     char_1D probname, int *numcols_p, int *numrows_p, 
                                     int *objsen_p, double_1D obj, double_1D rhs, char_1D sense, 
                                     int_1D matbeg, int_1D matcnt, int_1D matind, double_1D matval, 
                                     double_1D lb, double_1D ub, char_1D ctype,
                                     double_1D cprim, double_1D rprim, bival_1D b_primal,
                                     int_1D colindex, int_1D priority, int type_branchement, int_1D pr_indices,
                                     double_1D pr_value)
{
    /* ---------------------------------------------------------- */
	/* La derniere variable est la constante k*k : on la fixe a 1 */
    /* ---------------------------------------------------------- */

    char part1[LONGCHAINE];
	int  i, decalage, j, cpt, indice, minij, maxij;
	double val;

	pqe_common_intochar(part1,nbvar);
	strcpy(probname,"Lin2_stocha");
	strcat(probname,part1);

	*numcols_p = COLSPACE;  /* zn est inutile, toujours a zero */
	*numrows_p = ROWSPACE;  /* c1zn et c2zn sont inutiles      */
	*objsen_p = CPX_MAX;              /* Maximisation                    */

	/* ----------------- */
	/* Fonction objectif */
	/* ----------------- */

	for (i = 0; i < nbvar; i++)
	{
		indice = ordre_poids[i + 1].indice;
		obj[i] = q_lin[indice];
		if (i != nbvar - 1) obj[i] += sous_phi[indice];
	}
	for (i = nbvar; i <= 2 * nbvar - 2; i++) obj[i] = 1.0;
    obj[2*nbvar-1] = k * k; /* la constante */

	/* ----------------------- */
	/* Matrice des contraintes */
	/* ----------------------- */

    cpt = 0;
	for (j = 0; j <= nbvar - 2; j++)
	{
		indice = ordre_poids[j + 1].indice;

		matbeg[j] = cpt;
		matcnt[j] = 3 + j;

		matind[cpt] = 0;
		matval[cpt++] = cout[indice];                   /* Contrainte de capacite */

        matind[cpt] = 1;
        matval[cpt++] = moyenne[indice];                /* Contrainte 2 */

		matind[cpt] = j + 2;
		matval[cpt++] = - coef_c1z[indice];        /* Contrainte c1z */

		for (i = 2; i < matcnt[j]; i++)
		{
            /* -------------- */
			/* Contrainte c2z */
            /* -------------- */

			minij = (ordre_poids[i-1].indice < indice ? ordre_poids[i-1].indice : indice);
			maxij = (ordre_poids[i-1].indice > indice ? ordre_poids[i-1].indice : indice);

			matind[cpt] = nbvar + i - 1;
			matval[cpt++] = - q_qua[minij][maxij];
		}
	}

    /* --------------- */
	/* n-1 eme colonne */
    /* --------------- */

	j = nbvar - 1;
	indice = ordre_poids[nbvar].indice;
	matbeg[j] = cpt;
	matcnt[j] = 1 + nbvar;
	matind[cpt] = 0;
	matval[cpt++] = cout[indice];
    matind[cpt] = 1;
    matval[cpt++] = moyenne[indice];

	for (i = 2; i < matcnt[j]; i++)
	{
        /* -------------- */
		/* Contrainte c2z */
        /* -------------- */

		minij = (ordre_poids[i-1].indice < indice ? ordre_poids[i-1].indice : indice);
		maxij = (ordre_poids[i-1].indice > indice ? ordre_poids[i-1].indice : indice);

		matind[cpt] = nbvar + i - 1;
		matval[cpt++] = - q_qua[minij][maxij];
	}

    /* -------------------- */
	/* Colonnes de n a 2n-2 */
    /* -------------------- */

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		matbeg[j] = cpt;
		matcnt[j] = 2;

		matind[cpt] = j - nbvar + 2;
		matval[cpt++] = 1.0;

		matind[cpt] = j + 1;
		matval[cpt++] = 1.0;
	}

    /* --------------------------- */
    /* Colonne 2n-1 : la constante */
    /* --------------------------- */
    matbeg[2 * nbvar - 1] = cpt;
    matcnt[2 * nbvar - 1] = 1;
    matind[cpt]           = 2 * nbvar;
    matval[cpt]           = 1.0;

	/* ------ */
	/* Bornes */
	/* ------ */

	for (j = 0; j <= nbvar - 1; j++)
	{
		lb[j] = 0.0;
		ub[j] = 1.0;
		ctype[j] = 'B';
	}

	for (j = nbvar; j <= 2 * nbvar - 2; j++)
	{
		lb[j] = - CPX_INFBOUND;
		ub[j] = CPX_INFBOUND;
		ctype[j] = 'C';
	}
    lb[2*nbvar-1] = 0.0;
    ub[2*nbvar-1] = 1.0;
    ctype[2*nbvar-1] = 'B';

	/* ---------------- */
	/* Membre de droite */
	/* ---------------- */

	sense[0] = 'L';
	rhs[0]   = CAP;

    sense[1] = 'G';
    rhs[1]   = k;

	for (i = 2; i <= nbvar; i++)
	{
		sense[i] = 'L';
		rhs[i]   = 0.0;
	}
	for (i = nbvar + 1; i <= 2 * nbvar - 1; i++)
	{
		sense[i] = 'L';
		rhs[i]   = - sous_phi[ordre_poids[i - nbvar].indice];
	}
    sense[2*nbvar] = 'E';
    rhs[2*nbvar]   = 1.0;

	/* ----------------------------------------------------------------------------- */
	/* Preparation du vecteur contenant la meilleure solution admissible deja connue */
	/* Preparation des vecteurs colindex et priority pour l'ordre de branchement     */
	/* ----------------------------------------------------------------------------- */

	for (i = 1; i <= nbvar; i++)
	{
		indice = ordre_poids[i].indice;
		pr_indices[i - 1] = i - 1;
		pr_value[i - 1] = (double) b_primal[indice];
		if (i < nbvar) pr_indices[i + nbvar - 1] = i + nbvar - 1;
		val = 0.0;
		for (j = i + 1; j <= nbvar; j++)
		{
			minij = (ordre_poids[j].indice < indice ? ordre_poids[j].indice : indice);
			maxij = (ordre_poids[j].indice > indice ? ordre_poids[j].indice : indice); 
			val  += q_qua[minij][maxij] * b_primal[minij] * b_primal[maxij];	
		}
		/* if (i < nbvar) pr_value[i + nbvar - 1] = (val - sous_phi[indice]) * b_primal[indice]; */
		if (i < nbvar) pr_value[i + nbvar - 1] = CPX_INFBOUND;
	}
	pr_indices[2*nbvar-1]=2*nbvar-1;
	pr_value[2*nbvar-1]=1.0;

	for (i = 0; i <= nbvar - 1; i++)
	{
		indice = ordre_poids[i + 1].indice;
		/* ----------------- */
		/* Variables x1 a xn */
		/* ----------------- */
		colindex[i] = i;

		switch (type_branchement)
		{
			case BRAN_NAT :   priority[i] = nbvar - i + 60;
			case BRAN_PRI_1 : priority[i] = (b_primal[indice] == 1 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_PRI_0 : priority[i] = (b_primal[indice] == 0 ? 2 *nbvar - i + 60 : nbvar -i + 60);
			case BRAN_RAT   : priority[i] = nbvar - i + 60; /* Suppose ordre_var == VAR_RAT_DEC */
		}

		cprim[i] = (double) b_primal[indice];
		rprim[i] = 0.0;
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		/* ------------------- */
		/* Variables z1 a zn-1 */
		/* ------------------- */
		indice = ordre_poids[i - nbvar + 1].indice;
		rprim[i] = 0.0;
		cprim[i] = 0.0;
	}
    cprim[2*nbvar-1] = 1.0;
    rprim[2*nbvar-1] = 0.0;

} /* pqe_cplex_set_pb_data_lin2_stocha */
	
/*****************************************************************************/

void pqe_cplex_set_pb_names_demi(int nbvar, st_tri_1D ordre_poids, char_2D colname, char_1D colnamestore,
                              char_2D rowname, char_1D rownamestore, unsigned *p_colnamespace, 
                              unsigned *p_rownamespace)
{
	char chaine[16];
	int i, j, cpt;

	/* ---------------------------- */
	/* Nom de variables, le nom ;-) */
	/* ---------------------------- */

	cpt = 0;
	for (i = 0; i < nbvar; i++)
	{
		colname[i] = &colnamestore[cpt];
		colnamestore[cpt++] = 'x';
		pqe_common_intochar(chaine, ordre_poids[i + 1].indice);
		j = 0;
		while (chaine[j] != '\0') colnamestore[cpt++] = chaine[j++];
		colnamestore[cpt++] = '\0';
	}

	for (i = nbvar; i <= 2 * nbvar - 1; i++)
	{
		colname[i] = &colnamestore[cpt];
		colnamestore[cpt++] = 'z';
		pqe_common_intochar(chaine, ordre_poids[i - nbvar + 1].indice);
		j = 0;
		while (chaine[j] != '\0') colnamestore[cpt++] = chaine[j++];
		colnamestore[cpt++] = '\0';
	}

	*p_colnamespace = cpt;

	/* -------------------------- */
	/* Nom de contraintes, le nom */
	/* -------------------------- */

	strcpy(rownamestore, "capacite"); 
	rowname[0] = &rownamestore[0];
	cpt = 9;
	for (i = 1; i <= nbvar; i++)
	{
		/* -------------- */
		/* Contraintes c1 */
		/* -------------- */

		rowname[i] = &rownamestore[cpt];
		rownamestore[cpt++] = 'c'; 
		rownamestore[cpt++] = '1'; 
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, i);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, ordre_poids[i].indice);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '\0';
	}

	for (i = nbvar + 1; i <= 2 * nbvar; i++)
	{
		/* -------------- */
		/* Contraintes c2 */
		/* -------------- */

		rowname[i] = &rownamestore[cpt];
		rownamestore[cpt++] = 'c'; 
		rownamestore[cpt++] = '2'; 
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, i - nbvar);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, ordre_poids[i - nbvar].indice);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '\0';
	}

	*p_rownamespace = cpt;

} /* pqe_cplex_set_pb_names_demi */

/*****************************************************************************/

void pqe_cplex_set_pb_names(int nbvar, st_tri_1D ordre_poids, char_2D colname, char_1D colnamestore,
                              char_2D rowname, char_1D rownamestore, unsigned *p_colnamespace, 
                              unsigned *p_rownamespace)
{
	char chaine[16];
	int i, j, cpt;

	/* ---------------------------- */
	/* Nom de variables, le nom ;-) */
	/* ---------------------------- */

	cpt = 0;
	for (i = 0; i < nbvar; i++)
	{
		colname[i] = &colnamestore[cpt];
		colnamestore[cpt++] = 'x';
		pqe_common_intochar(chaine, ordre_poids[i + 1].indice);
		j = 0;
		while (chaine[j] != '\0') colnamestore[cpt++] = chaine[j++];
		colnamestore[cpt++] = '\0';
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		colname[i] = &colnamestore[cpt];
		colnamestore[cpt++] = 'z';
		pqe_common_intochar(chaine, ordre_poids[i - nbvar + 1].indice);
		j = 0;

		while (chaine[j] != '\0') colnamestore[cpt++] = chaine[j++];
		colnamestore[cpt++] = '\0';
	}

	*p_colnamespace = cpt;

	/* -------------------------- */
	/* Nom de contraintes, le nom */
	/* -------------------------- */

	strcpy(rownamestore, "capacite"); 
	rowname[0] = &rownamestore[0];
	cpt = 9;
	for (i = 1; i <= nbvar - 1; i++)
	{
		/* -------------- */
		/* Contraintes c1 */
		/* -------------- */

		rowname[i] = &rownamestore[cpt];
		rownamestore[cpt++] = 'c'; 
		rownamestore[cpt++] = '1'; 
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, i);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, ordre_poids[i].indice);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '\0';
	}

	for (i = nbvar; i <= 2 * nbvar - 2; i++)
	{
		/* -------------- */
		/* Contraintes c2 */
		/* -------------- */

		rowname[i] = &rownamestore[cpt];
		rownamestore[cpt++] = 'c'; 
		rownamestore[cpt++] = '2'; 
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, i - nbvar + 1);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '_'; 
		pqe_common_intochar(chaine, ordre_poids[i - nbvar + 1].indice);
		j = 0;
		while (chaine[j] != '\0') rownamestore[cpt++] = chaine[j++];
		rownamestore[cpt++] = '\0';
	}

	*p_rownamespace = cpt;

} /* pqe_cplex_set_pb_names */

/*****************************************************************************/

int pqe_cplex_convertit(int i, int j)
{
	/* ------------------------------------------------- */
	/* pre : i < j                                       */
	/* Retourne le numero de colonne de la variable wi_j */
	/* ------------------------------------------------- */

	return (nbvar * i - i * (i - 1) / 2 + j - i - 1);
} /* pqe_cplex_convertit */

/*****************************************************************************/

void pqe_cplex_prep_contraintes_produits (int nbvar, int_1D a, int capacite, int *p_ar_nzcount, double_1D ar_rhs,
                                            char_1D ar_sense, int_1D ar_matbeg, int_1D ar_matind, double_1D ar_matval, 
                                            char_2D newrowname)
{
	int i, j, cpt, lig;

	*p_ar_nzcount = nbvar * (3 * nbvar - 1);

	for (i = 0; i < nbvar; i++) 
	{
		ar_rhs[i] = 0.0;
		ar_sense[i] = 'L';
		ar_rhs[i + nbvar] = capacite;
		ar_sense[i + nbvar] = 'L';
	}

	/* ------------------------------------ */
	/* n premieres contraintes : produit xi */
	/* ------------------------------------ */

	cpt = 0; lig = 0;
	for (i = 1; i <= nbvar; i++) 
	{
		ar_matbeg[lig++] = cpt;

		ar_matind[cpt] = i - 1;
		ar_matval[cpt++] = a[i] - capacite;

		for (j = 1; j < i; j++)
		{
			ar_matind[cpt] = pqe_cplex_convertit(j,i);
			ar_matval[cpt++] = a[j];
		}
		for(j = i + 1; j <= nbvar; j++)
		{
			ar_matind[cpt] = pqe_cplex_convertit(i,j);
			ar_matval[cpt++] = a[j];
		}
	}
	
	/* ------------------------------------------------- */
	/* n contraintes suivantes : produit xi-barre (1-xi) */
	/* ------------------------------------------------- */

	for (i = 1; i <= nbvar; i++) 
	{
		ar_matbeg[lig++] = cpt;
	
		for (j = 1; j < i; j++)
		{
			ar_matind[cpt] = j - 1;
			ar_matval[cpt++] = a[j];
		}
		ar_matind[cpt] = i - 1;
		ar_matval[cpt++] = capacite;
		for(j = i + 1; j <= nbvar; j++)
		{
			ar_matind[cpt] = j - 1;
			ar_matval[cpt++] = a[j];
		}
		for (j = 1; j < i; j++)
		{
			ar_matind[cpt] = pqe_cplex_convertit(j,i);
			ar_matval[cpt++] = - a[j];
		}
		for(j = i + 1; j <= nbvar; j++)
		{
			ar_matind[cpt] = pqe_cplex_convertit(i,j);
			ar_matval[cpt++] = - a[j];
		}
	
	}

} /* pqe_cplex_prep_contraintes_produit */

/*****************************************************************************/

void pqe_cplex_prepare_nelle_contrainte(int nbvar, int *p_ar_nzcount, double_1D ar_rhs, char_1D ar_sense,
                                          int_1D ar_matbeg, int_1D ar_matind, double_1D ar_matval, 
                                          char_2D newrowname, double_1D coef_obj, float MeilleurPrimal)
{
	int i;
	
	*p_ar_nzcount = 0;

	ar_rhs[0] = (double) MeilleurPrimal;
	ar_sense[0] = 'G';
	ar_matbeg[0] = 0;
	
	for (i = 0; i <= 2 * nbvar - 2; i++) 
		if (coef_obj[i] != 0.0) 
		{
			ar_matind[*p_ar_nzcount] = i;
			ar_matval[*p_ar_nzcount] = (double) coef_obj[i];;
			(*p_ar_nzcount)++;
		}
		
	pqe_init_1D((void **) &newrowname[0], sizeof(char), 10);
	strcpy(newrowname[0], "c_cutoff");

} /* pqe_cplex_prepare_nelle_contrainte */


/*****************************************************************************/

float pqe_cplex_demi_lin (int nbvar, float_1D coef_c1z, st_tri_1D ordre_poids,
                          float_1D q_lin, float_2D q_qua, int capacite, double_1D v_sol, 
                          double *p_nb_noeuds, bival_1D b_primal, float MeilleurPrimal,
                          int type_branchement, int type_BB, int var_select,
                          bool cutoff, bool contrainte)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* ------------------------------------------------ */

	CPXLPptr    lp       = NULL;
	const int   COLSPACE = 2 * nbvar;
	const int   ROWSPACE = 2 * nbvar + 1;
	const int   NZSPACE  = nbvar * nbvar + 3 * nbvar;
	const int   AR_COLSPACE = 1;
	const int   AR_ROWSPACE = 1;
	const int   AR_NZSPACE  = COLSPACE;
	int         status, mipstat, numcols, numrows, objsen, ar_nzcount;
	unsigned    colnamespace;
	unsigned    rownamespace;

	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
	            x, pi, slack, dj;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
	
	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);
	
	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_demi_lin (nbvar, coef_c1z, ordre_poids, q_lin, q_qua, capacite,
                     probname, &numcols, &numrows, &objsen, coef_obj, rhs, sense, 
                     matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, b_primal,
                     colindex, priority, type_branchement);

	pqe_cplex_set_pb_names_demi(nbvar, ordre_poids, colname, colnamestore, rowname, rownamestore, &colnamespace, &rownamespace);
/* plus en usage avec cplex 9.0
 	lp = CPXloadlpwnames (env, probname, numcols, numrows, objsen, coef_obj, rhs,
                          sense, matbeg, matcnt, matind, matval,
                          lb, ub, NULL, colname, colnamestore, rowname, rownamestore,
                          COLSPACE, ROWSPACE, NZSPACE, colnamespace, rownamespace);
*/
	if ( lp == NULL ) pqe_common_problemo("Failed to load LP.\n");

	status = CPXcopyctype (env, lp, ctype);
	if ( status ) pqe_common_problemo("Failed to load ctype\n");


	/* ----------------------------------------------------------------------- */
	/* Chargement d'une bonne solution admissible et d'un ordre de branchement */
	/* ----------------------------------------------------------------------- */

	if (contrainte)
	{
		pqe_cplex_prepare_nelle_contrainte (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
		                                      ar_matval, newrowname, coef_obj, MeilleurPrimal);
		status = CPXaddrows (env, lp, 0, 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, 
	                     	 NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}

	CPXcopyorder(env, lp, nbvar, colindex, priority, NULL);

	if (cutoff) CPXsetdblparam (env, CPX_PARAM_CUTLO, (double) MeilleurPrimal - 0.5); 

	switch (type_BB)
	{
		case BB_PROF :       CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_DFS);
		                     break;
		case BB_BEST_BOUND : CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTBOUND);
		                     break;
		case BB_BEST_EST :   CPXsetintparam (env, CPX_PARAM_NODESEL, CPX_NODESEL_BESTEST);
		                     break;
	}

	switch (var_select)
	{
		case VAR_SEL_POUSSEE : CPXsetintparam (env, CPX_PARAM_VARSEL, CPX_VARSEL_STRONG);
		case VAR_SEL_DEFAUT :  break;
	}
	
    if (mode_debug) CPXlpwrite (env, lp, "voir.lp"); 

	/* ------------------------ */
	/* Optimisation du probleme */
	/* ------------------------ */

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut (status, "mipopt dans lin1");

	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */

	mipstat = CPXgetstat(env,lp);

	if (   ( mipstat != CPXMIP_OPTIMAL ) 
	    && ( mipstat != CPXMIP_OPTIMAL_TOL )
		&& ( mipstat != CPXMIP_NODE_LIM_FEAS) ) pqe_cplex_arret_sur_erreur(mipstat, "demi_lin");

	if ( mipstat == CPXMIP_NODE_LIM_FEAS) 
	{
		print ("\nNode limit exceeded, integer solution exists...");
	}

	CPXgetmipobjval(env,lp, &obj);
	CPXgetmipx (env, lp, v_sol, 0, CPXgetnumcols(env,lp) - 1);
	(*p_nb_noeuds) = CPXgetnodecnt(env, lp);

	/* if ((num_appel % 10) == 0) getbase(lp, cstat, rstat); */

	CPXfreeprob(env, &lp); 

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return ((float) obj);

} /* pqe_cplex_demi_lin */

/*****************************************************************************/

void pqe_cplex_preparer_forcer_egalite (int nbvar, int *p_ar_nzcount, double_1D ar_rhs, char_1D ar_sense,
                                          int_1D ar_matbeg, int_1D ar_matind, double_1D ar_matval, 
                                          char_2D newrowname, float_1D coef_c1z, st_tri_1D ordre_poids, float_2D q_qua)
{
	int i, j, indice, minij, maxij;
	
	*p_ar_nzcount = 0;

	for (i = 1; i < nbvar; i++) 
	{ 
		indice = ordre_poids[i].indice;
		ar_rhs[i-1] = coef_c1z[indice];
		ar_sense[i-1] = 'L';
		ar_matbeg[i-1] = *p_ar_nzcount;
		ar_matind[*p_ar_nzcount] = i - 1;
		ar_matval[*p_ar_nzcount] = coef_c1z[indice];
		(*p_ar_nzcount)++;
		for (j = i + 1; j <= nbvar; j++)
		{
			minij = (ordre_poids[j].indice < indice ? ordre_poids[j].indice : indice);
			maxij = (ordre_poids[j].indice > indice ? ordre_poids[j].indice : indice); 
			ar_matind[*p_ar_nzcount] = j - 1;
			ar_matval[*p_ar_nzcount] = q_qua[minij][maxij];
			(*p_ar_nzcount)++;
		}
		ar_matind[*p_ar_nzcount] = nbvar + i - 1;
		ar_matval[*p_ar_nzcount] = - 1.0;
		(*p_ar_nzcount)++;
	}	
		
/*
	pqe_init_1D((void **) &newrowname[0], sizeof(char), 10);
	sprintf
	strcpy(newrowname[0], "egalitef");
*/
}

/*****************************************************************************/

float pqe_cplex_lin_croisee (int nbvar, float_1D Lmax, float_1D sur_phi, float_1D Dmax, st_tri_1D ordre_poids,
                               float_1D q_lin, float_2D q_qua, int capacite, double_1D v_sol, 
                               double *p_nb_noeuds, bival_1D b_primal, float MeilleurPrimal,
                               int type_branchement, int type_BB, int var_select,
                               bool cutoff, bool contrainte)
{
/*
	/ * ------------------------ * /
	/ * Linearisation abandonnee * /
	/ * ------------------------ * /

	bool  forcer_egalite = FAUX;

	CPXLPptr    lp       = NULL;
	const int   COLSPACE = 4 * nbvar - 2;  / * Nombre de variables * / 
	const int   ROWSPACE = 6 * nbvar - 2;  / * Nb de contraintes * /
	const int   NZSPACE  = (3 * nbvar * nbvar + 23 * nbvar - 18) / 2; / * Nb d'entrees non nulles * /
	int         status, mipstat, numcols, numrows, objsen, cnt, ar_nzcount;
	unsigned    colnamespace;
	unsigned    rownamespace;

	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval,
	            x, pi, slack, dj;
	double      valeur;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
	
	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);
	
	/ * ---------------------------- * /
	/ * Ecriture et chargement du pb * /
	/ * ---------------------------- * /

	pqe_cplex_set_pb_data_lin_croisee (nbvar, Lmax, sur_phi, Dmax, ordre_poids, q_lin,
	             q_qua, capacite, probname, &numcols, &numrows, &objsen, coef_obj, rhs, 
	             sense, matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, 
	             b_primal, colindex, priority, type_branchement);

	/ * ----------------------------------------------------------------------- * /
	/ * Chargement d'une bonne solution admissible et d'un ordre de branchement * /
	/ * ----------------------------------------------------------------------- * /

	/ * ------------------------ * /
	/ * Optimisation du probleme * /
	/ * ------------------------ * /

	/ * --------------------------- * /
	/ * Recuperation de la solution * /
	/ * --------------------------- * /

	CPXfreeprob(env, &lp); 

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return ((float) obj);
*/
	return (0.0);
} /* pqe_cplex_lin_croisee */

/*****************************************************************************/

void pqe_cplex_changer_noms (CPXENVptr env, CPXLPptr lp, int nbvar)
{
	int i, j, cpt;
	char nom[20];

	/* Variables */

	cpt = nbvar;
	for (i = 1; i <= nbvar; i++)
		for (j= i + 1; j <= nbvar; j++)
		{
			sprintf(nom, "w%d_%d", i, j);
			CPXchgname(env, lp, 'c', cpt++, nom);
		}
	
	/* Contraintes */
	CPXchgname (env, lp, 'r', 0, "cap");
	
} /* pqe_cplex_changer_noms */

/*****************************************************************************/

void pqe_cplex_signaler_sol_entiere(int i, int nbvar)
{
	FILE *fich;

	fich = fopen("sol_entiere_trouvee","a");
	fprintf(fich, "Solution entiere trouvee pour sous_phi%d (%d var)", i, nbvar);
	fclose(fich);
}

/*****************************************************************************/

float pqe_cplex_LPi_stocha (bool prime, int i, int nbvar, float_1D sur_phi, st_tri_1D ordre_poids,
                         float_1D q_lin, float_2D q_qua, float_1D cout, float CAP, float_1D moyenne, float k,
                         float meilleur_primal, bival_1D b_primal, bool *p_faisable_continu, 
                         int *p_nbBB, int *p_nbsol01, int *p_nbopt01)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* prime indique si l'on resout LP'i ou LPi         */
	/* ------------------------------------------------ */

	CPXLPptr    lp          = NULL;
	const int   COLSPACE    = 2 * nbvar - 1;                   /* Nombre de variables */ 
	int   ROWSPACE    = 2 * nbvar + 2;                   /* Nb de contraintes */
	/* const int   NZSPACE     = nbvar * (nbvar + 11) / 2 + i - 7;  Nb d'entrees non nulles */
	const int   NZSPACE     = (nbvar * nbvar + 13 * nbvar - 6) / 2; /* Nb d'entrees non nulles */
	int         j, indice, indicej, status, mipstat, numcols = COLSPACE, 
	            numrows, objsen;
	double      temps_rel, val_sol_entiere;
	float       resu;
	const bool  sol_adm  = FAUX, politique_BB_1 = FAUX;
	bool        faisable_01 = VRAI;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices; 
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
	            x, pi, slack, dj; 
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
	char        nom[30];
    int l;

	if (!prime) ROWSPACE--;
	numrows = ROWSPACE;
	*p_faisable_continu = VRAI;
	indice = ordre_poids[i].indice;
	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_LPi_stocha (prime, i, nbvar, COLSPACE, ROWSPACE,
                 sur_phi, ordre_poids, q_lin, q_qua, cout, CAP, moyenne, k, &objsen,
                 coef_obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, 
                 ctype, meilleur_primal);

	lp = CPXcreateprob (env, &status, "LPi_stocha");
	if (lp == NULL) pqe_common_problemo("Failed to create 2LP.\n");
 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	sprintf(nom, "2LP%1d.lp",i);
	if (mode_debug) CPXlpwrite (env,lp, nom); 

	/* ----------------------------------------------------------- */
	/* Le probleme doit avoir au moins une solution entiere (en x) */
	/* Calcul de la relaxation continue -> resu                    */
	/* ----------------------------------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, &temps_rel, &resu, x, pi, slack, dj, p_faisable_continu);

	if (mode_debug) 
	{
		print ("i = %d, indice = %d", i, indice);
		/* if (FAUX) */
		if ((*p_faisable_continu) && (b_primal[indice] == 0)) 
			print("\n faisable en continu, b_primal[i] = 0 et resu = %f\n", resu);
	}



	if ((!(*p_faisable_continu)) || (!faisable_01)) resu = 0.0;
	else
	{
		/* --------------------------------------------------------------- */
		/* On retient comme borne inf la partie entiere superieure de resu */
		/* --------------------------------------------------------------- */
		resu = ceil(resu);
	}

	CPXfreeprob(env, &lp); 
	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);
	return (resu);

} /* pqe_cplex_LPi_stocha */

/*****************************************************************************/

float pqe_cplex_LPi (bool prime, int i, int nbvar, float_1D sur_phi, st_tri_1D ordre_poids,
                         float_1D q_lin, float_2D q_qua, int_1D a, int capacite, 
                         float meilleur_primal, bival_1D b_primal, bool *p_faisable_continu, 
                         int *p_nbBB, int *p_nbsol01, int *p_nbopt01, int sens)
{
    /* ------------------------------------------------ */
    /* Toutes les variables autres que env sont locales */
    /* prime indique si l'on resout LP'i ou LPi         */
    /* ------------------------------------------------ */

    CPXLPptr    lp          = NULL;
    const int   COLSPACE    = 2 * nbvar - 1;                   /* Nombre de variables */ 
    int         ROWSPACE    = 2 * nbvar + 1;                   /* Nb de contraintes */
    const int   NZSPACE     = 4 *(nbvar * nbvar); /* Nb d'entrees non nulles (au plus...) */
    int         j, indice, indicej, status, mipstat, numcols = COLSPACE, 
                numrows, objsen;
    double      temps_rel, val_sol_entiere;
    float       resu;
    const bool  sol_adm  = FAUX, politique_BB_1 = FAUX;
    bool        faisable_01 = VRAI;
    int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices; 
    double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
                x, pi, slack, dj; 
    char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
    char_2D     colname, rowname, newrowname;
    char        nom[30];

    if (!prime) ROWSPACE--;
    numrows = ROWSPACE;
    *p_faisable_continu = VRAI;
    indice = ordre_poids[i].indice;
    pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
                 &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
                 &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
                 &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
                 &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

    /* ---------------------------- */
    /* Ecriture et chargement du pb */
    /* ---------------------------- */

    pqe_cplex_set_pb_data_LPi (prime, i, nbvar, COLSPACE, ROWSPACE, sur_phi, ordre_poids,
                 q_lin, q_qua, a, capacite, &objsen, coef_obj, rhs, sense, matbeg, matcnt, 
                 matind, matval, lb, ub, ctype, meilleur_primal, sens);

    lp = CPXcreateprob (env, &status, "LPi");
    if (lp == NULL) pqe_common_problemo("Failed to create LP.\n");
     status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
    sprintf(nom, "LP%1d.lp",i);
    if (mode_debug) CPXlpwrite (env,lp, nom); 

    /* ----------------------------------------------------------- */
    /* Le probleme doit avoir au moins une solution entiere (en x) */
    /* Calcul de la relaxation continue -> resu                    */
    /* ----------------------------------------------------------- */

    pqe_cplex_relaxation_continue(env, lp, &temps_rel, &resu, x, pi, slack, dj, p_faisable_continu);

    if (mode_debug) 
    {
        print ("i = %d, indice = %d", i, indice);
        /* if (FAUX) */
        if ((*p_faisable_continu) && (b_primal[indice] == 0)) 
            print("\n faisable en continu, b_primal[i] = 0 et resu = %f\n", resu);
    }

    /* A SUPPRIMER plus proprement */
    if (FAUX && (*p_faisable_continu) && (resu > 0.0) && (b_primal[indice] == 0)) 
    {
        /* ------------------------------------------------------ */
        /* Il existe une solution continue. On doit s'assurer de  */
        /* la faisabilit\E9 en 0-1 avant de garder la borne obtenue */
        /* par relaxation continue. Le pb est toujours faisable   */
        /* si b_primal[indice] vaut 1                             */
        /* ------------------------------------------------------ */

        if (POLITIQUE_BB_0) 
        {
            /* ------------------------------------------------ */
    	    /* Il n'existe pas forcement de solutions en 0-1 :  */
            /* on lance le B&B sur quelques noeuds et si aucune */
            /* solution entiere admissible n'a ete trouvee, on  */
            /* considere que le pb est infaisable en 0-1 pour   */
            /* ne pas passer trop de temps dans ce B&B          */
            /* ------------------------------------------------ */

            (*p_nbBB)++;
            /* pqe_cplex_masquer_affichage(env); */
            pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, BB_BEST_BOUND, VAR_SEL_DEFAUT,
                                          nom, pr_indices, pr_value); 
            CPXsetintparam (env, CPX_PARAM_NODELIM, N_LIM_sous_phiI);
            status = CPXmipopt (env,lp);
            CPXsetintparam (env, CPX_PARAM_NODELIM, N_LIMITE);
            pqe_cplex_affichage_par_defaut(env, affichage);
            pqe_cplex_verif_statut(status, "mipopt dans sous_phii");
            mipstat = CPXgetstat(env,lp);
            faisable_01 = (   (mipstat == CPXMIP_OPTIMAL)       || (mipstat == CPXMIP_NODE_LIM_FEAS)
                           || (mipstat == CPXMIP_TIME_LIM_FEAS) || (mipstat == CPXMIP_FAIL_FEAS)
                           || (mipstat == CPXMIP_ABORT_FEAS));
            if (faisable_01) (*p_nbsol01)++;
            if (mipstat == CPXMIP_OPTIMAL)
            {
                (*p_nbopt01)++;
                CPXgetmipobjval(env, lp, &val_sol_entiere);
                resu = (float) val_sol_entiere;
                print("\n Solution entiere trouvee (%f) pour sous_phi%d\n", resu, i);
                pqe_cplex_signaler_sol_entiere(i, nbvar);
            }
            if (!faisable_01) print ("-> pas de sol entiere trouv\E9e (mais il en existe peut-etre)\n");
        }
        else
        { 
            faisable_01 = FAUX;
            resu = 0.0;
    	}
    }

	/* A supprimer aussi plus nettement */
	if (FAUX && (*p_faisable_continu) && (resu > 0.0) && (b_primal[indice] == 1) && (politique_BB_1)) 
	{
		/* -------------------------------------------------------------------------- */
		/* On essai d'ameliorer la borne inf en lancant un B&B limite en nb de noeuds */
		/* histoire de trouver l'optimum de LKPmin en 0-1. On connait deja une bonne  */
		/* solution admissible, b_primal.                                             */
		/* -------------------------------------------------------------------------- */
		(*p_nbBB)++;
		/* pqe_cplex_masquer_affichage(env); */
		pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, BB_BEST_BOUND, VAR_SEL_DEFAUT,
	                               nom, pr_indices, pr_value); 
		for (j = 1; j <= nbvar; j++)
		{
			indicej = ordre_poids[j].indice;
			pr_indices[j - 1] = j - 1;
			pr_value[j - 1] = (double) b_primal[indicej];
			if (j < nbvar-1) pr_indices[j + nbvar - 1] = j + nbvar - 1;
			if (j < nbvar-1) pr_value[j + nbvar - 1] = CPX_INFBOUND; 
		}
		//CPXsetintparam (env, CPX_PARAM_MIPSTART, CPX_ON); 
		CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON); 
		status = CPXcopymipstart (env, lp, COLSPACE, pr_indices, pr_value);
		if (status != 0) pqe_common_problemo("CPXcopymipstart");
		CPXsetintparam (env, CPX_PARAM_NODELIM, N_LIM_sous_phiI);
		status = CPXmipopt (env,lp);
		CPXsetintparam (env, CPX_PARAM_NODELIM, N_LIMITE);
		pqe_cplex_affichage_par_defaut(env, affichage);
		pqe_cplex_verif_statut(status, "mipopt dans sous_phii");
		mipstat = CPXgetstat(env,lp);
		(*p_nbsol01)++;
		if (mipstat == CPXMIP_OPTIMAL)
		{
			(*p_nbopt01)++;
			CPXgetmipobjval(env, lp, &val_sol_entiere);
			resu = (float) val_sol_entiere;
			print("\n Solution entiere trouvee (%f) pour sous_phi%d\n", resu, i);
			pqe_cplex_signaler_sol_entiere(i, nbvar);
		}
		else print("-> sol optimale en 0-1 non trouvee\n");
	}

	if ((!(*p_faisable_continu)) || (!faisable_01)) resu = 0.0;
	else
	{
		/* --------------------------------------------------------------- */
		/* On retient comme borne inf la partie entiere superieure de resu */
		/* --------------------------------------------------------------- */
		resu = ceil(resu);
	}

	CPXfreeprob(env, &lp); 
	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);
	return (resu);

} /* pqe_cplex_LPi */

/*****************************************************************************/

double pqe_cplex_max_LPi (bool prime, int i, int N, double_1D sur_phi, double_1D sous_phi,
                         double_1D rlin, double_2D rqua, double_1D weight, int CAPA, 
                         double meilleur_primal, bool *p_faisable_continu, 
                         int *p_nbBB, int *p_nbsol01, int *p_nbopt01, int sens)
{
    /* ------------------------------------------------ */
    /* Toutes les variables autres que env sont locales */
    /* prime indique si l'on resout LP'i ou LPi         */
    /* ------------------------------------------------ */

    CPXLPptr    lp          = NULL;
    const int   COLSPACE    = 2 * N - 1;                   /* Nombre de variables */ 
    int         ROWSPACE    = 2 * N + 1;                   /* Nb de contraintes */
    const int   NZSPACE     = 4 *(N * N); /* Nb d'entrees non nulles (au plus...) */
    int         j, status, mipstat, numcols = COLSPACE, 
                numrows, objsen;
    double      temps_rel, val_sol_entiere;
    float       resu;
    const bool  sol_adm  = FAUX, politique_BB_1 = FAUX;
    bool        faisable_01 = VRAI;
    int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices; 
    double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
                x, pi, slack, dj; 
    char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
    char_2D     colname, rowname, newrowname;
    char        nom[30];

    if (!prime) ROWSPACE--;
    numrows = ROWSPACE;
    *p_faisable_continu = VRAI;
    pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
                 &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
                 &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
                 &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
                 &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

    /* ---------------------------- */
    /* Ecriture et chargement du pb */
    /* ---------------------------- */

    pqe_cplex_set_pb_data_max_LPi (prime, i, N, COLSPACE, ROWSPACE, sur_phi, sous_phi,
                 rlin, rqua, weight, CAPA, &objsen, coef_obj, rhs, sense, matbeg, matcnt, 
                 matind, matval, lb, ub, ctype, meilleur_primal, sens);

    lp = CPXcreateprob (env, &status, "max_LPi");
    if (lp == NULL) pqe_common_problemo("Failed to create max_LPi\n");
    status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
    sprintf(nom, "maxLP%1d.lp",i);
    if (mode_debug) CPXlpwrite (env,lp, nom); 

    pqe_cplex_relaxation_continue(env, lp, &temps_rel, &resu, x, pi, slack, dj, p_faisable_continu);

    if (mode_debug) 
    {
        print ("i = %d", i);
        if (*p_faisable_continu) print(" faisable en continu - resu = %f\n", resu);
        else print(" il n'existe pas de solution admissible (en continu) \n");
    }

    if (!(*p_faisable_continu)) resu = 0.0;
    else {
        /* --------------------------------------------------------------- */
        /* On retient comme borne sup la partie entiere superieure de resu */
        /* --------------------------------------------------------------- */
	    resu = floor(resu);
    }

    CPXfreeprob(env, &lp); 
    pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);
    return (resu);

} /* pqe_cplex_max_LPi */

/*****************************************************************************/

double pqe_cplex_min_LPi (bool prime, int i, int N, double_1D sur_phi, double_1D sous_phi,
                         double_1D rlin, double_2D rqua, double_1D weight, int CAPA, 
                         double meilleur_primal, bool *p_faisable_continu, 
                         int *p_nbBB, int *p_nbsol01, int *p_nbopt01, int sens)
{
    /* ------------------------------------------------ */
    /* Toutes les variables autres que env sont locales */
    /* prime indique si l'on resout LP'i ou LPi         */
    /* ------------------------------------------------ */

    CPXLPptr    lp          = NULL;
    const int   COLSPACE    = 2 * N - 1;                   /* Nombre de variables */ 
    int         ROWSPACE    = 3 * N + 1;                   /* Nb de contraintes */
    const int   NZSPACE     = 4 *(N * N); /* Nb d'entrees non nulles (au plus...) */
    int         j, status, mipstat, numcols = COLSPACE, 
                numrows, objsen;
    double      temps_rel, val_sol_entiere;
    float       resu;
    const bool  sol_adm  = FAUX, politique_BB_1 = FAUX;
    bool        faisable_01 = VRAI;
    int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices; 
    double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
                x, pi, slack, dj; 
    char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
    char_2D     colname, rowname, newrowname;
    char        nom[30];

    if (!prime) ROWSPACE--;
    numrows = ROWSPACE;
    *p_faisable_continu = VRAI;
    pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
                 &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
                 &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
                 &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
                 &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

    /* ---------------------------- */
    /* Ecriture et chargement du pb */
    /* ---------------------------- */

    pqe_cplex_set_pb_data_min_LPi (prime, i, N, COLSPACE, ROWSPACE, sur_phi, sous_phi,
                 rlin, rqua, weight, CAPA, &objsen, coef_obj, rhs, sense, matbeg, matcnt, 
                 matind, matval, lb, ub, ctype, meilleur_primal, sens);

    lp = CPXcreateprob (env, &status, "max_LPi");
    if (lp == NULL) pqe_common_problemo("Failed to create min_LPi\n");
    status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
    sprintf(nom, "minLP%1d.lp",i);
    if (mode_debug) CPXlpwrite (env,lp, nom); 

    pqe_cplex_relaxation_continue(env, lp, &temps_rel, &resu, x, pi, slack, dj, p_faisable_continu);

    if (mode_debug) 
    {
        print ("i = %d", i);
        if (*p_faisable_continu) print(" faisable en continu - resu = %f\n", resu);
        else print(" il n'existe pas de solution admissible (en continu) \n");
    }

    if (!(*p_faisable_continu)) resu = 0.0;
    else
    {
        /* --------------------------------------------------------------- */
        /* On retient comme borne sup la partie entiere inf\8Erieure de resu */
        /* --------------------------------------------------------------- */
	    resu = ceil(resu);
    }

    CPXfreeprob(env, &lp); 
    pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);
    return (resu);

} /* pqe_cplex_min_LPi */

/*****************************************************************************/

float pqe_cplex_lin_classique (int nbvar, float_1D q_lin, float_2D q_qua, int capacite, int_1D a,
                                 double_1D v_sol, double *p_nb_noeuds, int type_branchement, int type_BB, 
                                 int var_select, float *p_meilleur_dual, double *p_temps_UB, 
                                 bool lin_produit, bool *p_temps_limite)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* ------------------------------------------------ */

	CPXLPptr    lp          = NULL;
	const bool  sol_adm     = FAUX;
	const int   COLSPACE    = nbvar * (nbvar + 1) / 2;   /* Nombre de variables */ 
	const int   ROWSPACE    = nbvar * nbvar - nbvar + 1; /* Nb de contraintes */
	const int   NZSPACE     = 2 * nbvar * nbvar - nbvar; /* Nb d'entrees non nulles */
	const int   AR_COLSPACE = 2 * nbvar;
	const int   AR_ROWSPACE = 2 * nbvar;
	const int   AR_NZSPACE  = nbvar * (3 * nbvar - 1);
	int         status, numcols, numrows, objsen, ar_nzcount;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
	double_1D   x, pi, slack, dj; 
	bool        faisable;


	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_lin_classique (nbvar, q_lin, q_qua, capacite,
                     probname, &numcols, &numrows, &objsen, coef_obj, rhs, sense, 
                     matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, b_primal,
                     colindex, priority, type_branchement, pr_indices, pr_value);

	lp = CPXcreateprob (env, &status, probname);
	if (lp == NULL) pqe_common_problemo("Failed to create LP.\n");
 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                            sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
	pqe_cplex_changer_noms (env, lp, nbvar);
	if (lin_produit)
	{	
		pqe_cplex_prep_contraintes_produits (nbvar, a, capacite, &ar_nzcount, ar_rhs, ar_sense,
		                                       ar_matbeg, ar_matind, ar_matval, newrowname);
		status = CPXaddrows (env, lp, 0, 2 * nbvar, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, 
		                     ar_matind, ar_matval, NULL, NULL);
		if (status != 0) pqe_common_problemo("addrows");
	}

	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, p_temps_UB, p_meilleur_dual, x, pi, slack, dj, &faisable);

	/* ----------------------------- */
	/* Resolution du probleme en 0-1 */
	/* ----------------------------- */

	pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, type_BB, var_select,
	                        (lin_produit ? "linprod.lp" : "class.lp"), pr_indices, pr_value); 
	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut (status, "mipopt dans lin_classique");
	pqe_cplex_recuperation_solution(env, lp, p_temps_limite, &obj, v_sol, p_nb_noeuds);

	CPXfreeprob(env, &lp); 
	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return (*p_temps_limite ? 0.0 : (float) obj);

} /* pqe_cplex_lin_classique */

/*****************************************************************************/

float pqe_cplex_lin1 (int nbvar, float_1D coef_c1z, st_tri_1D ordre_poids,
                          float_1D q_lin, float_2D q_qua, int capacite, double_1D v_sol, 
                          double *p_nb_noeuds, bival_1D b_primal, float MeilleurPrimal,
                          int type_branchement, int type_BB, int var_select,
                          bool cutoff, bool contrainte, float *p_meilleur_dual, 
                          double *p_temps_UB, bool *p_temps_limite)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* ------------------------------------------------ */

	const bool  forcer_egalite = FAUX;
	const bool  sol_adm        = VRAI;
	CPXLPptr    lp             = NULL;
	const int   COLSPACE       = 2 * nbvar - 1;  /* Nombre de variables */ 
	const int   ROWSPACE       = 3 * nbvar - 2;  /* Nb de contraintes : Anciennement 2 * nbvar - 1 */
	const int   NZSPACE        = nbvar * nbvar + 5 * nbvar - 5; /* Nb d'entrees non nulles                  */
	const int   AR_COLSPACE    = nbvar;
	const int   AR_ROWSPACE    = nbvar;
	const int   AR_NZSPACE     = (nbvar * nbvar + 3 * nbvar - 4) / 2;
	int         status, numcols, numrows, objsen, ar_nzcount;
	unsigned    colnamespace;
	unsigned    rownamespace;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value;
	double      obj;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
	double_1D   x, pi, slack, dj;
	bool        faisable;

	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_lin1(nbvar, coef_c1z, ordre_poids, q_lin, q_qua, capacite,
                     probname, &numcols, &numrows, &objsen, coef_obj, rhs, sense, 
                     matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, b_primal,
                     colindex, priority, type_branchement, pr_indices, pr_value);

	pqe_cplex_set_pb_names(nbvar, ordre_poids, colname, colnamestore, rowname,
	             rownamestore, &colnamespace, &rownamespace);

	lp = CPXcreateprob (env, &status, probname);

	if ( lp == NULL ) pqe_common_problemo("Failed to create LP.\n");

 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                        sense, matbeg, matcnt, matind, matval, lb, ub, NULL);

	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, p_temps_UB, p_meilleur_dual, x, pi, slack, dj, &faisable);

	/* ----------------------------- */
	/* Resolution du probleme en 0-1 */
	/* ----------------------------- */

	pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, type_BB, var_select, "voir.lp", pr_indices, pr_value);

	if (contrainte)
	{
		pqe_cplex_prepare_nelle_contrainte (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
	                                       	  ar_matval, newrowname, coef_obj, MeilleurPrimal);
		status = CPXaddrows (env, lp, 0, 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	if (forcer_egalite)
	{
		pqe_cplex_preparer_forcer_egalite (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
		                                     ar_matval, newrowname, coef_c1z, ordre_poids, q_qua);
		newrowname = NULL;
		status = CPXaddrows (env, lp, 0, nbvar - 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, 
	                     	 NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	CPXcopyorder(env, lp, nbvar, colindex, priority, NULL);
	if (cutoff) CPXsetdblparam (env, CPX_PARAM_CUTLO, (double) MeilleurPrimal - 0.5);

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans lin1");

	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */

	pqe_cplex_recuperation_solution(env, lp, p_temps_limite, &obj, v_sol, p_nb_noeuds);

	CPXfreeprob(env, &lp); 

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return (*p_temps_limite ? 0.0 : (float) obj);

} /* pqe_cplex_lin1 */

/*****************************************************************************/

void pqe_cplex_prepare_fixations(int nbvar, st_tri_1D ordre_poids, int nbfix, la_st_var_1D variables,
                                   double_1D ar_rhs, char_1D ar_sense, int_1D ar_matbeg, 
                                   int_1D ar_matind, double_1D ar_matval)
{
	int i, indice, cpt = 0;

	for (i = 0; i < nbfix; i++)
	{
		ar_rhs[i] = 0.0; /* Fixation \E0 0 */
		ar_sense[i] = 'E';
	}

	for (i = 1; i <= nbvar; i++)
	{
		indice = ordre_poids[i].indice;
		if (variables[indice].fixee)
		{
			ar_matbeg[cpt] = cpt;
			ar_matind[cpt] = i - 1; /* et non indice - 1, bien sur */
			ar_matval[cpt] = 1.0;
			cpt++;
		}
	}
	
} /* pqe_cplex_prepare_fixations */

/*****************************************************************************/

int pqe_cplex_stocha_callback (CPXENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{
    double best_integer, best_remaining;
    int    nodes_solved, nodes_remaining, mip_iterations, statut = 0;
	FILE   *fich_erreur, *BBcourant;

    if (wherefrom == CPX_CALLBACK_MIP)
    {
        statut = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                    CPX_CALLBACK_INFO_BEST_INTEGER, &best_integer);
        if (statut) 
        { 
				fich_erreur = fopen("pb_detecte", "a");
				fprintf(fich_erreur, "callback, stocha, best integer");
				fclose(fich_erreur);
		}

        statut = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                    CPX_CALLBACK_INFO_BEST_REMAINING, &best_remaining);
        if (statut)
        { 
				fich_erreur = fopen("pb_detecte", "a");
				fprintf(fich_erreur, "callback, stocha, best remaining");
				fclose(fich_erreur);
		}

        statut = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                    CPX_CALLBACK_INFO_NODE_COUNT, &nodes_solved);
        if (statut)
        { 
				fich_erreur = fopen("pb_detecte", "a");
				fprintf(fich_erreur, "callback, stocha, nodes solved");
				fclose(fich_erreur);
		}

        statut = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                    CPX_CALLBACK_INFO_NODES_LEFT, &nodes_remaining);
        if (statut)
        { 
				fich_erreur = fopen("pb_detecte", "a");
				fprintf(fich_erreur, "callback, stocha, nodes remaining");
				fclose(fich_erreur);
		}

        statut = CPXgetcallbackinfo (env, cbdata, wherefrom,
                                    CPX_CALLBACK_INFO_MIP_ITERATIONS, &mip_iterations);
        if (statut) 
        { 
				fich_erreur = fopen("pb_detecte", "a");
				fprintf(fich_erreur, "callback, stocha, mip iterations");
				fclose(fich_erreur);
		}

        /* --------- */
        /* Affichage */
        /* --------- */
        if ((mip_iterations % 200000 == 0))
        {
            BBcourant = fopen("BBcourant", "a");
            fprintf(BBcourant, "\n|Meilleur prim |Meilleur rest | Nb noeuds| Nds rest | MIP iter |");
            fclose(BBcourant);
        }
        if (mip_iterations % 20000 == 0)
        {
            BBcourant = fopen("BBcourant", "a");
            fprintf(BBcourant, "\n|%13.2lf |%13.2lf |%9d |%9d |%9d |", best_integer,
                    best_remaining, nodes_solved, nodes_remaining, mip_iterations);
            fclose(BBcourant);
        }
        if (best_integer > 0.001) return 1;
    }
    else 
    { 
		fich_erreur = fopen("pb_detecte", "a");
		fprintf(fich_erreur, "callback, stocha, cplex_stocha_callback");
		fclose(fich_erreur);
	}

	return (statut);
} /* pqe_cplex_stocha_callback */

/*****************************************************************************/

float pqe_cplex_lin2_stocha (int nbvar, float_1D coef_c1z, float_1D sous_phi, st_tri_1D ordre_poids,
                          float_1D q_lin, float_2D q_qua, float_1D cout, float CAP,
                          float_1D moyenne, float k, double_1D v_sol, 
                          double *p_nb_noeuds, bival_1D b_primal,
                          int type_branchement, int type_BB, int var_select, bool cutoff,
                          bool contrainte, float *p_meilleur_dual, double *p_temps_UB, 
                          bool *p_temps_limite, la_st_var_1D variables, int nbfix,
                          bool *p_sol_positive, int jeu, float valeur_coeff, float alpha)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
    /* L'objectif possede une constante : k*k           */
	/* ------------------------------------------------ */

	CPXLPptr    lp             = NULL;
	const int   COLSPACE       = 2 * nbvar,  /* Nombre de variables */
	            ROWSPACE       = 2 * nbvar + 1,  /* Nb de contraintes : Anciennement 2 * nbvar - 1 */
	            NZSPACE        = (nbvar * nbvar + 11 * nbvar - 4) / 2, /* Nb d'entrees non nulles                  */
				                          /* avant \E8 au lieu de 9 et dans ar_nzspace 3 aulieu de 5 */
	            AR_COLSPACE    = nbvar,
	            AR_ROWSPACE    = nbvar,
	            AR_NZSPACE     = (nbvar * nbvar + 7 * nbvar) / 2;
	const bool  forcer_egalite = FAUX;
	const bool  sol_adm        = VRAI;
	int         status, numcols, numrows, objsen, ar_nzcount;
	double      obj;
	bool        faisable;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
	            x, pi, slack, dj;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
    FILE        *BBcourant;

    (*p_sol_positive) = FAUX;
    BBcourant = fopen("BBcourant", "w");
    fprintf(BBcourant, "\n JEU n\B0 %d   -   k = %.0f %% * opt_moyen  -  alpha = %.0f %%\n",
            jeu, valeur_coeff * 100.0, alpha * 100.0);
    fprintf(BBcourant, "\n|Meilleur prim |Meilleur rest | Nb noeuds| Nds rest | MIP iter |");
    fclose(BBcourant);

	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_lin2_stocha(nbvar, COLSPACE, ROWSPACE, coef_c1z, sous_phi,
                 ordre_poids, q_lin, q_qua, cout, CAP, moyenne, k,
                 probname, &numcols, &numrows, &objsen, coef_obj, rhs, sense, 
                 matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, b_primal,
                 colindex, priority, type_branchement, pr_indices, pr_value);

/* 
	pqe_cplex_set_pb_names(nbvar, ordre_poids, colname, colnamestore, rowname,
	             rownamestore, &colnamespace, &rownamespace);
*/

	lp = CPXcreateprob (env, &status, probname);
	if ( lp == NULL ) pqe_common_problemo("Failed to create LP.\n");

/*
	status = CPXreadcopyprob (env, lp, "cover.lp", NULL);
	CPXsetintparam (env, CPX_PARAM_CLIQUES, -1); 
	CPXsetintparam (env, CPX_PARAM_COVERS, -1); 
	CPXsetintparam (env, CPX_PARAM_MIPSTART, CPX_ON); 
	status = CPXcopymipstart (env, lp, COLSPACE, pr_indices, pr_value);
	if (status != 0) pqe_common_problemo("CPXcopymipstart");
	status = CPXmipopt (env,lp);
*/

 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                        sense, matbeg, matcnt, matind, matval, lb, ub, NULL);

	if (nbfix > 0)
	{
		pqe_cplex_prepare_fixations(nbvar, ordre_poids, nbfix, variables, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval);
		status = CPXaddrows (env, lp, 0, nbfix, nbfix, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, NULL);
		if (status != 0) pqe_common_problemo("addrows");
	}

	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, p_temps_UB, p_meilleur_dual, x, pi, slack, dj, &faisable);

	/* ----------------------------- */
	/* Resolution du probleme en 0-1 */
	/* ----------------------------- */

	pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, type_BB, var_select, "voir.lp", pr_indices, pr_value);

	CPXsetintparam (env, CPX_PARAM_CLIQUES, 0); 
	CPXsetintparam (env, CPX_PARAM_COVERS, 0); 


	if (contrainte)
	{
		pqe_cplex_prepare_nelle_contrainte (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
	                                       	  ar_matval, newrowname, coef_obj, MeilleurPrimal);
		status = CPXaddrows (env, lp, 0, 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	if (forcer_egalite)
	{
		pqe_cplex_preparer_forcer_egalite (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
		                                     ar_matval, newrowname, coef_c1z, ordre_poids, q_qua);
		newrowname = NULL;
		status = CPXaddrows (env, lp, 0, nbvar - 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, 
	                     	 NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	CPXcopyorder(env, lp, nbvar, colindex, priority, NULL);
	if (cutoff) CPXsetdblparam (env, CPX_PARAM_CUTLO, (double) MeilleurPrimal - 0.5);

	if (mode_debug) CPXlpwrite (env,lp, "lin2stocha.lp");

/* Attention : ne compile plus avec cplex 9.0 !!!
    status = CPXsetmipcallbackfunc (env, &pqe_cplex_stocha_callback, NULL);
    if (status) pqe_common_problemo("setlpcallback dans lin2stocha");
*/

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans lin2stocha");

    status = CPXsetmipcallbackfunc (env, NULL, NULL);
    if (status) pqe_common_problemo("setlpcallback dans lin2stocha");

	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */

	/* pqe_cplex_recuperation_solution(env, lp, p_temps_limite, &obj, v_sol, p_nb_noeuds); */
	status = CPXgetstat(env,lp);

	if (   ( status != CPXMIP_OPTIMAL )
	    && ( status != CPXMIP_OPTIMAL_TOL )
		&& ( status != CPXMIP_NODE_LIM_FEAS)
		&& ( status != CPXMIP_TIME_LIM_FEAS)
		&& ( status != CPXMIP_TIME_LIM_INFEAS) )
    {
        if (status == CPXMIP_ABORT_FEAS) (*p_sol_positive) = VRAI;
        else
        {
            if (status == CPXMIP_ABORT_INFEAS) pqe_common_problemo("cplex abort infeas");
            else pqe_cplex_arret_sur_erreur(status, "recuperation_sol");
        }
    }

	if (status == CPXMIP_NODE_LIM_FEAS)
	{
		pqe_common_problemo("Limite du nb de noeuds atteinte, il existe une solution entiere...");
	}
	
	*p_temps_limite = ((status == CPXMIP_TIME_LIM_FEAS) || (status == CPXMIP_TIME_LIM_INFEAS));

	CPXgetmipobjval(env, lp, &obj);
	CPXgetmipx (env, lp, v_sol, 0, CPXgetnumcols(env,lp) - 1);
	(*p_nb_noeuds) = CPXgetnodecnt(env, lp);

	CPXfreeprob(env, &lp);

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return (*p_temps_limite ? 0.0 : (float) obj);

} /* pqe_cplex_lin2_stocha */


/*****************************************************************************/

float pqe_cplex_lin2 (int nbvar, float_1D coef_c1z, float_1D sous_phi, st_tri_1D ordre_poids,
                          float_1D q_lin, float_2D q_qua, int capacite, double_1D v_sol, 
                          double *p_nb_noeuds, bival_1D b_primal, float meilleur_primal,
                          int type_branchement, int type_BB, int var_select, bool cutoff, 
                          bool contrainte, float *p_meilleur_dual, double *p_temps_UB, 
                          bool *p_temps_limite, la_st_var_1D variables, int nbfix, int seed,
                          float densite, int ordre_var)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* ------------------------------------------------ */

	CPXLPptr    lp             = NULL;
	const int   COLSPACE       = 2 * nbvar - 1,  /* Nombre de variables */ 
	            ROWSPACE       = 2 * nbvar - 1,  /* Nb de contraintes : Anciennement 2 * nbvar - 1 */
	            NZSPACE        = (nbvar * nbvar + 9 * nbvar - 6) / 2, /* Nb d'entrees non nulles        */
				                          /* avant \E8 au lieu de 9 et dans ar_nzspace 3 aulieu de 5 */
	            AR_COLSPACE    = nbvar,
	            AR_ROWSPACE    = nbvar,
	            AR_NZSPACE     = (nbvar * nbvar + 5 * nbvar - 2) / 2;
	const bool  forcer_egalite = FAUX;
	const bool  sol_adm        = VRAI;
	int         status, numcols, numrows, objsen, ar_nzcount;
	double      obj;
	bool        faisable;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
	            x, pi, slack, dj;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
    char        part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE], nom_pb[LONGCHAINE];
    FILE        *fichlp;

	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

	pqe_cplex_set_pb_data_lin2(nbvar, COLSPACE, ROWSPACE, coef_c1z, sous_phi, ordre_poids, q_lin, q_qua, capacite,
                 probname, &numcols, &numrows, &objsen, coef_obj, rhs, sense, 
                 matbeg, matcnt, matind, matval, lb, ub, ctype, cprim, rprim, b_primal,
                 colindex, priority, type_branchement, pr_indices, pr_value);

/* 
	pqe_cplex_set_pb_names(nbvar, ordre_poids, colname, colnamestore, rowname,
	             rownamestore, &colnamespace, &rownamespace);
*/

	lp = CPXcreateprob (env, &status, probname);
	if ( lp == NULL ) pqe_common_problemo("Failed to create LP.\n");

/*
	status = CPXreadcopyprob (env, lp, "cover.lp", NULL);
	CPXsetintparam (env, CPX_PARAM_CLIQUES, -1); 
	CPXsetintparam (env, CPX_PARAM_COVERS, -1); 
	CPXsetintparam (env, CPX_PARAM_MIPSTART, CPX_ON); 
	status = CPXcopymipstart (env, lp, COLSPACE, pr_indices, pr_value);
	if (status != 0) pqe_common_problemo("CPXcopymipstart");
	status = CPXmipopt (env,lp);
*/

 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                        sense, matbeg, matcnt, matind, matval, lb, ub, NULL);

	if (nbfix > 0)
	{
		pqe_cplex_prepare_fixations(nbvar, ordre_poids, nbfix, variables, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval);
		status = CPXaddrows (env, lp, 0, nbfix, nbfix, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, NULL);
		if (status != 0) pqe_common_problemo("addrows");
	}

	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, p_temps_UB, p_meilleur_dual, x, pi, slack, dj, &faisable);

	/* ----------------------------- */
	/* Resolution du probleme en 0-1 */
	/* ----------------------------- */

	pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, type_BB, var_select, "voir.lp", pr_indices, pr_value);

	CPXsetintparam (env, CPX_PARAM_CLIQUES, 0); 
	CPXsetintparam (env, CPX_PARAM_COVERS, 0); 


	if (contrainte)
	{
		pqe_cplex_prepare_nelle_contrainte (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
	                                       	  ar_matval, newrowname, coef_obj, MeilleurPrimal);
		status = CPXaddrows (env, lp, 0, 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	if (forcer_egalite)
	{
		pqe_cplex_preparer_forcer_egalite (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
		                                     ar_matval, newrowname, coef_c1z, ordre_poids, q_qua);
		newrowname = NULL;
		status = CPXaddrows (env, lp, 0, nbvar - 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, 
	                     	 NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	CPXcopyorder(env, lp, nbvar, colindex, priority, NULL);
	if (cutoff) CPXsetdblparam (env, CPX_PARAM_CUTLO, (double) MeilleurPrimal - 0.5);

	if (mode_debug || ESRA)
    {
	    pqe_common_intochar(part1,nbvar);
        pqe_common_intochar(part2,(int) (densite*100));
        pqe_common_intochar(part3,seed);
	    if (MAUVAISE_SOL_INIT) strcpy(nom_pb,"mauvaise_");
        else switch (ordre_var)
	         {
		         case VAR_NAT : strcpy(nom_pb, "naturel_"); break;
                 case VAR_PDS_CR : strcpy(nom_pb, "croissant_"); break;
                 case VAR_PDS_DEC : strcpy(nom_pb, "decroissant_"); break;
             }
	    strcat(nom_pb,part1);
	    strcat(nom_pb,"_");
	    strcat(nom_pb,part2);
	    strcat(nom_pb,"_");
	    strcat(nom_pb,part3);
        strcat(nom_pb,".lp");
        CPXlpwrite (env,lp, nom_pb);
        fichlp = fopen(nom_pb, "a");
        fprintf(fichlp, "\\ nbvar : %d densite : %f germe : %d meilleur primal : %f",
                nbvar, densite, seed, meilleur_primal);
        fclose(fichlp);
    }

	if (!JUSTE_LP) status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans linMM");

	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */

	pqe_cplex_recuperation_solution(env, lp, p_temps_limite, &obj, v_sol, p_nb_noeuds);

	CPXfreeprob(env, &lp); 

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return (*p_temps_limite ? 0.0 : (float) obj);

} /* pqe_cplex_lin2 */


/*****************************************************************************/

double pqe_cplex_lin3 (int N, double_1D sur_phi, double_1D sous_phi,
                          double_1D rlin, double_2D rqua, int CAPA, double_1D v_sol, 
                          double *p_nb_noeuds, bival_1D b_primal, float meilleur_primal,
                          int type_branchement, int type_BB, int var_select, bool cutoff, 
                          bool contrainte, float *p_meilleur_dual, double *p_temps_UB, 
                          bool *p_temps_limite, int nbfix, int seed,
                          float densite, double_1D weight)
{
	/* ------------------------------------------------ */
	/* Toutes les variables autres que env sont locales */
	/* ------------------------------------------------ */

	CPXLPptr    lp             = NULL;
	const int   COLSPACE       = 2 * N - 1,  /* Nombre de variables */ 
	            ROWSPACE       = 2 * N - 1,  /* Nb de contraintes : Anciennement 2 * nbvar - 1 */
	            NZSPACE        = (N * N + 9 * N - 6) / 2, /* Nb d'entrees non nulles                  */
				                          /* avant \E8 au lieu de 9 et dans ar_nzspace 3 aulieu de 5 */
	            AR_COLSPACE    = N,
	            AR_ROWSPACE    = N,
	            AR_NZSPACE     = (N * N + 5 * N - 2) / 2;
	const bool  forcer_egalite = FAUX;
	const bool  sol_adm        = VRAI;
	int         status, numcols = COLSPACE, numrows = ROWSPACE, objsen = CPX_MAX, ar_nzcount, sens;
	double      obj;
	bool        faisable;
	int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices;
	double_1D   coef_obj, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
	            x, pi, slack, dj;
	char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
	char_2D     colname, rowname, newrowname;
    char        part1[LONGCHAINE], part2[LONGCHAINE], part3[LONGCHAINE], nom_pb[LONGCHAINE];
    FILE        *fichlp;

	pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, AR_COLSPACE, AR_ROWSPACE, AR_NZSPACE,
	             &coef_obj, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
	             &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
	             &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
	             &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);

	/* ---------------------------- */
	/* Ecriture et chargement du pb */
	/* ---------------------------- */

    pqe_cplex_set_pb_data_lin3 (N, COLSPACE, ROWSPACE, sur_phi, sous_phi,
                 rlin, rqua, weight, CAPA, &objsen, coef_obj, rhs, sense, matbeg, matcnt, 
                 matind, matval, lb, ub, ctype, meilleur_primal, sens = 1);


/* 
	pqe_cplex_set_pb_names(nbvar, ordre_poids, colname, colnamestore, rowname,
	             rownamestore, &colnamespace, &rownamespace);
*/

	lp = CPXcreateprob (env, &status, probname);
	if ( lp == NULL ) pqe_common_problemo("Failed to create LP in reso_lin3\n");

/*
	status = CPXreadcopyprob (env, lp, "cover.lp", NULL);
	CPXsetintparam (env, CPX_PARAM_CLIQUES, -1); 
	CPXsetintparam (env, CPX_PARAM_COVERS, -1); 
	CPXsetintparam (env, CPX_PARAM_MIPSTART, CPX_ON); 
	status = CPXcopymipstart (env, lp, COLSPACE, pr_indices, pr_value);
	if (status != 0) pqe_common_problemo("CPXcopymipstart");
	status = CPXmipopt (env,lp);
*/

 	status = CPXcopylp (env, lp, numcols, numrows, objsen, coef_obj, rhs,
                        sense, matbeg, matcnt, matind, matval, lb, ub, NULL);

    CPXlpwrite (env,lp, "LIN3.lp");


/*	if (nbfix > 0)
	{ 
		pqe_cplex_prepare_fixations(nbvar, ordre_poids, nbfix, variables, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval);
		status = CPXaddrows (env, lp, 0, nbfix, nbfix, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, NULL); 
 		if (status != 0) pqe_common_problemo("addrows");
 	}
*/
	/* -------------------------------- */
	/* Calcul de la relaxation continue */
	/* -------------------------------- */

	pqe_cplex_relaxation_continue(env, lp, p_temps_UB, p_meilleur_dual, x, pi, slack, dj, &faisable);

	/* ----------------------------- */
	/* Resolution du probleme en 0-1 */
	/* ----------------------------- */

	/* pqe_cplex_parametrage(env, lp, COLSPACE, ctype, sol_adm, type_BB, var_select, "voir.lp", pr_indices, pr_value); */
	pqe_cplex_parametrage(env, lp, COLSPACE, ctype, FAUX, type_BB, var_select, "voir.lp", pr_indices, pr_value);

	CPXsetintparam (env, CPX_PARAM_CLIQUES, 0); 
	CPXsetintparam (env, CPX_PARAM_COVERS, 0); 

/*
	if (contrainte)
	{
		pqe_cplex_prepare_nelle_contrainte (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
	                                       	  ar_matval, newrowname, coef_obj, MeilleurPrimal);
		status = CPXaddrows (env, lp, 0, 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
	if (forcer_egalite)
	{
		pqe_cplex_preparer_forcer_egalite (nbvar, &ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind,
		                                     ar_matval, newrowname, coef_c1z, ordre_poids, q_qua);
		newrowname = NULL;
		status = CPXaddrows (env, lp, 0, nbvar - 1, ar_nzcount, ar_rhs, ar_sense, ar_matbeg, ar_matind, ar_matval, 
	                     	 NULL, newrowname);
		if (status != 0) pqe_common_problemo("addrows");
	}
*/
	/* CPXcopyorder(env, lp, nbvar, colindex, priority, NULL); */
	if (cutoff) CPXsetdblparam (env, CPX_PARAM_CUTLO, (double) MeilleurPrimal - 0.5);

	if (FAUX && (mode_debug || ESRA))
    {
	    pqe_common_intochar(part1,nbvar);
        pqe_common_intochar(part2,(int) (densite*100));
        pqe_common_intochar(part3,seed);
	    if (MAUVAISE_SOL_INIT) strcpy(nom_pb,"mauvaise_");
	    strcat(nom_pb,part1);
	    strcat(nom_pb,"_");
	    strcat(nom_pb,part2);
	    strcat(nom_pb,"_");
	    strcat(nom_pb,part3);
        strcat(nom_pb,".lp");
        CPXlpwrite (env,lp, nom_pb);
        fichlp = fopen(nom_pb, "a");
        fprintf(fichlp, "\\ lin3 nbvar : %d densite : %f germe : %d meilleur primal : %f",
                nbvar, densite, seed, meilleur_primal);
        fclose(fichlp);
    }
    
    print("\n top 1 \n");

	if (!JUSTE_LP) status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans lin3");

    print("\n top 2 \n");

	/* --------------------------- */
	/* Recuperation de la solution */
	/* --------------------------- */

	pqe_cplex_recuperation_solution(env, lp, p_temps_limite, &obj, v_sol, p_nb_noeuds);

	CPXfreeprob(env, &lp); 

	pqe_cplex_free_tab(COLSPACE, ROWSPACE, coef_obj, lb, ub, rhs, ctype, sense,
	                     matbeg, matcnt, colindex, priority, matind, matval, 
	                     cprim, rprim, pr_indices, pr_value, ar_rhs,  ar_sense, 
	                     ar_matbeg, ar_matind, ar_matval, newrowname, colname, 
	                     rowname,  colnamestore, rownamestore, probname, x, pi, slack, dj);

	return (*p_temps_limite ? 0.0 : obj);

} /* pqe_cplex_lin3 */

/*****************************************************************************/

float pqe_cplex_reso_lkp (float_1D coef_obj, float_1D coef_contrainte, int n, int capa, double_1D v_sol)
{
	int        i, mipstat;
	static int num_appel = 1;
	char       lenom[LONGCHAINE], ctype[MACSZ];
	
	pqe_cplex_masquer_affichage(env);
	/* -------------------------------------------------------- */
	/* Part I - Loading, optimizing and obtaining a solution to */
	/* the original problem.                                    */
	/* -------------------------------------------------------- */
	
	pqe_cplex_set_pb_data(coef_obj, coef_contrainte, n, capa,
                            probname, &mac, &mar, &objsen, objx, rhsx, senx, 
	                        matbeg, matcnt, matind, matval, bdl, bdu, &dataname,
	                        &objname, &rhsname, &rngname, &bndname, cname, 
	                        cstore, rname, rstore, &macsz, &marsz, &matsz, 
	                        &cstorsz, &rstorsz);

	sprintf(lenom,"lkp%d.lp",num_appel);
	for (i = 0; i < MACSZ; i++) ctype[i] = 'I';
	
	lp = CPXcreateprob (env, &status, probname);
	if (lp == NULL) pqe_common_problemo("create dans lkp");
	status = CPXcopylp (env, lp, mac, mar, objsen, objx, rhsx, senx, matbeg, 
	                    matcnt, matind, matval, bdl, bdu, NULL);
	status = CPXcopyctype (env, lp, ctype);
	/* if (num_appel > 1) CPXloadbase(env,lp, cstat, rstat); */
	if (mode_debug) CPXlpwrite (env,lp, lenom); 
	
	if (lp == NULL) pqe_cplex_arret_sur_erreur (-1, "lp == NULL dans reso_lkp");

	/* ----------------------------------------- */
	/* Optimize the problem and obtain solution. */
	/* ----------------------------------------- */

	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	CPXsetintparam(env,CPX_PARAM_AGGIND,0);

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans reso_lkp");

	mipstat = CPXgetstat(env,lp);

	if (   ( mipstat != CPXMIP_OPTIMAL ) 
	    && ( mipstat != CPXMIP_OPTIMAL_TOL )
	    && ( mipstat != CPXMIP_NODE_LIM_FEAS) ) pqe_cplex_arret_sur_erreur (mipstat, "mipstat dans reso_lkp");

	if ( mipstat == CPXMIP_NODE_LIM_FEAS) 
	{
		print ("\nNode limit exceeded, integer solution exists...");
	}
	CPXgetmipobjval(env,lp, &obj);
	CPXgetmipx (env,lp, v_sol, 0, CPXgetnumcols(env,lp) - 1);

	num_appel++;
	/* if ((num_appel % 10) == 0) getbase(lp, cstat, rstat); */

	CPXfreeprob(env, &lp); 
	/* print("\r Valeur de la solution exacte : %f      ",(float) obj); */
	if (affichage) fflush(stdout);

	pqe_cplex_affichage_par_defaut(env, affichage);

	return ((float) obj);

} /* pqe_cplex_reso_lkp */


/*****************************************************************************/

double pqe_cplex_reso_souskpi (double_1D coef_obj, double_1D coef_contrainte, int n, int capa, int ii)
{
    CPXLPptr    lp          = NULL;
    const int   COLSPACE    = 2 * n;                   /* Nombre de variables */ 
    int         ROWSPACE    = n + 2;                   /* Nb de contraintes */
    const int   NZSPACE     = (n * n) + n + 1; /* Nb d'entrees non nulles (au plus...) */
    int         j, status, mipstat, numcols = COLSPACE, 
                numrows, objsen;
    double      temps_rel, val_sol_entiere;
    float       resu;
    const bool  sol_adm  = FAUX, politique_BB_1 = FAUX;
    bool        faisable_01 = VRAI;
    int_1D      matbeg, matcnt, matind, colindex, priority, ar_matbeg, ar_matind, pr_indices; 
    double_1D   coef_obj2, lb, ub, rhs, matval, cprim, rprim, ar_rhs, ar_matval, pr_value,
                x, pi, slack, dj; 
    char_1D     ctype, sense, colnamestore, rownamestore, probname, ar_sense; 
    char_2D     colname, rowname, newrowname;
    char        nom[30];

    static int num_appel = 1;
    char       lenom[LONGCHAINE];
    FILE       *verif;

    numrows = ROWSPACE;

    pqe_cplex_init_tab(COLSPACE, ROWSPACE, NZSPACE, 0, 0, 0,
                 &coef_obj2, &lb, &ub, &rhs, &ctype, &sense, &matbeg, &matcnt, &colindex, 
                 &priority, &matind, &matval, &cprim, &rprim, &pr_indices, &pr_value, 
                 &ar_rhs, &ar_sense, &ar_matbeg, &ar_matind, &ar_matval, &newrowname, &colname,  
                 &rowname, &colnamestore, &rownamestore, &probname, &x, &pi, &slack, &dj);


	
    pqe_cplex_masquer_affichage(env);
    /* -------------------------------------------------------- */
    /* Part I - Loading, optimizing and obtaining a solution to */
    /* the original problem.                                    */
    /* -------------------------------------------------------- */
	
    /* print("\n entree dans souskpi, ii = %d\n",ii); */
    
    pqe_cplex_set_souskpi_data(coef_obj, coef_contrainte, n, capa,
                            probname, &mac, &mar, &objsen, objx, rhsx, senx, 
	                        matbeg, matcnt, matind, matval, lb, ub, &dataname,
	                        &objname, &rhsname, &rngname, &bndname, cname, 
	                        cstore, rname, rstore, &macsz, &marsz, &matsz, 
	                        &cstorsz, &rstorsz, ii);

    /* print("\n top 1 \n"); */

    sprintf(lenom,"souskp%d.lp",ii);
    for (j = 0; j < n; j++) ctype[j] = 'B';
	
    lp = CPXcreateprob (env, &status, probname);
    if (lp == NULL) pqe_common_problemo("create dans reso_souskpi");
    status = CPXcopylp (env, lp, mac, mar, objsen, objx, rhsx, senx, matbeg, 
                        matcnt, matind, matval, lb, ub, NULL);
    status = CPXcopyctype (env, lp, ctype);
   
    /* print("\n top 2 \n");  */ 
   
    if (mode_debug) {
        CPXlpwrite (env,lp, lenom); 
         /* print("\n top 2.1 \n"); */
        /* verif = fopen("w", "verifsouskpi.txt");
        fprintf(verif, "ii = %d\n Coef objectif : \n",ii);
        for (j = 0; j < n; j++) {
            fprintf(verif,"%.2f ", coef_obj[j]);
        }
        fprintf(verif, "\n\n Coef contrainte de capacite : \n");
        for (j = 0; j < n; j++) {
            fprintf(verif, "%.2f ", coef_contrainte[j]); 
        }
        fprintf(verif, "\n");    
        fclose(verif); 
         print("\n top 2.2 \n"); */
    }
	
    if (lp == NULL) pqe_cplex_arret_sur_erreur (-1, "lp == NULL dans reso_souskpi");

    /* ----------------------------------------- */
    /* Optimize the problem and obtain solution. */
    /* ----------------------------------------- */

    CPXsetintparam(env,CPX_PARAM_PREIND,0);
    CPXsetintparam(env,CPX_PARAM_AGGIND,0);

    status = CPXmipopt (env,lp);
    pqe_cplex_verif_statut(status, "mipopt dans reso_souskpi");

    mipstat = CPXgetstat(env,lp);

    if (   ( mipstat != CPXMIP_OPTIMAL ) 
        && ( mipstat != CPXMIP_OPTIMAL_TOL )
        && ( mipstat != CPXMIP_NODE_LIM_FEAS) ) 
             pqe_cplex_arret_sur_erreur (mipstat, "mipstat dans reso_souskpi");

    if (mipstat == CPXMIP_NODE_LIM_FEAS) 
    {
        print ("\nNode limit exceeded, integer solution exists...");
    }
    CPXgetmipobjval(env, lp, &obj);

    /* CPXgetmipx (env, lp, v_sol, 0, CPXgetnumcols(env,lp) - 1); */

    num_appel++;

    CPXfreeprob(env, &lp); 
    /* print("\r Valeur de la solution exacte (souskpi) : %.2f      ", obj); */
    if (affichage) fflush(stdout);

    pqe_cplex_affichage_par_defaut(env, affichage);

    return (obj);

} /* pqe_cplex_reso_souskpi */

/*****************************************************************************/

float pqe_cplex_reso_2lkp (float_1D coef_obj, float_1D coef_contrainte1, float_1D coef_contrainte2, int n,
                             float capa1, float capa2 , double_1D v_sol)
{
	int        i, mipstat;
	static int num_appel = 1;
	char       lenom[LONGCHAINE], ctype[MACSZ];
	
	pqe_cplex_masquer_affichage(env);
	/* -------------------------------------------------------- */
	/* Part I - Loading, optimizing and obtaining a solution to */
	/* the original problem.                                    */
	/* -------------------------------------------------------- */
	
	pqe_cplex_set_pb_data_lkp2(coef_obj, coef_contrainte1, coef_contrainte2, n,
	                        capa1, capa2, probname, &mac, &mar, &objsen, objx, 
                                rhsx, senx, matbeg, matcnt, matind, matval, bdl, bdu, 
                                &dataname, &objname, &rhsname, &rngname, &bndname, 
                                cname, cstore, rname, rstore, &macsz, &marsz, &matsz, 
	                        &cstorsz, &rstorsz);

	sprintf(lenom,"2lkp%d.lp",num_appel);
	for (i = 0; i < MACSZ; i++) ctype[i] = 'I';
	
	lp = CPXcreateprob (env, &status, probname);
	if (lp == NULL) pqe_common_problemo("create dans 2lkp");
	status = CPXcopylp (env, lp, mac, mar, objsen, objx, rhsx, senx, matbeg, 
	                    matcnt, matind, matval, bdl, bdu, NULL);
	status = CPXcopyctype (env, lp, ctype);
	/* if (num_appel > 1) CPXloadbase(env,lp, cstat, rstat); */
	if (mode_debug) CPXlpwrite (env,lp, lenom); 
	
	if (lp == NULL) pqe_cplex_arret_sur_erreur (-1, "lp == NULL dans reso_lkp");

	/* ----------------------------------------- */
	/* Optimize the problem and obtain solution. */
	/* ----------------------------------------- */

	CPXsetintparam(env,CPX_PARAM_PREIND,0);
	CPXsetintparam(env,CPX_PARAM_AGGIND,0);

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans reso_2lkp");

	mipstat = CPXgetstat(env,lp);

	if (   ( mipstat != CPXMIP_OPTIMAL ) 
	    && ( mipstat != CPXMIP_OPTIMAL_TOL )
	    && ( mipstat != CPXMIP_NODE_LIM_FEAS) ) pqe_cplex_arret_sur_erreur (mipstat, "mipstat dans reso_2lkp");

	if ( mipstat == CPXMIP_NODE_LIM_FEAS) 
	{
		print ("\nNode limit exceeded, integer solution exists...");
	}
	CPXgetmipobjval(env,lp, &obj);
	CPXgetmipx (env,lp, v_sol, 0, CPXgetnumcols(env,lp) - 1);

	num_appel++;
	/* if ((num_appel % 10) == 0) getbase(lp, cstat, rstat); */

	CPXfreeprob(env, &lp); 
	print("\r Valeur de la solution exacte : %f      ",(float) obj); 
	if (affichage) fflush(stdout);

	pqe_cplex_affichage_par_defaut(env, affichage);

	return ((float) obj);

} /* pqe_cplex_reso_2lkp */

/*****************************************************************************/

void pqe_cplex_populate_lp01 (CPXCENVptr env, CPXLPptr lp, int nblignes,
                               int nbcols, float_1D coef_obj, float_2D mat_cont,
                               float_1D vec_rhs, char_1D sense, char_2D colname)
{
    /* Variables suppos\E9es toutes 0-1          */
    /* User's manual de cplex p. 107           */

    int        i, j, status = 0, nzcnt;
    double_1D  lb, ub, cost, rhs, rmatval;
    int_1D     rmatbeg, rmatind;
    char_1D    ctype;

    pqe_init_1D ((void **) &lb, sizeof(double), nbcols);
    pqe_init_1D ((void **) &ub, sizeof(double), nbcols);
    pqe_init_1D ((void **) &ctype, sizeof(char), nbcols);
    pqe_init_1D ((void **) &cost, sizeof(double), nbcols);
    pqe_init_1D ((void **) &rhs, sizeof(double), nblignes);
    pqe_init_1D ((void **) &rmatval, sizeof(double), nblignes * nbcols);
    pqe_init_1D ((void **) &rmatbeg, sizeof(int), nblignes);
    pqe_init_1D ((void **) &rmatind, sizeof(int), nblignes * nbcols);


    for (j = 0; j < nbcols; j++) {
        lb[j] = 0.0;
        ub[j] = 1.0;
        ctype[j] = 'B';
        cost[j] = (double) coef_obj[j];
    }

    /* --------------------- */
    /* Creation des colonnes */
    /* --------------------- */

    status = CPXnewcols(env, lp, nbcols, cost , lb, ub, ctype, colname);
    pqe_cplex_verif_statut(status, "populate_lp01_newcols");

    /* Calcul du nb de coefs non nuls (nzcnt) et remplissage de */
    /* rmatbeg, rmatind, rmatval, sense, rhs                    */

    nzcnt = 0;

    for (i = 0; i < nblignes; i++) {
        rhs[i]   = (double) vec_rhs[i];
        rmatbeg[i] = nzcnt;
        for (j = 0; j < nbcols; j++) {
            if (mat_cont[i][j] != 0.0) {
                rmatval[nzcnt] = (double) mat_cont[i][j];
                rmatind[nzcnt] = j;
                nzcnt++;
            }
        }
    }

    /* ------------------- */
    /* Creation des lignes */
    /* ------------------- */

    status = CPXaddrows (env, lp, 0, nblignes, nzcnt, rhs, sense, rmatbeg,
                         rmatind, rmatval, NULL, NULL);
    pqe_cplex_verif_statut(status, "populate_lp01_addrows");

    pqe_init_free_1D (lb);
    pqe_init_free_1D (ub);
    pqe_init_free_1D (ctype);
    pqe_init_free_1D (cost);
    pqe_init_free_1D (rhs);
    pqe_init_free_1D (rmatval);
    pqe_init_free_1D (rmatbeg);
    pqe_init_free_1D (rmatind);

} /* pqe_cplex_populate_lp01 */

/*****************************************************************************/

float pqe_cplex_reso_lp01(int nblignes, int nbcols, float_1D coef_obj,
                          float_2D mat_cont, float_1D vec_rhs,
                          int *p_statut, double_1D vec_sol, char_1D sense,
                          bool maximisation, char_2D colname,
                          bool sol_adm, double_1D primal, bool *p_temps_limite,
			  double *p_nb_noeuds)
{
    /* R\E9solution d'un PL en 0-1 quelconque. */
    /* p_statut \E0 la fin de la proc\E9dure vaut
            -> 0 si infaisabilit\E9e prouv\E9e,
            -> 1 si solution optimale existe et trouv\E9e
            -> 2 si infaisabilit\E9 non prouv\E9e mais limite de noeuds atteinte
    */
/*printf("on rentre dans pqe_cplex_reso_lp01\n");*/
    int        i, mipstat, status;
    static int num_appel = 1;
    char       lenom[LONGCHAINE];
    double     obj;
    int_1D     pr_indices;


    pqe_init_1D ((void **) &pr_indices, sizeof(int), nbcols);
    for (i = 0; i < nbcols; i++) pr_indices[i] = i;

    pqe_cplex_masquer_affichage(env);

    /* -------------------- */
    /* Cr\E9ation du probl\E8me */
    /* -------------------- */
	
    lp = CPXcreateprob (env, &status, probname);
    if (lp == NULL) pqe_common_problemo("create dans reso_lp01");

    pqe_cplex_populate_lp01(env, lp, nblignes, nbcols, coef_obj,
                            mat_cont, vec_rhs, sense, colname);

    if (maximisation) CPXchgobjsen(env, lp, CPX_MAX);
    else              CPXchgobjsen(env, lp, CPX_MIN);

    sprintf(lenom,"lp01_%d.lp",num_appel);
    
    if (mode_debug) CPXlpwrite (env,lp, lenom);
    

/*
     v\E9rifier les chargements de bases anciennes pour acc\E9l\E9ration
     if (num_appel > 1) CPXloadbase(env,lp, cstat, rstat);

*/

	/* ----------------------------------------- */
	/* Optimize the problem and obtain solution. */
	/* ----------------------------------------- */

/*
    Pre-processing de cplex : \E0 essayer \E0 1 pour voir
*/
    CPXsetintparam(env,CPX_PARAM_PREIND,0);
    CPXsetintparam(env,CPX_PARAM_AGGIND,0);


/*    status = CPXwriteprob (env, lp, "myprob.lp", "LP");
    exit(1);
*/
    if (sol_adm) {
        //CPXsetintparam(env, CPX_PARAM_MIPSTART, CPX_ON);
	CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON); 
        status = CPXcopymipstart (env, lp, nbcols, pr_indices, primal);
        if (status != 0) pqe_common_problemo("PB dans CPXcopymipstart");
    }

    CPXsetdblparam (env, CPX_PARAM_TILIM, T_LIMITE);

    status = CPXmipopt (env,lp);
    pqe_cplex_verif_statut(status, "mipopt dans reso_lp_01");

    if ((status == CPXMIP_TIME_LIM_FEAS) || (status == CPXMIP_TIME_LIM_INFEAS))
        *p_temps_limite = VRAI;

    mipstat = CPXgetstat(env,lp);

    if (mipstat == CPXMIP_INFEASIBLE) {
        *p_statut = 0;
    }
    else {
        if (   ( mipstat != CPXMIP_OPTIMAL )
	        && ( mipstat != CPXMIP_OPTIMAL_TOL )) *p_statut = 2;
        else {
            *p_statut = 1;
            CPXgetmipobjval(env, lp, &obj);
            CPXgetmipx (env, lp, vec_sol, 0, CPXgetnumcols(env, lp) - 1);
            print("\r Valeur de la solution exacte : %f\n ", (float) obj);
		
            for(i=0 ; i < CPXgetnumcols(env, lp); i++){
                print("\rcolonne %d:   vect_sol[%d] = %f\n",i, i, vec_sol[i]);
            }
		/* decalage d'indices sur vect_sol */
	/*	for(i = CPXgetnumcols(env, lp); i >=1; i--){
		    vec_sol[i] = vec_sol[i-1];
		    print("\r vec_sol[%d] = %f\n",i, vec_sol[i]);
		}*/
            if (affichage) fflush(stdout);
        }
    }

    *p_nb_noeuds = (double) CPXgetnodecnt (env, lp);

    num_appel++;
	/* if ((num_appel % 10) == 0) getbase(lp, cstat, rstat); */

    CPXfreeprob(env, &lp); 
    pqe_cplex_affichage_par_defaut(env, affichage);

    pqe_init_free_1D(pr_indices);

    return ((float) obj);

} /* pqe_cplex_reso_lp01 */

/*****************************************************************************/

void pqe_cplex_populate_lp_mixte (CPXCENVptr env, CPXLPptr lp, int nblignes,
                               int nbcols, float_1D coef_obj, float_2D mat_cont,
                               float_1D vec_rhs, char_1D sense, char_1D ctype,
                               double_1D lb, double_1D ub, char_2D colname)
{
    /* User's manual de cplex p. 107           */

    int        i, j, status = 0, nzcnt;
    double_1D  cost, rhs, rmatval;
    int_1D     rmatbeg, rmatind;

    pqe_init_1D ((void **) &cost, sizeof(double), nbcols);
    pqe_init_1D ((void **) &rhs, sizeof(double), nblignes);
    pqe_init_1D ((void **) &rmatval, sizeof(double), nblignes * nbcols);
    pqe_init_1D ((void **) &rmatbeg, sizeof(int), nblignes);
    pqe_init_1D ((void **) &rmatind, sizeof(int), nblignes * nbcols);


    for (j = 0; j < nbcols; j++) cost[j] = (double) coef_obj[j];

    /* --------------------- */
    /* Creation des colonnes */
    /* --------------------- */

    status = CPXnewcols(env, lp, nbcols, cost , lb, ub, ctype, colname);
    pqe_cplex_verif_statut(status, "populate_lp_mixte_newcols");

    /* Calcul du nb de coefs non nuls (nzcnt) et remplissage de */
    /* rmatbeg, rmatind, rmatval, sense, rhs                    */

    nzcnt = 0;

    for (i = 0; i < nblignes; i++) {
        rhs[i]   = (double) vec_rhs[i];
        rmatbeg[i] = nzcnt;
        for (j = 0; j < nbcols; j++) {
            if (mat_cont[i][j] != 0.0) {
                rmatval[nzcnt] = (double) mat_cont[i][j];
                rmatind[nzcnt] = j;
                nzcnt++;
            }
        }
    }

    /* ------------------- */
    /* Creation des lignes */
    /* ------------------- */

    status = CPXaddrows (env, lp, 0, nblignes, nzcnt, rhs, sense, rmatbeg,
                         rmatind, rmatval, NULL, NULL);
    pqe_cplex_verif_statut(status, "populate_lp_mixte_addrows");

    pqe_init_free_1D (cost);
    pqe_init_free_1D (rhs);
    pqe_init_free_1D (rmatval);
    pqe_init_free_1D (rmatbeg);
    pqe_init_free_1D (rmatind);

} /* pqe_cplex_populate_lp_mixte */

/*****************************************************************************/

float pqe_cplex_reso_lp_mixte(int nblignes, int nbcols, float_1D coef_obj,
                           float_2D mat_cont, float_1D vec_rhs,
                           int *p_statut, double_1D vec_sol, double_1D vec_dul,                            
			   char_1D sense, char_1D ctype, double_1D lb, double_1D ub,
                           bool resolution_mixte, bool maximisation,
                           char_2D colname, int_1D copycstat, int_1D copyrstat,
                           int_1D getcstat, int_1D getrstat)
{
    /* R\E9solution d'un PL en 0-1 quelconque. */
    /* On ne prend pas en compte le cas d'un optimum non borne :
       Le pb est suppose born\E9
       p_statut \E0 la fin de la proc\E9dure vaut
            -> 0 si infaisabilit\E9e en 0-1 prouv\E9e,
            -> 1 si solution optimale existe et trouv\E9e
            -> 2 si infaisabilit\E9 non prouv\E9e mais limite de noeuds atteinte
            -> 3 si infaisabilit\E9 en continu prouv\E9e
	    -> 4 si limite de temps depassee : abandon
       resolution_mixte indique si on doit lancer une reso en 0-1 mixte (mipopt)
           ou une r\E9solution en continu (ce qui suppose ctype = 'C' pour toutes
           les variables du pb
    */

    int        i, mipstat, lpstat, solstat, status;
    static int num_appel = 1;
	char       lenom[LONGCHAINE];
    double     obj;
    double_1D  slack, dj;



	pqe_cplex_masquer_affichage(env);
	CPXsetdblparam (env, CPX_PARAM_TILIM, T_LIMITE);

	/* -------------------- */
	/* Cr\E9ation du probl\E8me */
	/* -------------------- */
	
 	lp = CPXcreateprob (env, &status, probname);
	if (lp == NULL) pqe_common_problemo("create dans reso_lp_mixte");

	pqe_cplex_populate_lp_mixte(env, lp, nblignes, nbcols, coef_obj,
                            mat_cont, vec_rhs, sense, ctype, lb, ub, colname);

    if (maximisation) CPXchgobjsen (env, lp, CPX_MAX);
    else              CPXchgobjsen (env, lp, CPX_MIN);

	sprintf(lenom, "lp_mixte_%d.lp", num_appel);
	if (mode_debug) CPXlpwrite (env, lp, lenom);
	/*
	CPXlpwrite(env, lp, lenom);
	*/

/*
     v\E9rifier les chargements de bases anciennes pour acc\E9l\E9ration
     if (num_appel > 1) CPXloadbase(env,lp, cstat, rstat);

*/

	/* ----------------------------------------- */
	/* Optimize the problem and obtain solution. */
	/* ----------------------------------------- */

/*
    Pre-processing de cplex : \E0 essayer \E0 1 pour voir
*/

    CPXsetintparam(env,CPX_PARAM_PREIND,0);
    CPXsetintparam(env,CPX_PARAM_AGGIND,0);

/*
 * CPXlpwrite (env,lp, "nonag.lp");
 *
 */
   /* Now copy the basis * /

       status = CPXcopybase (env, lp, copycstat, copyrstat);
       if ( status ) {
        printf ("Failed to copy the basis.\n");
      }
  */ 

     /* CPXsetintparam(env,CPX_PARAM_ITLIM,100); * /
	CPXsetdblparam (env, CPX_PARAM_TILIM, 0.1);
	*/
    
    if (resolution_mixte) {

	status = CPXmipopt (env,lp);
	pqe_cplex_verif_statut(status, "mipopt dans reso_lp_mixte");

	mipstat = CPXgetstat(env,lp);
	if ((mipstat == CPXMIP_TIME_LIM_FEAS) || (mipstat == CPXMIP_TIME_LIM_INFEAS))
	    *p_statut = 4;

        if (mipstat == CPXMIP_INFEASIBLE) {
            *p_statut = 0;
        }
        else {
            if (*p_statut != 4) {
	        if (   ( mipstat != CPXMIP_OPTIMAL )
	            && ( mipstat != CPXMIP_OPTIMAL_TOL )) *p_statut = 2;
                else {
                *p_statut = 1;
	            CPXgetmipobjval(env, lp, &obj);
	            CPXgetmipx (env, lp, vec_sol, 0, CPXgetnumcols(env, lp) - 1);
         	    print("\r Valeur de la solution exacte : %f ", (float) obj);
                if (affichage) fflush(stdout);
            }
            }
        }
    }
    else {

        /* ------------------------------ */
	    /* Calcul de la solution continue */
	    /* ------------------------------ */

        lpstat = CPXchgprobtype (env, lp, CPXPROB_LP);

        print("\n Valeur de la relaxation continue : ");
    	lpstat = CPXlpopt (env, lp);
        print(" lpstat = %d ", lpstat);
	    if (lpstat) pqe_common_problemo("relaxation continue");
	    if ((lpstat = CPXgetstat(env, lp)) != CPX_STAT_INFEASIBLE) *p_statut = 1;
        else *p_statut = 3;
	    if (lpstat == CPX_STAT_INFEASIBLE) print("Infaisabilite en continu prouv\E9ee\n");

        pqe_init_1D ((void **) &slack, sizeof(double), nbcols + nblignes);
        pqe_init_1D ((void **) &dj, sizeof(double), nbcols + nblignes);

	    if ((lpstat = CPXsolution (env, lp, &solstat, &obj, vec_sol, vec_dul, slack, dj)))
		    pqe_common_problemo("Cpxsolution");
	    if (*p_statut = 1) print (" %lf      \n ", obj); if (affichage) fflush(stdout);

        pqe_init_free_1D(slack);
        pqe_init_free_1D(dj);
    }

     status = CPXgetbase (env, lp, getcstat, getrstat);
	num_appel++;
	/* if ((num_appel % 10) == 0) getbase(lp, cstat, rstat); */

	CPXfreeprob(env, &lp); 
	pqe_cplex_affichage_par_defaut(env, affichage);

    return ((float) obj);

} /* pqe_cplex_reso_lp_mixte */

/*****************************************************************************/

void pqe_cplex_reso_pq_fichier(double *p_val, double *p_tps, char_1D nom, bool entier,
                               double *p_nb_noeuds, bool *p_temps_limite, bool sol_adm, int n,
			       double_1D primalx, int_1D pr_indices)
{
    int solstat, mipstat;
    /* double x[300], pi[300], slack[300], dj[300]; */


    pqe_common_init_temps(p_tps);

    lp = CPXcreateprob(env, &status, nom);
	if (lp == NULL) pqe_common_problemo("create dans relaxation_continue");

    status = CPXreadcopyprob(env, lp, nom, NULL);
    pqe_cplex_verif_statut(status, "readcopy dans relaxation_continue");

    if (entier) {
	if (sol_adm) {
            //CPXsetintparam(env, CPX_PARAM_MIPSTART, CPX_ON);
		CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON); 
	    status = CPXcopymipstart (env, lp, n, pr_indices, primalx);
	    if (status != 0) pqe_common_problemo("PB dans CPXcopymipstart");
        }

	CPXsetdblparam (env, CPX_PARAM_TILIM, T_LIMITE);
	CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON); 
        status = CPXmipopt(env, lp);
	mipstat = CPXgetstat(env,lp);
	if ((mipstat == CPXMIP_TIME_LIM_FEAS) || (mipstat == CPXMIP_TIME_LIM_INFEAS))
	    *p_temps_limite = VRAI;
        pqe_cplex_verif_statut(status, "mipopt dans relaxation_continue");
        status = CPXgetmipobjval (env, lp, p_val);
        pqe_cplex_verif_statut(status, "getmipobjval dans relaxation_continue");
        *p_nb_noeuds = (double) CPXgetnodecnt (env, lp);

    }
    else {
        status = CPXqpopt(env, lp);
        pqe_cplex_verif_statut(status, "qpopt dans relaxation_continue");
        status = CPXsolution(env, lp, &solstat, p_val, NULL, NULL, NULL, NULL);
        pqe_cplex_verif_statut(status, "solution dans relaxation_continue");
    }

    /*    status = CPXmipopt(env, lp);
          CPXgetmipobjval(env, lp, &obj) ; */

    print(" Solution trouvee : %.3f \n", *p_val);

    pqe_common_calcule_temps(p_tps);

} /* pqe_cplex_relaxation_continue */

/*****************************************************************************/

void pqe_cplex_fin()
{
	status = CPXfreeprob (env, &lp);
	status = CPXcloseCPLEX (&env);
      
} /* pqe_cplex_fin */

