/* --------------------------------- */
/* Fichier pqe_prototypes.h         */
/* Prototype des fonctions utilisees */
/* --------------------------------- */

#include "pqe_types.h"


/* ************************** /
/ * Fichier : pqe_branch.c * /
/ *************************** /

void pqe_branch_trivar (int, int);
bool pqe_branch_amelioration_primal (float, int_1D, bool);
void pqe_branch_examen (void);
void pqe_branch_init (int, la_st_var_1D, bool *, bool *, int *, int *, double *);
void pqe_branch_fixation (bool, int *);
void pqe_branch_and_bound (void);
*/

/**********************/
/* Fichier : common.c */
/**********************/

void  pqe_common_init_temps(double *);
void  pqe_common_calcule_temps(double *);
float pqe_common_frand (void);
int   pqe_common_irand (int);
void  pqe_common_init_rand (int);
float pqe_common_f_economique (bival_1D, float, float_1D, float_2D);
void  pqe_common_intochar (char *, int);
void  pqe_common_init_cluster_fix (int, int);
void  pqe_common_problemo (char *);
void  pqe_common_probleme (int n, int m, int germe, char *);


/*************************/
/* Fichier : pqe_init.c */
/*************************/

void  pqe_init_1D (void **, int, int);
void  pqe_init_2D (void ***, int, int, int);
void  pqe_init_3D (void ****, int, int, int, int);
void  pqe_init_free_1D (void *);
void  pqe_init_free_2D (void **, int);
void  pqe_init_free_3D (void ***, int, int);
void  pqe_init_dynamique (void);
void  pqe_init_1er_cluster_2D (void ***, int, int, int, int);
void  pqe_init_1er_cluster_3D (void ****, int, int, int, int, int, int);
void  pqe_init_free_1er_cluster_3D (void ***, int, int, int);
void  pqe_init_alloc_var_relag (int, int, int, int);
void  pqe_init_free (void);
void  pqe_init_free_var_relag (int, int, int);
int   pqe_init_genere_contrainte (bool);
void  pqe_init_affiche_probl (void);
void  pqe_init_ecrit_probl (char *);
void  pqe_init_ecrit_probleme (char *, st_tri_1D, float_1D);
void  pqe_init_test_cplex (void);
void  pqe_init_JE_Sourour (int, FILE **, int, float, int, char *, bool, st_tri_1D);
void  pqe_init_construit_exemple (int);
void  pqe_init_lecture (void);
void  pqe_init_ecriture (bool);
float pqe_init_verif_densite (float_1D q_lin, float_2D q_qua);
void  pqe_init_transformation_pb (void);
float pqe_init_alpha (float, bival_1D, float *);

/* ************************** /
/ * Fichier : pqe_primal.c * /
/ *************************** /

void  pqe_primal_main (float *, bival_1D, float *, la_st_ord_1D);
int   pqe_primal_capacite_solution (int_1D, int_1D);
bool  pqe_primal_admissible(int_1D, int_1D, int);
void  pqe_primal_enum_recurs (int, int_1D, float *);
void  pqe_primal_calc_borne_inf (float *, int_1D, double_1D, float_1D, float_2D, int_1D, int, la_st_var *, bool);
void  pqe_primal_ameliore (float *, int_1D, float_1D, float_2D, int_1D, int, la_st_var *, bool);


/ *************************** /
/ * Fichier : pqe_recuit.c * /
/ *************************** /

void  pqe_recuit_augmente(float *, int, int_1D);
void  pqe_recuit_diminue (float *, int, int_1D);
void  pqe_recuit_echanger_articles (int *, int *, float *, int_1D);
void  pqe_recuit_mise_a_jour (int *, float *, int_1D, int_1D, float *, float *, int, int_1D, float, int *, int *);
void  pqe_recuit_sol_init (float *, int_1D);
float pqe_recuit_bonne_temp (float *, int_1D, float *, int_1D);
float pqe_recuit_reglage_temp (int, float);
bool  pqe_recuit (float *, int_1D);


/ ************************** /
/ * Fichier : pqe_relag.c * /
/ ************************** /

void  pqe_relag_verif_quad (void);
void  pqe_relag_heuristique (float_1D);
void  pqe_relag_decomp_rel (float_1D);
int   pqe_relag_capasurX (int, int **);
void  pqe_relag_echange (int, int, int);
void  pqe_relag_quick_sort_Y (int, int, int);
float pqe_relag_linear_pqe (int, int, int);
bool  pqe_relag_rentre (int, int, int, int);
void  pqe_relag_mettre (int, int, int, int *, float *);
void  pqe_relag_couper (int, int, int, float *);
void  pqe_relag_met_a_zero (int, int, int);
void  pqe_relag_NKR (int, int, int, int, float *);
float pqe_relag_partie_constante (int, int **);
void  pqe_relag_enum_X (int, int, float_1D);
bool  pqe_relag_verifie_contrainteA (int, int **);
bool  pqe_relag_verifie_contrainteB (int, int **);
void  pqe_relag_echange_tvI (int, int);
void  pqe_relag_tri_decroissant_tvI (int, int);
void  pqe_relag_tri_croissant_tvI (int, int);
void  pqe_relag_reso_relachee (int, int, int, float *);
float pqe_relag_respect_contrainte_AouB (int, int, int *, int **, int, int *);
void  pqe_relag_construction_tvI (int, int, int *, int **, int, bool *);
float pqe_relag_CKP_FA (int, int **, int, int);
float pqe_relag_CKP_FB (int, int **, int, int);
void  pqe_relag_enum_contrainteA (int, int, float_1D);
void  pqe_relag_enum_contrainteB (int, int, float_1D);
void  pqe_relag_calcul_dual (int, int, bool, int *, int *, int, float *);
void  pqe_relag_init (void);
int   pqe_relag_maj_copies (int, int *);
void  pqe_relag_nouveau_dual(float *, float, int *, int *, int, int *, float, float *, char *, clock_t *,
                              clock_t *, clock_t *, clock_t *, clock_t *, int, int, float);
int   pqe_relag (int *, int, int, int, int *, float *);
void  pqe_relag_main (int *, int *, float *);


/ *************************** /
/ * Fichier : pqe_clique.c * /
/ *************************** /

maillon * pqe_clique_new(void);
void      pqe_clique_free(maillon **);
void      pqe_clique_copy(maillon *, maillon **);
int       pqe_clique_comp(const void *, const void *);
int       pqe_clique(void);


/***************************/
/* Fichier : pqe_cplex.c */
/***************************/

void  pqe_cplex_init(void);
void  pqe_cplex_set_pb_data(float coef_obj[], float coef_contrainte[], int n, int capa, char probname[],
                              int *p_mac, int *p_mar, int *p_objsen, double objx[], double rhsx[],
                              char senx[], int matbeg[], int matcnt[], int matind[], double matval[],
                              double bdl[], double bdu[], char **p_dataname, char **p_objname, char **p_rhsname,
                              char **p_rngname, char **p_bndname, char *cname[], char cstore[], char *rname[],
                              char rstore[], int *p_macsz, int *p_marsz, int *p_matsz, unsigned *p_cstorsz,
                              unsigned *p_rstorsz);
void  pqe_cplex_set_pb_data_demi_lin (int, float_1D, st_tri_1D, float_1D, float_2D, int, char_1D, int *, int *,
                                      int *, double_1D, double_1D, char_1D, int_1D, int_1D, int_1D, double_1D, 
                                      double_1D, double_1D, char_1D, double_1D, double_1D, bival_1D, int_1D, int_1D, int);
void  pqe_cplex_set_pb_data_lin1 (int, float_1D, st_tri_1D, float_1D, float_2D, int, char_1D, int *, int *,
                                      int *, double_1D, double_1D, char_1D, int_1D, int_1D, int_1D, double_1D, 
                                      double_1D, double_1D, char_1D, double_1D, double_1D, bival_1D, int_1D, int_1D, int,
                                      int_1D, double_1D);
void  pqe_cplex_set_pb_data_lin_croisee (int, float_1D, float_1D, float_1D, int, int, st_tri_1D, float_1D, float_2D, int,
                                      char_1D, int *, double_1D, double_1D, char_1D, int_1D, int_1D, int_1D, double_1D, 
                                      double_1D, double_1D, char_1D, double_1D, double_1D, bival_1D, int_1D, int_1D, int);
float pqe_cplex_lin2 (int, float_1D, float_1D , st_tri_1D ,
                          float_1D , float_2D , int , double_1D , 
                          double *, bival_1D , float ,
                          int , int , int , bool , 
                          bool , float *, double *, 
                          bool *, la_st_var_1D , int , int ,
                          float , int );
double pqe_cplex_lin3 (int, double_1D, double_1D, double_1D, double_2D, int, double_1D,
                        double *, bival_1D, float, int, int, int, bool, bool, float *, double *,
                        bool *, int, int, float, double_1D);

void  pqe_cplex_set_pb_names_demi(int, st_tri_1D, char_2D, char_1D, char_2D, char_1D, unsigned *, unsigned *);
void  pqe_cplex_set_pb_names(int, st_tri_1D, char_2D, char_1D, char_2D, char_1D, unsigned *, unsigned *);
void  pqe_cplex_prepare_nelle_contrainte(int, int *, double_1D, char_1D, int_1D, int_1D,
                                           double_1D, char_2D, double_1D, float);
float pqe_cplex_demi_lin (int, float_1D, st_tri_1D, float_1D, float_2D, int, double_1D, double *, bival_1D, float,
                            int, int, int, bool, bool);
float pqe_cplex_LPi (bool, int, int, float_1D, st_tri_1D, float_1D, float_2D, int_1D, int, float, bival_1D, bool *,
                         int *, int *, int *, int);
double pqe_cplex_max_LPi (bool, int, int, double_1D, double_1D, double_1D, double_2D, double_1D, int,
                         double, bool *, int *, int *, int *, int);
double pqe_cplex_min_LPi (bool, int, int, double_1D, double_1D, double_1D, double_2D, double_1D, int,
                         double, bool *, int *, int *, int *, int);
float  pqe_cplex_lin1 (int, float_1D, st_tri_1D, float_1D, float_2D, int, double_1D, double *, bival_1D, float,
                          int, int, int, bool, bool, float *, double *, bool *);
float  pqe_cplex_lin_classique (int, float_1D, float_2D, int, int_1D, double_1D, double *, int, int,
                                 int, float *, double *, bool, bool *);
float  pqe_cplex_lin_croisee (int, float_1D, float_1D, float_1D, st_tri_1D, float_1D, float_2D, int,
                               double_1D, double *, 
                               bival_1D, float, int, int, int, bool, bool);
float  pqe_cplex_reso_lkp (float coef_obj[], float coef_contrainte[], int n, int capa, double v_sol[]);
double pqe_cplex_reso_souskpi(double coef_obj2[], double a[], int n, int capa, int i);
void   pqe_cplex_fin (void);
float  pqe_cplex_LPi_stocha (bool prime, int i, int nbvar, float_1D cmax, st_tri_1D ordre_poids,
                               float_1D q_lin, float_2D q_qua, float_1D cout, float CAP,
	                       float_1D moyenne, float k, float meilleur_primal, bival_1D b_primal,
                               bool *p_faisable_continu, int *p_nbBB, int *p_nbsol01,
                               int *p_nbopt01);
float pqe_cplex_lin2_stocha (int nbvar, float_1D coef_c1z, float_1D cmin, st_tri_1D ordre_poids,
                               float_1D q_lin, float_2D q_qua, float_1D cout, float CAP,
                               float_1D moyenne, float k, double_1D v_sol,
                               double *p_nb_noeuds, bival_1D b_primal,
                               int type_branchement, int type_BB, int var_select, bool cutoff,
                               bool contrainte, float *p_meilleur_dual, double *p_temps_UB,
                               bool *p_temps_limite, la_st_var_1D variables, int nbfix,
                               bool *p_sol_positive, int jeu, float valeur_coeff, float alpha);
float pqe_cplex_reso_2lkp (float_1D coef_obj, float_1D coef_contrainte1, float_1D coef_contrainte2, int n,
                             float capa1, float capa2 , double_1D v_sol);
float pqe_cplex_reso_lp01(int nblignes, int nbcols, float_1D coef_obj,
                           float_2D mat_cont, float_1D vec_rhs,
                           int *p_statut, double_1D vec_sol, char_1D sense,
                           bool maximisation, char_2D colname, bool sol_adm, double_1D primal, 
			   bool *p_temps_limite, double *p_nb_noeuds);
float pqe_cplex_reso_lp_mixte(int nblignes, int nbcols, float_1D coef_obj,
                           float_2D mat_cont, float_1D vec_rhs,
                           int *p_statut, double_1D vec_sol, double_1D vec_dul,
                          char_1D sense, char_1D ctype, double_1D lb, double_1D ub,
                           bool resolution_mixte, bool maximisation,
                           char_2D colname, int_1D copystat, int_1D copyrstat,
                           int_1D getcstat, int_1D getrstat);


void pqe_cplex_reso_pq_fichier(double *p_val, double *p_tps_opt, char_1D nom, bool entier, 
                               double *p_nb_noeuds, bool *p_temps_limite, bool sol_adm, int n, double_1D primalx,
                               int_1D pr_indices);


/**********************/
/* Fichier : pqe_bb.c */
/**********************/

void pqe_bb_trouver_var_frac(double_1D vec, int n, int *p_var, double *p_val);
void pqe_bb_branch_and_bound(int n, int m, int germe, double_1D c, double_1D d, double_1D b,
                             double_2D a,
                             int_1D u, double borne_inf, double_1D x_inf, double borne_racine,
                             int var_fract_racine, double val_fract_racine,
                             bool *p_opt_trouve, double *p_valopt, double_1D x_opt, double *p_temps,
			     int_1D copycstat, int_1D copyrstat, double_1D multiplicateur_racine,
			     double *p_nb_noeuds, double t_limite, bool *p_temps_lim_atteint);
void pqe_bb_branch_and_bound_borne1(int n, int m, int germe, double_1D c, double_1D d, double_1D b,
                             double_2D a, int_1D u, double borne_inf, double_1D x_inf,
                             double borne_racine, int var_fract_racine, double val_fract_racine,
                             bool *p_opt_trouve, double *p_valopt, double_1D x_opt, double *p_temps,
                             int_1D copycstat_racine, int_1D copyrstat_racine, double_1D w_racine,
                             double *p_nb_noeuds_borne1, double t_limite,
                             bool *p_temps_limite_atteint_Djerdjour);

void pqe_bb_eliminer_contraintes(int n, int m, int *p_mprime, double_2D a, double_1D b,
                                 int *p_nb_suppr, int_1D ind_cont_suppr);

void pqe_bb_fixations_paquets_0 (int n, int m, int germe, double_2D a, double_1D b, double_1D c, double_1D d,
                       int_1D u, double_2D s, double_1D x_inf, double borne_inf, double_1D vec_dul, double *p_temps_paq0,
                       st_minilib2 *minilib2);

void pqe_bb_fixations_01(int n, int m, int germe, double_2D a, double_1D b, double_1D c, double_1D d,
		                       int_1D u, double_1D x_inf, double borne_inf, double *p_temps_fix01,
				                              st_minilib2 *minilib2);
void pqe_bb_procedure_P(int n, int m, double_2D a_init, double_1D b_init, int_1D u_init,
		       st_minilib2 *minilib2, double *p_temps);
void test_rectif(void);


/**********************/
/* Fichier : pqe_main.c */
/**********************/
void pqe_main_tri_rapide (int g, int d, elt_1D tab);
bool pqe_main_appartient(int i, int_1D I, int taille);
void pqe_main_calcule_w_etoile(int n, int m, double_1D w0, double_2D s, double_2D a,
		               double_1D b, int_1D u, double_1D w_etoile, double *p_borne1,
	                       int *p_nbiter, double *p_temps, double_2D y_cont);
void pqe_main_linea_fction_quadra (int n, int_1D u, double_1D c, double_1D d, double_2D s);
void pqe_main_reso_kp_y(int n, int_1D u, double_1D A, double b, double_2D s,
	               double_2D y, int *p_ietoile, int *p_pietoile, int_1D p);

void pqe_main_calcul_borne3(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
                            double *p_borne3, double *p_temps, double_1D vec_dul,
                            double_1D solution_y_borne3, int_1D copycstat, int_1D copyrstat,
                            int_1D getcstat, int_1D getrstat);
void pqe_main_calcul_borne2(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
                   double_1D vec_agregation, double *p_borne2, double *p_temps, double_1D solution_y);
void pqe_main_calcul_borne2_I(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
		              double_1D vec_agregation, int_1D I, int taille, double *p_borne2,
		             double *p_temps, double_1D solution_y);
bool pqe_main_est_admissible_y1D(double_1D y, int n, int m, double_1D b, double_2D a, int_1D u);
void pqe_main_primalisation(int n, int m, int_1D u, double_2D s, double_2D a, double_1D b,
                   double_1D c, double_1D d, double constante, double_1D solution_y,
                   double *p_borne_inf, bool *p_optimum_trouve, double_1D primalx,
                   double_1D primaly, double borne2, double *p_temps);
int somme_ui(int n, int_1D u);
int pqe_main_max_ui(int n, int_1D u);
void pqe_main_heur_cube(double_1D solution_y_borne3, int n, int m, double_2D a,
                        double_1D b, double_1D c, double_1D d, int_1D u, double *p_borne_inf2,
                        double_1D sol_inf2, double_1D sol_x_b3, double *p_temps_inf2);





