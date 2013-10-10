/* ---------------------- */
/* Fichier : knap_const.h */
/* ---------------------- */
#if !defined(KNAP_CONST)

#define KNAP_CONST 0


/* ----------- */
/* Generalites */
/* ----------- */

/* #define     SYSTEME_ESRA       1     mettre cette ligne en commentaire si on n'est pas sur esra */
#define     SYSTEME_SIAMOIS    1    /* mettre cette ligne en commentaire si on n'est pas sur siamois */
#define     BENCHMARK          1    /* pour l'écriture du fichier correspondant au jeu */
#define     JUSTE_FIX          0    /* calcul du nb de variables fixees uniquement  sans B&B */ 
#define     JUSTE_LP           0    /* on ne fait qu'ecrire le fichier LP sans le resoudre */
#define     MAUVAISE_SOL_INIT  0    /* Test influence d'une mauvaise solution initiale */
#define     AFFICHAGE_PG       1    /* affichage du programme */
#define     ESRA               1    /* indique si on doit creer un fichier lp pour esra */
#define     STOCHA             0    /* indique si l'on resout des instances stochastiques */
#define     TPS_PRE_PROCESSING 0    /* juste pour calculer uniquement le temps de preprocessing */
#define     HELMBERG          0
#define     JEUX_COMPILATEUR  0     /* indique si l'on doit lire les instance */
                                    /* de Helmberg sur la conception de       */
									/* compilateurs ou bien generer           */
                                    /* aleatoirement les donnees
									 * */
#define     INFINI            10e20 /* avant 10e30 */
#define     INFINI_int        100000 /* Sur les HP, les int vont jusqu'a 2^31 */
#define     LONGCHAINE        50
#define     PETIT_PB          5
#define	    NB_PB	      20     /* Nb de pbms pour jeux d'essai */
#define     TAILLE_PAQUETS    20

#define     DEMI_LIN          0    /* Linearisation tenant compte des termes qij avec i>j - Formulation legerement differente */
#define     LIN1            1    /* Formulation AB (Glover amelioree) */
#define     LIN_CR            2    /* Formulation ES Linearisation croisee */
#define     LIN2            3    /* Formulation Min - Max                */
#define     LIN_CL            4    /* Linearisation classique              */
#define     LIN_PR            5    /* Linearisation produit                */
#define     DECOMPO           6
#define     LIN3              7

#define nb_col_max       39 /* Nombre max de colonnes du tableau des resultats */
#define nb_col_donnees    5 /* Nb de colonnes pour la presentation des donnees du pb */
#define nb_col_communes   7 /* Nb de colonnes communes a tous les types de resolution */
#define nb_col_decompo    2 /* Nb de colonnes propres a la decomposition lagrangienne */


/* ------------------------------------------ */
/* Constantes concernant les criteres d'arret */
/* ------------------------------------------ */

#define     FIN_1       1       /* Sous-gradient nul                            */
#define     FIN_2       2       /* Primal = Dual ou Primal > Dual pour fixation */
#define     FIN_3       3       /* Convergence specifique : sur-place           */
#define     FIN_4       4       /* Iteration maximum                            */
#define     FIN_5       5       /* Pas quasiment nul                            */


/* ------------------- */
/* Constantes du B & B */
/* ------------------- */

#define CHOIX 1   /* 0 si branchement sur valeur supposee, */
				  /* 1 si branchement sur l'opposee de la valeur supposee */
#define VITESSE_NORMALE     0 /* Fixe la vitesse du calul de la borne                                 */
#define VITESSE_MOYENNE     1
#define GRANDE_VITESSE1     2 /* Plus la vitesse est grande, moins on passe de temps sur chaque noeud */
#define GRANDE_VITESSE2     3 /* Plus la vitesse est grande, moins on passe de temps sur chaque noeud */
#define TRES_GRANDE_VITESSE 4 /* mais plus le nombre de noeud est eleve !                             */


/* --------------------------------------------------- */
/* Macros pour l'affichage a l'ecran, le min et le max */
/* --------------------------------------------------- */

#define    print         if (affichage) printf
#define    min(A, B)     ((A) < (B) ? (A) : (B))   /* Attention : A et B sont evalues deux fois : pas de ++ */
#define    max(A, B)     ((A) > (B) ? (A) : (B))

/* -------------------------------- */
/* Taille des fonctions economiques */
/* -------------------------------- */

#define 	COUTMAX		100
#define		COEFCAPMAX	50

/* ------------------------------------------------------------------- */
/* Precision fixant l'arret des iterations pour fixation des variables */
/* ------------------------------------------------------------------- */

#define		PRECI_STOP	     0.05 
#define		MARGE_ARRONDIS   0.5

/* ------------- */
/* Miscellaneous */
/* ------------- */

#define PAS_MAX          2.0
#define PAS_MIN          0.5
#define TAILLE_CLUSTER   5           /* Taille des clusters de variables non encore fixees */
#define TMAX             5000
#define NBITER           1000
#define ECARTMIN         0.00000009 
#define REFROIDI         0.99   /* Anciennement .99 */
#define NBPALIERS        500    
#define PRECIS           0.001
#define MIN_ENUM         14    /* Lorsqu'il reste MIN_ENUM var. libres dans le B&B, on enumere au lieu d'evaluaer */
                               /* Anciennement 14 */ 
#define RECONSTRUCTION   20    /* Nb de variables a fixer avant reconstruction des clusters */
#define NB_ALPHA         10    /* Nb de solutions generees aleatoirement pour le calcul de alpha */
#define SEUIL_VITESSE    0.6   /* Au dela de  SEUIL_VITESSE pourcents de var. fixees, on va vite */

/* --------------------------- */
/* Reglage de la decomposition */
/* --------------------------- */

#define     PRECISION         1.    /* Anciennement 1e-3                       */
#define     PRECIS_REL        1e-6  /* Precision relative                      */
#define     PRECIS_ABS        1.    /* Anciennement 1e-2 Precision absolue     */
#define     PRECIS_ABSCLIQUE  0.5   /* Anciennement 1e-2 Precision absolue     */
#define     PAS_QUASI_NUL     1e-4  /* Anciennement 1e-7                       */
#define     MAXSP             15    /* Anciennement 15 */
#define     MAXSPCLIQUE       15
#define     MAXITER_0         1000  /* ou 50 ou 125    -> correspond a VITESSE_NORMALE      */
#define     MAXITER_1         200   /* ou 50 ou 125    -> correspond a VITESSE_MOYENNE      */
#define     MAXITER_2         20    /* -> correspond a GRANDE_VITESSE1                      */
#define     MAXITER_3         10    /* -> correspond a GRANDE_VITESSE2                      */
#define     MAXITER_4         10    /* 75 ou 50 ou 125  -> correspond a TRES_GRANDE_VITESSE */
#define     MINITER           10    /* Anciennement 200 puis 10 : Nbr min d'iteration  */
#define     MAXCONV            4    /* Nbre d'iteration avant reduction du pas  Anciennement 5 */
#define     MAXCONVCLIQUE      3    /* Nbre d'iteration avant reduction du pas */
                                    /* Anciennement 4 puis 5                   */
/* ---------------- */
/* Linearisation AB */
/* ---------------- */

#define     NB_JEUX_SIP     20     /* Nb de jeux d'essai pour Stochastic Integer Programming */
#define     POLITIQUE_BB_0  FAUX   /* Test de faisabilite en 0-1 par B&B dans sous_phii (linearisation Min-Max) */
#define     N_LIMITE  2000000000   /* Limite du nb de noeuds des B&B */ 
#define     N_LIM_sous_phiI 125000  /* Limite du nb de noeuds ds sous_phii essai avec 4, 20, 200, 2000, 125000 */
#define     T_LIMITE    10800.0   /* temps limite en secondes pour resolution de la linearisation */
/*#define     T_LIMITE    60.0  */ /* temps limite en secondes pour resolution de la linearisation */

#define     BRAN_NAT          0    /* Pas de priorite dans le choix des variables de branchement */
#define     BRAN_PRI_1        1    /* Lors du B&B, les variables prioritaires pour le branchement sont celles a 1 dans la sol adm */
#define     BRAN_PRI_0        2
#define     BRAN_RAT          3    /* Les variables sont considerees dans l'ordre du ratio pour le choix de branchement */

#define     BB_PROF           0    /* B&B profondeur d'abord */
#define     BB_BEST_BOUND     1    /* B&B largeur d'abord - le noeud choisi est celui de cout minimum */
#define     BB_BEST_EST       2    /* Largeur d'abord mais le choix du noeud est fait par une boite noire de Cplex */

#define     VAR_SEL_DEFAUT    2    /* Quand Cplex choisit la variable de branchement, il les prend ds l'ordre naturel */
#define     VAR_SEL_POUSSEE   3    /* Cplex fait un calcul pour determiner la variable de branchement */

#define     VAR_NAT           0    /* Au depart, les variables ne sont pas renumerotees dans un ordre particulier */
#define     VAR_PDS_CR        1    /* Au depart, les variables sont renumerotees dans l'ordre de leur poids croissant */
#define     VAR_PDS_DEC       2    /* id mais dans l'ordre de leur poids decroissant */
#define     VAR_RAT_CR        3    /* id mais dans l'ordre du ratio d'interet croissant */
#define     VAR_RAT_DEC       4    /* id mais dans l'ordre du ratio d'interet decroissant */

#endif
