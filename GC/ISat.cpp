/* MIN_APPROCHEE_A_EXPANSION_GRAPH_CUT  algorithme des alpha-expansions
 *       calculées avec la construction de graphes de Kolmogorov (2004)
 *
 * Usage: Ir=min_approchee_a_expansion_graph_cut(I,labels,beta);
 
 */


#include "mex.h"
#include "math.h"

#include "graph.h"
#include "graph.cpp"
#include "maxflow.cpp"



// Terme d'attache aux données
double attache_aux_donnees(double noisy, double reg)
{
    // fonction à compléter...

    return (noisy-reg)*(noisy-reg);
}

// Terme de régularisation
double regularisation(double reg1, double reg2)
{
    // fonction à compléter...
  return abs(reg1-reg2);
}

// Message d'erreur
void error_msg(char * msg)
{
    mexWarnMsgTxt("Erreur!");
    mexWarnMsgTxt("Usage: [o_dMap, errMap, x0Map] = ISat(ILeft, IRight, i_dMap, gradI, sigma, lambda, Nlabels)");
    mexErrMsgTxt(msg);
}

// alpha expansion
bool alpha_expansion(double* ILeft, double* IRight, double* errMap, double* x0Map, double int W, int H, int iteration, float lambda)
{
    bool change=false, non_sousmodulaire=false;
    long Npix=W*H;
    
    Graph* g;              // graphe
    Graph::node_id *nodes; // tableau des noeuds du graphe
       
     /* Allocation du graphe */    
    
     g = new Graph(Npix,2*Npix,&error_msg);
     if(g==NULL)
        error_msg("Memory allocation for the graph failed!");
    
     nodes=new Graph::node_id[Npix];
     if(nodes==NULL)
        error_msg("Memory allocation failed! (array of nodes could not be created)");
        
     // Création des noeuds
     for(long p=0; p<Npix; p++)
         nodes[p]=g->add_node(); 
     
     
  /* Construction du graphe */
        
    // Création des arcs
    for(long i=0; i<W; i++)
    {
        long offset=i*H;
        for(long j=0; j<H; j++)
        {
            long pix=offset+j; // coordonnée linéaire du pixel courant
            
            // Termes d'attache aux données
            double E0=attache_aux_donnees(Id[pix],Ir[pix]);    // vraisemblance si l'on garde la même valeur
            double E1=attache_aux_donnees(Id[pix],alpha);      // vraisemblance si l'on prend la valeur alpha
            
            if(E0<E1)
                g->add_tweights(nodes[pix],E1-E0,0); // arc source -> noeud i
            else
                g->add_tweights(nodes[pix],0,E0-E1); // arc noeud i -> puits
            
            // Termes de régularisation
            if(j<H-1)
            {
                double E00=beta*regularisation(Ir[pix],Ir[pix+1]);
                double E01=beta*regularisation(Ir[pix],alpha);
                double E10=beta*regularisation(alpha,Ir[pix+1]);
                double E11=beta*regularisation(alpha,alpha);
                
                if(E10>E00)
                    g->add_tweights(nodes[pix],E10-E00,0); // arc source -> noeud i
                else
                    g->add_tweights(nodes[pix],0,E00-E10); // arc noeud i -> puits

                if(E11>E10)
                    g->add_tweights(nodes[pix+1],E11-E10,0); // arc source -> noeud i
                else
                    g->add_tweights(nodes[pix+1],0,E10-E11); // arc noeud i -> puits
            
                if(E01+E10-E00-E11>=0)
                    g->add_edge(nodes[pix],nodes[pix+1],E01+E10-E00-E11,0);
                else
                    non_sousmodulaire=true;
                    
            }
                
            if(i<W-1)
            {
                double E00=beta*regularisation(Ir[pix],Ir[pix+W]);
                double E01=beta*regularisation(Ir[pix],alpha);
                double E10=beta*regularisation(alpha,Ir[pix+W]);
                double E11=beta*regularisation(alpha,alpha);
                
                if(E10>E00)
                    g->add_tweights(nodes[pix],E10-E00,0); // arc source -> noeud i
                else
                    g->add_tweights(nodes[pix],0,E00-E10); // arc noeud i -> puits

                if(E11>E10)
                    g->add_tweights(nodes[pix+W],E11-E10,0); // arc source -> noeud i
                else
                    g->add_tweights(nodes[pix+W],0,E10-E11); // arc noeud i -> puits
            
                if(E01+E10-E00-E11>=0)
                    g->add_edge(nodes[pix],nodes[pix+W],E01+E10-E00-E11,0);
                else
                    non_sousmodulaire=true;
            }
        }
    }
    
    if(non_sousmodulaire)
        mexWarnMsgTxt("Energie non sous-modulaire, troncature réalisée");
        
  /* Calcul de la coupe minimale */
    double cout=g->maxflow();
    
    
    
  /* Déduction de l'image régularisée */
    long Npixchanged=0;
    for(long pix=0; pix<Npix; pix++)
    {
        // si le noeud n'est plus relié à la source, c'est qu'il faut associer
        // la valeur alpha au pixel correspondant dans l'image régularisée
        if(g->what_segment(nodes[pix])!=(Graph::SOURCE))
        {
            Ir[pix]=alpha;
            //mexPrintf("*");
            Npixchanged++;
            change=true;
        }
    }
    if(change)
        mexPrintf("*%d*",Npixchanged);
    
    /* Liberation de la mémoire (destruction du graphe de l'itération courante)*/
    delete g;
    delete [] nodes;
    
    return change;
}

// Algorithme des alpha-expansions
void minimize_energy(double * Id, double * Ir, long W, long H, double * label, long Nlabels, double beta)
{
    bool change=false;
    long iter=1;
    long MAXITER=5;
    
    // initialisation au maximum de vraisemblance
    for(long pix=0; pix<W*H; pix++)
    {
        double Emin=1E20;
        for(long i=0; i<Nlabels; i++)
        {
            double E=attache_aux_donnees(Id[pix],label[i]);
            if(E<Emin)
            {
                Emin=E;
                Ir[pix]=label[i];
            }
        }
    }
    
    
 // tant qu'il y a des changements, réaliser des alpha-expansions pour tous les labels alpha
    for(iter=1; iter<=MAXITER; iter++)
    {
        change=false;
        mexPrintf("Iteration %d:\n",iter);
        for(long index_alpha=0; index_alpha<Nlabels; index_alpha++)
        {
            // affichage
            if(((index_alpha+1)%10)==0)
                mexPrintf("\n");
            mexPrintf("%.1f",label[index_alpha]);
            mexPrintf("(%d)--",index_alpha+1);
            
            // calcul de l'expansion sur le niveau label[index_alpha]
            bool has_changed=alpha_expansion(Id,Ir,W,H,label[index_alpha],beta);
            change=change||has_changed;
        }
        if(!change)
            break;
        mexPrintf("\n\n");
    }
}

// weight function for errors computation
float phiX0(int i, int j){
	return float(exp(float(-(i^2+j^2))));
}

// compute error map of the input disparity map
void compute_error(double* errMap, double* i_dispMap, double* gradI, int w, int h, float sigma){
	int i, j, ii, jj;
	int m, M;
	int wi = 10, wj = 10;
	int wh = w*h;
	
	float num, den;
	float Nerr, Aerr;
	
	for(i=0; i<w; ++i){
		for(j=0; j<h; ++j){
			m = 255; M = 0;
			num = 0; den = 0;
			for(ii=-wi; ii<=wi; ++ii){
				for(jj=-wj; jj<=wj; ++jj){
					if(i+ii >= 0 && i+ii < w && j+jj >= 0 && j+jj < h && i_dispMap[i+ii +w*(j+jj)] == i_dispMap[i+ii +w*(j+jj)+2*wh]){
						m	=	(m < i_dispMap[i+ii +w*(j+jj)])? m:i_dispMap[i+ii +w*(j+jj)];
						M	=	(M > i_dispMap[i+ii +w*(j+jj)])? M:i_dispMap[i+ii +w*(j+jj)];
						den	=	den + gradI[i+ii + w*(j+jj)]*gradI[i+ii +w*(j+jj)]*phiX0(ii+1,jj+1);
						num	=	num + phiX0(ii+1,jj+1);
					}
				}
			}
			Aerr = (M-m > 0)? M-m:0;
			Nerr = sigma*sigma/(den/num+2*sigma*sigma);
			errMap[i+w*j] = double(Aerr + Nerr);
		}
	}
}

// Matlab liaison
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// check if argument number is correct
    if ( nrhs != 7 ) {
        error_msg("Wrong input number");
    } 
	else if ((nlhs != 3)) {
        error_msg("Wrong output number");
    }
    
	// compute pictures dimensions
    int W, H;
    W = int(mxGetM(prhs[0]));
    H = int(mxGetN(prhs[0])/3);
       
	// standard deviation of noise
	float sigma = float(*mxGetPr(prhs[4]));
	
	// smoothness coefficient
    float lambda = float(*mxGetPr(prhs[5]));
	
	// number of labels
	int Nlabels = int(*mxGetPr(prhs[6]));
    
	// output pictures initialization
	plhs[0] = mxCreateDoubleMatrix(W,H, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(W,H, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(W,H, mxREAL);  
	
	// transform matrices into arrays
    double *Il, *Ir, *gradI, *i_dispMap, *o_dispMap, *errMap, *x0Map;
    
	// input
    Il = mxGetPr(prhs[0]);
    Ir = mxGetPr(prhs[1]);
	i_dispMap = mxGetPr(prhs[2]);
	gradI = mxGetPr(prhs[3]);
    
	// output
	o_dispMap = mxGetPr(plhs[0]);
	errMap = mxGetPr(plhs[1]);
	x0Map = mxGetPr(plhs[2]);
    
	compute_error(errMap, x0Map, i_dispMap, gradI, W, H, sigma);
	//minimize_energy(Id,Ir,W,H,label,Nlabels,beta);
}

