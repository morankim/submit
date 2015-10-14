#include "crust.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"

/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void	Delaunay(void)
{


	tVertex  ptr_v;
	tVertex * all_v = NULL;
	tVertex * selec_v = NULL;
	int vsize = 0;
	int id = 0;

	//global varibles for qhull
	static char * options = (char *)"qhull QJ Pp";
	int curlong, totlong;
	coordT * pt = NULL;
	facetT *facet = NULL;
	facetT *facett = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	facetT *neighbor, **neighborp;
	int vid = 0;
	tTetra tetra;

	//count number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);


	//allocate memory
	pt = (coordT*)calloc(vsize * 4, sizeof(coordT)); //each point will have three coord
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	//copy points
	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		pt[id++] = (ptr_v->v[0])*(ptr_v->v[0]) + (ptr_v->v[1])*(ptr_v->v[1]) + (ptr_v->v[2])*(ptr_v->v[2]);
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//using qhull

	qh_init_A(stdin, stdout, stderr, 0, NULL);

	//qh DELAUNAY= True;     /* 'd'   */
	//qh SCALElast= True;    /* 'Qbb' */
	//qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();


	//loop through all faces
	FORALLfacets
	{
	
		if (facet->normal[3] < 0)
		{

			tetra = MakeNullTetra(); //make a face in 4D

			//get vertices of facet
			//loop through each vertex
			vid = 0;

			FOREACHvertex_(facet->vertices)
			{
				//get the id of the vertex
				tetra->vertex[vid++] = all_v[qh_pointid(vertex->point)];
			}

			FOREACHneighbor_(facet)
			{
				if (neighbor->normal[3]>0)
				{
					tVertex vertices[4];
					vid = 0;
					FOREACHvertex_(neighbor->vertices)
					{
						//get vertex
						vertices[vid++] = all_v[qh_pointid(vertex->point)];
					}

					tFace face = MakeNullFace();

					int k=0;
					for (int i = 0; i < 4; i++)
						for (int j = 0; j < 4; j++)
							if (tetra->vertex[i] == vertices[j])
							{
								face->vertex[k++] = vertices[j]; break;
							}
			
				}
			}//FOREACHneighbor_
		}
	}

		//not used
	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;

	//free mem
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}




