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
extern unsigned int alpha;
/*-------------------------------------------------------------------*/
//given a vertex to go round facets
struct vertexNode
{
	int cnt_facet;
	bool isOnqhull;
	int* nbd_facet;
	//also means direction
	double n[3];
	int id_pos;
	int id_nega;
	//double direction[3];
};
typedef struct vertexNode vertexNode;

//the 4th vertex is added to check whether the normal direction is right
void getNormal(tVertex v1, tVertex v2, tVertex v3, tVertex v4, double n[3])
{
	double v12[3], v13[3], v14[3];
	v12[0] = v1->v[0] - v2->v[0]; v13[0] = v1->v[0] - v3->v[0];
	v12[1] = v1->v[1] - v2->v[1]; v13[1] = v1->v[1] - v3->v[1];
	v12[2] = v1->v[2] - v2->v[2]; v13[2] = v1->v[2] - v3->v[2];
	n[0] = v12[1] * v13[2] - v12[2] * v13[1]; v14[0] = v1->v[0] - v4->v[0];
	n[1] = v12[2] * v13[0] - v12[0] * v13[2]; v14[1] = v1->v[1] - v4->v[1];
	n[2] = v12[0] * v13[1] - v12[1] * v13[0]; v14[2] = v1->v[2] - v4->v[2];

	double isOuterN = n[0] * v14[0] + n[1] * v14[1] + n[2] * v14[2];
	double normalize = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	if (normalize == 0){ n[0] = n[0]; n[1] = n[1]; n[2] = n[2]; }
	else{
		n[0] = n[0] / normalize; n[1] = n[1] / normalize; n[2] = n[2] / normalize;
		if (isOuterN < 0) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }
	}

}


void Crust(void)
{
	tVertex  ptr_v;
	tVertex * all_v = NULL;
	tVertex * selec_v = NULL;
	int vsize = 0;
	int id = 0;
	int cnt_pole = 0;

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

	//it returns to voroni center which means the center 
	//of the circumshpere of each facet.(in our case facet means a tetrahedron)
	qh_setvoronoi_all();


	/* **************************************ADDED************************************************ */
	vertexNode* vnodes = (vertexNode*)malloc(vsize*sizeof(vertexNode));
	facetT** allfacets = (facetT**)malloc(qh num_facets*sizeof(facetT*));
	//memory allocate to isPolefacets and at the same time initialize it as false.
	bool* isPolefacets = (bool*)calloc(qh num_facets, sizeof(bool));


	for (int i = 0; i < vsize; i++)
	{
		vnodes[i].cnt_facet = 0; vnodes[i].isOnqhull = false; 
		vnodes[i].n[0] = 0; vnodes[i].n[1] = 0; vnodes[i].n[2] = 0;
		//If positive(negative) pole does not exist.
		vnodes[i].id_pos = -1;
		//If we exclude the centers which are far from the data points
		//the facets will be excluded and therefore there does not exists a pole
		//in that case, not to include the poles, we set the index as negative. (since index starts from 0)
		vnodes[i].id_nega = -1;
	}

	int facet_indx = 0;
	FORALLfacets
	{
		if (facet->normal[3] < 0)
		{
			//get the number of neighboring facet for each vertex.
			FOREACHvertex_(facet->vertices)
			{
				int id = qh_pointid(vertex->point);
				vnodes[id].cnt_facet++;
			}
		}
		//store all down-side facets obtained from Qhull().
		allfacets[facet_indx++] = facet;
	}

		//malloc to neighbor facet for each vertex. 
	for (int i = 0; i < vsize; i++){
		vnodes[i].nbd_facet = (int*)malloc(vnodes[i].cnt_facet*sizeof(int));
		vnodes[i].cnt_facet = 0;
	}

	int id_facet = 0;
	FORALLfacets
	{
		if (facet->normal[3] < 0)
		{
			//give index for neighboring facets, loop through all the vertices.
			FOREACHvertex_(facet->vertices)
			{
				int id = qh_pointid(vertex->point);
				vnodes[id].nbd_facet[vnodes[id].cnt_facet++] = id_facet;//store neighbor facets for each vertex.
			}
		}
		id_facet++;
	}

	//loop through all faces
	FORALLfacets
	{
		if (facet->normal[3] < 0)
		{
			vid = 0;
			int facet_down[4];
			FOREACHvertex_(facet->vertices)
			{
				//get the id of the vertex
				facet_down[vid++] = qh_pointid(vertex->point);
			}

			FOREACHneighbor_(facet)
			{
				//if (neighbor->normal[3] > 0)
				{
					int vertices[4];

					vid = 0;
					FOREACHvertex_(neighbor->vertices)
					{
						vertices[vid++] = qh_pointid(vertex->point);
					}

					//compute normal vector 
					int vid = 0;
					int a[3], a4;
					for (int i = 0; i < 4; i++)
					{
						bool check = false;
						for (int j = 0; j < 4; j++)
						{
							if (facet_down[i] == vertices[j])
							{
								vnodes[vertices[j]].isOnqhull = true;
								a[vid++] = vertices[j];
								check = true;
							}
						}
						if (!check){ a4 = facet_down[i]; }
					}
					double N[3];
					getNormal(all_v[a[0]], all_v[a[1]], all_v[a[2]], all_v[a4], N);
					vnodes[a[0]].n[0] += N[0]; vnodes[a[0]].n[1] += N[1]; vnodes[a[0]].n[2] += N[2];
					vnodes[a[1]].n[0] += N[0]; vnodes[a[1]].n[1] += N[1]; vnodes[a[1]].n[2] += N[2];
					vnodes[a[2]].n[0] += N[0]; vnodes[a[2]].n[1] += N[1]; vnodes[a[2]].n[2] += N[2];

				}//FOREACHneighbor_
			}
		}
	}

	for (int i = 0; i < vsize; i++)
	{
		//if not on the convex hull
		if (!vnodes[i].isOnqhull)
		{
			double max = 0;
			for (int j = 0; j < vnodes[i].cnt_facet; j++)
			{
				//전체 facet들 중에서의 index
				int id_facet = vnodes[i].nbd_facet[j];
				double temp = sqrt((all_v[i]->v[0] - allfacets[id_facet]->center[0])*(all_v[i]->v[0] - allfacets[id_facet]->center[0]) +
					(all_v[i]->v[1] - allfacets[id_facet]->center[1])*(all_v[i]->v[1] - allfacets[id_facet]->center[1]) +
					(all_v[i]->v[2] - allfacets[id_facet]->center[2])*(all_v[i]->v[2] - allfacets[id_facet]->center[2]));
				//get the positive pole and store the index of that facet.
				if (max <= temp)
				{
					max = temp;
					vnodes[i].id_pos = id_facet;
				}
			}
			if (vnodes[i].id_pos >= 0)
			{
				vnodes[i].n[0] = allfacets[vnodes[i].id_pos]->center[0] - all_v[i]->v[0];
				vnodes[i].n[1] = allfacets[vnodes[i].id_pos]->center[1] - all_v[i]->v[1];
				vnodes[i].n[2] = allfacets[vnodes[i].id_pos]->center[2] - all_v[i]->v[2];
			}
		}
		double min = 10;
		for (int j = 0; j < vnodes[i].cnt_facet; j++)
		{
			//전체 facet들 중에서의 index
			int id_facet = vnodes[i].nbd_facet[j];
			double temp = sqrt(vnodes[i].n[0] * (all_v[i]->v[0] - allfacets[id_facet]->center[0]) +
				vnodes[i].n[1] * (all_v[i]->v[1] - allfacets[id_facet]->center[1]) +
				vnodes[i].n[2] * (all_v[i]->v[2] - allfacets[id_facet]->center[2]));
			if (min >= temp)
			{
				min = temp;
				vnodes[i].id_nega = id_facet;
			}
		}
		//bb경우 facet의 positive pole이 없는 경우가 발생
		for (int i = 0; i < vsize; i++)
		{
			if (vnodes[i].id_pos >= 0)
			{
				isPolefacets[vnodes[i].id_pos] = true;
			}
			if (vnodes[i].id_nega >= 0)
			{
				isPolefacets[vnodes[i].id_nega] = true;
			}
		}
	}

	for (int i = 0; i < qh num_facets; i++)
	{
		if (isPolefacets[i])
		{
			cnt_pole++;
		}
	}

	pt = (coordT*)realloc(pt, (vsize + cnt_pole) * 4 * sizeof(coordT)); //each point will have three coord
	all_v = realloc(all_v, (vsize + cnt_pole)*sizeof(tVertex));

	for (int i = 0; i < qh num_facets; i++)
	{
		//if the facet has a pole
		//add the poles to the point set. To construct Delauany Triangulation.
		if (isPolefacets[i])
		{
			pt[id++] = allfacets[i]->center[0];
			pt[id++] = allfacets[i]->center[1];
			pt[id++] = allfacets[i]->center[2];
			pt[id++] = (allfacets[i]->center[0]) * (allfacets[i]->center[0]) + (allfacets[i]->center[1]) * (allfacets[i]->center[1]) + (allfacets[i]->center[2]) * (allfacets[i]->center[2]);
		}
	}

	/*int vnum = vsize;
	for (int i = 0; i < qh num_facets; i++)
	{
	if (isPolefacets[i])
	{
	tVertex vertex = MakeNullVertex();
	for (int j = 0; j < 3; j++)
	vertex->v[j] = allfacets[i]->center[j];
	all_v[vnum] = vertex;
	vertex->vnum = vnum++;
	}
	}*/

	qh_freeqhull(!qh_ALL);


	//using qhull

	qh_init_A(stdin, stdout, stderr, 0, NULL);

	//qh DELAUNAY= True;     /* 'd'   */
	//qh SCALElast= True;    /* 'Qbb' */
	//qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

	qh_initflags(options);
	qh_init_B(pt, vsize + cnt_pole, 4, false);
	qh_qhull();
	qh_check_output();
	qh_setvoronoi_all();


	//디버깅
	/*int vnum = vsize;
	for (int v = 0; v < vsize; v++)
	{
	if (vnodes[v].isOnqhull)
	{
	printf("%d\n", v);
	tVertex vertex = MakeNullVertex();
	vertex->v[0] = all_v[v]->v[0] + vnodes[v].n[0];
	vertex->v[1] = all_v[v]->v[1] + vnodes[v].n[1];
	vertex->v[2] = all_v[v]->v[2] + vnodes[v].n[2];
	printf("%f %f %f\n", vertex->v[0], vertex->v[1], vertex->v[2]);
	vertex->vnum = vnum++;
	tEdge edge = MakeNullEdge();
	edge->endpts[0] = all_v[v];
	edge->endpts[1] = vertex;
	}
	}*/

	//loop through all faces
	FORALLfacets
	{
		double radius = 0;
		for (int i = 0; i < 3; i++)//3차원으로 계산 왜그런지는 .. 
		{
			//e is a set element which means one of four vertices of facet.
			radius += SQR(facet->center[i] - (*(vertexT**)&facet->vertices->e[3])->point[i]);
		}

		if (facet->normal[3] < 0 && radius < alpha)
		{
			//tTetra tetra = MakeNullTetra(); //make a face in 4D

			//get vertices of facet
			//loop through each vertex
			vid = 0;
			int temp[4]; int j = 0;
			for (int i = 0; i < 4; i++){ temp[i] = -1; }
			int cnt = 0;
			FOREACHvertex_(facet->vertices)
			{
				//if it is a point in a given point set
				if (qh_pointid(vertex->point) < vsize)
				{
					temp[j++] = qh_pointid(vertex->point);
					cnt++;
				}
			}
			//if three vertices are on the same plane
			if (cnt >= 3)
			{
				FOREACHneighbor_(facet)
				{
					double nbr_radius = 0;
					for (int i = 0; i < 3; i++)//3차원으로 계산 왜그런지는 .. 
					{
						//e is a set element which means one of four vertices of facet.
						nbr_radius += SQR(neighbor->center[i] - (*(vertexT**)&neighbor->vertices->e[3])->point[i]);
					}

					if (nbr_radius<alpha)
					{
						tVertex vertices[3];

						vid = 0;
						FOREACHvertex_(neighbor->vertices)
						{
							for (int i = 0; i < 4; i++)
							{
								if (temp[i] == qh_pointid(vertex->point))
								{
									vertices[vid++] = all_v[temp[i]];
								}
							}
						}
						if (vid >= 3)
						{
							tFace face = MakeNullFace();
							for (int k = 0; k < 3; k++)
							{
								face->vertex[k] = vertices[k];
							}
						}
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
