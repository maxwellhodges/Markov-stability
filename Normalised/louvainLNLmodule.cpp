
#if defined(macintosh) || defined(MACOS) || defined(_MACOS)  || defined(__APPLE__)  || defined(_MAC)   || defined(MAC)  || defined(mac)  || defined(MACINTOSH)
#define __MAC__
#endif

#if defined(_WIN32)  || defined(WIN32)  || defined(__WIN32__) || defined(_WIN64)  || defined(WIN64) || defined(__WIN64__) || defined(_WINDOWS) || defined(__WINDOWS__)
#define __WIN__
#endif

#if defined(linux)  || defined(__linux__)  || defined(_linux) || defined(LINUX)  || defined(_LINUX) || defined(_UNIX) || defined(__UNIX__) || defined(__gnu_linux__) || defined(__unix__) || defined(UNIX) || defined(unix) || defined(sparc)
#define __lin__
#endif

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>

#ifdef __WIN__
#include <time.h>
#include <sys/timeb.h> 
#include <process.h>
#include <windows.h>
#endif

#include "graph_binary.h"
#include "community.h"

#ifdef __WIN__
int gettimeofday (struct timeval *tp, void *tz)
{
	struct _timeb timebuffer;
	_ftime (&timebuffer);
	tp->tv_sec = timebuffer.time;
	tp->tv_usec = timebuffer.millitm * 1000;
	return 0;
}
#endif

#ifdef __lin__
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h> 
#endif

#ifdef __MAC__
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h> 
#endif

/*start of python wrapper */
extern "C" {
#include <Python.h>
#include "numpy/arrayobject.h"
#include <math.h>
}

using namespace std;

double *data = NULL;
int display_level = -1;
int length_data = -1;
int type = WEIGHTED;
int nb_pass = 0;
bool hierarchy = false; //can probably delete this


void display_time(const char *str) {
	time_t rawtime;
	time(&rawtime);
	cerr << str << " : " << ctime(&rawtime);
}

extern "C" {
static PyObject* stability(PyObject* self, PyObject* args) {

	 
	 //input variables
    PyArrayObject *input=NULL;
    double timet = 1.0;
    double precision = 0.000001;
    char weighted;


    //output variables
    PyArrayObject *community_id=NULL;
    double stability;
    int number_of_comms_out;

    if (!PyArg_ParseTuple(args, "Odd" , &input, &timet, &precision))
        return NULL;

    if (NULL == input)  return NULL;

    // double *data = PyArray_DATA(input);
	double *data;
	data = (double*)(input -> data);

    int length_data = PyArray_DIM(input, 1);  //number of edges - length of edge list


	int nd_columns = PyArray_DIM(input, 1);  //number of edges - length of edge list



	
	 /*This is the Louvain code */
	
	struct timeval tv;
	gettimeofday(&tv, NULL);

	// Initialize the random generator
	srand(((tv.tv_sec * 1000) + (tv.tv_usec / 1000)) * getpid());


		// Vector containing the hierarchy of the partitioning
		vector<vector<int> > output;

		// Creation of the Community object
		Community *c = new Community(data, length_data, -1, precision, timet,
				type);

		// Initial number of nodes
		int numberinitial = (*c).g.nb_nodes;

		// Calculation of the initial modularity
		double mod = (*c).modularity();

		// First partitioning
		bool improvement = (*c).one_level();

		// Calculation of its modularity
		double new_mod = (*c).modularity();

		output = (*c).display_partition2(output);

		Graph g = (*c).partition2graph_binary();

		int numberofcommunities = (*c).g.nb_nodes;
		int level = 0;
		while (improvement) {

			mod = new_mod;

			delete c;
			c = new Community(g, -1, precision, timet);

			numberofcommunities = (*c).g.nb_nodes;

			improvement = (*c).one_level();
			if ((*c).nb_pass != -10)
				new_mod = (*c).modularity();
			else
				new_mod = 0;

			output = (*c).display_partition2(output);

			g = (*c).partition2graph_binary();
			level++;

		}

		numberofcommunities = (*c).g.nb_nodes;

		free(g.weights);
		free(g.links);
		free(g.degrees);



		/* end */



	int no_of_nodes = (int)data[nd_columns-1] + 1;
	int dims[0];
    dims[0] = no_of_nodes;

    community_id = (PyArrayObject *) PyArray_FromDims(1, dims, NPY_DOUBLE);
    // double *pointer_to_data = PyArray_DATA(community_id); //pointer to actual array
	double *pointer_to_data;
	pointer_to_data = (double*)community_id -> data;

    number_of_comms_out = numberofcommunities; //probably don't even need
                                               //to copy this over...
    stability = new_mod; //again can probably output directly

    vector<int> n2c(numberinitial);

	for (unsigned int i = 0; i < numberinitial; i++)
				n2c[i] = i;

    for (int l = 0; l < output.size(); l++) {
				for (unsigned int node = 0; node < numberinitial; node++) {
					n2c[node] = output[l][n2c[node]];
                }
    }
	
	for (unsigned int node = 0; node < numberinitial; node++) {
				pointer_to_data[node] = (double)n2c[node];
    }


		

	delete c;


	return Py_BuildValue("diO", stability, number_of_comms_out, community_id);
	// return Py_BuildValue("d", out2);

}

static PyMethodDef methods[] = {
    {"stability", stability, METH_VARARGS, "This is the docstring for the function - fill this in"},
    {NULL, NULL, 0, NULL} //
};



PyMODINIT_FUNC PyInit_louvainLNL (void)
{
    PyObject *module;
	static struct PyModuleDef louvainLNLmoduledef = {
		PyModuleDef_HEAD_INIT,
		"louvainLNL", /*name of module */
		NULL, /*module documentation, need to fill this in*/
		-1,
		methods
};

	module =  PyModule_Create(&louvainLNLmoduledef);
    import_array();

	return module;
}

}