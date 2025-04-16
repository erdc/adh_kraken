#ifndef H_SPARAM_
#define H_SPARAM_

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++/
typedef struct {

	//simple structure
	//contains the map
	//from adh_def.param
	//to the supermodel params
	//and reverse map
	//as well as stores the prams themselves

	int n_int; //number of active params that are integers
	int n_dbl;
	int *param_int_map; //array of adh_def.n_param that maps adh_def.param stuff to active param #
	int *param_int_reverse_map; //array of n_int that maps the active param # to place in adh_def.param
	int *param_dbl_map; //array of adh_def.n_param that maps adh_def.param stuff to active param #
	int *param_dbl_reverse_map; //array of n_db that maps the active param # to place in adh_def.param

	//stores for data
	// params could be double or int
	int *int_data; //array of n_int containing all integer parameters
	double *dbl_data; //array of n_dbl contatining all double parameters



} SPARAM;



#endif
