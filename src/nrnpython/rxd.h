#pragma once
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include "ocmatrix.h"

#define SPECIES_ABSENT -1
#define PREFETCH       4

typedef void (*fptr)(void);

// @olupton 2022-09-16: deleted a declaration of OcPtrVector that did not match
// the one in ocptrvector.h

typedef struct {
    Reaction* reaction;
    int idx;
} ReactSet;

typedef struct {
    ReactSet* onset;
    ReactSet* offset;
} ReactGridData;


typedef struct {
    Grid_node* g;
    int onset, offset;
    double* val;
} CurrentData;


typedef struct SpeciesIndexList {
    int id;
    double atolscale;
    int* indices;
    int length;
    struct SpeciesIndexList* next;
} SpeciesIndexList;

typedef struct ICSReactions {
    ReactionRate reaction;
    int num_species;
    int num_regions;
    int num_params;
    int num_segments;
    int*** state_idx; /*[segment][species][region]*/
    int icsN;         /*total number species*regions per segment*/
    /*NOTE: icsN != num_species*num_regions as every species may not be defined
     *on every region - the missing elements of state_idx are SPECIES_ABSENT*/

    /*ECS for MultiCompartment reactions*/
    int num_ecs_species;
    int num_ecs_params;
    double*** ecs_state; /*[segment][ecs_species]*/
    int* ecs_offset_index;
    ECS_Grid_node** ecs_grid;
    int** ecs_index;
    int ecsN; /*total number of ecs species*regions per segment*/

    int num_mult;
    double** mc_multiplier;
    int* mc_flux_idx;
    double** vptrs;
    struct ICSReactions* next;

    // cache lower-upper decomposition and
    // old_state (?) to check for updates.
    // also track cache hits and misses
    std::unique_ptr<OcFullMatrix> cached_jacobian = nullptr;
    //std::unique_ptr<double*> cached_state = nullptr;

    std::unique_ptr<double[][]> = nullptr;


} ICSReactions;

struct ReactionStateCache {
    bool allocated = false;

    double **states_for_reaction = NULL;
    double **params_for_reaction = NULL;
    double *ecs_states_for_reaction = NULL;
    double *ecs_params_for_reaction = NULL;

    int num_params = 0;
    int num_species = 0;

    int num_ecs_params = 0;
    int num_ecs_species = 0;

    int num_regions = 0;

    void allocate(ICSReactions *react) {
      free_cache();

      num_params = react->num_params;
      num_species = react->num_species;
      num_ecs_params = react->num_ecs_params;
      num_ecs_params = react->num_ecs_params;
      num_regions = react->num_regions;

      // NB: can malloc(0) occur here? (Implementation defined behavior)
      states_for_reaction = (double **)malloc(num_species * sizeof(double *));
      for (int i = 0; i < num_species; i++) {
        states_for_reaction[i] = (double *)malloc(num_regions * sizeof(double));
      }

      params_for_reaction = (double **)malloc(num_params * sizeof(double *));
      for (int i = 0; i < num_params; i++) {
        params_for_reaction[i] = (double *)malloc(num_regions * sizeof(double));
      }

      if (num_ecs_species > 0) {
        ecs_states_for_reaction =
            (double *)malloc(num_ecs_species * sizeof(double));
      } else {
        free(ecs_states_for_reaction);
        ecs_states_for_reaction = NULL;
      }

      if (num_ecs_params > 0) {
        ecs_params_for_reaction =
            (double *)malloc(num_ecs_params * sizeof(double));
      } else {
        free(ecs_params_for_reaction);
        ecs_params_for_reaction = NULL;
      }

      allocated = true;
    }

    void save(double **states_for_reaction, double **params_for_reaction,
              double *ecs_states_for_reaction,
              double *ecs_params_for_reaction) {

        for (int i = 0; i < num_species; i++) {
            memcpy(this->states_for_reaction[i], states_for_reaction[i],
                   num_regions * sizeof(double));
        }
        for (int i = 0; i < num_species; i++) {
            memcpy(this->params_for_reaction[i], params_for_reaction[i],
                   num_regions * sizeof(double));
        }

        if (num_ecs_species > 0) {
            memcpy(this->ecs_states_for_reaction, ecs_states_for_reaction,
                   num_ecs_species * sizeof(double));
        }

        if (num_ecs_params > 0) {
            memcpy(this->ecs_params_for_reaction, ecs_params_for_reaction,
                   num_ecs_params * sizeof(double));
        }
    }

    bool state_changed(ICSReactions * react){        
        if(!allocated)
            return true;

        // continue on Monday
        // check if the states have changed, 

        // change to contiguous memory allocation?


    }

    void free_cache(){
        for (int i=0; i < num_species; i++){
            free(states_for_reaction[i]);
        }
        free(states_for_reaction);

        for (int i=0; i < num_params; i++){
            free(params_for_reaction[i]);
        }
        free(params_for_reaction);

        free(ecs_states_for_reaction);
        free(ecs_params_for_reaction);

        states_for_reaction = NULL;
        params_for_reaction = NULL;
        ecs_states_for_reaction = NULL;
        ecs_params_for_reaction = NULL;

        allocated = false;

    }

    ~ReactionStateCache(){
        free_cache();
    }
};

typedef struct TaskList {
    void* (*task)(void*);
    void* args;
    void* result;
    struct TaskList* next;
} TaskList;

typedef struct TaskQueue {
    std::condition_variable task_cond, waiting_cond;
    std::mutex task_mutex, waiting_mutex;
    std::vector<bool> exit;
    int length{};
    struct TaskList* first;
    struct TaskList* last;
} TaskQueue;

extern "C" void set_num_threads(const int);
void _fadvance(void);
void _fadvance_fixed_step_3D(void);

extern "C" int get_num_threads(void);
void ecs_set_adi_tort(ECS_Grid_node*);
void ecs_set_adi_vol(ECS_Grid_node*);
void ecs_set_adi_homogeneous(ECS_Grid_node*);

void dg_transfer_data(AdiLineData* const, double* const, int const, int const, int const);
void ecs_run_threaded_dg_adi(const int, const int, ECS_Grid_node*, ECSAdiDirection*, const int);
ReactGridData* create_threaded_reactions(const int);
void* do_reactions(void*);

void current_reaction(double* states);

void run_threaded_deltas(ICS_Grid_node* g, ICSAdiDirection* ics_adi_dir);
void run_threaded_ics_dg_adi(ICS_Grid_node* g, ICSAdiDirection* ics_adi_dir);
void ics_dg_adi_x(ICS_Grid_node* g,
                  int,
                  int,
                  int,
                  double,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*);
void ics_dg_adi_y(ICS_Grid_node* g,
                  int,
                  int,
                  int,
                  double,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*);
void ics_dg_adi_z(ICS_Grid_node* g,
                  int,
                  int,
                  int,
                  double,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*,
                  double*);
void ics_dg_adi_x_inhom(ICS_Grid_node* g,
                        int,
                        int,
                        int,
                        double,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*);
void ics_dg_adi_y_inhom(ICS_Grid_node* g,
                        int,
                        int,
                        int,
                        double,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*);
void ics_dg_adi_z_inhom(ICS_Grid_node* g,
                        int,
                        int,
                        int,
                        double,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*,
                        double*);

/*Variable step function declarations*/
void _rhs_variable_step(const double*, double*);

void _ode_reinit(double*);

int ode_count(const int);

extern "C" void scatter_concentrations(void);


int find(const int, const int, const int, const int, const int);

void _ics_hybrid_helper(ICS_Grid_node*);
void _ics_variable_hybrid_helper(ICS_Grid_node*,
                                 const double*,
                                 double* const,
                                 const double*,
                                 double* const);

void _ics_rhs_variable_step_helper(ICS_Grid_node*, double const* const, double*);
void _rhs_variable_step_helper(Grid_node*, double const* const, double*);

void ics_ode_solve(double, double*, const double*);
void ics_ode_solve_helper(ICS_Grid_node*, double, double*);

void _rhs_variable_step_helper_tort(Grid_node*, double const* const, double*);

void _rhs_variable_step_helper_vol(Grid_node*, double const* const, double*);

void set_num_threads_3D(int n);

void _rhs_variable_step_ecs(const double*, double*);

void clear_rates_ecs();
void do_ics_reactions(double*, double*, double*, double*);
void get_all_reaction_rates(double*, double*, double*);
void _ecs_ode_reinit(double*);
void do_currents(Grid_node*, double*, double, int);
void TaskQueue_add_task(TaskQueue*, void* (*task)(void* args), void*, void*);
void TaskQueue_exe_tasks(std::size_t, TaskQueue*);
void TaskQueue_sync(TaskQueue*);
void ecs_atolscale(double*);
void apply_node_flux3D(Grid_node*, double, double*);
