#ifndef GRIDCONSTANTS_H
#define GRIDCONSTANTS_H

#include <math.h>


// #define TESTS  // Uncomment to prevent all testing

#if defined(TESTS)
    // #define TEST_SW1            // Advection of a cosine bell curve
    // #define TEST_SW2               // Steady-state solution to nonlinear eqs
    #define TEST_GAUSS_HILLS
    // #define TEST_SW5
    // #define TEST_OPERATORS
#endif

constexpr int GRID_LVL = 8;

constexpr int get_node_num(int glvl)
{
  // geodesic node num expression (Lee and Macdonald, 2008)
  //   return 10 * pow(pow(2, glvl - 1), 2) + 2;

  int i=0;
  int x = 2;
  for (i=1; i<glvl-1; i++) {
    x *= 2; 
  }

  return 10 * x * x + 2;
  
}

constexpr int get_face_num(int node_num)
{
    return ((node_num-12)*6 + 12*5)/2;
}

constexpr int get_vertex_num(int node_num)
{
    return ((node_num-12)*6 + 12*5)/3;
}

const int NODE_NUM = get_node_num(GRID_LVL);
const int FACE_NUM = get_face_num(NODE_NUM);
const int VERTEX_NUM = get_vertex_num(NODE_NUM);

#endif

