/*
Copyright (c) 2023 Science and Technology Facilities Council

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef PLAN_H_COMPILE_GUARD
#define PLAN_H_COMPILE_GUARD

#include "define.h"
#include <map>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <vector>
#include <list>
#include <algorithm>
#include <vector>
#include <cmath>
#include <tuple>

using std::unique_ptr;

using std::cout;
using std::endl;
using std::flush;
using std::hash;
using std::unordered_map;





struct Op {
   int t, la, lb, lc, ld, m;
};

struct OpHash{
   std::size_t operator()(const Op& sh) const {
      return hash<int>()( ((((sh.t*32+sh.la)*32+sh.lb)*32+sh.lc)*32+sh.ld)*32+sh.m );
   }
};

struct OpEqual {
   bool operator()(const Op& lhs, const Op& rhs) const {
      return (lhs.t==rhs.t) and (lhs.la==rhs.la) and (lhs.lb==rhs.lb) and (lhs.lc==rhs.lc) and (lhs.ld==rhs.ld) and (lhs.m==rhs.m);
   }
};

typedef unordered_map<Op,int,OpHash,OpEqual> OpDict;
typedef std::array<int,OP_SIZE> OpArray;

bool ArrayOpLess( const OpArray& a, const std::array<int,OP_SIZE>& b );

inline bool ArrayOpGreater( const OpArray& a, const std::array<int,OP_SIZE>& b );

bool contains(const OpDict* G, Op op );
bool contains(const OpDict& G, Op op );

void rec_add( int h, int la, int lb, int lc, int ld, int m, OpDict *G );
int mem( const int la, const int lb, const int lc, const int ld);
int mem( const Op& sh );


std::vector<int> plan( int la, int lb, int lc, int ld, int* pr_size, int* sh_size, int* numV, int* numVC, int* numVCH );

typedef std::tuple< std::vector<int>, int, int, int, int, int > Plan_t;

bool contains(const std::unordered_map< int, Plan_t >& PC , const int la, const int lb, const int lc, const int ld );

class PlanCollection{
public:

   std::unordered_map< int, Plan_t > plans;

   void get(
          const int la, const int lb, const int lc, const int ld, 
          std::vector<int>** plan, unsigned int* pr_size, unsigned int* sh_size, unsigned int* numV, unsigned int* numVC, unsigned int* numVCH ){

      if ( not contains(plans,la,lb,lc,ld) ){
         add(la,lb,lc,ld);
      }
      Plan_t* tmp = &plans[SERIAL_L4(la,lb,lc,ld)];
      (*plan)     = &std::get<0>(*tmp);
      (*pr_size) = std::get<1>(*tmp);
      (*sh_size) = std::get<2>(*tmp);
      (*numV)    = std::get<3>(*tmp);
      (*numVC)   = std::get<4>(*tmp);
      (*numVCH)  = std::get<5>(*tmp);
   }

   void add( const int la, const int lb, const int lc, const int ld ){
      int t1,t2,numV,numVC,numVCH;
      std::vector<int> tmp = plan(la,lb,lc,ld,&t1,&t2,&numV,&numVC,&numVCH);
      Plan_t tmp_plan = Plan_t(tmp, t1, t2, numV, numVC, numVCH);
      plans[SERIAL_L4(la,lb,lc,ld)] = tmp_plan;
   }
};





#endif // PLAN_H_COMPILE_GUARD
