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

#include <array>
#include "plan.h"
#include "define.h"

bool ArrayOpLess( const OpArray& a, const std::array<int,OP_SIZE>& b ){
   return a[H__OFFSET] < b[H__OFFSET];
}

inline bool ArrayOpGreater( const OpArray& a, const std::array<int,OP_SIZE>& b ){
   if (a[H__OFFSET] > b[H__OFFSET]) return true;
   if (a[H__OFFSET] < b[H__OFFSET]) return false;
   if (a[M__OFFSET] < b[M__OFFSET]) return true;
   return false;
}

inline int NLco( const int L ){ return (L+1)*(L+2)/2; }


int mem( const int la, const int lb, const int lc, const int ld){
   return NLco(la)*NLco(lb)*NLco(lc)*NLco(ld);
}
int mem( const Op& sh ){
   return NLco(sh.la)*NLco(sh.lb)*NLco(sh.lc)*NLco(sh.ld);
}

bool contains(const OpDict* G, Op op ) { return G->count(op) > 0; }
bool contains(const OpDict& G, Op op ) { return G.count(op) > 0; }


bool contains(const std::unordered_map< int, Plan_t >& PC , const int la, const int lb, const int lc, const int ld ){
  return PC.count(SERIAL_L4(la,lb,lc,ld)) > 0;
}

void rec_add( int h, int la, int lb, int lc, int ld, int m, OpDict *G ){

   if ( la+lb+lc+ld == 0 ) return;

   Op tmp; tmp.la=la; tmp.lb=lb; tmp.lc=lc; tmp.ld=ld; tmp.m=m; 

   // Add Second Horizontal Recursive Relation (a,b|c,d) = (a+1,b-1,c,d) + AB (a,b-1,c,d)
   if ( lb > 0 ){
      tmp.t = HRR2;
      (*G)[tmp] = h;
      rec_add( h+1, la+1,lb-1,lc  ,ld  ,0  ,G);
      rec_add( h+1, la  ,lb-1,lc  ,ld  ,0  ,G);
      return;
   }
   // Add First Horizontal Recursive Relation (a,b|c,d) = (a,b,c+1,d-1) + CD (a,b,c,d-1)
   if ( ld > 0 ){
      tmp.t = HRR1;
      (*G)[tmp] = h;
      rec_add( h+1, la  ,lb  ,lc+1,ld-1,0  ,G);
      rec_add( h+1, la  ,lb  ,lc  ,ld-1,0  ,G);
      return;
   }
   // find which VRR should be used
   if (lc == 0){
      tmp.t = VRR1;
      if (la > 1) tmp.t = VRR2;
   } else {
      tmp.t = VRR3;
      if (lc > 1){
         tmp.t = VRR4;
         if (la > 0) tmp.t = VRR5;
      } else {
         if (la > 0) tmp.t = VRR6;
      }
   }
   // update the urgency of the op
   if ( not contains(G,tmp) ){
      (*G)[tmp] = h;
   } else {
      (*G)[tmp] = std::max( (*G)[tmp], h );
      return;
   }
   int t = tmp.t;
   // Apply the correct VRR
   if (t == VRR1){
      rec_add( h+1, la-1,lb,lc-0,ld, m+0, G );
      rec_add( h+1, la-1,lb,lc-0,ld, m+1, G );
      return;
   }
   if (t == VRR2){
      rec_add( h+1, la-1,lb,lc-0,ld, m+0, G );
      rec_add( h+1, la-1,lb,lc-0,ld, m+1, G );
      rec_add( h+1, la-2,lb,lc-0,ld, m+0, G );
      rec_add( h+1, la-2,lb,lc-0,ld, m+1, G );
      return;
   }
   if (t == VRR3){
      rec_add( h+1, la-0,lb,lc-1,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-1,ld, m+1, G );
       return;
   }
  if (t == VRR4){
      rec_add( h+1, la-0,lb,lc-1,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-1,ld, m+1, G );
      rec_add( h+1, la-0,lb,lc-2,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-2,ld, m+1, G );
      return;
   }
   if (t == VRR5){
      rec_add( h+1, la-0,lb,lc-1,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-1,ld, m+1, G );
      rec_add( h+1, la-0,lb,lc-2,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-2,ld, m+1, G );
      rec_add( h+1, la-1,lb,lc-1,ld, m+1, G );
      return;
   }
   if (t == VRR6){
      rec_add( h+1, la-0,lb,lc-1,ld, m+0, G );
      rec_add( h+1, la-0,lb,lc-1,ld, m+1, G );
      rec_add( h+1, la-1,lb,lc-1,ld, m+1, G );
      return;
   }
}


bool is_vrr(const int t){
   return t==VRR1 or t==VRR2 or t==VRR3 or t==VRR4 or t==VRR5 or t==VRR6;
}



bool CompareOp( const OpArray& op1, const OpArray& op2 ){
   // returns true if op1 goes before op2, false otherwise.
   int t1 = op1[T__OFFSET];
   int a1 = op1[LA_OFFSET]; 
   int b1 = op1[LB_OFFSET]; 
   int c1 = op1[LC_OFFSET]; 
   int d1 = op1[LD_OFFSET];
//   int m1 = op1[M__OFFSET];
   int h1 = op1[H__OFFSET];

   int t2 = op2[T__OFFSET];
   int a2 = op2[LA_OFFSET]; 
   int b2 = op2[LB_OFFSET]; 
   int c2 = op2[LC_OFFSET]; 
   int d2 = op2[LD_OFFSET];
//   int m2 = op2[M__OFFSET];
   int h2 = op2[H__OFFSET];

   int l1 = a1+b1+c1+d1;
   int l2 = a2+b2+c2+d2;

   const bool op1_first = true  ;
   const bool op2_first = false ;

   // vrr1 first, then vrr2 and vrr3. Any other op is always after them
   if (t1 == VRR1 and t2 != VRR1 ){ return op1_first; }
   if (t1 != VRR1 and t2 == VRR1 ){ return op2_first; }
   if (t1 == VRR1 and t2 == VRR1 ){ return  l1 < l2 ; }
   // no more vrr1. 
   if (t1 == VRR2 and t2 != VRR2 ){ return op1_first; }
   if (t1 != VRR2 and t2 == VRR2 ){ return op2_first; }
   if (t1 == VRR2 and t2 == VRR2 ){ return  l1 < l2 ; }

   if (t1 == VRR3 and t2 != VRR3 ){ return op1_first; }
   if (t1 != VRR3 and t2 == VRR3 ){ return op2_first; }
   if (t1 == VRR3 and t2 == VRR3 ){ return  l1 < l2 ; }

   bool is1_vrr = is_vrr(t1);
   bool is2_vrr = is_vrr(t2);

   // prefer vrrs over hrrs
   if ( is1_vrr and not is2_vrr) { return op1_first; }
   if ( is2_vrr and not is1_vrr) { return op2_first; }

   // if they are both hrr, return order by h. Higher h first
   if ( ( not is1_vrr ) and ( not is2_vrr ) ){ return h1 > h2; }

   // now we know we only have vrr4,vrr5,vrr6
   // order by moment, lower moment first
   return l1 < l2;
}

bool contains ( const std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int> > > >& shared_addr, int la, int lb, int lc, int ld ){
   return shared_addr.count(la) == 1 and shared_addr.at(la).count(lb) == 1 and shared_addr.at(la).at(lb).count(lc) == 1 and shared_addr.at(la).at(lb).at(lc).count(ld) == 1;
}


bool contains ( const std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int> > >& priv_addr, int la, int lc, int m ){
   return priv_addr.count(la) == 1 and priv_addr.at(la).count(lc) == 1 and priv_addr.at(la).at(lc).count(m) == 1;
}


std::vector<int> plan( const int LA, const int LB, const int LC, const int LD, int* pr_size, int* sh_size, int* numV, int* numVC, int* numVCH ){
   int L = LA + LB + LC + LD;

//   cout << "Planning " << LA << " " << LB << " " << LC << " " << LD << endl; cout.flush();

   // given the target angular moments la,lb,lc,ld prepare a plan of execution
   // the plan at the end will be a vector of ( array of int of size OP_SIZE )
   // each element describe how to compute one of the intermidiate product (or the final element)
   // each element of the plan contains 
   // 0    -> type of operation. Can be HRR, VRR or (Contract from Private To Shared memory) CP2D
   // 1-4  -> l of the target
   // 5    -> m of the target. Or maximum m of the target
   // 6    -> (order of priority) h. Or minimum m of the target
   // 7-12 -> memory offsets.
   // For example, if op = [6, 2 0 1 0 2, 7, 405 45 51 18 0 0]
   // We need to perform a VRR of type 6, with the target AL=2, CL=1, m=2 (BL and DL are always 0 for a vrr6 op)
   // this means dsps(2) = QC dsss(2) + WQ dsss(3) + 1/2z psss(3)
   // the output is at pr_mem+405, the inputs are at pr_mem+45,51 and 18

   // We use 3 main memory structure for the plan:
   // unordered map so that we can easily add and check if on op already exist 
   OpDict G = OpDict();
   // a list of OpDict so that we can easily insert new ops in any place 
   std::list<OpArray> G_list;
   std::list<OpArray> G_list2;
   // a vector of int (declared later, once we knwow its size)
   // std::vector<int> G_vec;
   //
   // We also use temporary object to store operations and memory addresses
   // temporary holders for an operation. Also used for memory addresses lookup
   Op tmp;
   OpArray tmp_op;
   // max and min m associated to a given la,lc
   std::unordered_map<int, std::unordered_map<int, int> > vrr_min_m;
   std::unordered_map<int, std::unordered_map<int, int> > vrr_max_m;
   // global/register/per thread memory
   std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int> > > priv_addr;
   // Addresses in shared memory for the HRR and output of CP2S
   std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, int> > > > shared_addr;

   int curr_pr_mem_loc = 0;
   int curr_sh_mem_loc = 0;

   // always insert the final results in the zero address
   shared_addr[LA][LB][LC][LD] = 0;
   curr_sh_mem_loc += mem(LA,LB,LC,LD);

   // make it so that the address for ssssm values starts at 0 and are consecutive
   for( int m=0; m < L+1; m++ ){
      priv_addr[0][0][m] = curr_pr_mem_loc;
      curr_pr_mem_loc += 1;
   }

   // create the tree of operations
   rec_add( 0, LA, LB, LC, LD, 0, &G );

   // change from unordered map to list
   for ( const auto& Op_h : G ) {
      tmp = Op_h.first;
      int h = Op_h.second;
      tmp_op[T__OFFSET] = tmp.t;
      tmp_op[LA_OFFSET] = tmp.la;
      tmp_op[LB_OFFSET] = tmp.lb;
      tmp_op[LC_OFFSET] = tmp.lc;
      tmp_op[LD_OFFSET] = tmp.ld;
      tmp_op[M__OFFSET] = tmp.m;
      tmp_op[H__OFFSET] = h;
      G_list.push_back(tmp_op);
   }

   // sort the order of operations by the value of h
   G_list.sort( CompareOp );

   // find the max and min m for all a0c0 operations
   int n_ops = int( G_list.size() );
   for ( std::list<OpArray>::const_iterator it = G_list.begin(); it != G_list.end(); it++ ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int m  = tmp_op[M__OFFSET];

      if ( is_vrr(t) ){
          if ( vrr_min_m.count(la) == 1 and vrr_min_m[la].count(lc) == 1 ){
             vrr_min_m[la][lc] = std::min(vrr_min_m[la][lc], m);
             vrr_max_m[la][lc] = std::max(vrr_max_m[la][lc], m);
          } else {
             vrr_min_m[la][lc] = m;
             vrr_max_m[la][lc] = m;
          }
//          printf(" la %d lc %d min %d max %d \n", la, lc, vrr_min_m[la][lc], vrr_max_m[la][lc] );
      }
   }

   // compress G by collapsing vrr operations with the same la and lc, but different m
   n_ops = int( G_list.size() );
   for ( std::list<OpArray>::const_iterator it = G_list.begin(); it != G_list.end(); it++ ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int m  = tmp_op[M__OFFSET];
      if ( is_vrr(t) ){
         if ( m == vrr_max_m[la][lc] ){
            tmp_op[H__OFFSET] = vrr_min_m[la][lc]; // note, we reuse h to store minimum m
            G_list2.push_back( tmp_op );
         } // else pass and DO NOT ADD THE OP
      } else {
         G_list2.push_back( tmp_op );
      }
   }

   G_list = G_list2 ;

   // assign memory locations to the operations
   // In the case of the HRRs, we add both input and output to the 
   // memory locations in shared memory. In the case of the VRR
   // this is not necessary, we only need to add the output to the memory
   n_ops = int( G_list.size() );
   for ( std::list<OpArray>::const_iterator it = G_list.begin(); it != G_list.end(); it++ ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int ld = tmp_op[LD_OFFSET];
      int m  = tmp_op[M__OFFSET];
      int h  = tmp_op[H__OFFSET];
      if (t == HRR1){
         if ( !contains(shared_addr, la,lb,lc  ,ld  ) ){ shared_addr[la][lb][lc  ][ld  ] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la,lb,lc  ,ld  ); }
         if ( !contains(shared_addr, la,lb,lc+1,ld-1) ){ shared_addr[la][lb][lc+1][ld-1] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la,lb,lc+1,ld-1); }
         if ( !contains(shared_addr, la,lb,lc  ,ld-1) ){ shared_addr[la][lb][lc  ][ld-1] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la,lb,lc  ,ld-1); }
      } else if (t == HRR2) {
         if ( !contains(shared_addr, la  ,lb  ,lc,ld) ){ shared_addr[la  ][lb  ][lc][ld] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la  ,lb  ,lc,ld); }
         if ( !contains(shared_addr, la+1,lb-1,lc,ld) ){ shared_addr[la+1][lb-1][lc][ld] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la+1,lb-1,lc,ld); }
         if ( !contains(shared_addr, la  ,lb-1,lc,ld) ){ shared_addr[la  ][lb-1][lc][ld] = curr_sh_mem_loc; curr_sh_mem_loc += mem(la  ,lb-1,lc,ld); }
      } else { // if is vrr
         int max_m = m; // added for clarity
         int min_m = h; // note, we reused h to store minimum m
         if ( !contains(priv_addr, la,lc,0) ){
            for ( int mm=min_m; mm <= max_m; mm++){ // long proof, but the m are always continuous from min_m to max_m
               priv_addr[la][lc][mm] = curr_pr_mem_loc;
               curr_pr_mem_loc += mem(la,0,lc,0); // when done this way, pr data with the same a c but different m are countiguous in memory
            }
         }
      }
   }
   // save the amount of memory needed in the /output/ variables pr_size and sh_size
   (*pr_size) = curr_pr_mem_loc;
   (*sh_size) = curr_sh_mem_loc;

   // if we only require the s integrals, the op list is empty.
   // We need to manually add the sum of the K contributions to the final integral
   tmp_op[T__OFFSET] = CP2S;
   tmp_op[LA_OFFSET] = 0;
   tmp_op[LB_OFFSET] = 0;
   tmp_op[LC_OFFSET] = 0;
   tmp_op[LD_OFFSET] = 0;
   tmp_op[M__OFFSET] = 0;
   tmp_op[H__OFFSET] = 0;
   if (n_ops == 0){
      G_list.insert( G_list.begin(), tmp_op );
   } else {
      // in the next step, we are going to add reductions from the VRR.
      // Since 0000 is never computed from a VRR, although it might be needed,
      // we check if it is needed (aka is on shared addr) and manually add the reduction
      if ( contains(shared_addr, 0, 0, 0, 0) ){
         G_list.insert( G_list.begin()++, tmp_op );
      }
   }

   // Add a CP2S operation after an element if ((is computed using a VRR) and (it is needed for a HRR))
   OpArray tmp_op2;
   for ( auto it = G_list.begin() ; it != G_list.end(); it++ ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int ld = tmp_op[LD_OFFSET];
//      int max_m = tmp_op[M__OFFSET];
      int min_m = tmp_op[H__OFFSET];
      bool is_op_vrr = is_vrr(t); 
      bool move_pr_to_sh = ( (min_m == 0) and (contains(shared_addr,la,lb,lc,ld)) );
      if ( is_op_vrr and move_pr_to_sh ){
         tmp_op2 = tmp_op;
         tmp_op2[T__OFFSET] = CP2S;
         G_list.insert(++it, tmp_op2 );// TODO ++it, followed by it-- ?
         it--;
      }
   }

   // add a sync team operation before a VRR operation if its total moment is different from the one of the vrr preceding it
   int curr_l = 1;
   for ( std::list<OpArray>::iterator it = G_list.begin() ; it != G_list.end(); ++it ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int ld = tmp_op[LD_OFFSET];
//      int m  = tmp_op[M__OFFSET];
//      int h  = tmp_op[H__OFFSET];
      int l  = la + lb + lc + ld ;
      if ( is_vrr(t) and l != curr_l ){
         tmp_op2 = tmp_op;
         tmp_op2[T__OFFSET] = SYTM;
         G_list.insert(it, tmp_op2 );
         curr_l = l;
      }
   }

   // add a sync block operation before a HRR operation if its priority is different from the one of the hrr preceding it
   // note: this is because (1) the hrr are sorted by already priority in G_list (2) hrr with the same priority can always be done in whatever order with no sync
   int curr_h = -1; // special value to mark first h
   for ( std::list<OpArray>::iterator it = G_list.begin() ; it != G_list.end(); ++it ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
//      int la = tmp_op[LA_OFFSET];
//      int lb = tmp_op[LB_OFFSET];
//      int lc = tmp_op[LC_OFFSET];
//      int ld = tmp_op[LD_OFFSET];
//      int m  = tmp_op[M__OFFSET];
      int h  = tmp_op[H__OFFSET];
      if ( t == HRR1 or t == HRR2 ){
         if ( curr_h == -1 ) { curr_h = h; }
         if ( h != curr_h ){
            tmp_op2 = tmp_op;
            tmp_op2[T__OFFSET] = SYBL;
            G_list.insert(it, tmp_op2 );
            curr_h = h;
         }
      }
   }

   // Convert the list to the final output vector
   // Reads the t,a,b,c,d,m,h from the list element.
   // Then does a lookup of the memory addresses for the given operation
   // and add them to the data in the vector
   n_ops = int( G_list.size() );
   std::vector<int> G_vec((n_ops+2)*OP_SIZE);
   int i = 0;
   // part1 : ALL VRRs and sync team
   for ( std::list<OpArray>::iterator it = G_list.begin() ; it != G_list.end(); ++it ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int ld = tmp_op[LD_OFFSET];
      int M  = tmp_op[M__OFFSET]; // max_m
      int m  = tmp_op[H__OFFSET]; // min_m
      if ( is_vrr(t) or t == SYTM ){
         G_vec[i*OP_SIZE + T__OFFSET] = t  ; // tmp_op[T__OFFSET];
         G_vec[i*OP_SIZE + LA_OFFSET] = la ; // tmp_op[LA_OFFSET];
         G_vec[i*OP_SIZE + LB_OFFSET] = lb ; // tmp_op[LB_OFFSET];
         G_vec[i*OP_SIZE + LC_OFFSET] = lc ; // tmp_op[LC_OFFSET];
         G_vec[i*OP_SIZE + LD_OFFSET] = ld ; // tmp_op[LD_OFFSET];
         G_vec[i*OP_SIZE + M__OFFSET] = M  ; // tmp_op[M__OFFSET];
         G_vec[i*OP_SIZE + H__OFFSET] = m  ; // tmp_op[H__OFFSET];
         // sorry for the block of code. For each op type, find the memory addresses needed
         // and add them in the G vector in the correct positions
         // Note that some memory in G will remain unused. We do so to keep the size of op uniform
         if (t == VRR1){
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la-1).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la-1).at(lc  ).at(m+1);
         } else if ( t == VRR2 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la-1).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la-1).at(lc  ).at(m+1);
            G_vec[i*OP_SIZE + M4_OFFSET] = priv_addr.at(la-2).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M5_OFFSET] = priv_addr.at(la-2).at(lc  ).at(m+1);
         } else if ( t == VRR3 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m+1);
         } else if ( t == VRR4 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m+1);
            G_vec[i*OP_SIZE + M4_OFFSET] = priv_addr.at(la  ).at(lc-2).at(m  );
            G_vec[i*OP_SIZE + M5_OFFSET] = priv_addr.at(la  ).at(lc-2).at(m+1);
         // VRR5 requires 6 memory addresses. It (currently) sets the value of OP_SIZE
         } else if ( t == VRR5 ){         
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m+1);
            G_vec[i*OP_SIZE + M4_OFFSET] = priv_addr.at(la  ).at(lc-2).at(m  );
            G_vec[i*OP_SIZE + M5_OFFSET] = priv_addr.at(la  ).at(lc-2).at(m+1);
            G_vec[i*OP_SIZE + M6_OFFSET] = priv_addr.at(la-1).at(lc-1).at(m+1);
         } else if ( t == VRR6 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la  ).at(lc  ).at(m  );
            G_vec[i*OP_SIZE + M2_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m  );
            G_vec[i*OP_SIZE + M3_OFFSET] = priv_addr.at(la  ).at(lc-1).at(m+1);
            G_vec[i*OP_SIZE + M4_OFFSET] = priv_addr.at(la-1).at(lc-1).at(m+1);
         } else if ( t == SYTM ){
            // We don't need any other info
         }
      i++;
      }
   }
   // saves number of vvr and team sync to output variable numV
   (*numV) = i;

   // add a sync team?block? after the vrr, before the cp2s // almost sure is not needed
   G_vec[i*OP_SIZE + T__OFFSET] = SYBL; i++;


   // now all cp2s
   for ( std::list<OpArray>::iterator it = G_list.begin() ; it != G_list.end(); ++it ){
      tmp_op = (*it);
      int t = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
//      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
//      int ld = tmp_op[LD_OFFSET];
//      int m  = tmp_op[M__OFFSET];
//      int h  = tmp_op[H__OFFSET];
      if ( t == CP2S ){
         G_vec[i*OP_SIZE + T__OFFSET] = tmp_op[T__OFFSET];
         G_vec[i*OP_SIZE + LA_OFFSET] = tmp_op[LA_OFFSET];
         G_vec[i*OP_SIZE + LB_OFFSET] = tmp_op[LB_OFFSET];
         G_vec[i*OP_SIZE + LC_OFFSET] = tmp_op[LC_OFFSET];
         G_vec[i*OP_SIZE + LD_OFFSET] = tmp_op[LD_OFFSET];
         G_vec[i*OP_SIZE + M__OFFSET] = tmp_op[M__OFFSET];
         G_vec[i*OP_SIZE + H__OFFSET] = tmp_op[H__OFFSET];
         G_vec[i*OP_SIZE + M1_OFFSET] = priv_addr.at(la).at(lc).at(0);
         G_vec[i*OP_SIZE + M2_OFFSET] = shared_addr.at(la).at(0).at(lc).at(0);
         i++;
      }
   }
   // saves number of vvr, team sync and cp2s to output variable numVC
   (*numVC) = i;

   // add a sync block after the cp2s, before the hrr
   G_vec[i*OP_SIZE + T__OFFSET] = SYBL; i++;

   for ( std::list<OpArray>::iterator it = G_list.begin() ; it != G_list.end(); ++it ){
      tmp_op = (*it);
      int t  = tmp_op[T__OFFSET];
      int la = tmp_op[LA_OFFSET];
      int lb = tmp_op[LB_OFFSET];
      int lc = tmp_op[LC_OFFSET];
      int ld = tmp_op[LD_OFFSET];
//      int m  = tmp_op[M__OFFSET];
//      int h  = tmp_op[H__OFFSET];
      if ( t == HRR1 or t == HRR2 or t == SYBL ){
         G_vec[i*OP_SIZE + T__OFFSET] = tmp_op[T__OFFSET];
         G_vec[i*OP_SIZE + LA_OFFSET] = tmp_op[LA_OFFSET];
         G_vec[i*OP_SIZE + LB_OFFSET] = tmp_op[LB_OFFSET];
         G_vec[i*OP_SIZE + LC_OFFSET] = tmp_op[LC_OFFSET];
         G_vec[i*OP_SIZE + LD_OFFSET] = tmp_op[LD_OFFSET];
         G_vec[i*OP_SIZE + M__OFFSET] = tmp_op[M__OFFSET];
         G_vec[i*OP_SIZE + H__OFFSET] = tmp_op[H__OFFSET];
         if ( t == HRR1 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = shared_addr.at(la).at(lb).at(lc  ).at(ld  );
            G_vec[i*OP_SIZE + M2_OFFSET] = shared_addr.at(la).at(lb).at(lc+1).at(ld-1);
            G_vec[i*OP_SIZE + M3_OFFSET] = shared_addr.at(la).at(lb).at(lc  ).at(ld-1);
         } else if ( t == HRR2 ){
            G_vec[i*OP_SIZE + M1_OFFSET] = shared_addr.at(la  ).at(lb  ).at(lc).at(ld);
            G_vec[i*OP_SIZE + M2_OFFSET] = shared_addr.at(la+1).at(lb-1).at(lc).at(ld);
            G_vec[i*OP_SIZE + M3_OFFSET] = shared_addr.at(la  ).at(lb-1).at(lc).at(ld);
         }
         i++;
      }
   }
   (*numVCH) = i;
   //
   return G_vec;
}
