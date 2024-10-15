#ifndef BSUMTREE_H
#define BSUMTREE_H

#include <valarray>
#include <assert.h>
#include <iostream>

struct s_bsumtree {
  void init(long nrelements_) { // Must be called before use
    // treesize is the smallest power of two above nrelements minus 1
    nrelements = nrelements_;
    if (nrelements == 0)
      treesize = 0;
    else
      treesize = (long) pow(2,ceil(log((double) nrelements)/log((double) 2)))-1;
    
    //std::cout << "nrelements = " << nrelements << "   treesize = " << treesize << std::endl;

    dirty_array.resize(treesize);
	  if ((long) element_array.size()<nrelements) { element_array.resize(nrelements); }
    partsum_array.resize(treesize);
    
    // Initialize arrays
    for (long i=0;i<treesize;i++) {
      dirty_array[i] = false;
      partsum_array[i] = 0;
    }
    for (long i=0;i<nrelements;i++) {
      element_array[i] = 0;
    }
  }
  void setelement(long i, double value) { // 0 <= i < nrelements
    element_array[i] = value;
    long j = i+treesize;
    j = div(j-1,(long) 2).quot; // Parent node

    while (!dirty_array[j]) { // Mark this node and all parents dirty if not already
      dirty_array[j] = true;
      if (j != 0) { // Make sure we stop at the root node
        j = div(j-1,(long) 2).quot; // Parent node
      }
    }
  }
  double getelement(long i) {
    return element_array[i];
  }
  double compute_sum() { // Returns total sum of all elements

    //if (treesize == 0)
    //  return 0;

    // recursively recompute all dirty nodes
    long i = 0; // Start at root node
    while (dirty(i)) {
      if (dirty(2*i + 1)) { // Is left subtree dirty ?
        i = 2*i + 1;
      }
      else {
        if (dirty(2*i + 2)) { // Is right subtree dirty?
          i = 2*i + 2;
        }
        else { // Both subtrees are clean, update this node
          partsum_array[i] = partsum(2*i+1) + partsum(2*i+2);
          dirty_array[i] = false;
          if (i != 0) { // Make sure we stop at the root node
            i = div(i-1,(long) 2).quot; // Parent node
          }
        }
      }
    }
    return partsum_array[0];
  }

  // Search returns index to element i: sum(0..i) <= searchkey < sum(0..i+1),
  // where the sum is taken over the succesive elements.
  long search(double searchkey) { // Returns index to element
    long maxindex = treesize + nrelements;
    long i = 0; // value must be located in subtree denoted by index i

    //while (2*i+2 < maxindex) {      BUG: if 2*i+1 is the last element of partsum then the while
                                  //    loop should be entered and the first if statement will be true
    while (2*i+1 < maxindex) {
      if (searchkey <= partsum(2*i+1)) { // value is located in left subtree
        i = 2*i+1;
      }
      else { // value is located in right subtree
        searchkey -= partsum(2*i+1); // values are relative
        i = 2*i+2;
      }
    }
    i -= treesize;
    return i;
  }

  void resize(long newsize) { // Resize arrays. Expensive, so use with care!
    /*
     *  When newsize >= oldsize: all elements are copied, new elements are 0.
     *  When newsize < oldsize: excess elements are thrown away.
     */
    std::valarray<double> temp_element_array(double(0),newsize); // Temporary storage
    for (long i=0;i<nrelements && i<newsize;i++) {
      temp_element_array[i] = element_array[i];
    }
    init(newsize);
    for (long i=0;i<newsize;i++) {
      setelement(i, temp_element_array[i]);
    }
  }

  long getnrelements() {
    return nrelements;
  }

private:
  bool dirty(long i) {
    if (i<treesize) {
      return dirty_array[i];
    }
    else {
      return false; // Nodes outside the partial rate sum tree are always clean
    }
  }
  double partsum(long i) {
    if (i < treesize) {
      return partsum_array[i];
    }
    else {
      if (i<treesize + nrelements) {
        return element_array[i-treesize];
      }
      else {
        return 0; // Non-existent nodes have partial rate sum equal to 0
      }
    }
  }
  std::valarray<bool> dirty_array; // Are the subtrees dirty?
  std::valarray<double> element_array; // The elements (summands)
  std::valarray<double> partsum_array; // Array of partial sums
  long treesize;
  long nrelements;
};

#endif