"use strict";

function Octree(height) {
  this.height = height;
  this.children = new Array(8);
  Object.freeze(this);
}

Octree.prototype.get = function(i, j, k) {
  var bit = (1 << this.height),
      idx = ((i & bit) >> this.height) +
            ((j & bit) >> (this.height-1)) +
            ((k & bit) >> (this.height-2));
  if(idx in this.children) {
    if(this.height > 0) {
      return this.children[idx].get(i, j, k);
    }
    return this.children[idx];
  }
  return 0;
}

Octree.prototype.set = function(i,j,k,v) {
  var bit = (1 << this.height),
      idx = ((i & bit) >> this.height) +
            ((j & bit) >> (this.height-1)) +
            ((k & bit) >> (this.height-2));

  if(this.height > 0) {
    if(idx in this.children) {
        return this.children[idx].set(i, j, k, v);
    }
    else {
      this.children[idx] = new Octree(this.height - 1);
      return this.children[idx].set(i,j,k,v);
    }
  }
  else {
    return this.children[idx] = v;
  }
}

exports.ChunkSet = Octree;
