"use strict";

function ChunkSet(dims) {
  this.voxels = new Array(dims[0] * dims[1] * dims[2]);
  for(var l=dims[0]*dims[1]*dims[2]-1; l>=0; --l) {
    this.voxels[l]=0;
  }
  this.dims = dims;
  Object.freeze(this);
}

ChunkSet.prototype.get = function(x, y, z) {
  return this.voxels[ x + this.dims[0] * (y + this.dims[1] * z)];
}

ChunkSet.prototype.set = function(x, y, z, v) {
  return this.voxels[ x + this.dims[0] * (y + this.dims[1] * z)] = v;
}

exports.ChunkSet = ChunkSet;
