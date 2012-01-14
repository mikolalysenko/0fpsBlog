"use strict";

var CHUNK_SHIFT_X = 5,
    CHUNK_SHIFT_Y = 5,
    CHUNK_SHIFT_Z = 5,
    CHUNK_X       = (1 << CHUNK_SHIFT_X),
    CHUNK_Y       = (1 << CHUNK_SHIFT_Y),
    CHUNK_Z       = (1 << CHUNK_SHIFT_Z),
    CHUNK_SIZE    = CHUNK_X * CHUNK_Y * CHUNK_Z,
    CHUNK_MASK_X  = CHUNK_X - 1,
    CHUNK_MASK_Y  = CHUNK_Y - 1,
    CHUNK_MASK_Z  = CHUNK_Z - 1;

function expand(x) {
  x &= 0x3FF;
  x  = (x | (x<<16)) & 4278190335;
  x  = (x | (x<<8))  & 251719695;
  x  = (x | (x<<4))  & 3272356035;
  x  = (x | (x<<2))  & 1227133513;
  return x;
};

function hashCode(i,j,k) {
  return expand(i)+(expand(j)<<1)+(expand(k)<<2);
};


function ChunkSet() {
  this.chunks = [];
  Object.freeze(this);
}

ChunkSet.prototype.get = function(i,j,k) {
  var cx = i >> CHUNK_SHIFT_X,
      cy = j >> CHUNK_SHIFT_Y,
      cz = k >> CHUNK_SHIFT_Z,
      chunk = this.chunks[hashCode(cx, cy, cz)];
  if(chunk) {
    return chunk[ (i&CHUNK_MASK_X) + ((((k&CHUNK_MASK_Z)<< CHUNK_SHIFT_Y)+(j&CHUNK_MASK_Y)) << CHUNK_SHIFT_X)];
  }
  return 0;
}

ChunkSet.prototype.set = function(i,j,k,v) {
  var cx = i >> CHUNK_SHIFT_X,
      cy = j >> CHUNK_SHIFT_Y,
      cz = k >> CHUNK_SHIFT_Z,
      key = hashCode(cx, cy, cz),
      chunk = this.chunks[key];
  if(!chunk) {
    chunk = this.chunks[key] = new Array(CHUNK_SIZE);
    for(var l=CHUNK_SIZE-1; l>=0; --l) {
      chunk[l] = 0;
    }
  }
  return chunk[ (i&CHUNK_MASK_X) + ((((k&CHUNK_MASK_Z)<< CHUNK_SHIFT_Y)+(j&CHUNK_MASK_Y)) << CHUNK_SHIFT_X)] = v;
}

exports.ChunkSet = ChunkSet;

