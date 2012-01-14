"use strict";

var FlatArray = require("./flatarray.js").ChunkSet,
    Octree = require("./octree.js").ChunkSet,
    VirtualArray = require("./virtual.js").ChunkSet,
    IntervalTree = require("./intervaltree.js").ChunkSet;

//Since flatarray, virtual and octrees don't support ranges, we add a generic naive foreach to each of them
function naiveForeach(lo, hi, r, cb) {
  var nmod = 2*r + 1,
      size = nmod * nmod * nmod,
      nbhd = new Array(size),
      i, j, k, l, o, p, q;
      
  for(i=0; i<size; ++i) {
    nbhd[i] = 0;
  }
  
  for(i=lo[0]; i<hi[0]; ++i) {
    for(j=lo[1]; j<hi[1]; ++j) {
      for(k=lo[2]; k<hi[2]; ++k) {
        l = 0;
        for(o=0; o<nmod; ++o) {
          for(p=0; p<nmod; ++p) {
            for(q=0; q<nmod; ++q) {
              nbhd[l++] = this.get(i + o - r, j + p - r, k + q - r);
            }
          }
        }
        //Process neighborhood
        cb(i,j,k,nbhd,1);
      }
    }
  }
}

FlatArray.prototype.rangeForeach = naiveForeach;
Octree.prototype.rangeForeach = naiveForeach;
VirtualArray.prototype.rangeForeach = naiveForeach;

//Layered world
function initializeLayeredTerrain(voxels, radius, depth, layers) {
  for(var l=0; l<layers; ++l) {
    for(var z=l*depth; z<(l+1)*depth; ++z) {
      for(var x=0; x<radius; ++x) {
        for(var y=0; y<radius; ++y) {
          voxels.set(x,y,z, l+1);
        }
      }
    }
  }
}

//Adds speckly noise to the terrain
function randomWrite(voxels, radius, count, max_val) {
  for(var i=0; i<count; ++i) {
    voxels.set(Math.floor(Math.random() * radius),
      Math.floor(Math.random() * radius), 
      Math.floor(Math.random() * radius), 
      Math.floor(Math.random() * max_val));
  }
}

//Random read test
function testRandom(voxels, radius, count) {
  //Need some dummy variable to keep optimizer from killing loop
  var s = 0;
  for(var i=0; i<count; ++i) {
    s += voxels.get(Math.floor(Math.random() * radius),
      Math.floor(Math.random() * radius), 
      Math.floor(Math.random() * radius));
  }
  return s;
}

//Sequential read test
function testSequential(voxels, rangelo, rangehi, n, cases) {
  var count = 0;
  for(var i=0; i<cases; ++i) {
    voxels.rangeForeach(rangelo, rangehi, n, function(i,j,k,nbhd,l) {
      //do nothing
      ++count;
    });
  }
  return count;
}

//Generic benchmark
function runBenchmark(voxels, nrandom, nsequential, max_r) {

  //Initialize terrain
  console.log("Initializing terrain...");
  initializeLayeredTerrain(voxels, [0, 0, 0], [256, 256, 256], 4, 32);
  randomWrite(voxels, 256, 65536, 32);
  
  //Do random read test
  console.log("Starting random test...");
  var start = (new Date).getTime();
  testRandom(voxels, 256, nrandom);
  var elapsed = (new Date).getTime() - start;
  console.log("Total time = " + elapsed + "ms, average time = " + (elapsed / nrandom) + "ms/read");

  //Do sequential test
  for(var r=0; r<=max_r; ++r) {
    console.log("Starting sequential read test, Moore radius = " + r);
    start = (new Date).getTime();
    testSequential(voxels, [r+1,r+1,r+1], [255-r, 255-r, 255-r], r, nsequential)
    elapsed = (new Date).getTime() - start;
    var nreads = (2*r+1)*(2*r+1)*(2*r+1)*nsequential*(254-2*r)*(254-2*r)*(254-2*r);
    console.log("Total time = " + elapsed + "ms, average time = " + (elapsed / nsequential) + "ms/iter, amortized read time = " + (elapsed / nreads) + "ms/read");
  }
}


//Test configuration parameters
var nrandom = (1<<23),
    nsequential = 10,
    max_r = 1;

console.log("Testing flat array");
runBenchmark(new FlatArray([256,256,256]), nrandom, nsequential, max_r);

console.log("Testing octree");
runBenchmark(new Octree(8), nrandom, nsequential, max_r);

console.log("Testing virtual array");
runBenchmark(new VirtualArray(), nrandom, nsequential, max_r);

console.log("Testing interval tree");
runBenchmark(new IntervalTree(), nrandom, nsequential, max_r);


