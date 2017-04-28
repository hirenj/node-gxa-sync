"use strict";

const stream = require('stream');
const util = require('util');

const Transform = stream.Transform;
const PassThrough = stream.PassThrough;

function CheapJSON(meta,options) {
  if (!(this instanceof CheapJSON)) {
    return new CheapJSON(meta,options);
  }

  if (!options) options = {};
  options.objectMode = true;
  this.meta = meta;
  this.first = false;
  Transform.call(this, options);
}

util.inherits(CheapJSON, Transform);

CheapJSON.prototype._transform = function (obj,enc,cb) {
  let sep = ',';
  if (! this.first) {
    sep = '{\n"data":{';
  }
  this.first = this.first || true;
  this.push(sep+'\n\t"'+obj[0]+'":'+JSON.stringify(obj[1]));
  cb();
};

CheapJSON.prototype._flush = function(cb) {
  this.push('\n},\n"metadata":'+JSON.stringify(this.meta)+'\n}\n');
  cb();
};

exports.CheapJSON = CheapJSON;