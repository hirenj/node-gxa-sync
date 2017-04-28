"use strict";
// ES6

const request = require('request');

const experiments_url = 'http://www.ebi.ac.uk/gxa/json/experiments';

const ftp = require('ftp');
const parse_url = require('url').parse;
const csv = require('csv-parse');
const xmlNodes = require('xml-nodes');
const xmlObjects = require('xml-objects');
const fs = require('fs');
const converter = require('node-uberon-mappings');
const CheapJSON = require('./output').CheapJSON;
const temp = require('temp');
const path = require('path');

const baseline_whitelist = [ 'E-MTAB-2836', 'E-MTAB-5214','E-MTAB-4344' ];

const stream = require('stream');
const util = require('util');
const Transform = stream.Transform;

const nconf = require('nconf');

nconf.argv();

function EntryTransform(ontology,gene,config) {
  if (!(this instanceof EntryTransform)) {
    return new EntryTransform(mappings);
  }
  Transform.call(this, {objectMode: true});
  this.ontology = ontology;
  this.gene = gene;
  this.config = config;
};

util.inherits(EntryTransform, Transform);

temp.track();

EntryTransform.prototype._transform = function(dat,enc,cb) {
  let gene_id = this.gene[dat['Gene ID']];
  delete dat['Gene ID'];
  delete dat['Gene Name'];
  if ( ! gene_id ) {
    cb();
    return;
  }
  let entries = Object.keys(dat).map( sample => {
    let term_info = this.ontology[this.config[sample].sample] || {};
    let vals = dat[sample].split(',').map( val => +val );
    if ( ! term_info.uberon ) {
      console.log(sample,this.config[sample]);
    }
    return {
      gene_id: gene_id,
      uberon: term_info.uberon,
      simple_tissue: term_info.simple_tissue,
      simple_uberon: term_info.simple_uberon,
      exp: vals[2],
      annotation: { exp: vals } };
  });
  this.push( [ gene_id, entries ]);
  cb();
};

let getExperiments = () => {
  return new Promise(function(resolve,reject) {
    request(experiments_url,function(err,response,body) {
      if ( err ) {
        return resolve(err);
      }
      let parsed = JSON.parse(body);
      if ( ! parsed.aaData) {
        return resolve(new Error('Invalid data format'));
      }
      resolve(parsed.aaData);
    });
  });
};

let filterSpecies = (species,experiments) => {
  let species_set = new Set(species);
  return experiments.filter((exp) => {
    return [].concat((exp.species || [])).filter((tax) => species_set.has(tax)).length > 0;
  });
};

let filterBaseline = (experiments) => {
  return experiments.filter((exp) => {
    return exp.experimentType === 'RNASEQ_MRNA_BASELINE';
  }).filter( exp => {
    return (baseline_whitelist.length == 0) || baseline_whitelist.indexOf(exp.experimentAccession) >= 0;
  });
};

let getExperimentInfo = (exp) => {
  return {
    accession: exp.experimentAccession,
    date: exp.lastUpdate,
    sources: exp.experimentalFactors,
    taxonomy: [].concat(exp.species || []).join(','),
    description: exp.experimentDescription
  };
};

let getBaselineForSpecies = (species) => {
  return getExperiments()
  .then(filterSpecies.bind(null,species))
  .then(filterBaseline)
  .then((exps) => exps.map(getExperimentInfo));
};

let make_ontology_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}.condensed-sdrf.tsv`;
};

let make_configuration_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}-configuration.xml`;
};

let make_data_url = function(experiment_id) {
  return `ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/${experiment_id}/${experiment_id}.tsv`;
};

const read_configuration = function(stream) {
  let result = {};
  let data_stream = stream.pipe(xmlNodes('assay_group'))
  .pipe(xmlObjects({explicitRoot: false, explicitArray: false, mergeAttrs: true}))
  .on('data', dat => {
    result[dat.id] = { label: dat.label, sample: Array.isArray(dat.assay) ? dat.assay.map( assay_val => assay_val['_'] || assay_val )[0] : (dat.assay['_'] || dat.assay) };
  });
  return new Promise( (resolve,reject) => {
    data_stream.on('end',resolve);
    stream.on('error',reject);
    data_stream.on('error',reject);
  }).then( () => result );
};


const read_ontology = function(stream) {
  return new Promise( (resolve,reject) => {
    let parser = csv({delimiter: "\t", relax_column_count: true}, function(err, data){
      if (err) {
        reject(err);
        return;
      }
      let ids = data.filter( row => row[3] == 'characteristic' && row[4] == 'organism part')
      .map( row => { return { ontology: row[6], sample_id: row[2]  } })
      .filter( purl => purl && purl.ontology )
      .map( purl => {
        purl.ontology = purl.ontology.split('/').reverse().shift().replace('_',':');
        return purl;
      });
      resolve(ids);
    });
    stream.pipe(parser);
  });
};

const read_data = function(stream) {
  return new Promise( (resolve,reject) => {
    resolve(stream.pipe(csv({delimiter: "\t", relax_column_count: true, columns: true})));
  });
};

const read_geneids = (function() {
  return new Promise( (resolve,reject) => {
    let stream = fs.createReadStream('gene2ensembl');
    stream.on('error',reject);
    let parser = csv({delimiter: "\t", relax_column_count: true, columns: true}, function(err, data){
      if (err) {
        reject(err);
        return;
      }
      let results = {};
      data.forEach( row => {
        results[row.Ensembl_gene_identifier] = row.GeneID;
      });
      resolve(results);
    });
    stream.pipe(parser);
  });
})();

const connect_ftp = function(url) {
  let req = parse_url(url);
  let ftp_site  = new ftp();
  let result = new Promise(function(resolve,reject) {
    ftp_site.on('ready',function() {
      resolve(ftp_site);
    });
    ftp_site.on('error',reject);
  });
  req.connTimeout = 3000;
  ftp_site.connect(req);
  return result;
};

const download_file = function(url,ftp) {
  let outstream = temp.createWriteStream();
  let out_path = outstream.path;
  let path = parse_url(url).path;
  return new Promise(function(resolve,reject) {
    ftp.get(path, function(err, stream) {
      if (err) {
        reject(err);
      } else {
        stream.pipe(outstream);
        outstream.once('close', () => {
          ftp.end();
          resolve(fs.createReadStream(out_path));
        })
        stream.on('error', err => reject(err));
      }
    });
  });
};

const get_ebi_file = function(url) {
  let filename = path.join(nconf.get('workdir') || '',url.split('/').reverse()[0]);
  console.log(filename);
  if (fs.existsSync(filename)) {
    return Promise.resolve(fs.createReadStream(filename));
  }
  console.log('Downloading',url,'from ftp');
  return connect_ftp(url).then(download_file.bind(null,url));
};

let downloadExpData = function(experiment_id,description) {
  let ontology_url = make_ontology_url(experiment_id);
  console.log(ontology_url);
  let ontology_ids = get_ebi_file(ontology_url)
  .then( read_ontology )
  .then( ontology => {
    let ontology_map = {};
    return Promise.all(ontology.map( map => {
      return converter.convert(map.ontology).then( converted => {
        ontology_map[ map.sample_id ] = {
          'uberon' : map.ontology,
          'simple_tissue' : converted.name,
          'simple_uberon' : converted.root,
        }
      });
    })).then( () => ontology_map );
  });

  let data_url = make_data_url(experiment_id);
  console.log(data_url);
  let data = get_ebi_file(data_url)
  .then( read_data );

  let gene_ids = read_geneids;

  let configuration_url = make_configuration_url(experiment_id);
  console.log(configuration_url);
  let configuration = get_ebi_file(configuration_url)
  .then( read_configuration );

  let writer = new CheapJSON({
    'mimetype' : 'application/json+expression',
    'title' : `${experiment_id} ${description}`,
    'sample' : { 'species': 9606, 'tissue' : 'bto:0001489' },
    'data-version' : ''+nconf.get('version'),
    "software" : {"ARRAY": "true", "0" : { "name" : "hirenj/node-gxa-sync", "version" : nconf.get('git') , "run-date" : nconf.get('timestamp') }},

  });
  let output = fs.createWriteStream(`gxa_${experiment_id}.json`);
  writer.pipe(output);

  return Promise.all([data,ontology_ids,gene_ids,configuration]).then( ids => {
    let transformer = new EntryTransform( ids[1],ids[2],ids[3] );
    return ids[0].pipe(transformer);
  }).then( entry_stream => {
    entry_stream.pipe(writer);
    return new Promise( (resolve,reject) => {
      output.on('finish',resolve);
      output.on('close',resolve);
      output.on('error',reject);
      entry_stream.on('error',reject);
    });
  });
};

const handle_exps = function(exps) {
  if ( exps.length < 1 ) {
    return;
  }
  let curr_exp = exps.shift();
  return downloadExpData(curr_exp.accession,curr_exp.description).then( () => exps ).then( handle_exps );
};

if (nconf.get('print-sets')) {
  getBaselineForSpecies(['Homo sapiens']).then( exps => {
    exps.forEach(exp => {
      console.log(exp.accession);
    });
  })
  .then( () => process.exit(0) )
  .catch((err) => {
    console.error("Error is ",err);
    process.exit(1);
  });
} else {
  let foo = getBaselineForSpecies(['Homo sapiens']).then(handle_exps)
  .then( () => process.exit(0) )
  .catch((err) => {
    console.error("Error is ",err);
    process.exit(1);
  });
}

